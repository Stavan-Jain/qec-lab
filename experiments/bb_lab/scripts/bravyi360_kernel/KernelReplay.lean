module

import all Lean.Replay
import Lean.Elab.Import
import Lean.Util.CollectAxioms

open Lean

/-- Read every compiled constant from all available object-file parts. -/
unsafe def readConstants (moduleName : Name) : IO (Array ConstantInfo) := do
  let file ← findOLean moduleName
  let mut files := #[file]
  let server := OLeanLevel.server.adjustFileName file
  if ← server.pathExists then
    files := files.push server
    let privateFile := OLeanLevel.private.adjustFileName file
    if ← privateFile.pathExists then files := files.push privateFile
  let parts ← readModuleDataParts files
  if h : parts.size = 0 then
    throw <| IO.userError s!"No compiled module data for {moduleName}"
  else
    if parts[0].1.isModule && files.size < 3 then
      throw <| IO.userError s!"Private object-file data is missing for {moduleName}"
    return parts[parts.size - 1].1.constants

/-- Put each new constant after its available dependencies for granular progress reports. -/
partial def visitReplayOrder (constants : Std.HashMap Name ConstantInfo) (name : Name) :
    StateM (NameSet × Array Name) Unit := do
  let (seen, _) ← get
  if seen.contains name then return
  let some ci := constants[name]? | return
  modify fun (seen, order) => (seen.insert name, order)
  for dependency in ci.getUsedConstantsAsSet do
    visitReplayOrder constants dependency
  modify fun (seen, order) => (seen, order.push name)

/-- Produce a dependency-first order without trusting any stored declaration. -/
def replayOrder (constants : Std.HashMap Name ConstantInfo) : Array Name :=
  ((do
    for (name, _) in constants.toList do
      visitReplayOrder constants name) : StateM (NameSet × Array Name) Unit).run (∅, #[]) |>.2.2

/-- Profile the standard replay routine one dependency-ordered declaration at a time. -/
def replayWithProgress (constants : Std.HashMap Name ConstantInfo)
    (env : Environment) : IO Environment := do
  let order := replayOrder constants
  let mut remaining : NameSet := ∅
  for (name, _) in constants.toList do
    remaining := remaining.insert name
  let program : Environment.Replay.M Unit := do
    let mut completed := 0
    for name in order do
      if (← get).remaining.contains name then
        IO.println s!"CHECK {completed + 1}/{constants.size} {name}"
        (← IO.getStdout).flush
        let started ← IO.monoMsNow
        Environment.Replay.replayConstant name
        IO.println s!"TIME {(← IO.monoMsNow) - started} ms {name}"
        (← IO.getStdout).flush
        completed := constants.size - (← get).remaining.size
    Environment.Replay.checkPostponedConstructors
    Environment.Replay.checkPostponedRecursors
  let (_, state) ← StateRefT'.run
    (s := { env := env.toKernelEnv, remaining : Environment.Replay.State })
    (ReaderT.run program { newConstants := constants })
  return .ofKernelEnv state.env

/-- The only axiom dependencies permitted in the new proof development. -/
def permittedAxiom (name : Name) : Bool :=
  name == ``propext || name == ``Classical.choice || name == ``Quot.sound

/-- Audit actual new declaration expressions against trusted baseline axiom metadata.
Every new declaration is inspected, so following each edge into the baseline suffices:
edges to another new declaration are audited when that declaration is visited. New
object-file axiom metadata is never used. Baseline results are cached by constant name. -/
def checkAxiomDependencies (constants : Std.HashMap Name ConstantInfo)
    (baseline : Environment) : IO Nat := do
  let mut cache : Std.HashMap Name (Array Name) := {}
  for (name, ci) in constants.toList do
    for dependency in ci.getUsedConstantsAsSet do
      if constants.contains dependency then continue
      if (baseline.toKernelEnv.find? dependency).isNone then
        throw <| IO.userError s!"Unresolved dependency {dependency} in new declaration {name}"
      let axioms ← match cache[dependency]? with
        | some axioms => pure axioms
        | none =>
          (collectAxioms dependency : CoreM (Array Name)).toIO'
            { fileName := "kernel-replay-axiom-audit", fileMap := FileMap.ofString "" }
            { env := baseline }
      cache := cache.insert dependency axioms
      for axiomName in axioms do
        if !permittedAxiom axiomName then
          throw <| IO.userError
            s!"Forbidden axiom {axiomName} in new declaration {name} through baseline dependency {dependency}"
  return cache.size

/-- Only identical safe theorem exports may be shared by separately compiled modules.
Strict expression equality includes binder names and annotations; proof values are
compared as well as types, universe parameters, and mutual-block metadata. -/
def identicalSafeTheorems (first second : ConstantInfo) : Bool :=
  match first, second with
  | .thmInfo a, .thmInfo b =>
    a.name == b.name && a.levelParams == b.levelParams &&
      a.type.equal b.type && a.value.equal b.value && a.all == b.all
  | _, _ => false

/-- Recheck the entire new development against only its preexisting library imports. -/
public unsafe def main (args : List String) : IO UInt32 := do
  initSearchPath (← findSysroot)
  enableInitializersExecution
  let (baseFile, newFile, capstoneFile, profile) ← match args with
    | [baseFile, newFile, capstoneFile] => pure (baseFile, newFile, capstoneFile, false)
    | [baseFile, newFile, capstoneFile, "--profile"] => pure (baseFile, newFile, capstoneFile, true)
    | _ => throw <| IO.userError "Usage: QEC360Replay baseline-modules.txt new-modules.txt capstones.txt [--profile]"
  let baseNames := (← IO.FS.lines baseFile).filter (· != "")
  let newNames := (← IO.FS.lines newFile).filter (· != "")
  let capstones := (← IO.FS.lines capstoneFile).filter (· != "")
  if capstones.isEmpty then throw <| IO.userError "At least one required capstone theorem is necessary"
  if newNames.isEmpty then throw <| IO.userError "At least one new module is necessary"
  let started ← IO.monoMsNow
  let imports : Array Import := baseNames.map fun s => { module := s.toName }
  let env ← importModules imports {} 0 (loadExts := true)
  let mut constants : Std.HashMap Name ConstantInfo := {}
  let mut runtimeDeclarations := 0
  let mut duplicateTheorems := 0
  let mut origins : Std.HashMap Name String := {}
  for name in newNames do
    let entries ← readConstants name.toName
    for ci in entries do
      match ci with
      | .axiomInfo _ => throw <| IO.userError s!"New axiom or erased theorem body: {ci.name}"
      | _ => pure ()
      if ci.isUnsafe || ci.isPartial then
        runtimeDeclarations := runtimeDeclarations + 1
        continue
      if (env.toKernelEnv.find? ci.name).isSome then
        throw <| IO.userError s!"A new declaration was already imported: {ci.name}"
      if let some first := constants[ci.name]? then
        if identicalSafeTheorems first ci then
          duplicateTheorems := duplicateTheorems + 1
          continue
        throw <| IO.userError
          s!"Differing or non-theorem duplicate declaration {ci.name} in {origins[ci.name]!} and {name}"
      constants := constants.insert ci.name ci
      origins := origins.insert ci.name name
  for theoremName in capstones do
    match constants[theoremName.toName]? with
    | some (.thmInfo _) => pure ()
    | _ => throw <| IO.userError s!"Required safe theorem is absent: {theoremName}"
  IO.println s!"Coalesced {duplicateTheorems} structurally identical safe theorem exports"
  IO.println s!"Required {capstones.size} safe capstone theorems"
  IO.println s!"Excluded {runtimeDeclarations} unsafe or partial runtime declarations"
  IO.println s!"Rechecking {constants.size} declarations from {newNames.size} modules"
  (← IO.getStdout).flush
  let dependencyCount ← checkAxiomDependencies constants env
  IO.println s!"Axiom audit accepted {dependencyCount} distinct baseline dependencies"
  (← IO.getStdout).flush
  let replayStarted ← IO.monoMsNow
  let checked ← if profile then replayWithProgress constants env else Environment.replay constants env
  for (name, _) in constants.toList do
    if (checked.toKernelEnv.find? name).isNone then
      throw <| IO.userError s!"Kernel did not retain new declaration {name}"
  for theoremName in capstones do
    match checked.toKernelEnv.find? theoremName.toName with
    | some (.thmInfo _) => pure ()
    | _ => throw <| IO.userError s!"Kernel did not retain required theorem {theoremName}"
  IO.println s!"Kernel accepted all declarations: replay {(← IO.monoMsNow) - replayStarted} ms; total {(← IO.monoMsNow) - started} ms"
  return 0
