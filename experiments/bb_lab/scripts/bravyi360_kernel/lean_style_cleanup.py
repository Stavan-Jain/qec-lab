"""Keep generated certificate formatting and heartbeat budgets reproducible.

This conservative formatter preserves mathematical tokens and numeral values.
It rejects unsupported declaration forms instead of guessing their scopes.
"""
from __future__ import annotations
from bisect import bisect_right
import re


REASON = (
    "-- Elaborating this packed certificate can exceed the default heartbeat budget.\n"
)
LOCAL_BUDGET = "set_option maxHeartbeats 0 in\n" + REASON
SPECIAL = re.compile(r'/\-|--|"|«|\'(?:\\.|[^\'\\\n])\'')
DECL = (
    r"(?:(?:private|protected|noncomputable|unsafe|partial)\s+)*"
    r"(?:theorem|lemma|def|abbrev|opaque|instance|structure|inductive|class|example)\b"
)
COMMAND = re.compile(
    rf"(?m)^(?:{DECL}|(?:noncomputable\s+)?section\b|namespace\b|end\b|"
    r"set_option\s+maxHeartbeats\b|mutual\b)"
)
BULLET = re.compile(r"(?m)^([ \t]*)·[ \t]*\n\1  (?!(?:--|/-))(?=\S)")


def protected_spans(text: str) -> list[tuple[int, int, str]]:
    """Locate comments and quoted text without parsing or evaluating Lean."""
    result = []
    position = 0
    while match := SPECIAL.search(text, position):
        start = match.start()
        token = match.group()
        if token == "/-":
            depth, position = 1, match.end()
            while depth:
                left = text.find("/-", position)
                right = text.find("-/", position)
                if right < 0:
                    raise ValueError("Unclosed Lean block comment")
                if 0 <= left < right:
                    depth += 1
                    position = left + 2
                else:
                    depth -= 1
                    position = right + 2
            kind = "doc" if text.startswith("/--", start) else "comment"
        elif token == "--":
            position = text.find("\n", match.end())
            if position < 0:
                position = len(text)
            kind = "comment"
        elif token == '"':
            position = match.end()
            while position < len(text):
                if text[position] == "\\":
                    position += 2
                elif text[position] == '"':
                    position += 1
                    break
                else:
                    position += 1
            else:
                raise ValueError("Unclosed Lean string")
            kind = "quoted"
        elif token == "«":
            position = text.find("»", match.end())
            if position < 0:
                raise ValueError("Unclosed Lean quoted identifier")
            position += 1
            kind = "quoted"
        else:
            position = match.end()
            kind = "quoted"
        result.append((start, position, kind))
    return result


def _in_span(position: int, spans, starts) -> bool:
    index = bisect_right(starts, position) - 1
    return index >= 0 and position < spans[index][1]


def scope_heartbeats(text: str) -> str:
    """Localize every previously active global zero budget, including data defs.

    Namespace/section scopes and docstrings retain their original meaning.
    Already-clean files are unchanged, making regeneration idempotent.
    """
    if not re.search(r"(?m)^set_option maxHeartbeats 0[ \t]*$", text):
        return text
    spans = protected_spans(text)
    starts = [start for start, _, _ in spans]
    doc_by_end = {end: start for start, end, kind in spans if kind == "doc"}
    # These forms are absent from the target modules. Reject them rather than
    # silently miss a declaration or alter an attribute/command wrapper.
    unsupported = re.compile(
        rf"(?m)^(?:[ \t]+(?:{DECL}|namespace\b|(?:noncomputable\s+)?section\b|end\b)"
        r"|[ \t]*@\["
        r"|(?:local|scoped)\b"
        r"|(?:axiom|constant|initialize|builtin_initialize|macro|macro_rules|syntax|"
        r"elab|elab_rules|notation|infix|infixl|infixr|prefix|postfix|attribute)\b"
        r"|set_option[ \t]+[^\n]+[ \t]+in[ \t]+\S[^\n]*)"
    )
    for match in unsupported.finditer(text):
        if not _in_span(match.start(), spans, starts):
            line = text.count("\n", 0, match.start()) + 1
            raise ValueError(f"Declaration indentation/attributes/wrapper at line {line} "
                             "requires manual heartbeat scoping")
    scope = [None]
    edits = []
    for match in COMMAND.finditer(text):
        start = match.start()
        if _in_span(start, spans, starts):
            continue
        end = text.find("\n", start)
        if end < 0:
            end = len(text)
        line = text[start:end]
        command = match.group()
        if command.startswith("set_option"):
            if not re.fullmatch(r"set_option maxHeartbeats 0[ \t]*", line):
                raise ValueError("Mixed or nonzero heartbeat settings require manual review")
            scope[-1] = 0
            edits.append((start, min(end + 1, len(text)), ""))
        elif command == "namespace" or command.endswith("section"):
            scope.append(scope[-1])
        elif command == "end":
            if len(scope) == 1:
                raise ValueError("Unmatched namespace/section end")
            scope.pop()
        elif command == "mutual":
            raise ValueError("Mutual declarations require manual heartbeat scoping")
        elif scope[-1] == 0:
            prefix = start
            # Keep the option outside the declaration's attached docstring.
            before = prefix
            while before and text[before - 1].isspace():
                before -= 1
            if before in doc_by_end:
                prefix = doc_by_end[before]
                if prefix and text[prefix - 1] != "\n":
                    raise ValueError("Inline declaration docstring requires manual review")
            previous_end = prefix
            while previous_end and text[previous_end - 1].isspace():
                previous_end -= 1
            previous_start = text.rfind("\n", 0, previous_end) + 1
            previous = text[previous_start:previous_end]
            if previous.startswith("@["):
                raise ValueError("Standalone attributes require manual heartbeat scoping")
            edits.append((prefix, prefix, LOCAL_BUDGET))
    cursor, pieces = 0, []
    for start, end, replacement in sorted(edits):
        if start < cursor:
            raise ValueError("Overlapping declaration-prefix edits")
        pieces.extend((text[cursor:start], replacement))
        cursor = end
    pieces.append(text[cursor:])
    return "".join(pieces)


def join_standalone_bullets(text: str) -> str:
    """Join a bullet to its following tactic, preserving tactic indentation."""
    spans = protected_spans(text)
    starts = [start for start, _, _ in spans]
    return BULLET.sub(
        lambda match: match.group() if _in_span(match.start(), spans, starts)
        else match.group(1) + "· ", text
    )


def clean_certificate_source(text: str) -> str:
    """Apply only heartbeat scoping and standalone-bullet cleanup."""
    return join_standalone_bullets(scope_heartbeats(text))
