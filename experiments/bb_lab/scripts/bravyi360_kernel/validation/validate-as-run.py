#!/usr/bin/env python3
"""Freshly compile and kernel-replay the complete Bravyi360 proof development."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import time

PACKAGE = Path(__file__).resolve().parent
CAPSTONES = (
    'Quantum.Stabilizer.Homological.BB.Bravyi360.distance_lower_bound',
    'Quantum.Stabilizer.Homological.BB.Bravyi360.stabilizerCode_hasDistance24',
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def imports(path: Path) -> list[str]:
    return [name for row in re.findall(r'^import (.+)$', path.read_text(), re.M)
            for name in row.split() if name != 'all']


def verify_unchanged_inputs(project: Path, output: Path, package: Path,
                            report: dict, *, replay_artifacts: bool = False) -> None:
    """Bind the validation report to unchanged sources and replay inputs."""
    def require_hash(path: Path, expected: str) -> None:
        if not path.is_file() or digest(path) != expected:
            raise RuntimeError(f'Validation input changed or disappeared: {path}')

    for key in ('source_sha256', 'baseline_source_sha256'):
        for module, expected in report[key].items():
            require_hash(project / (module.replace('.', '/') + '.lean'), expected)
    for path, key in (
        (project / 'lean-toolchain', 'lean_toolchain_sha256'),
        (project / 'lake-manifest.json', 'lake_manifest_sha256'),
        (project / 'lakefile.toml', 'lakefile_sha256'),
        (package / 'module_map.json', 'module_map_sha256'),
    ):
        require_hash(path, report[key])
    if not replay_artifacts:
        return

    require_hash(package / 'KernelReplay.lean', report['checker_source_sha256'])
    artifacts = {str(path.relative_to(output))
                 for path in (output / 'QEC').rglob('*.olean*')}
    if artifacts != set(report['compiled_sha256']):
        raise RuntimeError('The copied object-file inventory changed during validation.')
    for relative, expected in report['compiled_sha256'].items():
        require_hash(output / relative, expected)
    for relative, expected in report['replay_support_sha256'].items():
        require_hash(output / relative, expected)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qec-root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path,
                        help='A fresh directory for compiled modules, logs, and measurements.')
    parser.add_argument('--base-revision', default='HEAD',
                        help='QECLean revision before this development was added.')
    parser.add_argument('--max-replay-seconds', type=float, default=900)
    parser.add_argument('--profile', action='store_true',
                        help='Report individual declarations using the same kernel replay routine.')
    args = parser.parse_args()
    project, output = args.qec_root.resolve(), args.output.resolve()
    if output.exists():
        parser.error('Refusing to reuse an output directory; choose a fresh path.')
    if args.max_replay_seconds <= 0:
        parser.error('The replay time limit must be positive.')
    output.mkdir(parents=True)
    env = os.environ.copy()
    env['LEAN_NUM_THREADS'] = '1'
    env['LAKE_ARTIFACT_CACHE'] = 'false'

    def capture(command: list[str]) -> str:
        return subprocess.check_output(command, cwd=project, env=env, text=True).strip()

    base_revision = capture(['git', 'rev-parse', args.base_revision])
    mapping = json.loads((PACKAGE / 'module_map.json').read_text())['modules']
    new = {item['target'] for item in mapping.values() if not item.get('existing')}
    if not new:
        raise RuntimeError('The proof module map is empty.')
    old_paths = set(capture(['git', 'ls-tree', '-r', '--name-only', base_revision]).splitlines())
    for module in new:
        if module.replace('.', '/') + '.lean' in old_paths:
            raise RuntimeError(f'{module} already exists at the baseline revision; '
                               'select the revision before this development.')
    baseline, visiting, visited, order = set(), set(), set(), []

    def source(module: str) -> Path:
        return project / (module.replace('.', '/') + '.lean')

    def visit(module: str) -> None:
        if module not in new:
            baseline.add(module)
            return
        if module in visited:
            return
        if module in visiting:
            raise RuntimeError(f'Import cycle at {module}')
        visiting.add(module)
        for dependency in imports(source(module)):
            visit(dependency)
        visiting.remove(module)
        visited.add(module)
        order.append(module)

    for module in sorted(new):
        visit(module)

    # Check the whole QEC baseline source closure against the selected old revision.
    # Thus a changed umbrella cannot quietly introduce the new proof as an import.
    baseline_sources = {}

    def verify_baseline(module: str) -> None:
        if module in new:
            raise RuntimeError(f'A baseline import depends on a new module: {module}')
        if not module.startswith('QEC.') or module in baseline_sources:
            return
        path = source(module)
        relative = path.relative_to(project).as_posix()
        old = subprocess.check_output(['git', 'show', f'{base_revision}:{relative}'], cwd=project)
        if old != path.read_bytes():
            raise RuntimeError(f'Baseline source differs from the old revision: {relative}')
        baseline_sources[module] = digest(path)
        for dependency in imports(path):
            verify_baseline(dependency)

    for module in sorted(baseline):
        verify_baseline(module)
    source_hashes = {module: digest(source(module)) for module in order}
    report = {'status': 'running', 'base_revision': base_revision,
              'platform': platform.platform(), 'machine': platform.machine(),
              'toolchain': (project / 'lean-toolchain').read_text().strip(),
              'lean_toolchain_sha256': digest(project / 'lean-toolchain'),
              'lake_manifest_sha256': digest(project / 'lake-manifest.json'),
              'lakefile_sha256': digest(project / 'lakefile.toml'),
              'module_map_sha256': digest(PACKAGE / 'module_map.json'),
              'lake_artifact_cache': False, 'lean_num_threads': 1,
              'replay_profile': args.profile,
              'source_sha256': source_hashes, 'baseline_source_sha256': baseline_sources,
              'capstones': CAPSTONES, 'steps': [], 'started_at': time.time()}
    if platform.system() == 'Darwin':
        try:
            report['cpu'] = capture(['sysctl', '-n', 'machdep.cpu.brand_string'])
            report['memory_bytes'] = int(capture(['sysctl', '-n', 'hw.memsize']))
        except (subprocess.CalledProcessError, ValueError):
            report['hardware_details'] = 'unavailable'

    def save() -> None:
        (output / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')

    def verify_inputs(phase: str, *, replay_artifacts: bool = False) -> None:
        started = time.monotonic()
        verify_unchanged_inputs(project, output, PACKAGE, report,
                               replay_artifacts=replay_artifacts)
        report.setdefault('integrity_checks', []).append({
            'phase': phase, 'seconds': round(time.monotonic() - started, 3)})
        save()

    def run(label: str, command: list[str], timeout: float) -> float:
        log = output / (label.replace('.', '_') + '.log')
        started = time.monotonic()
        with log.open('w') as handle:
            process = subprocess.Popen(command, cwd=project, env=env, stdout=handle,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            timed_out = False
            try:
                process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        seconds = time.monotonic() - started
        report['steps'].append({'label': label, 'seconds': round(seconds, 3),
                                'exit_code': process.returncode, 'timeout': timed_out})
        save()
        print(f'{label}: {seconds:.3f}s, exit={process.returncode}', flush=True)
        if process.returncode or timed_out:
            raise RuntimeError(f'{label} failed; see {log}')
        return seconds

    try:
        save()
        targets = sorted(module for module in baseline if module.startswith('QEC.'))
        if targets:
            run('baseline-build', ['lake', 'build', *targets], 1800)
        # Preserve and remove only artifacts belonging to the new development.
        # With traces absent and the cache disabled, Lake must check their sources.
        # This is also the required whole-library build, with its actual options.
        library = project / '.lake/build/lib/lean'
        for module in order:
            target = library / (module.replace('.', '/') + '.olean')
            artifacts = list(target.parent.glob(target.name + '*'))
            trace = target.with_suffix('.trace')
            if trace.exists():
                artifacts.append(trace)
            for artifact in artifacts:
                backup = output / 'previous-artifacts' / artifact.relative_to(library)
                backup.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(artifact), backup)
        started = time.monotonic()
        run('fresh-library-build', ['lake', '--rehash', '--no-cache', 'build'], 3600)
        report['source_check_seconds'] = round(time.monotonic() - started, 3)
        verify_inputs('after-source-build')
        # Lean resolves the root module directory once, so an isolated QEC tree
        # must also contain its preexisting QEC import closure.
        for module in sorted(baseline_sources):
            artifact = library / (module.replace('.', '/') + '.olean')
            if not artifact.exists():
                raise RuntimeError(f'A baseline artifact is missing: {module}')
            for part in artifact.parent.glob(artifact.name + '*'):
                target = output / part.relative_to(library)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(part, target)
        for module in order:
            artifact = library / (module.replace('.', '/') + '.olean')
            if not artifact.exists():
                raise RuntimeError(f'The full build omitted a new module: {module}')
            for part in artifact.parent.glob(artifact.name + '*'):
                target = output / part.relative_to(library)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(part, target)
        report['compiled_sha256'] = {
            str(path.relative_to(output)): digest(path)
            for path in sorted((output / 'QEC').rglob('*.olean*'))}
        base_path = capture(['lake', 'env', 'printenv', 'LEAN_PATH'])
        lean = capture(['lake', 'env', 'which', 'lean'])
        env['LEAN_PATH'] = str(output) + os.pathsep + base_path

        base_file, new_file, capstone_file = [output / name for name in
                                            ('baseline.txt', 'new-modules.txt', 'capstones.txt')]
        base_file.write_text('\n'.join(sorted(baseline)) + '\n')
        new_file.write_text('\n'.join(order) + '\n')
        capstone_file.write_text('\n'.join(CAPSTONES) + '\n')
        checker = PACKAGE / 'KernelReplay.lean'
        report['checker_source_sha256'] = digest(checker)
        run('checker-compile', [lean, '--root=' + str(PACKAGE),
                               '-o', str(output / 'KernelReplay.olean'),
                               '-c', str(output / 'KernelReplay.c'), str(checker)], 60)
        executable = output / 'kernel-replay'
        run('checker-link', ['lake', 'env', 'leanc', '-rdynamic', '-o', str(executable),
                            str(output / 'KernelReplay.c')], 60)
        # These files determine the executed checker and its exact module inputs.
        # Hashing happens outside the measured kernel-replay interval.
        support = [base_file, new_file, capstone_file, output / 'KernelReplay.c', executable]
        support.extend(output.glob('KernelReplay.olean*'))
        report['replay_support_sha256'] = {
            str(path.relative_to(output)): digest(path) for path in sorted(support)}
        command = [str(executable), str(base_file), str(new_file), str(capstone_file)]
        if args.profile:
            command.append('--profile')
        verify_inputs('before-replay', replay_artifacts=True)
        seconds = run('kernel-replay', command, args.max_replay_seconds)
        report['replay_wall_seconds'] = round(seconds, 3)
        report['replay_limit_seconds'] = args.max_replay_seconds
        if seconds >= args.max_replay_seconds:
            raise RuntimeError('Kernel replay exceeded the requested wall-clock bound.')
        verify_inputs('after-replay', replay_artifacts=True)
        report['status'] = 'passed'
        save()
        print(f'All {len(new)} new modules passed; kernel replay {seconds:.3f}s.', flush=True)
    except BaseException:
        report['status'] = 'failed'
        save()
        raise


if __name__ == '__main__':
    main()
