#!/usr/bin/env python3
"""Reproduce Bravyi360 Lean certificates from the checked-in candidate data."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

from lean_style_cleanup import clean_certificate_source

PACKAGE = Path(__file__).resolve().parent
REPO = PACKAGE.parents[3]
MANIFEST_NAMES = ('root_manifest.json', 'audit_manifest.json',
                  'research_manifest.json', 'compression_manifest.json')
GENERATOR = 'experiments/bb_lab/scripts/bravyi360_kernel/generate.py'


def rewrite_imports(text: str, mapping: dict) -> str:
    """Retarget only module imports; preserve definitions and proofs verbatim."""
    def rewrite(match):
        return 'import ' + ' '.join(mapping.get(name, {}).get('target', name)
                                   for name in match.group(1).split())
    return re.sub(r'^import (.+)$', rewrite, text, flags=re.MULTILINE)



def kernel_only_source(text: str) -> str:
    """Retain the logical values while omitting unneeded runtime code for concrete data."""
    text = re.sub(r'(?m)^(private )?def\b', lambda m: (m.group(1) or '') + 'noncomputable def', text)
    text = re.sub(r'(?m)^(private )?instance\b', lambda m: (m.group(1) or '') + 'noncomputable instance', text)
    if 'noncomputable section' not in text:
        matches = list(re.finditer(r'^import .+$', text, re.MULTILINE))
        position = matches[-1].end() + 1 if matches else 0
        text = text[:position] + '\nnoncomputable section\n' + text[position:]
    return clean_certificate_source(text)


def banner(manifest: str, emitter: str) -> str:
    return ('/-\nGENERATED FILE — DO NOT HAND-EDIT\n'
            f'Generator: qec-lab:{GENERATOR}\n'
            'Emitter, relative to generator directory:\n'
            f'  {emitter}\n'
            f'Data and SHA-256 hashes: {manifest}\n'
            'Regenerate:\n'
            '  QECLEAN_ROOT=/path/to/QECLean \\\n'
            f'    python3 {GENERATOR} --force\n'
            '-/\n\n')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qec-root', type=Path,
                        default=Path(os.environ.get('QECLEAN_ROOT', REPO.parent / 'QECLean')))
    parser.add_argument('--force', action='store_true', help='Allow replacing generated Lean files.')
    parser.add_argument('--verify', action='store_true', help='Compare regenerated files without writing QECLean.')
    parser.add_argument('--stage', type=Path, help='Use a fresh persistent staging directory.')
    args = parser.parse_args()
    qec_root = args.qec_root.resolve()
    if not (qec_root / 'QEC').is_dir():
        parser.error(f'QECLean checkout not found: {qec_root}')
    manifests = [(name, json.loads((PACKAGE / name).read_text())) for name in MANIFEST_NAMES]
    mapping = json.loads((PACKAGE / 'module_map.json').read_text())['modules']
    outputs = {}
    for name, manifest in manifests:
        for step in manifest['steps']:
            for output in step['outputs']:
                module = Path(output).stem
                if module not in mapping:
                    continue
                item = mapping[module]
                if item.get('existing', False):
                    continue
                outputs[module] = (output, name, step['script'], item['target'])
    if not args.verify and not args.force:
        existing = [target for _, _, _, target in outputs.values()
                    if (qec_root / (target.replace('.', '/') + '.lean')).exists()]
        if existing:
            parser.error(f'{len(existing)} generated files already exist; use --force to replace them.')
    temporary = None
    if args.stage:
        stage = args.stage.resolve()
        if stage.exists():
            parser.error('Staging directory already exists; choose a fresh directory.')
        stage.mkdir(parents=True)
    else:
        temporary = tempfile.TemporaryDirectory(prefix='bravyi360-certificates-')
        stage = Path(temporary.name)
    try:
        seen = {}
        for _, manifest in manifests:
            for item in manifest['inputs']:
                data = gzip.decompress((PACKAGE / item['compressed']).read_bytes())
                digest = hashlib.sha256(data).hexdigest()
                if digest != item['sha256']:
                    raise RuntimeError(f'Certificate digest mismatch: {item["path"]}')
                if item['path'] in seen and seen[item['path']] != digest:
                    raise RuntimeError(f'Conflicting shared certificate: {item["path"]}')
                seen[item['path']] = digest
                target = stage / item['path']
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(data)
            for item in manifest.get('copies', []):
                target = stage / item['stage']
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(PACKAGE / item['source'], target)
        env = os.environ.copy()
        env['QEC360_STAGE'] = str(stage)
        for _, manifest in manifests:
            for step in manifest['steps']:
                arguments = [arg.format(stage=stage, package=PACKAGE) for arg in step.get('args', [])]
                subprocess.run([sys.executable, str(PACKAGE / step['script']), *arguments],
                               env=env, check=True, cwd=PACKAGE)
        mismatches = []
        for module, (source, manifest, emitter, target_module) in sorted(outputs.items()):
            text = (stage / source).read_text()
            text = banner(manifest, emitter) + kernel_only_source(rewrite_imports(text, mapping))
            target = qec_root / (target_module.replace('.', '/') + '.lean')
            if args.verify:
                if not target.exists() or target.read_text() != text:
                    mismatches.append(str(target.relative_to(qec_root)))
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(text)
        if mismatches:
            raise RuntimeError('Regeneration mismatch:\n' + '\n'.join(mismatches))
        action = 'Verified' if args.verify else 'Generated'
        print(f'{action} {len(outputs)} Lean modules from {len(seen)} hashed certificate inputs.')
    finally:
        if temporary is not None:
            temporary.cleanup()


if __name__ == '__main__':
    main()
