"""Emit kernel-checkable certificates into QEC360_STAGE."""
import os
from pathlib import Path
STAGE=Path(os.environ["QEC360_STAGE"]).resolve()
from pathlib import Path
p = Path(str(STAGE / 'OrbitBatchMetadata.lean'))
s = p.read_text()
if 'import OrbitBatchMetadataChecks' in s:
    raise SystemExit('already split')
k = s.index('/-- Orbit-stabilizer counts sum')
Path(str(STAGE / 'OrbitBatchMetadataChecks.lean')).write_text(s[:k] + '\nend OrbitBatchData\n')
p.write_text('import OrbitBatchMetadataChecks\n\nnamespace OrbitBatchData\n\n' + s[k:])
