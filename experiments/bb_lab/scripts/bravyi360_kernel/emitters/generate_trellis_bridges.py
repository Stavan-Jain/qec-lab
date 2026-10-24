"""Insert explicit finite certificate data into the proved trellis bridge templates."""
import argparse,json
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data',type=Path,required=True)
parser.add_argument('--templates',type=Path,default=Path(__file__).with_name('bridge_templates'))
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
data=json.loads(args.data.read_text())
for module,entries in data.items():
    source=args.templates.joinpath(module+'.lean.in').read_text()
    for key,value in entries.items():
        placeholder='@@'+key+'@@'
        assert source.count(placeholder)==1,(module,key)
        source=source.replace(placeholder,value)
    assert '@@' not in source,module
    args.out.joinpath(module+'.lean').write_text(source)

# Certificate data is evaluated by the kernel; no runtime duplicate is needed.
from kernel_only_data import mark_kernel_data
mark_kernel_data(args.out / name for name in ['TrellisBasisBridge.lean', 'ConcreteKernelSpan.lean', 'ReachableBasisChecks.lean', 'BaseCycleCensus.lean', 'ConcretePermutations.lean'])
