"""Stand-in for external docking programs in tests.

Records its command line, CUDA_VISIBLE_DEVICES and the files in its working directory to
`<cwd>/call.json`, then behaves like:

- smina/gnina (`--ligand ... --out ...`): writes `--num_modes` copies of the ligand with
  increasing minimizedAffinity;
- rbdock (`-i ... -o ...`): writes poses with unsorted SCORE values to `<o>.sd`;
- rbcavity / obabel: does nothing.

Set FAKE_ENGINE_FAIL=1 to exit with an error instead.
"""

import json
import os
import sys
from pathlib import Path

from rdkit import Chem

args = sys.argv[1:]
log = Path("call.json")
calls = json.loads(log.read_text()) if log.exists() else []
calls.append(
    {
        "args": args,
        "cuda": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "files": sorted(p.name for p in Path().iterdir()),
    }
)
log.write_text(json.dumps(calls))
# Keep a copy outside the engine's temporary directory for assertions.
if "FAKE_ENGINE_LOG" in os.environ:
    Path(os.environ["FAKE_ENGINE_LOG"]).write_text(json.dumps(calls))

if os.environ.get("FAKE_ENGINE_FAIL"):
    print("simulated failure", file=sys.stderr)
    sys.exit(3)


def value(flag):
    return args[args.index(flag) + 1]


if "--ligand" in args:
    mol = next(Chem.SDMolSupplier(value("--ligand"), removeHs=False))
    with Chem.SDWriter(value("--out")) as w:
        for i in range(int(value("--num_modes"))):
            mol.SetProp("minimizedAffinity", f"{-9.0 + i:.2f}")
            w.write(mol)
elif "-i" in args:
    mol = next(Chem.SDMolSupplier(value("-i"), removeHs=False))
    with Chem.SDWriter(value("-o") + ".sd") as w:
        for score in (-5.0, -12.0, -8.0):
            mol.SetProp("SCORE", str(score))
            w.write(mol)
