# solutionAnalysis

Some code to find ion pairs in solution. WIP.

## Installation from Source

```bash
micromamba create -f environment.yml # For the first time only
micromamba activate soluenv
rm -rf subprojects 
git restore subprojects
meson setup build --wipe
pip install -e . --no-build-isolation
```

## Usage

`solvlib` is the C++ library. The library `graph_lib` contains classes for the Network class.

```bash
import soluanalysis as solu
system = solu.solvlib.System()
```

### Testing

We test with `pytest`. Run it like so (in verbose mode): 

```bash
pytest -v
```