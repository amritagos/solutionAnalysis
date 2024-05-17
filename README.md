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

For now, I've just implemented a simple add function in C++
`solvlib` is the C++ library.

```bash
import soluanalysis as solu
solu.solvlib.add(1,2)
```

### Testing

We test with `pytest`. Run it like so (in verbose mode): 

```bash
pytest -v
```