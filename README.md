# solutionAnalysis

Some code to find ion pairs in solution. WIP.

## Installation from Source

```bash
micromamba create -f environment.yml # For the first time only
micromamba activate soluenv
rm -rf subprojects 
git restore subprojects
meson setup build --wipe
pip install -e . #--no-build-isolation
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

### Features

<p float="left">
    <img src="https://github.com/amritagos/james/blob/develop/resources/oct_with_hydrogens.png?raw=true" width="300" />
    <img src="https://github.com/amritagos/james/blob/develop/resources/oct_no_hydrogens.png?raw=true" width="300" />
</p>

Can be used to determine ion pairs. Has bindings to the `C++` library [James](https://github.com/amritagos/james). 