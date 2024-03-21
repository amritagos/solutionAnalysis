# solutionAnalysis

Some code to find ion pairs in solution. WIP.

## Installation from Source

```bash
micromamba create -f environment.yml # For the first time only
micromamba activate soluenv
pdm install
```

## Usage

For now, I've just implemented a simple add function in C++
`solvlib` is the C++ library.

```python
import soluanalysis as solu
solu.solvlib.add(1,2)
```