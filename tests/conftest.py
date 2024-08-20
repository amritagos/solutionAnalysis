from pathlib import Path
import pytest
import soluanalysis as solu


@pytest.fixture
def octahedral_system():
    """
    Reads a list of System objects (containing a single frame) with 1 Fe3+ ion, 3 Cl- ions, and 6 water molecules (22 atoms in total).
    """
    test_dir = Path(__file__).resolve().parent
    infilename = test_dir / "../resources/oct.lammpstrj"
    # Read in the trajectory
    systems, timesteps = solu.io.read_lammps_dump(infilename)
    return systems, timesteps
