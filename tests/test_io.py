import pytest
from pathlib import Path
import soluanalysis as solu
from soluanalysis.io import read_lammps_dump


def test_read_lammps_dump():
    """
    Tests whether you can read in a trajectory into a list of System objects
    """

    test_dir = Path(__file__).resolve().parent
    infilename = test_dir / "../resources/oct.lammpstrj"
    # Read in the trajectory
    # Since this is a single frame, there is just one System object
    atoms, timesteps = read_lammps_dump(infilename)

    # There is only one timestep
    assert timesteps == 0

    # The number of atoms should be correct
    assert atoms.n_atoms() == 22

    # Check the box size and lower box limits
    assert atoms.box == [49.4752565268, 49.4752565268, 49.4752565268]
    assert atoms.boxLo == [-12.457628299, -12.457628299, -12.457628299]

    # Check that the positions were properly added
    assert atoms.atoms[-1].position == [25.1689, 20.8364, 22.0004]
