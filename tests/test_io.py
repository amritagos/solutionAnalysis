import pytest
from pathlib import Path
import soluanalysis as solu
from soluanalysis.io import read_lammps_dump, write_lammps_dump


@pytest.fixture
def small_system():
    """
    Creates a System object with two water molecules (without reading in a file)
    """
    system = solu.james.System()
    system.box = [49.35, 49.35, 49.35]  # Box dimensions
    system.boxLo = [0, 0, 0]  # Lower limits of the box
    mol_ids = [31, 1, 1, 2, 2, 2]
    # id, type, mol_id, position
    system.push_back(
        solu.james.Atom(1, 1, mol_ids[0], [33.701145, 31.147538, 48.163726])
    )  # O
    system.push_back(
        solu.james.Atom(2, 2, mol_ids[1], [33.608039, 32.001989, 48.014306])
    )  # H1
    system.push_back(
        solu.james.Atom(3, 2, mol_ids[2], [32.880148, 30.854267, 48.195911])
    )  # H2
    system.push_back(
        solu.james.Atom(4, 1, mol_ids[3], [37.219210, 42.159763, 29.721422])
    )  # O
    system.push_back(
        solu.james.Atom(5, 2, mol_ids[4], [36.836746, 41.424653, 29.994215])
    )  # H1
    system.push_back(
        solu.james.Atom(6, 2, mol_ids[5], [36.614956, 42.771812, 29.867550])
    )  # H2

    return system


def test_write_read_single_dump(small_system):
    """
    Tests that you can write out and then read back in a System object containing two molecules.
    """
    test_dir = Path(__file__).resolve().parent
    file_path = test_dir / "test_single_frame.lammpstrj"
    timestep = 100
    # Write this out to a file
    write_lammps_dump(file_path, small_system, timestep)

    # Read the file in to a new Systems object
    system, timestep_read = read_lammps_dump(file_path)

    assert system.atoms[0].mol_id == 31
    assert system.n_atoms() == small_system.n_atoms()


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

    # Check that the molecule IDs are correct
    mol_ids_desired = [
        333,
        333,
        333,
        521,
        687,
        687,
        687,
        1627,
        1683,
        1683,
        1683,
        1924,
        1924,
        1924,
        2319,
        2319,
        2319,
        3273,
        3273,
        3273,
        3968,
        4099,
    ]
    for atom, mol_id_expected in zip(atoms.atoms, mol_ids_desired):
        assert atom.mol_id == mol_id_expected
