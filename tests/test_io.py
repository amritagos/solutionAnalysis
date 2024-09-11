import pytest
from pathlib import Path
import soluanalysis as solu
from soluanalysis.io import read_lammps_dump, write_lammps_dump
from soluanalysis.hdf5_io import read_ion_pairs_from_hdf5, save_ion_pairs_to_hdf5


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
    systems, timesteps = read_lammps_dump(file_path)

    assert systems[0].atoms[0].mol_id == 31
    assert systems[0].n_atoms() == small_system.n_atoms()
    assert timesteps[0] == timestep


def test_read_lammps_dump(octahedral_system):
    """
    Tests whether you can read in a trajectory into a list of System objects

    octahedral_system (tuple[List[System], List[int]]): a list of System objects and a list of timesteps, corresponding to the octahedral system
    """
    # Read in the trajectory
    # Since this is a single frame, there is just one System object in the list
    atoms, timesteps = octahedral_system

    # There is only one timestep
    assert timesteps[0] == 0

    # The number of atoms should be correct
    assert atoms[0].n_atoms() == 22

    # Check the box size and lower box limits
    assert atoms[0].box == [49.4752565268, 49.4752565268, 49.4752565268]
    assert atoms[0].boxLo == [-12.457628299, -12.457628299, -12.457628299]

    # Check that the positions were properly added
    assert atoms[0].atoms[-1].position == [28.8714, 21.5182, 20.8169]

    # Check that the molecule IDs are correct
    mol_ids_desired = [
        521,
        1627,
        3968,
        4099,
        333,
        333,
        333,
        687,
        687,
        687,
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
    ]
    for atom, mol_id_expected in zip(atoms[0].atoms, mol_ids_desired):
        assert atom.mol_id == mol_id_expected


def test_hdf5_files(octahedral_system):
    """Tests whether you can output the time series information of the ion pairs to an HDF5 file format. The keys are timesteps, and the values are are a dictionary corresponding to ion pais (sorted by length, which are the keys at that timestep"""
    # Sample output dictionary for a time series of ion pairs
    time_series_data = {
        100: {2: [[5, 1]], 3: [[2, 11, 1], [3, 5, 1], [4, 14, 1]]},
        200: {2: [[5, 1]]},
    }
    timesteps = [100, 200]
    systems, times = octahedral_system

    test_dir = Path(__file__).resolve().parent
    file_path = test_dir / "output_ion_pairs.h5"
    max_depth = 3
    write_identifier = solu.james.WriteIdentifier.AtomID

    save_ion_pairs_to_hdf5(
        file_path, time_series_data, systems[0], max_depth, write_identifier
    )

    time_series_read, timesteps_read, system_read, max_depth_read, identifier_read = (
        read_ion_pairs_from_hdf5(file_path)
    )

    assert time_series_read == time_series_data
    assert timesteps_read == timesteps
    assert max_depth_read == max_depth
    assert identifier_read == write_identifier
    assert system_read.collect_ids() == systems[0].collect_ids()
    assert system_read.boxLo == systems[0].boxLo
