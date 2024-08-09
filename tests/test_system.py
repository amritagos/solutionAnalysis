import pytest
import soluanalysis as solu


@pytest.fixture
def simple_system():
    """
    Creates a System object with two water molecules (without reading in a file)
    """
    system = solu.james.System()
    system.box = [49.35, 49.35, 49.35]  # Box dimensions
    system.boxLo = [0, 0, 0]  # Lower limits of the box
    mol_ids = [1, 1, 1, 2, 2, 2]
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


def test_atom():
    """
    Test that you can create an Atom object
    """
    id = 1
    type = 8
    mol_id = 31
    atom = solu.james.Atom(id, type, mol_id, [33.701145, 31.147538, 48.163726])

    assert atom.id == id
    assert atom.type == type
    assert atom.mol_id == mol_id
    assert atom.position == [33.701145, 31.147538, 48.163726]


def test_system_creation():
    """
    Test that you can create a System object from vectors of Atom member parameters (such as positions).
    """
    ids = [1, 2, 3]
    types = [1, 2, 2]
    pos = [
        [33.701145, 31.147538, 48.163726],
        [33.608039, 32.001989, 48.014306],
        [32.880148, 30.854267, 48.195911],
    ]
    mol_ids = [1, 31, 1]
    box = [49.35, 49.35, 49.35]
    boxLo = [0, 0, 0]
    system = solu.james.System(ids, types, pos, mol_ids, box, boxLo)

    # Test that the IDs are correct
    assert system.collect_ids() == ids
    # Test the molecule IDs
    assert system.atoms[1].mol_id == 31


def test_system_object(simple_system):
    """
    Test that you can get the number of atoms and modify the System object etc.
    """
    n_atoms = simple_system.n_atoms()
    # There should be 6 atoms in the System at first
    assert n_atoms == 6

    # Delete the last atom
    simple_system.delete(5)
    assert simple_system.n_atoms() == 5

    # Delete the intermediate H atoms
    simple_system.delete(1, 3)
    assert simple_system.n_atoms() == 3

    # Check that the correct Atom objects have been deleted
    ids = simple_system.collect_ids()
    assert ids == [1, 4, 5]

    # # Check the molecule IDs
    mol_ids = [1, 2, 2]
    # Check the molecule IDs
    for i, (atom, mol_id) in enumerate(zip(simple_system.atoms, mol_ids)):
        print(f"Atom i={i}, id = {atom.id}, mol_id = {atom.mol_id}:\n")
        assert atom.mol_id == mol_id
