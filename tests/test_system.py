import pytest
import soluanalysis as solu


@pytest.fixture
def simple_system():
    """
    Creates a System object with two water molecules (without reading in a file)
    """
    system = solu.solvlib.System()
    system.box = [49.35, 49.35, 49.35]  # Box dimensions
    system.boxLo = [0, 0, 0]  # Lower limits of the box
    # id, type, mol_id, position
    system.push_back(solu.solvlib.Atom(1, 1, 1, [33.701145, 31.147538, 48.163726]))  # O
    system.push_back(
        solu.solvlib.Atom(2, 2, 1, [33.608039, 32.001989, 48.014306])
    )  # H1
    system.push_back(
        solu.solvlib.Atom(3, 2, 1, [32.880148, 30.854267, 48.195911])
    )  # H2
    system.push_back(solu.solvlib.Atom(4, 1, 2, [37.219210, 42.159763, 29.721422]))  # O
    system.push_back(
        solu.solvlib.Atom(5, 2, 2, [36.836746, 41.424653, 29.994215])
    )  # H1
    system.push_back(
        solu.solvlib.Atom(6, 2, 2, [36.614956, 42.771812, 29.867550])
    )  # H2

    return system


def test_system(simple_system):
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
