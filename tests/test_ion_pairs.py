import pytest
import soluanalysis as solu
from pathlib import Path


def test_ion_pairs(octahedral_system):
    """Tests that ion pairs in an octahedral system can be found.

    Args:
        octahedral_system (tuple[List[System], List[int]]): a list of System objects and a list of timesteps, corresponding to the octahedral system
    """
    fe_type = 3
    o_type = 1
    h_type = 2
    cl_type = 4
    output_dict = {}

    # Distance-based bonds
    pair_feo = solu.james.Pair(fe_type, o_type)
    pairs = [pair_feo]
    cutoffs = [2.6]  # in Angstrom
    donor_atom_types = [o_type]
    acceptor_atom_types = [cl_type, o_type]
    h_atom_types = [h_type]
    donor_acceptor_cutoff = 3.2
    max_angle_deg = 30  # in degrees
    # Information for the ion pairs
    destination_atom_types = [cl_type]
    intermediate_atom_types = [
        o_type,
        h_type,
    ]  # Technically just o_type would have sufficed since hydrogens are being ignored
    max_depth = 3
    identifier = solu.james.WriteIdentifier.AtomID

    # List of System objects and timestep from a LAMMPS trajectory
    systems, timesteps = octahedral_system

    for system, timestep in zip(systems, timesteps):
        ion_pair_list = []
        n_atoms = system.n_atoms()
        # Construct the network and fill it with bonds (distance-based and hydrogen bonds)
        network = solu.graphlib.UndirectedNetwork(n_atoms)
        # distance-based bonds
        solu.james.add_distance_based_bonds(network, system, pairs, cutoffs)
        # hydrogen bonds
        solu.james.add_hbonds(
            network,
            system,
            donor_atom_types,
            acceptor_atom_types,
            h_atom_types,
            donor_acceptor_cutoff,
            max_angle_deg,
            True,
        )
        # Find the source indices (Fe in this case)
        source_indices = [
            system.atoms.index(atom) for atom in system.atoms if atom.type in [fe_type]
        ]
        assert source_indices == [0]
        # Get the ion pairs for each source index
        # The ion pairs will be saved with the atom IDs as elements
        for source in source_indices:
            ion_pairs = solu.james.find_ion_pairs(
                source,
                network,
                system,
                destination_atom_types,
                intermediate_atom_types,
                max_depth,
                identifier,
            )
            ion_pair_list += ion_pairs
        output_dict[timestep] = ion_pair_list
    assert output_dict == {timesteps[0]: [[2, 11, 1], [3, 5, 1], [4, 14, 1]]}
