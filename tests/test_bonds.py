import pytest
import soluanalysis as solu


def test_bond_formation(octahedral_system):
    """Tests that you can form bonds (distance-based and hydrogen bonds)

    Args:
        octahedral_system (tuple[List[System], List[int]]): a list of System objects and a list of timesteps, corresponding to the octahedral system
    """
    systems, timesteps = octahedral_system
    fe_type = 3
    o_type = 1
    h_type = 2
    cl_type = 4

    # There is only one frame (and only one system)
    for system in systems:
        n_atoms = system.n_atoms()
        # Distance-based bonds
        pair_feo = solu.james.Pair(fe_type, o_type)
        pairs = [pair_feo]
        cutoffs = [2.6]  # in Angstrom
        # Create the network
        network = solu.graphlib.UndirectedNetwork(n_atoms)
        solu.james.add_distance_based_bonds(network, system, pairs, cutoffs)
        # 6 bonds should have been created emanating from the Fe3+ center
        # We know the Fe ion has an index of 0
        assert network.n_edges(0) == 6
        # There should also only be a total of 6 bonds
        assert network.n_edges() == 6

        # Add intramolecular water molecule bonds
        # We will ignore hydrogens so there is no need for intramolecular hydrogen bonds
        donor_atom_types = [o_type]
        acceptor_atom_types = [cl_type, o_type]
        h_atom_types = [h_type]
        donor_acceptor_cutoff = 3.2
        max_angle_deg = 30  # in degrees
        # Continuous bonds
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
        # There will be three hydrogen bonds formed between the donor (O) and acceptors (Cl)
        # Therefore there should be a total of 9 bonds
        assert network.n_edges() == 9
