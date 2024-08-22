import pytest
import soluanalysis as solu


def create_network_list():
    """Generate a list of networks for the test case"""
    network_list = []
    n_atoms = 4
    weight_edge = 1

    # First time step
    network = solu.graphlib.UndirectedNetwork(n_atoms)
    network.push_back_neighbour_and_weight(0, 1, weight_edge)
    network.push_back_neighbour_and_weight(0, 2, weight_edge)
    network.push_back_neighbour_and_weight(0, 3, weight_edge)
    network.push_back_neighbour_and_weight(1, 2, weight_edge)
    network.push_back_neighbour_and_weight(1, 3, weight_edge)
    network.push_back_neighbour_and_weight(2, 3, weight_edge)
    network_list.append(network)

    # Second time step
    network = solu.graphlib.UndirectedNetwork(n_atoms)
    network.push_back_neighbour_and_weight(0, 1, weight_edge)
    network.push_back_neighbour_and_weight(0, 3, weight_edge)
    network.push_back_neighbour_and_weight(1, 2, weight_edge)
    network.push_back_neighbour_and_weight(1, 3, weight_edge)
    network.push_back_neighbour_and_weight(2, 3, weight_edge)
    network_list.append(network)
    # Third time step
    network = solu.graphlib.UndirectedNetwork(n_atoms)
    network.push_back_neighbour_and_weight(0, 1, weight_edge)
    network.push_back_neighbour_and_weight(1, 2, weight_edge)
    network.push_back_neighbour_and_weight(1, 3, weight_edge)
    network.push_back_neighbour_and_weight(2, 3, weight_edge)
    network_list.append(network)
    # Fourth time step
    network = solu.graphlib.UndirectedNetwork(n_atoms)
    network.push_back_neighbour_and_weight(1, 2, weight_edge)
    network.push_back_neighbour_and_weight(1, 3, weight_edge)
    network.push_back_neighbour_and_weight(2, 3, weight_edge)
    network_list.append(network)
    # 5th time step
    network = solu.graphlib.UndirectedNetwork(n_atoms)
    network.push_back_neighbour_and_weight(1, 2, weight_edge)
    network.push_back_neighbour_and_weight(1, 3, weight_edge)
    network.push_back_neighbour_and_weight(2, 3, weight_edge)
    network_list.append(network)
    # 6th time step
    network = solu.graphlib.UndirectedNetwork(n_atoms)
    network.push_back_neighbour_and_weight(1, 2, weight_edge)
    network.push_back_neighbour_and_weight(1, 3, weight_edge)
    network.push_back_neighbour_and_weight(2, 3, weight_edge)
    network_list.append(network)

    return network_list


def test_hbond_correlation():
    """Test to obtain the correlation function given a list of UndirectedNetwork objects"""
    times = [0, 10, 20]  # Times at timesteps required for the correlation function
    networks = create_network_list()
    # Check that the network list was created properly
    assert solu.graphlib.get_neighbours(networks[0], 0) == [1, 2, 3]
    assert solu.graphlib.get_neighbours(networks[-1], 0) == []

    # Get the time correlation function
    tau_values_expected = [0, 10, 20]
    tcf_avg_expected = [1.0, 0.7944444444444444, 0.6722222222222222]
    tau_values, tcf_avg, tcf_error = solu.james.time_correlation_function(
        networks, times, 0, 1, 1, None
    )
    assert tau_values == tau_values_expected
    assert tcf_avg == tcf_avg_expected
