import pytest
import soluanalysis as solu


@pytest.fixture
def undirected_network():
    """
    Creates an UndirectedNetwork
    """
    n_nodes = 7
    weight_edge = 1
    network = solu.graphlib.UndirectedNetwork(n_nodes)
    network.push_back_neighbour_and_weight(0, 1, weight_edge)
    network.push_back_neighbour_and_weight(0, 2, weight_edge)
    network.push_back_neighbour_and_weight(1, 2, weight_edge)
    network.push_back_neighbour_and_weight(1, 3, weight_edge)
    network.push_back_neighbour_and_weight(1, 6, weight_edge)
    network.push_back_neighbour_and_weight(2, 6, weight_edge)
    network.push_back_neighbour_and_weight(3, 4, weight_edge)
    network.push_back_neighbour_and_weight(3, 6, weight_edge)
    network.push_back_neighbour_and_weight(4, 5, weight_edge)
    network.push_back_neighbour_and_weight(4, 6, weight_edge)

    return network


def test_undirected_network(undirected_network):
    """
    Test that you can get the number of nodes, edges etc.
    """
    n_nodes_required = 7
    # There should be 7 nodes in the undirected network
    assert undirected_network.n_nodes() == n_nodes_required

    # Test the number of edges
    n_edges_required = 10
    assert undirected_network.n_edges() == n_edges_required
    # There are two edges connected to 0
    assert undirected_network.n_edges(0) == 2
