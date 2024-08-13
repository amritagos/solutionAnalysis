#pragma once
#include "network_base.hpp"
#include "undirected_network.hpp"
#include <cstddef>
#include <optional>
#include <vector>
#define PYBIND11_DETAILED_ERROR_MESSAGES
// Standard libraries
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

// Basics
#include "bondfinder.hpp"
#include "directed_network.hpp"
#include "network_base.hpp"
#include "pairtypes.hpp"
#include "system.hpp"
#include "undirected_network.hpp"
// Bindings
#include <pybind11/pybind11.h>
// Additional
#include <pybind11/operators.h>
#include <pybind11/stl.h>
// Namespaces
using namespace std::string_literals; // For ""s
using namespace pybind11::literals;   // For ""_a
namespace py = pybind11;              // Convention

PYBIND11_MODULE(james, m) {
  m.doc() = "Python bindings for james"; // optional module docstring

  py::class_<James::Atoms::Atom>(
      m, "Atom", "A struct for holding the ID, type, molecule ID and position")
      .def(py::init<>())
      .def(py::init<int, int, std::optional<int>, std::vector<double> &>())
      .def_readwrite("id", &James::Atoms::Atom::id)
      .def_readwrite("mol_id", &James::Atoms::Atom::mol_id)
      .def_readwrite("type", &James::Atoms::Atom::type)
      .def_readwrite("position", &James::Atoms::Atom::position);

  py::class_<James::Atoms::System>(m, "System",
                                   "A class for a collection of Atom objects")
      .def(py::init<>())
      .def(py::init<const std::vector<James::Atoms::Atom> &,
                    std::optional<std::vector<double>> &,
                    std::optional<std::vector<double>> &>())
      .def(py::init<const std::vector<int> &, const std::vector<int> &,
                    const std::vector<std::vector<double>> &,
                    std::optional<std::vector<int>>,
                    std::optional<std::vector<double>>,
                    std::optional<std::vector<double>>>(),
           "Constructor for System that takes in the ids, types, positions, "
           "optionally the molecule IDs, box size and lower box limits")
      .def_readwrite("atoms", &James::Atoms::System::atoms)
      .def_readwrite("box", &James::Atoms::System::box)
      .def_readwrite("boxLo", &James::Atoms::System::boxLo)
      .def("n_atoms", &James::Atoms::System::n_atoms)
      .def("collect_ids", &James::Atoms::System::collect_ids)
      .def("delete",
           static_cast<void (James::Atoms::System::*)(int)>(
               &James::Atoms::System::del),
           "Delete the (n+1)^th Atom from the System")
      .def("delete",
           static_cast<void (James::Atoms::System::*)(int, int)>(
               &James::Atoms::System::del),
           "Delete a range of Atom objects, in the range [first, last)")
      .def("push_back", &James::Atoms::System::push_back)
      .def("distance", &James::Atoms::System::distance,
           "Gets the distance between two Atom objects in the System, using "
           "the minimum image convention if the box has been defined.")
      .def("find_atoms_in_molecule",
           &James::Atoms::System::find_atoms_in_molecule,
           "Finds all indices in atoms such that the molecule ID is the same");

  // Bindings to commutative pair types and lambda binding to the function for
  // getting distance based bonds
  py::class_<James::Bond::Pair>(m, "Pair", "Commutative pair class")
      .def(py::init<int, int>(), py::arg("typeA"), py::arg("typeB"))
      .def_readwrite("typeA", &James::Bond::Pair::typeA)
      .def_readwrite("typeB", &James::Bond::Pair::typeB)
      .def(py::self == py::self);

  m.def("add_distance_based_bonds", [](Graph::NetworkBase<double> &network,
                                       const James::Atoms::System &system,
                                       std::vector<James::Bond::Pair> &pairs,
                                       std::vector<double> &cutoffs) {
    return James::Bond::add_distance_based_bonds<double>(network, system, pairs,
                                                         cutoffs);
  });
  // Binding for the templated function add_hbonds. This is templated on the
  // WeightType of the network
  m.def("add_hbonds",
        [](Graph::NetworkBase<double> &network,
           const James::Atoms::System &system,
           const std::vector<int> &donor_atom_types,
           const std::vector<int> &acceptor_atom_types,
           const std::vector<int> &h_atom_types, double cutoff_distance = 3.2,
           double max_angle_deg = 30, bool ignore_hydrogens = true) {
          return James::Bond::add_hbonds<double>(
              network, system, donor_atom_types, acceptor_atom_types,
              h_atom_types, cutoff_distance, max_angle_deg, ignore_hydrogens);
        });
}

PYBIND11_MODULE(graphlib, m) {
  m.doc() = "Python bindings for graph_lib"; // optional module docstring

  py::class_<Graph::NetworkBase<double>> give_me_a_name(
      m, "NetworkBase",
      "An abstract base class for undirected and directed networks");

  py::class_<Graph::UndirectedNetwork<double>, Graph::NetworkBase<double>>(
      m, "UndirectedNetwork",
      "A class that represents an undirected graph using adjacency lists")
      .def(py::init<>())
      .def(py::init<size_t>())
      .def(py::init<std::vector<std::vector<size_t>> &&,
                    std::vector<std::vector<double>> &&>())
      .def("n_nodes", &Graph::UndirectedNetwork<double>::n_agents,
           "A function that gives the number of nodes, or agents or atoms in "
           "the network")
      .def("n_edges", &Graph::UndirectedNetwork<double>::n_edges,
           "Gives the number of edges connected to the node. If node is None, "
           "gives the total number of edges",
           py::arg("agent_idx") = std::nullopt)
      .def("clear", &Graph::UndirectedNetwork<double>::clear,
           "Clears the network")
      .def("get_neighbours", &Graph::UndirectedNetwork<double>::get_neighbours,
           "Gives a view into the neighbour indices connected to node index")
      .def("get_weights", &Graph::UndirectedNetwork<double>::get_weights,
           "Gives a view into the edge weights connected to node index")
      .def("set_edge_weight",
           &Graph::UndirectedNetwork<double>::set_edge_weight,
           "Sets the weight for a node index, for an existing neighbour index")
      .def("get_edge_weight",
           &Graph::UndirectedNetwork<double>::get_edge_weight,
           "Gets the weight for a node index, for an existing neighbour index")
      .def("push_back_neighbour_and_weight",
           &Graph::UndirectedNetwork<double>::push_back_neighbour_and_weight,
           "Adds an edge between node_i and node_j with weight w. This could "
           "cause double counting, if not carefully called.")
      .def("remove_double_counting",
           &Graph::UndirectedNetwork<double>::remove_double_counting,
           "Sorts the neighbours by index and removes doubly counted edges by "
           "summing the weights");

  py::enum_<Graph::DirectedNetwork<double>::EdgeDirection>(
      m, "EdgeDirection",
      "Enum class for handling the direction of the edges in the directed "
      "network")
      .value("Incoming",
             Graph::DirectedNetwork<double>::EdgeDirection::Incoming)
      .value("Outgoing",
             Graph::DirectedNetwork<double>::EdgeDirection::Outgoing);

  py::class_<Graph::DirectedNetwork<double>, Graph::NetworkBase<double>>(
      m, "DirectedNetwork",
      "A class that represents a directed graph using adjacency lists")
      .def(py::init<>())
      .def(py::init<size_t>())
      .def(py::init<std::vector<std::vector<size_t>> &&,
                    std::vector<std::vector<double>> &&,
                    Graph::DirectedNetwork<double>::EdgeDirection>())
      .def("n_nodes", &Graph::DirectedNetwork<double>::n_agents,
           "A function that gives the number of nodes, or agents or atoms in "
           "the network")
      .def("n_edges", &Graph::DirectedNetwork<double>::n_edges,
           "Gives the number of edges connected to the node. If node is None, "
           "gives the total number of edges",
           py::arg("agent_idx") = std::nullopt)
      .def("clear", &Graph::DirectedNetwork<double>::clear,
           "Clears the network")
      .def("get_neighbours", &Graph::DirectedNetwork<double>::get_neighbours,
           "Gives a view into the neighbour indices connected to node index")
      .def("get_weights", &Graph::DirectedNetwork<double>::get_weights,
           "Gives a view into the edge weights connected to node index")
      .def("set_edge_weight", &Graph::DirectedNetwork<double>::set_edge_weight,
           "Sets the weight for a node index, for an existing neighbour index")
      .def("get_edge_weight", &Graph::DirectedNetwork<double>::get_edge_weight,
           "Gets the weight for a node index, for an existing neighbour index")
      .def("push_back_neighbour_and_weight",
           &Graph::DirectedNetwork<double>::push_back_neighbour_and_weight,
           "Adds an edge between node_i and node_j with weight w. This could "
           "cause double counting, if not carefully called.")
      .def("remove_double_counting",
           &Graph::DirectedNetwork<double>::remove_double_counting,
           "Sorts the neighbours by index and removes doubly counted edges by "
           "summing the weights");
}