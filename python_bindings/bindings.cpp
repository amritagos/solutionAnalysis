#pragma once
#include <optional>
#include <vector>
#define PYBIND11_DETAILED_ERROR_MESSAGES
// Standard libraries
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

// Basics
#include "system.hpp"
// Bindings
#include <pybind11/pybind11.h>
// Additional
#include <pybind11/operators.h>
#include <pybind11/stl.h>
// Namespaces
using namespace std::string_literals; // For ""s
using namespace pybind11::literals;   // For ""_a
namespace py = pybind11;              // Convention

PYBIND11_MODULE(solvlib, m) {
  m.doc() = "Python bindings for solvlib"; // optional module docstring

  py::class_<SolvLib::Atom>(
      m, "Atom", "A struct for holding the ID, type, molecule ID and position")
      .def(py::init<>())
      .def(py::init<int, int, std::optional<int>, std::vector<double> &>())
      .def_readwrite("id", &SolvLib::Atom::id)
      .def_readwrite("mol_id", &SolvLib::Atom::mol_id)
      .def_readwrite("type", &SolvLib::Atom::type)
      .def_readwrite("position", &SolvLib::Atom::position);

  py::class_<SolvLib::System>(m, "System",
                              "A class for a collection of Atom objects")
      .def(py::init<>())
      .def(py::init<const std::vector<SolvLib::Atom> &,
                    std::optional<std::vector<double>> &,
                    std::optional<std::vector<double>> &>())
      .def_readwrite("atoms", &SolvLib::System::atoms)
      .def_readwrite("box", &SolvLib::System::box)
      .def_readwrite("boxLo", &SolvLib::System::boxLo)
      .def("n_atoms", &SolvLib::System::n_atoms)
      .def("delete",
           static_cast<void (SolvLib::System::*)(int)>(&SolvLib::System::del),
           "Delete the (n+1)^th Atom from the System")
      .def("delete",
           static_cast<void (SolvLib::System::*)(int, int)>(
               &SolvLib::System::del),
           "Delete a range of Atom objects, in the range [first, last)")
      .def("push_back", &SolvLib::System::push_back);
}