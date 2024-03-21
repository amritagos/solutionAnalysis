#pragma once
#define PYBIND11_DETAILED_ERROR_MESSAGES
// Standard libraries
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

// Basics
#include "add.hpp"
// Bindings
#include <pybind11/pybind11.h>
// Additional
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
// Namespaces
using namespace std::string_literals; // For ""s
using namespace pybind11::literals;   // For ""_a
namespace py = pybind11;              // Convention

PYBIND11_MODULE(solvlib, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &my_add, "A function that adds two numbers");
}