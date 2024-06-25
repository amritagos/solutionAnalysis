#pragma once
#include <optional>
#include <system_error>
#include <vector>

namespace SolvLib {

// Holds information about an individual atom
struct Atom {
  int id{};                                 // Identifier number
  std::optional<int> mol_id = std::nullopt; // Molecule identifier
  int type{}; // Type number of the atom (could be an atomic number, or type
              // number as in LAMMPS)
};

class System {
public:
  std::vector<Atom> atoms{};
  std::optional<std::vector<double>> box = std::nullopt;

  System(const std::vector<Atom> &atoms, std::optional<std::vector<double>> box)
      : atoms(atoms), box(box) {}

  System() = default;

  // Delete the (n+1)^th Atom in the System object
  void del(int n) { atoms.erase(atoms.begin() + n); }

  // Delete a range of Atom objects, in the range [first, last)
  // For instance, if first=0 and last=3, the first three Atom objects will be
  // deleted
  void del(int first, int last) {
    atoms.erase(atoms.begin() + first, atoms.begin() + last);
  }

  // Add an atom to the System object
  void push_back(const Atom &atom) { atoms.push_back(atom); }

  // Get the number of atoms
  int n_atoms() { return atoms.size(); }
};
} // namespace SolvLib