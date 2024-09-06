import numpy as np
from typing import Dict, List, Optional, TextIO, Tuple, Union
from pathlib import Path
from collections import deque
from os.path import splitext
from soluanalysis.james import Atom, System
import h5py
import soluanalysis as solu
import numbers


def save_system_to_hdf5(system: solu.james.System, hdf5_group: h5py.Group):
    """Save a System object to an HDF5 group.

    Args:
        system (solu.james.System): System object to serialize
        hdf5_group (h5py.Group): The HDF5 group to be saved to
    """

    atom_ids = np.array([atom.id for atom in system.atoms], dtype=np.int32)
    atom_types = np.array([atom.type for atom in system.atoms], dtype=np.int32)
    mol_ids = np.array(
        [atom.id if atom.mol_id is None else atom.mol_id for atom in system.atoms],
        dtype=np.int32,
    )
    positions = np.array([atom.position for atom in system.atoms], dtype=np.float64)

    hdf5_group.create_dataset("atom_ids", data=atom_ids)
    hdf5_group.create_dataset("atom_types", data=atom_types)
    hdf5_group.create_dataset("mol_ids", data=mol_ids)
    hdf5_group.create_dataset("positions", data=positions)

    if system.box is not None:
        hdf5_group.create_dataset("box", data=np.array(system.box, dtype=np.float64))
    if system.boxLo is not None:
        hdf5_group.create_dataset(
            "boxLo", data=np.array(system.boxLo, dtype=np.float64)
        )


def read_ion_pairs_from_hdf5(
    file_path: Path,
) -> Tuple[
    Dict[int, Dict[int, List[List[int]]]],
    List[int],
    solu.james.System,
    int,
    solu.james.WriteIdentifier,
]:
    """Reads the HDF5 file and reconstructs a dictionary with the time series information about the ion pairs

    Args:
        file_path (Path): The HDF5 file to read from

    Returns:
        Tuple[Dict[int, Dict[int, List[List[int]]]], List[int], int, solu.james.WriteIdentifier]: A tuple containing
        1) the dictionary with the ion pairs,
        2) timesteps,
        3) System object
        4) max_depth,
        5) writeIdentifier
    """
    time_series_dict = {}

    enum_mapping = {
        "WriteIdentifier.AtomID": solu.james.WriteIdentifier.AtomID,
        "WriteIdentifier.Index": solu.james.WriteIdentifier.Index,
    }

    with h5py.File(file_path, "r") as file:
        # Read the metadata
        max_depth = file.attrs["max_depth"]
        identifier_str = file.attrs["writeIdentifier"]

        # Read the timesteps
        timesteps = file["timesteps"][:].tolist()

        # Read the representative System object
        system_group = file["system"]
        system = read_system_from_hdf5(system_group)

        # Iterate over the timesteps
        for timestep in timesteps:
            timestep_group = file[str(timestep)]
            groups = {}

            # Iterate over the lengths within each timestep
            for length in timestep_group.keys():
                length_group = timestep_group[length]

                # Read the numpy array and convert it back to a list of lists
                data_array = length_group["ion_pairs"][:]
                lists = data_array.tolist()

                groups[int(length)] = lists

            time_series_dict[int(timestep)] = groups

        # Return the time series, timesteps, max_depth, the enum, and the number of atoms
        return (
            time_series_dict,
            timesteps,
            system,
            max_depth,
            enum_mapping.get(identifier_str),
        )


def read_system_from_hdf5(hdf5_group: h5py.Group) -> solu.james.System:
    """Read a System object from an HDF5 group.

    Args:
        hdf5_group (h5py.Group): HDF5 group, from which the System object will be reconstructed

    Returns:
        solu.james.System: Reconstructed System object
    """
    atom_ids = hdf5_group["atom_ids"][:]
    atom_types = hdf5_group["atom_types"][:]
    mol_ids = hdf5_group["mol_ids"][:]
    positions = hdf5_group["positions"][:]

    # Reconstruct atoms list
    atoms = [
        solu.james.Atom(atom_id, atom_type, mol_id, position)
        for atom_id, atom_type, mol_id, position in zip(
            atom_ids, atom_types, mol_ids, positions
        )
    ]

    # Read optional attributes
    box = hdf5_group["box"][:] if "box" in hdf5_group else None
    boxLo = hdf5_group["boxLo"][:] if "boxLo" in hdf5_group else None

    # Reconstruct the System object
    system = solu.james.System(atoms, box, boxLo)

    return system


def save_ion_pairs_to_hdf5(
    file_path: Path,
    time_series_dict: Dict[int, Dict[int, List[List[int]]]],
    system: solu.james.System,
    max_depth: int,
    write_identifier: solu.james.WriteIdentifier,
    **compression_kwargs: Union[str, int],
) -> None:
    """Save the ion pairs per time step, sorted according to length into an HDF5 file.

    Args:
        file_path (Path): File path of the HDF5 file to write to
        time_series_dict (Dict[int, Dict[int, List[List[int]]]]): Dictionary containing timesteps and ion pairs.
        The keys of the outer dictionary are timesteps, and the keys of the inner dictionary are ion pair lengths
        system (solu.james.System): Representative System object, containing indices, atom IDs, atom types, molecular IDs
        max_depth (int): Maximum length of the ion pair
        write_identifier (solu.james.WriteIdentifier): enum class which describes whether the elements correspond to
        atom IDs or indices in the System object.
        compression_kwargs(Union[str, int]): additional compression options for the create_dataset command in h5py.
        For instance, compression="gzip" and compression_opts=4
    """
    # Extract timesteps from the keys of time_series_dict and sort
    timesteps = sorted(time_series_dict.keys())

    with h5py.File(file_path, "w") as file:
        # Save metadata
        file.attrs["max_depth"] = max_depth
        file.attrs["writeIdentifier"] = str(write_identifier)  # convert enum to string

        # Save the timesteps as a separate dataset
        file.create_dataset(
            "timesteps", data=np.array(timesteps, dtype=np.int32), **compression_kwargs
        )

        # Save the System object
        system_group = file.create_group("system")
        save_system_to_hdf5(system, system_group)

        # Now save the ion pairs per timestep into separate groups (each length would be in a different group)
        for timestep in timesteps:
            groups = time_series_dict[timestep]

            # Create a group for each timestep (timesteps are unique)
            timestep_group = file.create_group(str(timestep))

            for length, data in groups.items():
                # Create a subgroup for each length (can go upto max_length)
                length_group = timestep_group.create_group(str(length))

                # Convert the list of lists to a numpy array
                ion_pair_data = np.array(data, dtype=np.int32)

                # Save the numpy array to the HDF5 file
                length_group.create_dataset(
                    "ion_pairs", data=ion_pair_data, **compression_kwargs
                )
