import numpy as np
from typing import Dict, List, Optional, TextIO, Tuple, Union
from pathlib import Path
from collections import deque
from os.path import splitext
from soluanalysis.james import Atom, System
import h5py
import soluanalysis as solu


def lammps_data_to_system(
    data,
    colnames,
    box,
    boxLo,
    order=True,
):
    """Extract positions and other per-atom parameters and create System
    Similar to lammps_data_to_ase_atoms in ASE

    :param data: per atom data
    :param colnames: index for data
    :param box: cell dimensions
    :param boxLo: lower limits of the box
    :param order: sort atoms by id. Might be faster to turn off.
    :returns: System object
    :rtype: System

    """
    if len(data.shape) == 1:
        data = data[np.newaxis, :]

    # read IDs if given and order if needed
    if "id" in colnames:
        ids = data[:, colnames.index("id")].astype(int)
        if order:
            sort_order = np.argsort(ids)
            data = data[sort_order, :]

    # determine LAMMPS types
    if "type" in colnames:
        # fall back to `types` otherwise
        types = data[:, colnames.index("type")].astype(int)
    else:
        # Error if type has not been provided
        raise ValueError("Cannot determine atom types from LAMMPS dump file")

    # Get the molecule IDs
    if "mol" in colnames:
        # fall back to `types` otherwise
        mol_ids = data[:, colnames.index("mol")].astype(int)
    else:
        # # Use the type IDs as molecule IDs
        # mol_ids = data[:, colnames.index("id")].astype(int)
        # Set the molecule IDs to None
        mol_ids = None

    def get_quantity(labels):
        try:
            cols = [colnames.index(label) for label in labels]

            return data[:, cols].astype(float)
        except ValueError:
            return None

    # Positions
    positions = None
    scaled_positions = None
    if "x" in colnames:
        # doc: x, y, z = unscaled atom coordinates
        positions = get_quantity(["x", "y", "z"])
    elif "xs" in colnames:
        # doc: xs,ys,zs = scaled atom coordinates
        scaled_positions = get_quantity(["xs", "ys", "zs"])
    elif "xu" in colnames:
        # doc: xu,yu,zu = unwrapped atom coordinates
        positions = get_quantity(["xu", "yu", "zu"])
    elif "xsu" in colnames:
        # xsu,ysu,zsu = scaled unwrapped atom coordinates
        scaled_positions = get_quantity(["xsu", "ysu", "zsu"])
    else:
        raise ValueError("No atomic positions found in LAMMPS output")

    # Convert everything to lists
    ids = ids.tolist()
    types = types.tolist()
    positions = positions.tolist()
    mol_ids = mol_ids.tolist()
    box = box.tolist()
    boxLo = boxLo.tolist()
    if positions is not None:
        out_atoms = System(ids, types, positions, mol_ids, box, boxLo)
    elif scaled_positions is not None:
        out_atoms = System(ids, types, positions, mol_ids, box, boxLo)

    return out_atoms


def get_max_index(index):
    if np.isscalar(index):
        return index
    elif isinstance(index, slice):
        return index.stop if (index.stop is not None) else float("inf")


def write_single_frame(fd: TextIO, system: System, timestep: int) -> None:
    """Writes out a single frame to a file handle. Inspired by ASE `write_lammps_data` function.

    Args:
        fd (TextIO): File being written to
        system (System): The System object whose data is being written out
        timestep (int): Current timestep corresponding to the System object
    """
    # The output format is:
    # ITEM: TIMESTEP
    # 0
    # ITEM: NUMBER OF ATOMS
    # 22
    # ITEM: BOX BOUNDS pp pp pp
    # -12.457628299 37.0176282278
    # -12.457628299 37.0176282278
    # -12.457628299 37.0176282278
    # ITEM: ATOMS id type mol x y z
    n_atoms = system.n_atoms()

    # Write out the single frame
    fd.write("ITEM: TIMESTEP\n")
    fd.write(f"{timestep}\n")
    fd.write("ITEM: NUMBER OF ATOMS\n")
    fd.write(f"{n_atoms}\n")
    fd.write("ITEM: BOX BOUNDS pp pp pp\n")
    for i_boxLow, i_boxsize in zip(system.boxLo, system.box):
        box_hi = i_boxLow + i_boxsize
        fd.write(f"{i_boxLow:23.17g} {box_hi:23.17g}\n")
    fd.write("ITEM: ATOMS id type mol x y z\n")
    # Loop through all the atoms
    for atom in system.atoms:
        if atom.mol_id == None:
            mol_id = atom.id
        else:
            mol_id = atom.mol_id
        line = f"{atom.id:>6} {atom.type:>3} {mol_id:>6}"
        # Add the positions
        for pos in atom.position:
            line += f" {pos:23.17g}"
        line += "\n"
        fd.write(line)  # actually write the line


def write_lammps_dump(
    file_path: Path,
    systems: Union[List[System], System],
    timesteps: Union[List[int], int],
) -> None:
    """Writes out a LAMMPS trajectory file using a list of System objects (or a single System object). No support for triclinic boxes

    Args:
        file_path (Path): File path
        systems (Union[List[System], System]): System object(s) whose information will be written out
        timesteps (Union[List[int], int]) : Timesteps corresponding to each system
    """

    # The output format is:
    # ITEM: TIMESTEP
    # 0
    # ITEM: NUMBER OF ATOMS
    # 22
    # ITEM: BOX BOUNDS pp pp pp
    # -12.457628299 37.0176282278
    # -12.457628299 37.0176282278
    # -12.457628299 37.0176282278
    # ITEM: ATOMS id type mol x y z

    if isinstance(systems, list):
        # There are multiple steps to be written out
        # Raise exception if the systems and timesteps are not of the same length
        try:
            if len(systems) != len(timesteps):
                raise ValueError("Systems and timesteps do not have the same size.")
        except ValueError as err:
            print(err)
            if len(systems) > len(timesteps):
                systems = systems[: len(timesteps)]
            else:
                timesteps = timesteps[: len(systems)]
        with open(file_path, "w") as writer:
            # Write out the multiple steps
            for i, (timestep, system) in enumerate(zip(timesteps, systems)):
                write_single_frame(writer, system, timestep)
    else:  # for a single frame
        # if timesteps is a list somehow, then take the first value only
        if isinstance(timesteps, list):
            timestep = timesteps[0]
        else:
            timestep = timesteps
        with open(file_path, "w") as writer:
            write_single_frame(writer, systems, timestep)


def read_lammps_dump(file_path: Path, index=-1) -> tuple[List[System], List[int]]:
    """
    Reads a LAMMPS trajectory file and returns a list of System objects for the specified range of timesteps.
    Similar to the ASE read function. Does not have support for triclinic boxes.

    Args:
        file_path (pathlib.Path): Path to the LAMMPS trajectory file.
        index: integer or slice object (by default, gets the last timestep)
        index=':' or index=slice(None) : all

    Returns:
        List[System], List[int]: A list of System objects representing each frame in the trajectory within the specified range, and a list of timesteps read from the LAMMPS dump file.
    """
    n_atoms = 0
    images = []
    box_bounds = None
    lower_box_limits = None
    timesteps = []

    if index is None or index == ":":
        index = slice(None, None, None)

    try:
        with open(file_path, "r") as file:
            lines = deque(file.readlines())
            index_end = get_max_index(index)
    except IOError as e:
        print(f"Error reading file {file_path}: {e}")
        return images

    while len(lines) > n_atoms:
        line = lines.popleft()

        if "ITEM: TIMESTEP" in line:
            n_atoms = 0
            line = lines.popleft()
            timestep = int(line.split()[0])
            timesteps.append(timestep)

        if "ITEM: NUMBER OF ATOMS" in line:
            line = lines.popleft()
            n_atoms = int(line.split()[0])

        # No support for triclinic boxes
        if "ITEM: BOX BOUNDS" in line:
            try:
                celldatarows = [lines.popleft() for _ in range(3)]
                box_bounds = np.array(
                    [list(map(float, l.split())) for l in celldatarows]
                )
                lower_box_limits = box_bounds[:, 0]
            except (ValueError, IndexError):
                print(f"Error parsing BOX BOUNDS")
                file.close()
                return images

        if "ITEM: ATOMS" in line:
            colnames = line.split()[2:]
            datarows = [lines.popleft() for _ in range(n_atoms)]
            data = np.loadtxt(datarows, dtype=str)
            box = box_bounds[:, 1] - box_bounds[:, 0]
            out_atoms = lammps_data_to_system(
                data=data,
                colnames=colnames,
                box=box,
                boxLo=lower_box_limits,
                order=True,
            )
            images.append(out_atoms)

        if len(images) > index_end >= 0:
            break

    if isinstance(images[index], list):
        return images[index], timesteps[index]
    else:
        return [images[index]], [timesteps[index]]


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


def read_ion_paird_from_hdf5(
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
