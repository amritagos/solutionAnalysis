__all__ = ["read_lammps_dump"]

import numpy as np
from typing import List, Optional
from pathlib import Path
from collections import deque
from os.path import splitext
from soluanalysis.solvlib import Atom, System


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


def read_lammps_dump(file_path: Path, index=-1) -> List[System]:
    """
    Reads a LAMMPS trajectory file and returns a list of System objects for the specified range of timesteps.
    Similar to the ASE read function. Does not have support for triclinic boxes.

    Args:
        file_path (pathlib.Path): Path to the LAMMPS trajectory file.
        index: integer or slice object (by default, gets the last timestep)

    Returns:
        List[System], List[int]: A list of System objects representing each frame in the trajectory within the specified range, and a list of timesteps read from the LAMMPS dump file.
    """
    n_atoms = 0
    images = []
    box_bounds = None
    lower_box_limits = None
    timesteps = []

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

    return images[index], timesteps[index]
