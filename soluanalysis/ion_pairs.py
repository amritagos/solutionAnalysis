from typing import Dict, List, Optional, Tuple
import soluanalysis as solu
import numpy as np


def get_ion_pairs_time_series(
    systems: List[solu.james.System],
    timesteps: List[int],
    donor_atom_types: List[int],
    acceptor_atom_types: List[int],
    h_atom_types: List[int],
    source_atom_types: List[int],
    destination_atom_types: List[int],
    intermediate_atom_types: List[int],
    max_depth: Optional[int],
    identifier: solu.james.WriteIdentifier = solu.james.WriteIdentifier.AtomID,
    donor_acceptor_cutoff: float = 3.2,
    max_angle_deg: float = 30,
    pairs: Optional[List[solu.james.Pair]] = None,
    cutoffs: Optional[List[float]] = None,
) -> Dict[int, Dict[int, List[List[int]]]]:
    """Gets ion pairs for each frame (System object) in a trajectory, sorting them according to length as well.
    Hydrogens are "ignored", in the sense that they are not saved as connections, although they are used in the
    hydrogen bond criterion. This reduces the number of connections in the UndirectedNetwork object.
    First, (optionally), bonds/connections can be created using distance-based cutoffs. Subsequently, hydrogen bonds
    are found. Then, ion pairs can be found using the connectivity information in the UndirectedNetwork object.

    Args:
        systems (List[solu.james.System]): List of System objects in the trajectory with atom ID, atom type, molecule ID, etc. information.
        timesteps (List[int]): A list of timesteps in the trajectory
        donor_atom_types (List[int]): The atom types for donors (required for hydrogen bonding)
        acceptor_atom_types (List[int]): The atom types for acceptors (required for hydrogen bonding)
        h_atom_types (List[int]): Atom types for the hydrogen atoms (also required for hydrogen bonding)
        source_atom_types (List[int]): The atom type(s) which correspond to the source atoms (required for the ion pair search)
        destination_atom_types (List[int]): The atom type(s) corresponding to the destination atoms. Currently does not
        differentiate between different types of destination atoms. (required for the ion pair search)
        intermediate_atom_types (List[int]): Atom type(s) corresponding to atom types allowed for the intermediate atoms
        in between the source and destination. For instance, these would correspond to atom types in water molecules.
        max_depth (Optional[int]): The maximum depth (or number of successive neighbours) upto which the bread-first-search
        will be conducted in the BFS.
        identifier (solu.james.WriteIdentifier, optional): Decides whether the elements in the ion pairs correspond
        to atom IDs in the System object, or to indices in the System object. Defaults to solu.james.WriteIdentifier.AtomID.
        donor_acceptor_cutoff (float, optional): Cutoff for the donor-acceptor distance, in the geometric hydrogen bond criterion. Defaults to 3.2.
        max_angle_deg (float, optional): Maximum angle cutoff (in degrees) for the HDA (hydrogen-donor-acceptor) angle. Defaults to 30.
        pairs (Optional[List[solu.james.Pair]], optional): Pair objects corresponding to a pair for which a distance
        based cutoff can be defined, for finding distance-based bonds. Defaults to None.
        cutoffs (Optional[List[float]], optional): Cutoffs for finding distance-based bonds. Defaults to None.

    Returns:
        Dict[int, Dict[int, List[List[int]]]]: The dictionary contains the time series information about ion pairs, such that
        the keys of the outer dictionary are the timesteps, and the inner dictionaries contain information about
        the ion pairs per timestep. The keys of the inner dictionary are the lengths of the ion pairs, and the
        values are lists of lists corresponding to the ion pairs.
    """
    output_dict = {}  # contains the ion pairs

    for i, (system, timestep) in enumerate(zip(systems, timesteps)):
        ion_pair_list = []
        n_atoms = system.n_atoms()
        # Construct the network and fill it with bonds (distance-based and hydrogen bonds)
        network = solu.graphlib.UndirectedNetwork(n_atoms)

        # Optionally create distance-based bonds
        if pairs is not None and cutoffs is not None:
            solu.james.add_distance_based_bonds(network, system, pairs, cutoffs)

        # Create and find the hydrogen bonds (ignoring the hydrogens)
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
        # Find the source indices (Fe in this case)
        source_indices = [
            system.atoms.index(atom)
            for atom in system.atoms
            if atom.type in source_atom_types
        ]

        # Get the ion pairs for each source index
        # The ion pairs will be saved with elements decided by the identifier
        for source in source_indices:
            ion_pairs = solu.james.find_ion_pairs(
                source,
                network,
                system,
                destination_atom_types,
                intermediate_atom_types,
                max_depth,
                identifier,
            )
            ion_pair_list += ion_pairs
        groups = {}
        for ion_pair in ion_pair_list:
            groups.setdefault(len(ion_pair), []).append(ion_pair)
        output_dict[timestep] = groups

    return output_dict


def get_end_point_indices_ion_pair(
    ion_pair: List[int],
    identifier: solu.james.WriteIdentifier,
    system: solu.james.System,
) -> Tuple[int, int]:
    """Obtain the indices (in the System object, and therefore in the UndirectedNetwork object) of the end points of the ion pair (which are connected)

    Args:
        ion_pair (List[int]): The ion pair in question whose end points are required
        identifier (solu.james.WriteIdentifier): Enum with a value that reveals whether the ion pair elements are indices or atom IDs
        system (solu.james.System): Representative System object

    Returns:
        Tuple[int, int]: First index, last index
    """
    if len(ion_pair) < 2:
        raise AttributeError("The length of the ion pair is less than 2")
    ele_first = ion_pair[0]
    ele_last = ion_pair[-1]

    if identifier == solu.james.WriteIdentifier.AtomID:
        index_first = system.index_from_id(ele_first)
        index_last = system.index_from_id(ele_last)
        if index_first is None or index_last is None:
            raise Exception(f"The atom IDs {ele_first} and {ele_last} do not exist.\n")
        else:
            return index_first, index_last
    # Atom indices are in the ion pairs
    else:
        return ele_first, ele_last


def network_from_ion_pairs(
    ion_pair_info: Dict[int, List[List[int]]],
    ion_pair_length: int,
    system: solu.james.System,
    identifier: solu.james.WriteIdentifier,
) -> solu.graphlib.UndirectedNetwork:
    """Generates a network with the connectivity information of the two end points of ion pairs of a
      prescribed length, at a particular timestep

    Args:
        ion_pair_info (Dict[int, List[List[int]]]): Keys are the lengths of ion pairs,
        and the values are lists of lists (corresponding to ion pairs)
        ion_pair_length (int): The desired length for the ion pairs (other ion pairs will be ignored)
        system (solu.james.System): Representative System object for the frame
        identifier (solu.james.WriteIdentifier): Enum with a value that reveals whether the ion pair elements are indices or atom IDs

    Returns:
        solu.graphlib.UndirectedNetwork: UndirectedNetwork object into
        which the ion pair connectivity will be stored
    """
    # Initialize the UndirectedNetwork object
    network = solu.graphlib.UndirectedNetwork(system.n_atoms())
    weight_edge = 1

    # Get the ion pairs corresponding to the desired ion pair length
    for ion_pair in ion_pair_info.get(ion_pair_length) or []:
        index_first, index_last = get_end_point_indices_ion_pair(
            ion_pair, identifier, system
        )
        network.push_back_neighbour_and_weight(index_first, index_last, weight_edge)
    # Network corresponding to desired ion pair length
    return network


def networks_from_ion_pair_series(
    ion_pair_series: Dict[int, Dict[int, List[List[int]]]],
    ion_pair_length: int,
    system: solu.james.System,
    identifier: solu.james.WriteIdentifier,
) -> List[solu.graphlib.UndirectedNetwork]:
    """Generates a list of networks with the connectivity information of the two end points of ion pairs of a
      prescribed length, for a time series

    Args:
        ion_pair_series (Dict[int, Dict[int, List[List[int]]]]): Contains the time series information about ion pairs, such that
        the keys of the outer dictionary are the timesteps, and the inner dictionaries contain information about
        the ion pairs per timestep. The keys of the inner dictionary are the lengths of the ion pairs, and the
        values are lists of lists corresponding to the ion pairs.
        ion_pair_length (int): The desired length for the ion pairs (other ion pairs will be ignored)
        system (solu.james.System): Representative System object containing the number of atoms, atom ID information etc.
        identifier (solu.james.WriteIdentifier): Enum with a value that reveals whether the ion pair elements are indices or atom IDs

    Returns:
        List[solu.graphlib.UndirectedNetwork]: List of UndirectedNetwork objects containing
        ion pair connectivity information
    """
    # Extract timesteps from the keys of time_series_dict and sort
    timesteps = sorted(ion_pair_series.keys())
    # Output list with UndirectedNetwork objects
    network_series = []

    for timestep in timesteps:
        network_series.append(
            network_from_ion_pairs(
                ion_pair_series[timestep], ion_pair_length, system, identifier
            )
        )

    return network_series
