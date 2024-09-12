import pytest
from pathlib import Path
import soluanalysis as solu
from soluanalysis.io import read_lammps_dump, write_lammps_dump
from soluanalysis.hdf5_io import read_ion_pairs_from_hdf5, save_ion_pairs_to_hdf5
import numpy as np


def test_average_frames():
    """
    Tests that you can read in a trajectory in "chunks", and also average frames, thereby creating a
    new list of System objects.
    """
    read_chunk_size = 10

    n_avg_frames = 5

    n_frames = 11  # known here, but there could be a function to get this?

    test_dir = Path(__file__).resolve().parent
    infilename = test_dir / "../resources/Fe3_water_cluster.lammpstrj"

    averaged_systems = []
    timesteps_avg = []

    for ichunk in range(0, n_frames, read_chunk_size):
        chunk_start = ichunk
        chunk_end = ichunk + read_chunk_size
        slice_str = str(chunk_start) + ":" + str(chunk_end)
        # Read in the trajectory (the last index in the slice is non inclusive)
        systems, timesteps = solu.io.read_lammps_dump(infilename, slice_str)

        n_frames_chunk = len(systems)

        # Now you have chunks of size 10
        # average every 5 frames
        for jchunk in range(0, n_frames_chunk, n_avg_frames):
            avg_chunk_start = jchunk
            avg_chunk_end = avg_chunk_start + n_avg_frames
            if avg_chunk_end > n_frames_chunk:
                avg_chunk_end = n_frames_chunk
            # Get the system (assume that all atoms in the correct order TODO: handle incorrect order)
            current_system = systems[avg_chunk_start]
            avg_pos = np.array(current_system.collect_positions())
            # Number of frames in averaging chunk
            n_chunk_avg = avg_chunk_end - avg_chunk_start

            for frame in range(avg_chunk_start + 1, avg_chunk_end, 1):
                avg_pos = avg_pos + np.array(systems[frame].collect_positions())

            # Get the mean by dividing the number
            avg_pos = avg_pos / n_chunk_avg
            # Reset the positions
            current_system.reset_positions(avg_pos)
            # save the timesteps and the averaged system
            timesteps_avg.append(timesteps[avg_chunk_start])
            averaged_systems.append(current_system)

    # Check that the timesteps are what are expected
    assert timesteps_avg == [0, 25, 50]
    assert averaged_systems[0].atoms[0].position == pytest.approx(
        [22.50546, 31.19614, 29.12514]
    )
    assert averaged_systems[1].atoms[0].position == pytest.approx(
        [22.58648, 31.1157, 29.0308]
    )
    assert averaged_systems[2].atoms[0].position == pytest.approx(
        [22.64, 31.0671, 28.9615]
    )
