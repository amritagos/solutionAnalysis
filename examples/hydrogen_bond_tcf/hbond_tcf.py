from pathlib import Path
import gdown
import soluanalysis as solu
from soluanalysis.io import read_lammps_dump
import numpy as np

hbond_example_dir = Path(__file__).resolve().parent
input_folder = hbond_example_dir / "input"
# Create the folder if it does not exist 
input_folder.mkdir(parents=True, exist_ok=True)
# Create an output folder if it does not exist
output_folder = hbond_example_dir / "output"
output_folder.mkdir(parents=True, exist_ok=True)
output_file =  output_folder / 'tcf.txt'

# URL for the trajectory file (571 MB)
file_id = '15sv3cc_Mx281C1RHDTzUiqMoVBjwcWNu'  
url = f"https://drive.google.com/uc?id={file_id}"
# File path for the trajectory 
input_file = input_folder / 'dump-tcf.lammpstrj'
# Download the file into the 'input' folder (if it does not already exist)
if not input_file.exists():
    gdown.download(url, str(input_file), quiet=False)

# Read in the trajectory
systems, timesteps = read_lammps_dump(input_file, ':')
print("Trajectory read.\n")
# General system information
fe_type = 3
o_type = 1
h_type = 2
cl_type = 4
# Cutoffs, types etc needed for hydrogen bonds
donor_atom_types = [o_type]
acceptor_atom_types = [cl_type, o_type]
h_atom_types = [h_type]
donor_acceptor_cutoff = 3.2
max_angle_deg = 30  # in degrees 
# List that will hold UndirectedNetwork objects corresponding to hydrogen bonds
networks = [] 

# Loop through the trajectory 
for system in systems:
    n_atoms = system.n_atoms()
    network = solu.graphlib.UndirectedNetwork(n_atoms) # Create the network
    # We will ignore hydrogens so there is no need for intramolecular hydrogen bonds
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
    networks.append(network)
print("Network list generated\n")
# Now that you have a list of UndirectedNetwork objects and the timesteps, 
# the time autocorrelation function can be calculated

# Calculate the time correlation function (using default values of start_t0 etc)
tau_values, tcf_avg, tcf_error = solu.james.time_correlation_function(
        networks, timesteps, 0, 1, 1, None
    )
print("Calculated the time correlation function\n")
# Write these out to a CSV 
# Write out to file  
header_string = 'tau\ttcf\ttcf_stderr'
np.savetxt(output_file, np.column_stack((tau_values, tcf_avg, tcf_error)), delimiter=' ', header = header_string)