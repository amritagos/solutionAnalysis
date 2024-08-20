from pathlib import Path
import gdown

hbond_example_dir = Path(__file__).resolve().parent
input_folder = hbond_example_dir / "input"
# Create the folder if it does not exist 
input_folder.mkdir(parents=True, exist_ok=True)

# URL for the trajectory file (571 MB)
file_id = '15sv3cc_Mx281C1RHDTzUiqMoVBjwcWNu'  
url = f"https://drive.google.com/uc?id={file_id}"
# Output file path 
output_file = input_folder / 'dump-tcf.lammpstrj'
# Download the file into the 'input' folder (if it does not already exist)
if not output_file.exists():
    gdown.download(url, str(output_file), quiet=False)