from spirit_extras.plotting import Paper_Plot
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from PIL import Image
import os
from os.path import isfile, join
import matplotlib.patheffects as path_effects
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerTuple

def insert_inset(image_path, axis, rel_height, 
    rel_width, margin_x, margin_y, x_align, y_align):
    image = pplot.open_image(image_path)  # Read the image as a numpy array
    # image = pplot.crop_to_content(image)  # crop the image to content
    # image = pplot.crop(image, 950,950)

    image = pplot.replace_background_color(image, [0,0,0,0], background_color=None)

    # Create an axis to hold the inset
    # The inset shall be located in the top right corner
    ax_inset = pplot.create_inset_axis(
        containing_ax=axis,
        rel_height=rel_height,
        rel_width=rel_width,
        margin_x=margin_x,
        margin_y=margin_y,
        x_align=x_align,
        y_align=y_align,
    )
    return ax_inset, image

# Figure attributes 
hbond_example_dir = Path(__file__).resolve().parent
input_dir = hbond_example_dir / "output"
fname = "tcf.txt"

# Lengths in inches in general 
CM = Paper_Plot.cm  # Use this to give lengths in cm...

# Params 
params = {
"font.size": 8,
"font.family": ("Arial","sans-serif"),
"mathtext.fontset": "dejavuserif",
"xtick.labelsize": 7,
"ytick.labelsize": 7,
"axes.labelsize": 8
}

# For legend labels
font = {'family': ("Arial","sans-serif"),
        'weight': 'normal',
        'size': 7.5,
        }

# Either set the absolute height here or set the aspect ratio
# inside apply_absolute_margins

pplot = Paper_Plot(
    width=3.25, height=2.75,nrows=1, ncols=1, rcParams=params
) 

# Vertical margin: bottom and then top 
# horizontal margin: left and then right
# Golden ratio aspect ratio=1.618 
# hspace -> height space 

pplot.apply_absolute_margins(
    aspect_ratio=None,
    abs_horizontal_margins=[1.3 * CM, 0.12 * CM],
    abs_vertical_margins=[1.0 * CM, 0.5 * CM],
    abs_wspace=0.0 * CM,
    abs_hspace=0.0 * CM,
)

print(pplot.info_string())

# Get the figure and gridpsec objects
fig = pplot.fig()
gs = pplot.gs()

# DATA 
# Get the data (with counterions)
infile = os.path.join(input_dir,fname)
data = np.loadtxt(infile,skiprows=1) # Load the text file 
tau_val = np.array(data[:,0]) # Tau values in fs
tau_val = tau_val/1000 # Tau values in ps 
tcf_val =  np.array(data[:,1]) # Time correlation function values  
tcf_val_err = np.array(data[:,2]) # Error bars for the TCF  

# ---------------------------------------------
ax1 = fig.add_subplot(gs[0,0])
ax1.set_xlabel(r'$\tau$ (ps)')
ax1.set_ylabel(r'$\mathrm{C_{HB}} (\tau)$',labelpad=7)  # we already handled the x-label with ax1

# Data points and line for non-octahedral state 
(tcf_1,) = ax1.plot(tau_val, tcf_val, marker=".", markersize=6,
    color="orangered",markeredgecolor='black',markeredgewidth=0.6,zorder=4,
    label=r'continuous bond')
# Shaded error region for non-octahedral states 
ax1.fill_between(tau_val, tcf_val-tcf_val_err, tcf_val+tcf_val_err,linewidth=0.1,color="peachpuff", alpha=0.8,zorder=0)

# Horizontal line through 0.0 
plt.axhline(y = 0.0, color = 'black', linestyle = '--', linewidth=1) 

# ax1.set_xlim([0,2.25]) # in percent 
# xtick_vec = np.arange(0.0,2.5,0.5)
# xtick_vec = np.append(xtick_vec, 2.25)
# ax1.set_xticks(xtick_vec)
# ax1.set_ylim([0.0,25])
ax1.set_xlim([0.0,10])

# ---------------------------------------------------------------------

# LEGEND

ax1.legend(handles=[tcf_1], fontsize=7.4)

# ---------------------------------------------------------------------

# plt.show()
fig.savefig(input_dir/'tcf.png', dpi=300)
