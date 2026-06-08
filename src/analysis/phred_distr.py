import numpy as np
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt

in_dir = Path("in/")
out_dir = Path("out/") 
fig_dir = out_dir / "figs"
fig_dir.mkdir(exist_ok=True)

phreds = np.loadtxt(out_dir / "phreds.txt", dtype = int)

fig, ax = plt.subplots(figsize = (6,4))
sns.kdeplot(phreds, ax = ax, bw_adjust= 5)
ax.set_title("Phred Scores KDE")
plt.tight_layout()
plt.savefig( fig_dir / "phred_dirstr.png")

