# linwalker/plotting.py

import matplotlib.pyplot as plt
from .palette import SOURCE_COLOURS

def plot_diversification(div_df, title=None):
    """Plot diversification curves from linwalker.diversification.lin_diversification output."""
    plt.figure(figsize=(8, 6))
    for grp, sub in div_df.groupby("group"):
        plt.plot(
            sub["LIN_level"],
            sub["n_unique_LINs"],
            label=grp,
            color=SOURCE_COLOURS.get(str(grp).lower(), "#000000"),
            linewidth=2.8
        )
    plt.xlabel("LIN threshold")
    plt.ylabel("Number of unique LIN IDs")
    if title:
        plt.title(title)
    plt.legend(title="Source", frameon=False)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.show()
