from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
import pandas as pd
import unique_count
from pathlib import Path
import numpy as np

if __name__ == "__main__":
  data = unique_count.read_data(Path("outfiles/correlation_file_colon.csv"),
                                delimiter=",")

  for strength in ["strong", "weak", "garbage"]:
    for percentage in range(0, 100, 10):
      new_data = data[data["strength"] == strength]
      ydata = (
        [list(new_data[f"m1 Pearson Length: {percentage} %"])[2]] +
        [list(new_data[f"m1 Pearson Length: {percentage} %"])[4]]
      )
      yerr = np.array([
        [y - x for x, y in zip([list(new_data[f"m1 Pearson Length Conf Low: {percentage} %"])[2]],
                                [list(new_data[f"m1 Pearson Length: {percentage} %"])[2]])] +
        [y - x for x, y in zip([list(new_data[f"m1 Pearson Length Conf Low: {percentage} %"])[4]],
                                [list(new_data[f"m1 Pearson Length: {percentage} %"])[4]])]
      ] + [
        [x - y for x, y in zip([list(new_data[f"m1 Pearson Length Conf High: {percentage} %"])[2]],
                                [list(new_data[f"m1 Pearson Length: {percentage} %"])[2]])] +
        [x - y for x, y in zip([list(new_data[f"m1 Pearson Length Conf High: {percentage} %"])[4]],
                                [list(new_data[f"m1 Pearson Length: {percentage} %"])[4]])]
      ])

      fig = plt.figure(figsize=(5, 5))
      ax = SubplotHost(fig, 111)
      fig.add_subplot(ax)

      ax.axhline(0.0, color="grey", ls="--")
      ax.errorbar([1, 2],
                  ydata,
                  yerr=0,
                  fmt="none",
                  capsize=20,
                  color="black")
      ax.errorbar([1, 2],
                  ydata,
                  yerr=yerr,
                  fmt="none",
                  capsize=10,
                  color="black")
      ax.set_xticks([1, 2])
      ax.set_xticklabels([
          "b2m mutant",
          "b2m wildtype"
      ])
      ax.set_xlim(0, 3)

      ax.axis["bottom"].label.set_fontweight("bold")
      plt.tight_layout()
      plt.savefig(f"outputimages/correlation_box_{strength}_{percentage}_new.svg")
