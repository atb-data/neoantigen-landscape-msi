import unique_count
import count_epitopes
import excel_processing
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from pathlib import Path
from scipy import stats
from scipy.optimize import curve_fit
import numpy as np

from unique_count import parse_renamings, rename_candidate

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

def textonly(ax, txt, fontsize = 14, loc = 2, *args, **kwargs):
  at = AnchoredText(txt,
                    prop=dict(size=fontsize), 
                    frameon=True,
                    loc=loc)
  at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
  ax.add_artist(at)
  return at

def scatter_show(name, row_1, row_2, c=None, names=None, renamings=None, topn=9, output="outputimages/scatter.svg"):
  fig, ax = plt.subplots(ncols=2,
                          figsize=(10, 5),
                          gridspec_kw = {
                              'width_ratios': [8, 1]
                          })
  if c != None and len(c) > 0:
    max_top = len(c)
    true_top = min(topn, max_top)
    top_c = sorted(c)[-true_top]
    count_id = 1
    leg_string = "\n"
    tops = []
    s = []
    for idx, name in enumerate(names):
      true_name = rename_candidate(name, renamings=renamings)
      if c[idx] >= top_c:
        tops.append((c[idx], idx, true_name))
        s.append(20)
      else:
        s.append(10)
    sc = ax[0].scatter(row_1, row_2, c=c, cmap="Reds", edgecolors="black")
    tops.sort(key=lambda x: -x[0])
    for val, idx, true_name in tops:
      leg_string += f"{count_id} {true_name}\n"
      count_id += 1
    textonly(ax[1], leg_string, fontsize=12)

  else:
    sc = plt.scatter(row_1, row_2)
  if len(row_1) == 0:
    plt.close()
    return
  ax[1].set_axis_off()
  ax[0].set_xlim(left=-0.1, right=1.1)
  ax[0].set_ylim(bottom=-0.1, top=0.6)
  ax[0].set_xticks([0.1 * x for x in range(11)])
  ax[0].set_yticks([0.1 * x for x in range(6)])
  ax[0].set_xticklabels([x for x in range(0, 110, 10)])
  ax[0].set_yticklabels([x for x in range(0, 60, 10)])

  for label in ax[0].get_xticklabels():
    label.set_fontsize(16)

  for label in ax[0].get_yticklabels():
    label.set_fontsize(16)

  divider = make_axes_locatable(ax[0])
  cax = divider.append_axes("right", size = "1%", pad = 0.2)
  cbar = fig.colorbar(sc, cax = cax, orientation = "vertical")
  cbar.ax.tick_params(labelsize=10)
  cbar.set_ticks([idx * 0.1 for idx in range(11)])
  cbar.set_ticklabels([f"{idx}" for idx in range(0, 110, 10)])

  plt.tight_layout()
  plt.savefig(output)
  plt.close()

def scatter_show_e(name, row_1, row_2, c=None, names=None, renamings=None, topn=9, output="outputimages/scatter.svg"):
  fig, ax = plt.subplots(ncols=2,
                          figsize=(10, 5),
                          gridspec_kw = {
                              'width_ratios': [8, 1]
                          })
  if c != None and len(c) > 0:
    max_top = len(c)
    true_top = min(topn, max_top)
    top_c = sorted(c)[-true_top]
    sc = ax[0].scatter(row_1, row_2, c=c, cmap="Reds", edgecolors="black")
    count_id = 1
    leg_string = "\n"
    tops = []
    for idx, name in enumerate(names):
        true_name = rename_candidate(name, renamings=renamings)
        if c[idx] >= top_c:
            tops.append((c[idx], idx, true_name))
    tops.sort(key=lambda x: -x[0])
    for val, idx, true_name in tops:
        leg_string += f"{count_id} {true_name}\n"
        count_id += 1
    textonly(ax[1], leg_string, fontsize=12)
  else:
    sc = plt.scatter(row_1, row_2)
  if len(row_1) == 0:
    plt.close()
    return
  ax[1].set_axis_off()
  ax[0].set_xlim(left=0, right=200)
  ax[0].set_ylim(bottom=-0.1, top=0.6)
  ax[0].set_xticks([x for x in range(0, 200, 20)])
  ax[0].set_yticks([0.1 * x for x in range(6)])
  ax[0].set_xticklabels([x for x in range(0, 200, 20)])
  ax[0].set_yticklabels([x for x in range(0, 60, 10)])

  divider = make_axes_locatable(ax[0])
  cax = divider.append_axes("right", size = "1%", pad = 0.2)
  cbar = fig.colorbar(sc, cax = cax, orientation = "vertical")
  cbar.ax.tick_params(labelsize=10)
  cbar.set_ticks([idx * 0.1 for idx in range(11)])
  cbar.set_ticklabels([f"{idx}" for idx in range(0, 110, 10)])

  plt.tight_layout()
  plt.savefig(output)
  plt.close()

def scatter_show_compare(name, row_1, row_2, row_11, row_22, c=None, cc=None, names=None, renamings=None, topn=9, output="outputimages/scatter.svg"):
  fig, ax = plt.subplots(ncols=2,
                          figsize=(10, 5),
                          gridspec_kw = {
                              'width_ratios': [8, 1]
                          })
  if c != None and len(c) > 0:
    max_top = len(c)
    true_top = min(topn, max_top)
    top_c = sorted(c)[-true_top]
    sc = ax[0].scatter(row_1, row_2, c=c, cmap="Reds", edgecolors="black")
    sc1 = ax[0].scatter(row_11, row_22, c=cc, cmap="Blues", marker='s', edgecolors="black")
    count_id = 1
    leg_string = "\n"
    tops = []
    for idx, name in enumerate(names):
      true_name = rename_candidate(name, renamings=renamings)
      if c[idx] >= top_c:
        tops.append((c[idx], idx, true_name))
    tops.sort(key=lambda x: -x[0])
    for val, idx, true_name in tops:
      ax[0].arrow(row_1[idx], row_2[idx], row_11[idx] - row_1[idx], row_22[idx] - row_2[idx])
      leg_string += f"{count_id} {true_name}\n"
      count_id += 1
    textonly(ax[1], leg_string, fontsize=12)
  else:
    sc = plt.scatter(row_1, row_2)
  if len(row_1) == 0:
    plt.close()
    return
  ax[1].set_axis_off()
  ax[0].set_xlim(left=-0.1, right=1.1)
  ax[0].set_ylim(bottom=-0.1, top=0.6)
  ax[0].set_xticks([0.1 * x for x in range(11)])
  ax[0].set_yticks([0.1 * x for x in range(6)])
  ax[0].set_xticklabels([x for x in range(0, 110, 10)])
  ax[0].set_yticklabels([x for x in range(0, 60, 10)])

  for label in ax[0].get_xticklabels():
    label.set_fontsize(16)

  for label in ax[0].get_yticklabels():
    label.set_fontsize(16)

  divider = make_axes_locatable(ax[0])
  cax = divider.append_axes("right", size = "1%", pad = 0.2)
  cbar = fig.colorbar(sc, cax = cax, orientation = "vertical")
  cbar.ax.tick_params(labelsize=10)
  cbar.set_ticks([idx * 0.1 for idx in range(11)])
  cbar.set_ticklabels([f"{idx}" for idx in range(0, 110, 10)])

  cax1 = divider.append_axes("right", size = "1%", pad = 0.6)
  cbar1 = fig.colorbar(sc1, cax = cax1, orientation = "vertical")
  cbar1.ax.tick_params(labelsize=10)
  cbar1.set_ticks([idx * 0.1 for idx in range(11)])
  cbar1.set_ticklabels([f"{idx}" for idx in range(0, 110, 10)])

  plt.tight_layout()
  plt.savefig(output)
  plt.close()

def scatter_show_l(name, row_1, row_2, c=None, names=None, renamings=None, topn=9, output="outputimages/scatter.svg"):
  fig, ax = plt.subplots(
    ncols=2,
    figsize=(10, 5),
    gridspec_kw = {'width_ratios': [8, 1]}
  )
  if c != None and len(c) > 0:
    max_top = len(c)
    true_top = min(topn, max_top)
    top_c = sorted(c)[-true_top]
    sc = ax[0].scatter(row_1, row_2, c=c, cmap="Reds", edgecolors="black")
    count_id = 1
    leg_string = "\n"
    tops = []
    for idx, name in enumerate(names):
      true_name = rename_candidate(name, renamings=renamings)
      if c[idx] >= top_c:
        tops.append((c[idx], idx, true_name))
    tops.sort(key=lambda x: -x[0])
    for val, idx, true_name in tops:
      leg_string += f"{count_id} {true_name}\n"
      count_id += 1
    textonly(ax[1], leg_string, fontsize=12)
  else:
    sc = plt.scatter(row_1, row_2)
  if len(row_1) == 0:
    plt.close()
    return
  ax[1].set_axis_off()
  ax[0].set_xlim(left=0, right=200)
  ax[0].set_ylim(bottom=-0.1, top=1.1)
  ax[0].set_xticks([x for x in range(0, 200, 20)])
  ax[0].set_yticks([0.1 * x for x in range(11)])
  ax[0].set_xticklabels([x for x in range(0, 200, 20)])
  ax[0].set_yticklabels([x for x in range(0, 110, 10)])

  divider = make_axes_locatable(ax[0])
  cax = divider.append_axes("right", size = "1%", pad = 0.2)
  cbar = fig.colorbar(sc, cax = cax, orientation = "vertical")
  cbar.ax.tick_params(labelsize=10)
  cbar.set_ticks([idx * 0.1 for idx in range(11)])
  cbar.set_ticklabels([f"{idx}" for idx in range(0, 110, 10)])

  plt.tight_layout()
  plt.savefig(output)
  plt.close()

def scatter_show_q(name, row_1, row_2, c=None, names=None, renamings=None, topn=4, output="outputimages/scatter.svg"):
  fig, ax = plt.subplots()
  if c != None:
    max_top = len(c)
    true_top = min(topn, max_top)
    top_c = sorted(c)[-true_top]
    sc = plt.scatter(row_1, row_2, c=c, cmap="Reds", edgecolors="black")
    count_id = 1
    leg_string = "\n"
    tops = []
    for idx, name in enumerate(names):
      true_name = rename_candidate(name, renamings=renamings)
      if c[idx] >= top_c:
          tops.append((c[idx], idx, true_name))
    tops.sort(key=lambda x: -x[0])
    for val, idx, true_name in tops:
      text = ax.annotate(str(count_id), xy=(row_1[idx], row_2[idx]),
                         xytext=(row_1[idx] + 0.01, row_2[idx] + 0.01),
                         rotation=45, ha="left", va="bottom")
      leg_string += f"{count_id} {true_name}\n"
      count_id += 1
    textonly(ax, leg_string, fontsize=12)
  else:
    sc = plt.scatter(row_1, row_2)
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size = "1%", pad = 0.2)
  cbar = fig.colorbar(sc, cax = cax, orientation = "vertical")
  cbar.ax.tick_params(labelsize=10)
  cbar.set_ticks([idx * 0.1 for idx in range(11)])
  cbar.set_ticklabels([f"{idx}" for idx in range(0, 110, 10)])

  plt.tight_layout()
  plt.savefig(output)
  plt.close()

def scatter_show_m(name, row_1, row_2, c=None, irs=None, output="outputimages/scatter.svg"):
  fig, ax = plt.subplots()
  row_n_1 = []
  row_n_2 = []
  c_n = []
  if irs != None:
    for idx, elem in enumerate(irs):
      if elem > 0.02:
        row_n_1.append(row_1[idx])
        row_n_2.append(row_2[idx])
        c_n.append(c[idx])
  row_n_1 = np.array(row_n_1)
  row_n_2 = np.array(row_n_2)
  c_n = np.array(c_n)
  if c != None:
    sc = plt.scatter(row_n_1, row_n_2, c=c_n, cmap="seismic", edgecolors="black")
  else:
    sc = plt.scatter(row_n_1, row_n_2)
  ax.set_xlim(left=-0.1, right=1.1)
  ax.set_ylim(bottom=-0.1, top=0.6)
  ax.set_xticks([0.1 * x for x in range(11)])
  ax.set_yticks([0.1 * x for x in range(6)])
  ax.set_xticklabels([x for x in range(0, 110, 10)])
  ax.set_yticklabels([x for x in range(0, 60, 10)])

  plt.tight_layout()
  plt.savefig(output)
  plt.close()

def fit_show(name, x, y, func):
  fig, ax = plt.subplots()
  plt.scatter(x, y)
  try:
    popt, pcov = curve_fit(func, x, y, bounds=(0.0, [10.0, 10.0, 10.0]))
    minx = min(x)
    maxx = max(x)
    xspace = np.linspace(1e-6, maxx, 1000)
    plt.plot(xspace, func(xspace, *popt))
  except:
    print("no fit found!")
  plt.show()
  plt.close()

if __name__ == "__main__":
    import math

    def r_to_z(r):
        return math.log((1 + r) / (1 - r)) / 2.0

    def z_to_r(z):
        e = math.exp(2 * z)
        return((e - 1) / (e + 1))

    def r_confidence_interval(r, alpha, n):
        z = r_to_z(r)
        se = 1.0 / math.sqrt(n - 3)
        z_crit = stats.norm.ppf(1 - alpha/2)  # 2-tailed z critical value

        lo = z - z_crit * se
        hi = z + z_crit * se

        # Return a sequence
        return (z_to_r(lo), z_to_r(hi))

    renamings = parse_renamings("Renamings.csv", delimiter=",")

    for location in ["colon", "endometrium"]:
        raw_data = unique_count.read_data(Path(f"table_dump_{location}.csv"))
        raw_data = raw_data[raw_data["cohort"] == "European Caucasian"]
        strengths = ["strong", "weak", "garbage"]
        output_data = f"outfiles/correlation_file_{location}.csv"

        peptide_lengths = unique_count.read_peptide_lengths("NetMHCPanInput/cMNR_peptides.csv")

        render = True
        full_data = {}
        significant_anti_correlation = {}

        rows = []
        header = [
            "strength", "filter", "relevance"
        ]
        for p in range(0, 100, 10):
            header += [
                f"Quotient Spearman: {p} %", f"Quotient Spearman p-value: {p} %",
                f"Quotient Pearson: {p} %", f"Quotient Pearson p-value: {p} %"
            ]
        for p in range(0, 100, 10):
            header += [
                f"m1 Spearman: {p} %", f"m1 Spearman p-value: {p} %",
                f"m2 Spearman: {p} %", f"m2 Spearman p-value: {p} %",
                f"All Spearman: {p} %", f"All Spearman p-value: {p} %",
                f"m1 Pearson: {p} %", f"m1 Pearson p-value: {p} %",
                f"m1 Pearson Conf Low: {p} %",  f"m1 Pearson Conf High: {p} %",
                f"m2 Pearson: {p} %", f"m2 Pearson p-value: {p} %",
                f"m2 Pearson Conf Low: {p} %",  f"m2 Pearson Conf High: {p} %",
                f"All Pearson: {p} %", f"All Pearson p-value: {p} %",
                f"All Pearson Conf Low: {p} %",  f"All Pearson Conf High: {p} %",
                f"m1 Pearson Length: {p} %", f"m1 Pearson Length p-value: {p} %",
                f"m1 Pearson Length Conf Low: {p} %",  f"m1 Pearson Length Conf High: {p} %",
                f"m2 Pearson Length: {p} %", f"m2 Pearson Length p-value: {p} %",
                f"m2 Pearson Length Conf Low: {p} %",  f"m2 Pearson Length Conf High: {p} %",
                f"All Pearson Length: {p} %", f"All Pearson Length p-value: {p} %",
                f"All Pearson Length Conf Low: {p} %",  f"All Pearson Length Conf High: {p} %",
            ]
        rows.append(header)

        for strength in range(3):
            subdata = []
            for idx in range(strength + 1):
                subdata.append(raw_data[raw_data["strength"] == strengths[strength]])
            data = pd.concat(subdata)
            data = data.drop_duplicates(subset="candidate")
            data_m1 = data[data["mm"] == 1]
            data_m2 = data[data["mm"] == 2]
            # Characteristics data:
            characteristics = excel_processing.tumor_characteristics("testfiles/CRCQMR.csv")
            # Characteristics filters:
            filters = {
                "all": lambda x: True,
                "hereditary": lambda x: characteristics[x]["hereditary"],
                "b2m_mutated": lambda x: characteristics[x]["b2m"],
                "sporadic": lambda x: not characteristics[x]["hereditary"],
                "b2m_wildtype": lambda x: not characteristics[x]["b2m"],
                "hereditary_b2m_mutated": lambda x: characteristics[x]["hereditary"] and characteristics[x]["b2m"],
                "hereditary_b2m_wildtype": lambda x: characteristics[x]["hereditary"] and (not characteristics[x]["b2m"]),
                "sporadic_b2m_mutated": lambda x: (not characteristics[x]["hereditary"]) and characteristics[x]["b2m"],
                "sporadic_b2m_wildtype": lambda x: not (characteristics[x]["hereditary"] or characteristics[x]["b2m"]),
            }

            compkeys = {
                "all": "all",
                "hereditary" : "sporadic",
                "b2m_mutated" : "b2m_wildtype",
                "sporadic" : "hereditary",
                "b2m_wildtype" : "b2m_mutated",
                "hereditary_b2m_mutated" : "sporadic_b2m_wildtype",
                "hereditary_b2m_wildtype" : "sporadic_b2m_mutated",
                "sporadic_b2m_mutated" : "hereditary_b2m_wildtype",
                "sporadic_b2m_wildtype" : "hereditary_b2m_mutated"
            }

            if location == "endometrium":
                filters = {
                    "all": lambda x: True
                }
            
            strength_level_data = [strengths[strength]]

            for filter_type, filter in filters.items():
                # Mutation data:
                mutation_counts = excel_processing.read_mutation_files("testfiles/QMR1.xlsx",
                                                                    "testfiles/QMR2.xlsx",
                                                                    filter=lambda x: (x in characteristics) and filter(x))
                mutation_comparison = excel_processing.read_mutation_files("testfiles/QMR1.xlsx",
                                                                    "testfiles/QMR2.xlsx",
                                                                    filter=lambda x: (x in characteristics) and filters[compkeys[filter_type]](x))
                if location == "endometrium":
                    mutation_counts = excel_processing.parse_endometrium_mutations(Path("testfiles/ECQMR.csv"))
                    mutation_comparison = mutation_counts
                # HLAs:
                hlas = list(data["hla"].sort_values().unique())
                candidates_m2 = list(data_m2["candidate"].sort_values().unique())
                candidates_m1 = list(map(lambda x: x[:-3] + "_m1", candidates_m2))

                mutations_available = list(mutation_counts.keys())
                mutations_available_comparison = list(mutation_comparison.keys())
                mut_m1 = sorted([elem for elem in mutations_available if elem.endswith("m1")])
                mut_m2 = sorted([elem for elem in mutations_available if elem.endswith("m2")])
                mut_m1_comp = sorted([elem for elem in mutations_available_comparison if elem.endswith("m1")])
                mut_m2_comp = sorted([elem for elem in mutations_available_comparison if elem.endswith("m2")])
                mut_q = []
                for nm1, nm2 in zip(mut_m1, mut_m2):
                    pm1 = mutation_counts[nm1]
                    pm2 = mutation_counts[nm2]
                    num = pm1
                    den = pm2
                    if den > 0.0:
                        mut_q.append(num / den)
                    else:
                        mut_q.append(float("nan"))

                data_by_candidate = {}
                for candidate in mutations_available:
                    data_by_candidate[candidate] = data[data["candidate"] == candidate]
                
                def powerlaw(x, a, b, c):
                    return a * np.power(x, -c) + b

                filter_level_data = strength_level_data + [filter_type]

                for relevance_threshold in range(0, 1, 5):
                    relevance_level_data = filter_level_data + [relevance_threshold]

                    percentage_level_data = relevance_level_data

                    for percentage in range(0, 100, 10):
                        geds_q = []
                        mut_no_nan = []
                        cand_no_nan = []
                        geds_no_nan = []
                        for idx, mut in enumerate(mut_m1):
                            local_data_m1 = data_by_candidate[mut]
                            local_data_m2 = data_by_candidate[mut[:-3] + "_m2"]
                            geds_m1 = list(local_data_m1[f"GEDS: {percentage} %"])[0] if len(list(local_data_m1[f"GEDS: {percentage} %"])) > 0 else 0.0
                            geds_m2 = list(local_data_m2[f"GEDS: {percentage} %"])[0] if len(list(local_data_m2[f"GEDS: {percentage} %"])) > 0 else 0.0
                            irs_m1 = list(local_data_m1[f"IRS: {percentage} %"])[0] if len(list(local_data_m1[f"IRS: {percentage} %"])) > 0 else 0.0
                            irs_m2 = list(local_data_m2[f"IRS: {percentage} %"])[0] if len(list(local_data_m2[f"IRS: {percentage} %"])) > 0 else 0.0
                            
                            num = geds_m1
                            den = geds_m2
                            geds_val = (num / den) if den > 0.0 else float("nan")
                            mut_val = mut_q[idx]
                            geds_q.append(geds_val)
                            if mut_val > 0.0 and geds_val > 0.0 and irs_m1 > relevance_threshold * 0.01:
                                mut_no_nan.append(mut_val)
                                cand_no_nan.append(mut[:-3])
                                geds_no_nan.append(geds_val)
                        spearman = stats.spearmanr(geds_no_nan, mut_no_nan)
                        pearson = stats.pearsonr(geds_no_nan, mut_no_nan)
                        percentage_level_data += [spearman[0], spearman [1],
                                                pearson[0], pearson[1]]
                        if render:
                            try:
                                scatter_show_q("quotients", geds_no_nan, mut_no_nan, c=geds_no_nan, names=cand_no_nan, renamings=renamings,
                                            output=f"outputimages/quotients_location_{location}_filter_{filter_type}_relevance_{relevance_threshold}_{strengths[strength]}_{percentage}.svg")
                            except:
                                pass
                    for percentage in range(0, 100, 10):
                        mutations_m1 = []
                        mutations_m2 = []
                        geds_m1 = []
                        geds_m2 = []
                        mut_data_m1 = data_m1[data_m1[f"IRS: {percentage} %"] >= relevance_threshold * 0.01]
                        mut_data_m2 = data_m2[data_m2[f"IRS: {percentage} %"] >= relevance_threshold * 0.01]
                        mutations_m1 = list(map(lambda x: mutation_counts[x], list(mut_data_m1["candidate"])))
                        mutations_m2 = list(map(lambda x: mutation_counts[x], list(mut_data_m2["candidate"])))
                        mutations_m1_comparison = list(map(lambda x: mutation_comparison[x], list(mut_data_m1["candidate"])))
                        mutations_m2_comparison = list(map(lambda x: mutation_comparison[x], list(mut_data_m2["candidate"])))
                        mutations = mutations_m1 + mutations_m2
                        mutations_comparison = mutations_m1_comparison + mutations_m2_comparison
                        geds_m1 = list(mut_data_m1[f"GEDS: {percentage} %"])
                        geds_m2 = list(mut_data_m2[f"GEDS: {percentage} %"])
                        irs_m1 = [x * y for x, y in zip(geds_m1, mutations_m1)]
                        irs_m2 = [x * y for x, y in zip(geds_m2, mutations_m2)]
                        irs_m1_comparison = [x * y for x, y in zip(geds_m1, mutations_m1_comparison)]
                        irs_m2_comparison = [x * y for x, y in zip(geds_m2, mutations_m2_comparison)]
                        epi_m1 = list(map(lambda x: peptide_lengths[x], list(mut_data_m1["candidate"])))
                        epi_m2 = list(map(lambda x: peptide_lengths[x], list(mut_data_m2["candidate"])))
                        epi = epi_m1 + epi_m2

                        irs = irs_m1 + irs_m2
                        geds = geds_m1 + geds_m2
                        mness = [-1] * len(irs_m1) + [1] * len(irs_m2)
                        if render:
                            # EXPERIMENTAL:
                            scatter_show_compare("m1 vs geds m1", geds_m1, mutations_m1, geds_m1, mutations_m1_comparison, c=irs_m1, cc=irs_m1_comparison, names=list(mut_data_m1["candidate"]), renamings=renamings,
                                        output=f"outputimages/m1_v_geds_comparison_location_{location}_filter_{filter_type}_relevance_{relevance_threshold}_{strengths[strength]}_{percentage}.svg")
                            scatter_show_compare("m2 vs geds m2", geds_m2, mutations_m2, geds_m2, mutations_m2_comparison, c=irs_m2, cc=irs_m2_comparison, names=list(mut_data_m2["candidate"]), renamings=renamings,
                                        output=f"outputimages/m2_v_geds_comparison_location_{location}_filter_{filter_type}_relevance_{relevance_threshold}_{strengths[strength]}_{percentage}.svg")
                            
                            scatter_show("m1 vs geds m1", geds_m1, mutations_m1, c=irs_m1, names=list(mut_data_m1["candidate"]), renamings=renamings,
                                        output=f"outputimages/m1_v_geds_location_{location}_filter_{filter_type}_relevance_{relevance_threshold}_{strengths[strength]}_{percentage}.svg")
                            scatter_show("m2 vs geds m2", geds_m2, mutations_m2, c=irs_m2, names=list(mut_data_m2["candidate"]), renamings=renamings,
                                        output=f"outputimages/m2_v_geds_location_{location}_filter_{filter_type}_relevance_{relevance_threshold}_{strengths[strength]}_{percentage}.svg")
                            scatter_show_m("mut vs geds", geds, mutations, c=mness, irs=irs,
                                        output=f"outputimages/m_v_geds_location_{location}_filter_{filter_type}_relevance_{relevance_threshold}_{strengths[strength]}_{percentage}.svg")
                            scatter_show_e("m1 vs len", epi_m1, mutations_m1, c=irs_m1, names=list(mut_data_m1["candidate"]), renamings=renamings,
                                        output=f"outputimages/m1_v_len_location_{location}_filter_{filter_type}_relevance_{relevance_threshold}_{strengths[strength]}_{percentage}.svg")
                            scatter_show_e("m2 vs len", epi_m2, mutations_m2, c=irs_m2, names=list(mut_data_m2["candidate"]), renamings=renamings,
                                        output=f"outputimages/m2_v_len_location_{location}_filter_{filter_type}_relevance_{relevance_threshold}_{strengths[strength]}_{percentage}.svg")
                            scatter_show_l("m1 vs len", epi_m1, geds_m1, c=irs_m1, names=list(mut_data_m1["candidate"]), renamings=renamings,
                                        output=f"outputimages/geds_m1_v_len_location_{location}_filter_{filter_type}_relevance_{relevance_threshold}_{strengths[strength]}_{percentage}.svg")
                            scatter_show_l("m1 vs len", epi_m2, geds_m2, c=irs_m2, names=list(mut_data_m2["candidate"]), renamings=renamings,
                                        output=f"outputimages/geds_m2_v_len_location_{location}_filter_{filter_type}_relevance_{relevance_threshold}_{strengths[strength]}_{percentage}.svg")

                        pcc_e1 = stats.pearsonr(mutations_m1, epi_m1)
                        pcc_conf_i_e1 = r_confidence_interval(pcc_e1[0], 0.05, len(mutations_m1))
                        pcc_e2 = stats.pearsonr(mutations_m2, epi_m2)
                        pcc_conf_i_e2 = r_confidence_interval(pcc_e2[0], 0.05, len(mutations_m2))
                        pcc_eA = stats.pearsonr(mutations, epi)
                        pcc_conf_i_eA = r_confidence_interval(pcc_eA[0], 0.05, len(mutations))
                        pcc1 = stats.pearsonr(mutations_m1, geds_m1)
                        pcc_conf_i_1 = r_confidence_interval(pcc1[0], 0.05, len(mutations_m1))
                        pcc2 = stats.pearsonr(mutations_m2, geds_m2)
                        pcc_conf_i_2 = r_confidence_interval(pcc2[0], 0.05, len(mutations_m2))
                        pccA = stats.pearsonr(mutations, geds)
                        pcc_conf_i_A = r_confidence_interval(pccA[0], 0.05, len(mutations))
                        spr1 = stats.spearmanr(mutations_m1, geds_m1)
                        spr2 = stats.spearmanr(mutations_m2, geds_m2)
                        sprA = stats.spearmanr(mutations, geds)

                        percentage_level_data += [spr1[0], spr1[1], spr2[0], spr2[1],
                                                sprA[0], sprA[1],
                                                pcc1[0], pcc1[1],
                                                pcc_conf_i_1[0], pcc_conf_i_1[1],
                                                pcc2[0], pcc2[1],
                                                pcc_conf_i_2[0], pcc_conf_i_2[1],
                                                pccA[0], pccA[1],
                                                pcc_conf_i_A[0], pcc_conf_i_A[1],
                                                pcc_e1[0], pcc_e1[1],
                                                pcc_conf_i_e1[0], pcc_conf_i_e1[1],
                                                pcc_e2[0], pcc_e2[1],
                                                pcc_conf_i_e2[0], pcc_conf_i_e2[1],
                                                pcc_eA[0], pcc_eA[1],
                                                pcc_conf_i_eA[0], pcc_conf_i_eA[1]]
                
                    rows.append(percentage_level_data)

        with open(output_data, "w") as csv_file:
            for row in rows:
                joined = ",".join(list(map(str, row)))
                csv_file.write(f"{joined}\n")