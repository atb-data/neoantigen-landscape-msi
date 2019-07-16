import csv
from pathlib import Path
import numpy as np
import matplotlib
from matplotlib import colors as col
from matplotlib import cm
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import plottools

_newmap = col.LinearSegmentedColormap.from_list("Magentas", [
    (1.0, 1.0, 1.0, 1.0),
    (226/255, 0, 116/225, 1.0)])
cm.register_cmap(name="Magentas", cmap=_newmap)

def prefilter(data):
    """Remove artifact candidates."""
    new_data = data[data["ID"] != "CASP5_A8_m1"]
    new_data = new_data[new_data["ID"] != "CASP5_A8_m2"]
    return new_data

def read_data(csv_path, delimiter = ";", filter=prefilter):
    """Read data."""
    with csv_path.open() as csv_file:
        data = pd.read_csv(csv_file, delimiter=delimiter)
    return data

def parse_renamings(renamings, delimiter=";"):
    """Reads all given internal to HUGO renamings from a csv file."""
    data = pd.read_csv(renamings, delimiter=delimiter)
    result = {}
    for row in data.itertuples():
        result[row[1]] = row[2]
    return result

def rename_candidate_hugo(candidate, renamings):
    """Renames a candidate name according to a renaming map."""
    candidate_split = candidate.split("_")
    name_expr = candidate_split[0].split(".")
    base_name = name_expr[0]
    if base_name in renamings:
        base_name = renamings[base_name]
    name_expr[0] = base_name
    name_start = ".".join(name_expr)
    candidate_split[0] = name_start
    result = "_".join(candidate_split)
    return result

def rename_candidate(candidate, renamings=None):
    """Renames a candidate name according to a set of renamings."""
    result = candidate
    if renamings != None:
        result = rename_candidate_hugo(candidate, renamings)
    if result == "RNF43_C6_m1":
        result = "RNF43.A2_m1"
    elif result == "RNF43_C6_m2":
        result = "RNF43.A2_m2"
    elif result == "RNF43_G7_m1":
        result = "RNF43.A3_m1"
    elif result == "RNF43_G7_m2":
        result = "RNF43.A3_m2"
    return result

def filter_same(data, min_distance=0, by_hla=None, length_total=200, renamings=None):
    """Filters a data frame removing duplicated epitopes."""
    candidate_rows = {}
    if by_hla != None:
        data = data[data["ID"] == by_hla]
        data = data.sort_values(["HLA", "Pos"])
        candidates = data["HLA"].unique()
    else:
        data = data.sort_values(["ID", "Pos"])
        candidates = data["ID"].unique()
    for cand in candidates:
        candidate_rows[rename_candidate(cand)] = np.zeros(length_total)
    candidate_dicts = []
    current_candidate = None
    current_position = 0
    current_length = 0
    for row in data.itertuples():
        if by_hla != None:
            candidate = getattr(row, "HLA")
        else:
            candidate = rename_candidate(getattr(row, "ID"))
        position = getattr(row, "Pos")
        length = len(getattr(row, "Peptide"))
        distance = position - current_position
        if current_candidate == None:
            current_candidate = candidate
            current_position = position
            current_length = length
            continue
        elif current_candidate != candidate:
            candidate_dicts.append({
                "ID": current_candidate, "Pos": current_position,
                "length": current_length, "HLA": current_candidate
            })
            current_candidate = candidate
            current_position = position
            current_length = length
            continue
        if distance < min_distance:
            current_length = position - current_position + length
        else:
            candidate_dicts.append({
                "ID": current_candidate, "Pos": current_position,
                "length": current_length, "HLA": current_candidate
            })
            current_candidate = candidate
            current_position = position
            current_length = length
    for candidate_dict in candidate_dicts:
        candidate = candidate_dict["ID"]
        position = candidate_dict["Pos"]
        length = candidate_dict["length"]
        for idx in range(position, position + length):
            candidate_rows[candidate][idx] += 1.0
    return candidate_dicts, candidate_rows

def maximum_overlap(posarray):
    maximum = 0.0
    for key in posarray.keys():
        for elem in posarray[key]:
            if elem > maximum:
                maximum = elem
    return int(maximum)

def read_peptide_lengths(path, renamings=None):
    """Reads a set of peptide lengths from file."""
    result = {}
    with open(path) as csv_file:
        rd = csv.reader(csv_file, delimiter=",")
        for row in rd:
            if not(rename_candidate(row[1], renamings=renamings) in result):
                result[rename_candidate(row[1], renamings=renamings)] = 0
            result[rename_candidate(row[1], renamings=renamings)] += len(row[0])
    return result

if __name__ == "__main__":
    def plot_many(posarray, hlas, candidates, cmap_code="Magentas", output_path="outputimages/out_epitopes.png", width=20, height=140,
                  x_label="foo", y_label="bar", plot_title="baz", cbar_label="quux"):
        length = len(posarray)
        fig, ax = plt.subplots(nrows=length, figsize=(width, height), sharex=True, sharey=True)
        fig.subplots_adjust(hspace=0, wspace=0)
        xticks = [idx * 10 for idx in range(20)] + [200]
        yticks = [idx for idx in range(len(hlas) + 1)]
        ax[0].set_xticks(xticks)
        ax[0].set_yticks(yticks)
        ax[0].set_yticklabels(["all"] + hlas)
        xl = ax[-1].set_xlabel(x_label)
        fig.text(0.06, 0.5, y_label, ha='center', va='center', rotation='vertical')
        xl.set_fontsize(40)
        ax[0].xaxis.labelpad = 20
        ax[0].yaxis.labelpad = 20
        for idx in range(length):
            pt = ax[idx].set_title(candidates[idx])
            pt.set_fontsize(24)
            pt.set_fontweight("bold")
            ax[idx].set_aspect("auto")
            ax[idx].axhline(0.5, color="black")
            cmap = matplotlib.cm.get_cmap(cmap_code)
            cmap.set_bad(color="black")
            allhlas = np.zeros(200)
            for row in hlas:
                allhlas += posarray[idx][row] if row in posarray[idx].keys() else np.zeros(200)
            array = np.asarray([allhlas] + [posarray[idx][row] if row in posarray[idx].keys() else np.zeros(200) for row in hlas])
            for label in ax[idx].get_yticklabels():
                label.set_fontsize(24)
                label.set_fontweight("bold")
            for label in ax[idx].get_xticklabels():
                label.set_rotation(90)
                label.set_horizontalalignment("center")
                label.set_fontsize(24)
                label.set_fontweight("bold")
            ax[idx].patch.set(hatch="xxxxxxx", edgecolor="black")
            im = ax[idx].imshow(array, cmap=cmap, vmax=4, vmin=0, aspect="auto")

        plt.savefig(output_path, dpi=300)
        plt.close()

        ## Plot all those things as separate figures also:
        for idx in range(length):
            fig, ax = plt.subplots(figsize=(20, 20))
            pt = ax.set_title(candidates[idx])
            pt.set_fontsize(24)
            pt.set_fontweight("bold")
            ax.set_aspect("auto")
            ax.axhline(0.5, color="black")
            cmap = matplotlib.cm.get_cmap(cmap_code)
            cmap.set_bad(color="black")
            allhlas = np.zeros(200)
            for row in hlas:
                allhlas += posarray[idx][row] if row in posarray[idx].keys() else np.zeros(200)
            array = np.asarray([allhlas] + [posarray[idx][row] if row in posarray[idx].keys() else np.zeros(200) for row in hlas])
            for label in ax.get_yticklabels():
                label.set_fontsize(24)
                label.set_fontweight("bold")
            for label in ax.get_xticklabels():
                label.set_rotation(90)
                label.set_horizontalalignment("center")
                label.set_fontsize(24)
                label.set_fontweight("bold")
            ax.patch.set(hatch="xxxxxxx", edgecolor="black")
            im = ax.imshow(array, cmap=cmap, vmax=4, vmin=0, aspect="auto")
            plt.savefig(output_path + ".subplot." + candidates[idx] + ".png", dpi=300)
            plt.close()

    def plot_many_complete(posarray_strong, posarray_weak, posarray_garbage, hlas, candidates, lengths, cmap_code="Magentas", output_path="outputimages/out_epitopes.png", width=20, height=140,
                  x_label="foo", y_label="bar", plot_title="baz", cbar_label="quux", renamings=None):
        length = len(posarray_garbage)
        indices = list(range(len(candidates)))
        indices.sort(key=lambda x: -lengths[rename_candidate(candidates[x])])
        candidates.sort(key=lambda x: -lengths[rename_candidate(x)])
        posarray_strong = [posarray_strong[idx] for idx in indices]
        posarray_weak = [posarray_weak[idx] for idx in indices]
        posarray_garbage = [posarray_garbage[idx] for idx in indices]
        max_weak = 14 + 1#max([maximum_overlap(elem) + 1 for elem in posarray_weak])
        max_strong = 7 + 1#max([maximum_overlap(elem) + 1 for elem in posarray_strong])
        max_garbage = 29 + 1#max([maximum_overlap(elem) + 1 for elem in posarray_garbage])
        length_total = max(map(lambda x: lengths[rename_candidate(x)], candidates))
        fig, ax = plt.subplots(nrows=length, ncols=3, figsize=(width, height), sharex=True, sharey=True)
        fig.subplots_adjust(hspace=0.05, wspace=0)
        xticks = [idx for idx in range(0, length_total - 5, 10)]
        yticks = [idx for idx in range(len(hlas) + 1)]
        if length > 1:
            ax[0, 0].set_xticks(xticks)
            ax[0, 0].set_yticks(yticks)
            ax[0, 0].set_yticklabels(["all"] + hlas)
        else:
            ax[0].set_xticks(xticks)
            ax[0].set_yticks(yticks)
            ax[0].set_yticklabels(["all"] + hlas)

        maxval = 4
        cmap_strong = matplotlib.cm.get_cmap(cmap_code)
        cmap_strong = plottools.cmap_discretize(cmap_strong, max_strong)
        cmap_strong.set_bad(color="#ccccccff")
        cmap_weak = plottools.desaturate(matplotlib.cm.get_cmap(cmap_code))
        cmap_garbage = plottools.desaturate(cmap_weak)
        cmap_weak = plottools.cmap_discretize(cmap_weak, max_weak)
        cmap_weak.set_bad(color="#ccccccff")
        cmap_garbage = plottools.cmap_discretize(cmap_garbage, max_garbage)
        cmap_garbage.set_bad(color="#ccccccff")
        for idx in range(length):
            for idy in range(3):
                # if length > 1:
                #     pt = plt.text(0.95, 0.05, rename_candidate(candidates[idx], renamings=renamings)[:-3],
                #                 transform=ax[idx, idy].transAxes,
                #                 horizontalalignment="right",
                #                 verticalalignment="bottom")
                # else:
                #     pt = plt.text(0.95, 0.05, rename_candidate(candidates[idx], renamings=renamings)[:-3],
                #                 transform=ax[idy].transAxes,
                #                 horizontalalignment="right",
                #                 verticalalignment="bottom")
                # pt.set_fontsize(24)
                # pt.set_fontweight("bold")
                if length > 1:
                    ax[idx, idy].set_aspect("auto")
                    ax[idx, idy].axhline(0.5, color="black")
                    if idy == 0:
                        for label in ax[idx, idy].get_yticklabels():
                            label.set_fontsize(24)
                            label.set_fontweight("bold")
                    if idx == length - 1:
                        for label in ax[idx, idy].get_xticklabels():
                            label.set_rotation(90)
                            label.set_horizontalalignment("center")
                            label.set_fontsize(24)
                            label.set_fontweight("bold")
                else:
                    ax[idy].set_aspect("auto")
                    ax[idy].axhline(0.5, color="black")
                    if idy == 0:
                        for label in ax[idy].get_yticklabels():
                            label.set_fontsize(24)
                            label.set_fontweight("bold")
                    if idx == length - 1:
                        for label in ax[idy].get_xticklabels():
                            label.set_rotation(90)
                            label.set_horizontalalignment("center")
                            label.set_fontsize(24)
                            label.set_fontweight("bold")
            allhlas_strong = np.zeros(length_total)
            allhlas_weak = np.zeros(length_total)
            allhlas_garbage = np.zeros(length_total)
            for row in hlas:
                if row in posarray_strong[idx].keys():
                    posarray_strong[idx][row][lengths[rename_candidate(candidates[idx])]:] = float("nan")
                else:
                    posarray_strong[idx][row] = np.zeros(length_total)
                    posarray_strong[idx][row][lengths[rename_candidate(candidates[idx])]:] = float("nan")
                if row in posarray_weak[idx].keys():
                    posarray_weak[idx][row][lengths[rename_candidate(candidates[idx])]:] = float("nan")
                else:
                    posarray_weak[idx][row] = np.zeros(length_total)
                    posarray_weak[idx][row][lengths[rename_candidate(candidates[idx])]:] = float("nan")
                if row in posarray_garbage[idx].keys():
                    posarray_garbage[idx][row][lengths[rename_candidate(candidates[idx])]:] = float("nan")
                else:
                    posarray_garbage[idx][row] = np.zeros(length_total)
                    posarray_garbage[idx][row][lengths[rename_candidate(candidates[idx])]:] = float("nan")
                allhlas_strong += posarray_strong[idx][row] if row in posarray_strong[idx].keys() else np.zeros(length_total)
                allhlas_weak += posarray_weak[idx][row] if row in posarray_weak[idx].keys() else np.zeros(length_total)
                allhlas_garbage += posarray_garbage[idx][row] if row in posarray_garbage[idx].keys() else np.zeros(length_total)
            array_strong = np.asarray([allhlas_strong] + [posarray_strong[idx][row] if row in posarray_strong[idx].keys() else np.zeros(length_total) for row in hlas])
            array_weak = np.asarray([allhlas_weak] + [posarray_weak[idx][row] if row in posarray_weak[idx].keys() else np.zeros(length_total) for row in hlas])
            array_garbage = np.asarray([allhlas_garbage] + [posarray_garbage[idx][row] if row in posarray_garbage[idx].keys() else np.zeros(length_total) for row in hlas])

            if length > 1:
                im_strong = ax[idx, 0].imshow(array_strong, cmap=cmap_strong, vmax=max_strong, vmin=0, aspect="auto")
                im_weak = ax[idx, 1].imshow(array_weak, cmap=cmap_weak, vmax=max_weak, vmin=0, aspect="auto")
                im_garbage = ax[idx, 2].imshow(array_garbage, cmap=cmap_garbage, vmax=max_garbage, vmin=0, aspect="auto")
            else:
                im_strong = ax[0].imshow(array_strong, cmap=cmap_strong, vmax=max_strong, vmin=0, aspect="auto")
                im_weak = ax[1].imshow(array_weak, cmap=cmap_weak, vmax=max_weak, vmin=0, aspect="auto")
                im_garbage = ax[2].imshow(array_garbage, cmap=cmap_garbage, vmax=max_garbage, vmin=0, aspect="auto")
            
        plt.savefig(output_path, dpi=300)
        plt.close()

    def plot_complete_h(strong_posarray_m1, weak_posarray_m1, garbage_posarray_m1,
                        strong_posarray_m2, weak_posarray_m2, garbage_posarray_m2,
                        mutation_counts, length_total = 200,
                        width=40, height=40, output_path="outputimages/plot_complete.png",
                        x_label="foo", y_label="bar", renamings=None):
        max_weak = max(maximum_overlap(weak_posarray_m1), maximum_overlap(weak_posarray_m2)) + 1
        max_strong = max(maximum_overlap(strong_posarray_m1), maximum_overlap(strong_posarray_m2)) + 1
        max_garbage = max(maximum_overlap(garbage_posarray_m1), maximum_overlap(garbage_posarray_m2)) + 1
        
        cmap_m1 = matplotlib.cm.get_cmap("Magentas")
        cmap_m1 = plottools.cmap_discretize(cmap_m1, max_strong)
        cmap_m1.set_bad(color="#ccccccff")
        cmap_weak_m1 = plottools.desaturate(matplotlib.cm.get_cmap("Magentas"))
        cmap_garbage_m1 = plottools.desaturate(cmap_weak_m1)
        cmap_weak_m1 = plottools.cmap_discretize(cmap_weak_m1, max_weak)
        cmap_weak_m1.set_bad(color="#ccccccff")
        cmap_garbage_m1 = plottools.cmap_discretize(cmap_garbage_m1, max_garbage)
        cmap_garbage_m1.set_bad(color="#ccccccff")
        cmap_m2 = matplotlib.cm.get_cmap("Greens")
        cmap_m2 = plottools.cmap_discretize(cmap_m2, max_strong)
        cmap_m2.set_bad(color="#ccccccff")
        cmap_weak_m2 = plottools.desaturate(matplotlib.cm.get_cmap("Greens"))
        cmap_garbage_m2 = plottools.desaturate(cmap_weak_m2)
        cmap_weak_m2 = plottools.cmap_discretize(cmap_weak_m2, max_weak)
        cmap_weak_m2.set_bad(color="#ccccccff")
        cmap_garbage_m2 = plottools.cmap_discretize(cmap_garbage_m2, max_garbage)
        cmap_garbage_m2.set_bad(color="#ccccccff")

        candidates_garbage_m1 = [key for key in mutation_counts.keys() if key.endswith("_m1")]
        candidates_garbage_m1.sort(key=lambda x: rename_candidate(x, renamings=renamings))
        candidates_garbage_m2 = [key for key in mutation_counts.keys() if key.endswith("_m2")]
        candidates_garbage_m2.sort(key=lambda x: rename_candidate(x, renamings=renamings))

        array_strong_m1 = np.asarray([strong_posarray_m1[row] if row in strong_posarray_m1.keys() else np.zeros(length_total) for row in candidates_garbage_m1])
        array_weak_m1 = np.asarray([weak_posarray_m1[row] if row in weak_posarray_m1.keys() else np.zeros(length_total) for row in candidates_garbage_m1])
        array_garbage_m1 = np.asarray([garbage_posarray_m1[row] if row in garbage_posarray_m1.keys() else np.zeros(length_total) for row in candidates_garbage_m1])
        array_strong_m2 = np.asarray([strong_posarray_m2[row] if row in strong_posarray_m2.keys() else np.zeros(length_total) for row in candidates_garbage_m2])
        array_weak_m2 = np.asarray([weak_posarray_m2[row] if row in weak_posarray_m2.keys() else np.zeros(length_total) for row in candidates_garbage_m2])
        array_garbage_m2 = np.asarray([garbage_posarray_m2[row] if row in garbage_posarray_m2.keys() else np.zeros(length_total) for row in candidates_garbage_m2])

        xticks = [idx * 10 for idx in range(21)]
        yticks = [idx for idx in range(len(candidates_garbage_m1))]

        fig, ax = plt.subplots(nrows=2, ncols=8,
                               gridspec_kw = {
                                   'height_ratios': [1, 80],
                                   'width_ratios': [3, 3, 3, 1,
                                                    3, 3, 3, 1]
                               }, sharex=False, sharey=False,
                               figsize=(width, height))
        fig.subplots_adjust(hspace=0, wspace=0)

        # axis sharing and deletion
        ax[1, 1].get_shared_y_axes().join(ax[1, 1], ax[1, 2])
        ax[1, 2].get_shared_y_axes().join(ax[1, 2], ax[1, 3])
        ax[1, 5].get_shared_y_axes().join(ax[1, 5], ax[1, 6])
        ax[1, 6].get_shared_y_axes().join(ax[1, 6], ax[1, 7])
        ax[0, 3].set_axis_off()
        ax[0, 7].set_axis_off()
        single_ax = ax[1, 0]
        single_ax.set_aspect("auto")
        single_ax.set_yticks(yticks)
        for elem in candidates_garbage_m1:
            print(rename_candidate(elem, renamings=renamings), renamings, elem)
        single_ax.set_yticklabels(map(lambda x: x[:-3],
            map(lambda x: rename_candidate(x, renamings=renamings), candidates_garbage_m1)))
        for label in single_ax.get_yticklabels():
            label.set_fontsize(24)
            label.set_fontweight("bold")
        for ax_idy in (list(range(3)) + list(range(4, 7))):
            single_ax = ax[1, ax_idy]
            single_ax.grid(True, which='major', axis='x', linestyle='--')
            single_ax.set_aspect("auto")
            if ax_idy > 0:
                single_ax.set_yticks(yticks)
                single_ax.set_yticklabels([])
            single_ax.set_adjustable("box-forced")
            single_ax.set_xticks(xticks)
            for label in single_ax.get_xticklabels():
                label.set_rotation(90)
                label.set_horizontalalignment("center")
                label.set_fontsize(24)
                label.set_fontweight("bold")

        # Mutation counts:
        ax[1, 3].barh(list(range(len(candidates_garbage_m1))), [200 * mutation_counts[key] if key in mutation_counts.keys() else 300.0 for key in candidates_garbage_m1], height=1.0, color="black")
        ax[1, 7].barh(list(range(len(candidates_garbage_m2))), [200 * mutation_counts[key] if key in mutation_counts.keys() else 300.0 for key in candidates_garbage_m2], height=1.0, color="black")
        ax[1, 3].set_aspect("auto")
        ax[1, 3].set_xticks([20 * x for x in range(11)])
        ax[1, 3].set_xticklabels([f"{x * 0.1:.1f}" for x in range(10)])
        ax[1, 3].set_yticks(yticks)
        ax[1, 3].set_yticklabels([])
        for label in ax[1, 3].get_xticklabels():
            label.set_rotation(90)
            label.set_horizontalalignment("center")
            label.set_fontsize(16)
            label.set_fontweight("bold")

        ax[1, 7].set_aspect("auto")
        ax[1, 7].set_xticks([20 * x for x in range(11)])
        ax[1, 7].set_xticklabels([f"{x * 0.1:.1f}" for x in range(11)])
        ax[1, 7].set_yticks(yticks)
        ax[1, 7].set_yticklabels([])
        for label in ax[1, 7].get_xticklabels():
            label.set_rotation(90)
            label.set_horizontalalignment("center")
            label.set_fontsize(16)
            label.set_fontweight("bold")

        im_strong_m1 = ax[1, 0].imshow(array_strong_m1, cmap=cmap_m1, vmin=0.0, vmax=max_strong, aspect="auto")
        im_weak_m1 = ax[1, 1].imshow(array_weak_m1, cmap=cmap_weak_m1, vmin=0.0, vmax=max_weak, aspect="auto")
        im_garbage_m1 = ax[1, 2].imshow(array_garbage_m1, cmap=cmap_garbage_m1, vmin=0.0, vmax=max_garbage, aspect="auto")
        im_strong_m2 = ax[1, 4].imshow(array_strong_m2, cmap=cmap_m2, vmin=0.0, vmax=max_strong, aspect="auto")
        im_weak_m2 = ax[1, 5].imshow(array_weak_m2, cmap=cmap_weak_m2, vmin=0.0, vmax=max_weak, aspect="auto")
        im_garbage_m2 = ax[1, 6].imshow(array_garbage_m2, cmap=cmap_garbage_m2, vmin=0.0, vmax=max_garbage, aspect="auto")

        cbar_strong_m1 = fig.colorbar(im_strong_m1, cax=ax[0, 0], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        cbar_weak_m1 = fig.colorbar(im_weak_m1, cax=ax[0, 1], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        cbar_garbage_m1 = fig.colorbar(im_garbage_m1, cax=ax[0, 2], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        cbar_strong_m2 = fig.colorbar(im_strong_m2, cax=ax[0, 4], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        cbar_weak_m2 = fig.colorbar(im_weak_m2, cax=ax[0, 5], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        cbar_garbage_m2 = fig.colorbar(im_garbage_m2, cax=ax[0, 6], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])

        cbar_strong_m1.set_ticks([float(idx) for idx in range(max_strong)])
        cbar_weak_m1.set_ticks([float(idx) for idx in range(0, max_weak, 2)])
        cbar_garbage_m1.set_ticks([float(idx) for idx in range(0, max_garbage, 5)])
        cbar_strong_m1.ax.set_xticklabels([idx for idx in range(max_strong)])
        cbar_weak_m1.ax.set_xticklabels([idx for idx in range(0, max_weak, 2)])
        cbar_garbage_m1.ax.set_xticklabels([idx for idx in range(0, max_garbage, 5)])
        cbar_strong_m1.ax.xaxis.set_ticks_position("top")
        cbar_weak_m1.ax.xaxis.set_ticks_position("top")
        cbar_garbage_m1.ax.xaxis.set_ticks_position("top")
        cbar_strong_m1.ax.xaxis.set_label_position("top")
        cbar_weak_m1.ax.xaxis.set_label_position("top")
        cbar_garbage_m1.ax.xaxis.set_label_position("top")

        cbar_strong_m2.set_ticks([float(idx) for idx in range(max_strong)])
        cbar_weak_m2.set_ticks([float(idx) for idx in range(0, max_weak, 2)])
        cbar_garbage_m2.set_ticks([float(idx) for idx in range(0, max_garbage, 5)])
        cbar_strong_m2.ax.set_xticklabels([idx for idx in range(max_strong)])
        cbar_weak_m2.ax.set_xticklabels([idx for idx in range(0, max_weak, 2)])
        cbar_garbage_m2.ax.set_xticklabels([idx for idx in range(0, max_garbage, 5)])
        cbar_strong_m2.ax.xaxis.set_ticks_position("top")
        cbar_weak_m2.ax.xaxis.set_ticks_position("top")
        cbar_garbage_m2.ax.xaxis.set_ticks_position("top")
        cbar_strong_m2.ax.xaxis.set_label_position("top")
        cbar_weak_m2.ax.xaxis.set_label_position("top")
        cbar_garbage_m2.ax.xaxis.set_label_position("top")

        for label in cbar_weak_m1.ax.get_xticklabels():
            label.set_fontsize(16)
            label.set_fontweight("bold")
        for label in cbar_strong_m1.ax.get_xticklabels():
            label.set_fontsize(16)
            label.set_fontweight("bold")
        for label in cbar_garbage_m1.ax.get_xticklabels():
            label.set_fontsize(16)
            label.set_fontweight("bold")

        for label in cbar_weak_m2.ax.get_xticklabels():
            label.set_fontsize(16)
            label.set_fontweight("bold")
        for label in cbar_strong_m2.ax.get_xticklabels():
            label.set_fontsize(16)
            label.set_fontweight("bold")
        for label in cbar_garbage_m2.ax.get_xticklabels():
            label.set_fontsize(16)
            label.set_fontweight("bold")

        plt.savefig(output_path, dpi=300)
        plt.close("all")

    def plot_complete_strong_h(strong_posarray_m1, strong_posarray_m2,
                        mutation_counts, length_total = 200,
                        width=40, height=40, output_path="outputimages/plot_complete.png",
                        x_label="foo", y_label="bar", renamings=None):
        max_strong = max(maximum_overlap(strong_posarray_m1), maximum_overlap(strong_posarray_m2)) + 1
        
        cmap_m1 = matplotlib.cm.get_cmap("Magentas")
        cmap_m1 = plottools.cmap_discretize(cmap_m1, max_strong)
        cmap_m1.set_bad(color="#ccccccff")
        cmap_m2 = matplotlib.cm.get_cmap("Greens")
        cmap_m2 = plottools.cmap_discretize(cmap_m2, max_strong)
        cmap_m2.set_bad(color="#ccccccff")

        candidates_garbage_m1 = [key for key in mutation_counts.keys() if key.endswith("_m1")]
        candidates_garbage_m1.sort(key=lambda x: rename_candidate(x, renamings=renamings))
        candidates_garbage_m2 = [key for key in mutation_counts.keys() if key.endswith("_m2")]
        candidates_garbage_m2.sort(key=lambda x: rename_candidate(x, renamings=renamings))
        
        array_strong_m1 = np.asarray([strong_posarray_m1[row][:length_total] if row in strong_posarray_m1.keys() else np.zeros(length_total) for row in candidates_garbage_m1])
        array_strong_m2 = np.asarray([strong_posarray_m2[row][:length_total] if row in strong_posarray_m2.keys() else np.zeros(length_total) for row in candidates_garbage_m2])

        xticks = [idx * 10 for idx in range(21)]
        yticks = [idx for idx in range(len(candidates_garbage_m1))]

        fig, ax = plt.subplots(nrows=2, ncols=4,
                               gridspec_kw = {
                                   'height_ratios': [1, 80],
                                   'width_ratios': [3, 1,
                                                    3, 1]
                               }, sharex=False, sharey=False,
                               figsize=(width, height))
        fig.subplots_adjust(hspace=0, wspace=0)

        # axis sharing and deletion
        ax[1, 1].get_shared_y_axes().join(ax[1, 1], ax[1, 2])
        ax[1, 2].get_shared_y_axes().join(ax[1, 2], ax[1, 3])
        ax[0, 1].set_axis_off()
        ax[0, 3].set_axis_off()
        single_ax = ax[1, 0]
        single_ax.set_aspect("auto")
        single_ax.set_yticks(yticks)
        for elem in candidates_garbage_m1:
            print(rename_candidate(elem, renamings=renamings), renamings, elem)
        single_ax.set_yticklabels(map(lambda x: x[:-3],
            map(lambda x: rename_candidate(x, renamings=renamings), candidates_garbage_m1)))
        for label in single_ax.get_yticklabels():
            label.set_fontsize(24)
            label.set_fontweight("bold")
        for ax_idy in [0, 2]:
            single_ax = ax[1, ax_idy]
            single_ax.grid(True, which='major', axis='x', linestyle='--')
            single_ax.set_aspect("auto")
            if ax_idy > 0:
                single_ax.set_yticks(yticks)
                single_ax.set_yticklabels([])
            single_ax.set_adjustable("box-forced")
            single_ax.set_xticks(xticks)
            for label in single_ax.get_xticklabels():
                label.set_rotation(90)
                label.set_horizontalalignment("center")
                label.set_fontsize(24)
                label.set_fontweight("bold")

        # Mutation counts:
        ax[1, 1].barh(list(range(len(candidates_garbage_m1))), [200 * mutation_counts[key] if key in mutation_counts.keys() else 300.0 for key in candidates_garbage_m1], height=1.0, color="black")
        ax[1, 3].barh(list(range(len(candidates_garbage_m2))), [200 * mutation_counts[key] if key in mutation_counts.keys() else 300.0 for key in candidates_garbage_m2], height=1.0, color="black")
        ax[1, 1].set_xticks([20 * x for x in range(11)])
        ax[1, 1].set_xticklabels([f"{x * 0.1:.1f}" for x in range(10)])
        ax[1, 1].set_yticks(yticks)
        ax[1, 1].set_yticklabels([])
        for label in ax[1, 1].get_xticklabels():
            label.set_rotation(90)
            label.set_horizontalalignment("center")
            label.set_fontsize(16)
            label.set_fontweight("bold")

        ax[1, 3].set_aspect("auto")
        ax[1, 3].set_xticks([20 * x for x in range(11)])
        ax[1, 3].set_xticklabels([f"{x * 0.1:.1f}" for x in range(10)])
        ax[1, 3].set_yticks(yticks)
        ax[1, 3].set_yticklabels([])
        for label in ax[1, 3].get_xticklabels():
            label.set_rotation(90)
            label.set_horizontalalignment("center")
            label.set_fontsize(16)
            label.set_fontweight("bold")

        im_strong_m1 = ax[1, 0].imshow(array_strong_m1, cmap=cmap_m1, vmin=0.0, vmax=max_strong, aspect="auto")
        im_strong_m2 = ax[1, 2].imshow(array_strong_m2, cmap=cmap_m2, vmin=0.0, vmax=max_strong, aspect="auto")
        
        cbar_strong_m1 = fig.colorbar(im_strong_m1, cax=ax[0, 0], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        cbar_strong_m2 = fig.colorbar(im_strong_m2, cax=ax[0, 2], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        
        cbar_strong_m1.set_ticks([float(idx) for idx in range(max_strong)])
        cbar_strong_m1.ax.set_xticklabels([idx for idx in range(max_strong)])
        cbar_strong_m1.ax.xaxis.set_ticks_position("top")
        cbar_strong_m1.ax.xaxis.set_label_position("top")
        
        cbar_strong_m2.set_ticks([float(idx) for idx in range(max_strong)])
        cbar_strong_m2.ax.set_xticklabels([idx for idx in range(max_strong)])
        cbar_strong_m2.ax.xaxis.set_ticks_position("top")
        cbar_strong_m2.ax.xaxis.set_label_position("top")


        for label in cbar_strong_m1.ax.get_xticklabels():
            label.set_fontsize(16)
            label.set_fontweight("bold")
        for label in cbar_strong_m2.ax.get_xticklabels():
            label.set_fontsize(16)
            label.set_fontweight("bold")
        
        plt.savefig(output_path, dpi=300)
        plt.close("all")


    def plot_complete(strong_posarray_m1, weak_posarray_m1, garbage_posarray_m1,
                      strong_posarray_m2, weak_posarray_m2, garbage_posarray_m2,
                      mutation_counts, length_total = 200,
                      width=40, height=40, output_path="outputimages/plot_complete.png",
                      x_label="foo", y_label="bar"):
        max_weak = max(maximum_overlap(weak_posarray_m1), maximum_overlap(weak_posarray_m2)) + 1
        max_strong = max(maximum_overlap(strong_posarray_m1), maximum_overlap(strong_posarray_m2)) + 1
        max_garbage = max(maximum_overlap(garbage_posarray_m1), maximum_overlap(garbage_posarray_m2)) + 1
        
        cmap_m1 = matplotlib.cm.get_cmap("Magentas")
        cmap_m1 = plottools.cmap_discretize(cmap_m1, max_strong)
        cmap_m1.set_bad(color="#ccccccff")
        cmap_weak_m1 = plottools.desaturate(matplotlib.cm.get_cmap("Magentas"))
        cmap_garbage_m1 = plottools.desaturate(cmap_weak_m1)
        cmap_weak_m1 = plottools.cmap_discretize(cmap_weak_m1, max_weak)
        cmap_weak_m1.set_bad(color="#ccccccff")
        cmap_garbage_m1 = plottools.cmap_discretize(cmap_garbage_m1, max_garbage)
        cmap_garbage_m1.set_bad(color="#ccccccff")
        cmap_m2 = matplotlib.cm.get_cmap("Greens")
        cmap_m2 = plottools.cmap_discretize(cmap_m2, max_strong)
        cmap_m2.set_bad(color="#ccccccff")
        cmap_weak_m2 = plottools.desaturate(matplotlib.cm.get_cmap("Greens"))
        cmap_garbage_m2 = plottools.desaturate(cmap_weak_m2)
        cmap_weak_m2 = plottools.cmap_discretize(cmap_weak_m2, max_weak)
        cmap_weak_m2.set_bad(color="#ccccccff")
        cmap_garbage_m2 = plottools.cmap_discretize(cmap_garbage_m2, max_garbage)
        cmap_garbage_m2.set_bad(color="#ccccccff")

        candidates_garbage_m1 = [key for key in mutation_counts.keys() if key.endswith("_m1")]
        candidates_garbage_m1.sort(key=lambda x: -sum([1.0 for elem in garbage_posarray_m1[x] if elem > -1.0]) if x in garbage_posarray_m1.keys() else 0.0)
        candidates_garbage_m2 = [key for key in mutation_counts.keys() if key.endswith("_m2")]
        candidates_garbage_m2.sort(key=lambda x: -sum([1.0 for elem in garbage_posarray_m1[x[:-3] + "_m1"] if elem > -1.0]) if x[:-3] + "_m1" in garbage_posarray_m1.keys() else 0.0)

        array_strong_m1 = np.asarray([strong_posarray_m1[row] if row in strong_posarray_m1.keys() else np.zeros(length_total) for row in candidates_garbage_m1])
        array_weak_m1 = np.asarray([weak_posarray_m1[row] if row in weak_posarray_m1.keys() else np.zeros(length_total) for row in candidates_garbage_m1])
        array_garbage_m1 = np.asarray([garbage_posarray_m1[row] if row in garbage_posarray_m1.keys() else np.zeros(length_total) for row in candidates_garbage_m1])
        array_strong_m2 = np.asarray([strong_posarray_m2[row] if row in strong_posarray_m2.keys() else np.zeros(length_total) for row in candidates_garbage_m2])
        array_weak_m2 = np.asarray([weak_posarray_m2[row] if row in weak_posarray_m2.keys() else np.zeros(length_total) for row in candidates_garbage_m2])
        array_garbage_m2 = np.asarray([garbage_posarray_m2[row] if row in garbage_posarray_m2.keys() else np.zeros(length_total) for row in candidates_garbage_m2])

        xticks = [idx * 10 for idx in range(21)]
        yticks = [idx for idx in range(len(candidates_garbage_m1))]

        fig, ax = plt.subplots(nrows=4, ncols=4, gridspec_kw = {'height_ratios': [1, 1, 80, 80], 'width_ratios': [3, 3, 3, 1]}, sharex=False, sharey=False, figsize=(width, height))
        fig.subplots_adjust(hspace=0, wspace=0)

        # axis sharing and deletion
        ax[2, 1].get_shared_y_axes().join(ax[2, 1], ax[2, 2])
        ax[2, 2].get_shared_y_axes().join(ax[2, 2], ax[2, 3])
        ax[3, 1].get_shared_y_axes().join(ax[3, 1], ax[3, 2])
        ax[3, 2].get_shared_y_axes().join(ax[3, 2], ax[3, 3])
        ax[0, 3].set_axis_off()
        ax[1, 3].set_axis_off()
        for ax_idx in range(2, 4):
            single_ax = ax[ax_idx, 0]
            single_ax.set_aspect("auto")
            single_ax.set_yticks(yticks)
            single_ax.set_yticklabels(map(lambda x: x[:-3], candidates_garbage_m1))
            for label in single_ax.get_yticklabels():
                label.set_fontsize(24)
                label.set_fontweight("bold")
            for ax_idy in range(3):
                single_ax = ax[ax_idx, ax_idy]
                single_ax.set_aspect("auto")
                if ax_idy > 0:
                    single_ax.set_yticks(yticks)
                    single_ax.set_yticklabels([])
                single_ax.set_adjustable("box-forced")
                if ax_idx == 3:
                    single_ax.set_xticks(xticks)
                    for label in single_ax.get_xticklabels():
                        label.set_rotation(90)
                        label.set_horizontalalignment("center")
                        label.set_fontsize(24)
                        label.set_fontweight("bold")
                else:
                    plt.setp(single_ax.get_xticklabels(), visible=False)

            # Mutation counts:
            if ax_idx == 2:
                ax[ax_idx, 3].barh(list(range(len(candidates_garbage_m1))), [200 * mutation_counts[key] if key in mutation_counts.keys() else 300.0 for key in candidates_garbage_m1], height=1.0, color="black")
            else:
                ax[ax_idx, 3].barh(list(range(len(candidates_garbage_m2))), [200 * mutation_counts[key] if key in mutation_counts.keys() else 300.0 for key in candidates_garbage_m2], height=1.0, color="black")
            ax[ax_idx, 3].set_aspect("auto")
            ax[ax_idx, 3].set_xticks([20 * x for x in range(11)])
            ax[ax_idx, 3].set_xticklabels([f"{x * 0.1:.1f}" for x in range(11)])
            ax[ax_idx, 3].set_yticks(yticks)
            ax[ax_idx, 3].set_yticklabels([])
            for label in ax[ax_idx, 3].get_xticklabels():
                label.set_rotation(90)
                label.set_horizontalalignment("center")
                label.set_fontsize(24)
                label.set_fontweight("bold")

        im_strong_m1 = ax[2, 0].imshow(array_strong_m1, cmap=cmap_m1, vmin=0.0, vmax=max_strong, aspect="auto")
        im_weak_m1 = ax[2, 1].imshow(array_weak_m1, cmap=cmap_weak_m1, vmin=0.0, vmax=max_weak, aspect="auto")
        im_garbage_m1 = ax[2, 2].imshow(array_garbage_m1, cmap=cmap_garbage_m1, vmin=0.0, vmax=max_garbage, aspect="auto")
        im_strong_m2 = ax[3, 0].imshow(array_strong_m2, cmap=cmap_m2, vmin=0.0, vmax=max_strong, aspect="auto")
        im_weak_m2 = ax[3, 1].imshow(array_weak_m2, cmap=cmap_weak_m2, vmin=0.0, vmax=max_weak, aspect="auto")
        im_garbage_m2 = ax[3, 2].imshow(array_garbage_m2, cmap=cmap_garbage_m2, vmin=0.0, vmax=max_garbage, aspect="auto")

        cbar_strong_m1 = fig.colorbar(im_strong_m1, cax=ax[0, 0], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        cbar_weak_m1 = fig.colorbar(im_weak_m1, cax=ax[0, 1], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        cbar_garbage_m1 = fig.colorbar(im_garbage_m1, cax=ax[0, 2], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        cbar_strong_m2 = fig.colorbar(im_strong_m2, cax=ax[1, 0], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        cbar_weak_m2 = fig.colorbar(im_weak_m2, cax=ax[1, 1], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])
        cbar_garbage_m2 = fig.colorbar(im_garbage_m2, cax=ax[1, 2], orientation = "horizontal", ticks=[0, 1, 2, 3, 4])

        # colorbars
        cbar_strong_m2.ax.set_xticklabels([])
        cbar_weak_m2.ax.set_xticklabels([])
        cbar_garbage_m2.ax.set_xticklabels([])

        cbar_strong_m1.set_ticks([float(idx) for idx in range(max_strong)])
        cbar_weak_m1.set_ticks([float(idx) for idx in range(0, max_weak, 2)])
        cbar_garbage_m1.set_ticks([float(idx) for idx in range(0, max_garbage, 5)])
        cbar_strong_m1.ax.set_xticklabels([idx for idx in range(max_strong)])
        cbar_weak_m1.ax.set_xticklabels([idx for idx in range(0, max_weak, 2)])
        cbar_garbage_m1.ax.set_xticklabels([idx for idx in range(0, max_garbage, 5)])
        cbar_strong_m1.ax.xaxis.set_ticks_position("top")
        cbar_weak_m1.ax.xaxis.set_ticks_position("top")
        cbar_garbage_m1.ax.xaxis.set_ticks_position("top")
        cbar_strong_m1.ax.xaxis.set_label_position("top")
        cbar_weak_m1.ax.xaxis.set_label_position("top")
        cbar_garbage_m1.ax.xaxis.set_label_position("top")

        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes("right", size = "3%", pad = 0.2)
        # divider.set_aspect(False)
        # cax.set_aspect("auto")
        # cbar = fig.colorbar(im, cax = cax, orientation = "vertical")
        # cbar.ax.tick_params(labelsize=32)
        # cbar.set_label(cbar_label, labelpad=20, fontsize=40)
        # plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close("all")

    def nan_post_length(posarray, lengths):
        result = {}
        for name in lengths.keys():
            if name in posarray.keys():
                row = posarray[name]
            else:
                row = np.zeros(200)
            row[lengths[name]:] = float("nan")
            result[name] = row
        return result

    def plot_v2():
        from excel_processing import read_mutation_files

        renamings = parse_renamings("Renamings.csv", delimiter=",")
        data_strong = read_data(Path("./NetMHCPanOutput/cMNR_peptides.csv.strong.csv"), delimiter=",")
        data_weak = read_data(Path("./NetMHCPanOutput/cMNR_peptides.csv.weak.csv"), delimiter=",")
        data_weak = pd.concat([data_strong, data_weak])
        data_garbage = read_data(Path("./NetMHCPanOutput/cMNR_peptides.csv.garbage.csv"), delimiter=",")
        data_garbage = pd.concat([data_weak, data_garbage])
        peptide_lengths = read_peptide_lengths("./NetMHCPanInput/cMNR_peptides.csv")
        candidates = list(data_garbage[list(map(lambda x: x.endswith("m1"), data_garbage["ID"]))]["ID"].unique())

        single_candidates = candidates + [cand[:-3] + "_m2" for cand in candidates]
        posarrays_strong = []
        posarrays_weak = []
        posarrays_garbage = []
        hlas = list(data_garbage["HLA"].sort_values().unique())
        length_total = max(map(lambda x: peptide_lengths[rename_candidate(x)], single_candidates))
        for candidate in single_candidates:
            try:
                _, posarray_strong = filter_same(data_strong, min_distance=0, by_hla=candidate, length_total=length_total)
                _, posarray_weak = filter_same(data_weak, min_distance=0, by_hla=candidate, length_total=length_total)
                _, posarray_garbage = filter_same(data_garbage, min_distance=0, by_hla=candidate, length_total=length_total)
                _, single_posarray_strong = filter_same(data_strong, min_distance=0, by_hla=candidate, length_total=peptide_lengths[rename_candidate(candidate)])
                _, single_posarray_weak = filter_same(data_weak, min_distance=0, by_hla=candidate, length_total=peptide_lengths[rename_candidate(candidate)])
                _, single_posarray_garbage = filter_same(data_garbage, min_distance=0, by_hla=candidate, length_total=peptide_lengths[rename_candidate(candidate)])
                posarrays_strong.append(posarray_strong)
                posarrays_weak.append(posarray_weak)
                posarrays_garbage.append(posarray_garbage)
                cmap_color = "Magentas" if candidate.endswith("m1") else "Greens"
                plot_many_complete([single_posarray_strong], [single_posarray_weak], [single_posarray_garbage], hlas, [candidate], peptide_lengths, output_path=f"outputimages/multiple_hlas.{candidate}.svg", width=20, height=10,
                    x_label="amino acid position [1]", cbar_label="# overlapping epitopes [1]",
                    cmap_code=cmap_color,
                    plot_title="predicted binders by HLA-type, candidate peptide and position",
                    renamings=renamings)
            except Exception as e:
                print("BROKEN CANDIDATE:", candidate, "with exception: ", repr(e))
        
        fig, ax = plt.subplots(nrows=6, ncols=1)
        fig.subplots_adjust(hspace=1.5, wspace=0)

        def max_of_all(ps):
            max_single = max([maximum_overlap(posarray) for posarray in ps])
            max_sum = max([max(sum([posarray[key] for key in posarray])) for posarray in ps if len(posarray.keys()) > 0])
            return max_single, max_sum

        max_m1_strong = max_of_all(posarrays_strong)[0]
        max_m1_weak = max_of_all(posarrays_weak)[0]
        max_m1_garbage = max_of_all(posarrays_garbage)[0]
        max_m2_strong = max_of_all(posarrays_strong)[0]
        max_m2_weak = max_of_all(posarrays_weak)[0]
        max_m2_garbage = max_of_all(posarrays_garbage)[0]

        dummy_data_m1_strong = np.array([list(range(max_m1_strong))])
        dummy_data_m1_weak = np.array([list(range(max_m1_weak))])
        dummy_data_m1_garbage = np.array([list(range(max_m1_garbage))])
        dummy_data_m2_strong = np.array([list(range(max_m2_strong))])
        dummy_data_m2_weak = np.array([list(range(max_m2_weak))])
        dummy_data_m2_garbage = np.array([list(range(max_m2_garbage))])

        cmap_m1_strong = matplotlib.cm.get_cmap("Magentas")
        cmap_m1_strong = plottools.cmap_discretize(cmap_m1_strong, max_m1_strong)
        cmap_m1_strong.set_bad(color="#ccccccff")
        cmap_m1_weak = plottools.desaturate(matplotlib.cm.get_cmap("Magentas"))
        cmap_m1_garbage = plottools.desaturate(cmap_m1_weak)
        cmap_m1_weak = plottools.cmap_discretize(cmap_m1_weak, max_m1_weak)
        cmap_m1_weak.set_bad(color="#ccccccff")
        cmap_m1_garbage = plottools.cmap_discretize(cmap_m1_garbage, max_m1_garbage)
        cmap_m1_garbage.set_bad(color="#ccccccff")

        cmap_m2_strong = matplotlib.cm.get_cmap("Greens")
        cmap_m2_strong = plottools.cmap_discretize(cmap_m2_strong, max_m2_strong)
        cmap_m2_strong.set_bad(color="#ccccccff")
        cmap_m2_weak = plottools.desaturate(matplotlib.cm.get_cmap("Greens"))
        cmap_m2_garbage = plottools.desaturate(cmap_m2_weak)
        cmap_m2_weak = plottools.cmap_discretize(cmap_m2_weak, max_m2_weak)
        cmap_m2_weak.set_bad(color="#ccccccff")
        cmap_m2_garbage = plottools.cmap_discretize(cmap_m2_garbage, max_m2_garbage)
        cmap_m2_garbage.set_bad(color="#ccccccff")

        for idx in range(6):
            ax[idx].set_yticks([])

        ax[0].set_xticks([idx for idx in range(0, max_m1_strong, max(max_m1_strong // 10, 1))])
        ax[0].set_xticklabels([idx for idx in range(0, max_m1_strong, max(max_m1_strong // 10, 1))])
        ax[1].set_xticks([idx for idx in range(0, max_m1_weak, max(max_m1_weak // 10, 1))])
        ax[1].set_xticklabels([idx for idx in range(0, max_m1_weak, max(max_m1_weak // 10, 1))])
        ax[2].set_xticks([idx for idx in range(0, max_m1_garbage, max(max_m1_garbage // 10, 1))])
        ax[2].set_xticklabels([idx for idx in range(0, max_m1_garbage, max(max_m1_garbage // 10, 1))])
        ax[3].set_xticks([idx for idx in range(0, max_m2_strong, max(max_m2_strong // 10, 1))])
        ax[3].set_xticklabels([idx for idx in range(0, max_m2_strong, max(max_m2_strong // 10, 1))])
        ax[4].set_xticks([idx for idx in range(0, max_m2_weak, max(max_m2_weak // 10, 1))])
        ax[4].set_xticklabels([idx for idx in range(0, max_m2_weak, max(max_m2_weak // 10, 1))])
        ax[5].set_xticks([idx for idx in range(0, max_m2_garbage, max(max_m2_garbage // 10, 1))])
        ax[5].set_xticklabels([idx for idx in range(0, max_m2_garbage, max(max_m2_garbage // 10, 1))])

        ax[0].imshow(dummy_data_m1_strong, cmap=cmap_m1_strong, aspect="auto")
        ax[1].imshow(dummy_data_m1_weak, cmap=cmap_m1_weak, aspect="auto")
        ax[2].imshow(dummy_data_m1_garbage, cmap=cmap_m1_garbage, aspect="auto")
        ax[3].imshow(dummy_data_m2_strong, cmap=cmap_m2_strong, aspect="auto")
        ax[4].imshow(dummy_data_m2_weak, cmap=cmap_m2_weak, aspect="auto")
        ax[5].imshow(dummy_data_m2_garbage, cmap=cmap_m2_garbage, aspect="auto")

        plt.show()

        plot_many_complete(posarrays_strong, posarrays_weak, posarrays_garbage, hlas, single_candidates, peptide_lengths, output_path="outputimages/multiple_hlas.svg", width=40, height=80,
                  x_label="amino acid position [1]", cbar_label="# overlapping epitopes [1]",
                  plot_title="predicted binders by HLA-type, candidate peptide and position")

        posarrays = []
        strong_posarray_m1 = nan_post_length(filter_same(data_strong[list(map(lambda x: x.endswith("m1"), data_strong["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)
        weak_posarray_m1 = nan_post_length(filter_same(data_weak[list(map(lambda x: x.endswith("m1"), data_weak["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)
        garbage_posarray_m1 = nan_post_length(filter_same(data_garbage[list(map(lambda x: x.endswith("m1"), data_garbage["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)
        strong_posarray_m2 = nan_post_length(filter_same(data_strong[list(map(lambda x: x.endswith("m2"), data_strong["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)
        weak_posarray_m2 = nan_post_length(filter_same(data_weak[list(map(lambda x: x.endswith("m2"), data_weak["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)
        garbage_posarray_m2 = nan_post_length(filter_same(data_garbage[list(map(lambda x: x.endswith("m2"), data_garbage["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)

        top10 = ["TGFBR2", "ZNF294", "SLC22A9", "SLC35F5", "CASP5", "TTK", "TCF7L2", "MARCKS", "MYH11", "BANP"]
        def is_top10_m1(x):
            return x.endswith("m1") and (x[:-3] in top10)
        def is_top10_m2(x):
            return x.endswith("m2") and (x[:-3] in top10)
        strong_posarray_top10_m1 = nan_post_length(filter_same(data_strong[list(map(lambda x: is_top10_m1(x), data_strong["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)
        weak_posarray_top10_m1 = nan_post_length(filter_same(data_weak[list(map(lambda x: is_top10_m1(x), data_weak["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)
        garbage_posarray_top10_m1 = nan_post_length(filter_same(data_garbage[list(map(lambda x: is_top10_m1(x), data_garbage["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)
        strong_posarray_top10_m2 = nan_post_length(filter_same(data_strong[list(map(lambda x: is_top10_m2(x), data_strong["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)
        weak_posarray_top10_m2 = nan_post_length(filter_same(data_weak[list(map(lambda x: is_top10_m2(x), data_weak["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)
        garbage_posarray_top10_m2 = nan_post_length(filter_same(data_garbage[list(map(lambda x: is_top10_m2(x), data_garbage["ID"]))], min_distance=0, renamings=renamings)[1], peptide_lengths)


        mutation_counts = read_mutation_files("testfiles/QMR1.xlsx",
                                              "testfiles/QMR2.xlsx")
        mutation_counts_top10 = {
            (candidate + mness): mutation_counts[candidate + mness]
            for candidate in top10
            for mness in ["_m1", "_m2"]
        }

        # candidates.sort(key=lambda x: -sum([elem for elem in posarray[x]]))
        plot_complete_h(strong_posarray_m1, weak_posarray_m1, garbage_posarray_m1,
                        strong_posarray_m2, weak_posarray_m2, garbage_posarray_m2,
                        mutation_counts, width=60, height=20, renamings=renamings,
                        output_path="outputimages/minimum_distance_0_complete_result_horizontal.svg")
        
        plot_complete_strong_h(strong_posarray_top10_m1, strong_posarray_top10_m2,
                        mutation_counts_top10, length_total=100, width=40, height=5, renamings=renamings,
                        output_path="outputimages/minimum_distance_0_strong_top10_result_horizontal.svg")
        plot_complete_strong_h(weak_posarray_top10_m1, weak_posarray_top10_m2,
                        mutation_counts_top10, length_total=100, width=40, height=5, renamings=renamings,
                        output_path="outputimages/minimum_distance_0_weak_top10_result_horizontal.svg")
        plot_complete_strong_h(garbage_posarray_top10_m1, garbage_posarray_top10_m2,
                        mutation_counts_top10, length_total=100, width=40, height=5, renamings=renamings,
                        output_path="outputimages/minimum_distance_0_garbage_top10_result_horizontal.svg")

        plot_complete_strong_h(strong_posarray_m1, strong_posarray_m2,
                        mutation_counts, width=40, height=20, renamings=renamings,
                        output_path="outputimages/minimum_distance_0_strong_result_horizontal.svg")

        # Plot by HLA versions:
        for hla in hlas:
            filtered_data_strong = data_strong[data_strong["HLA"] == hla]
            filtered_data_weak = data_weak[data_weak["HLA"] == hla]
            filtered_data_garbage = data_garbage[data_garbage["HLA"] == hla]
            filtered_strong_posarray_m1 = nan_post_length(filter_same(filtered_data_strong[list(map(lambda x: x.endswith("m1"), filtered_data_strong["ID"]))], min_distance=0)[1], peptide_lengths)
            filtered_weak_posarray_m1 = nan_post_length(filter_same(filtered_data_weak[list(map(lambda x: x.endswith("m1"), filtered_data_weak["ID"]))], min_distance=0)[1], peptide_lengths)
            filtered_garbage_posarray_m1 = nan_post_length(filter_same(filtered_data_garbage[list(map(lambda x: x.endswith("m1"), filtered_data_garbage["ID"]))], min_distance=0)[1], peptide_lengths)
            filtered_strong_posarray_m2 = nan_post_length(filter_same(filtered_data_strong[list(map(lambda x: x.endswith("m2"), filtered_data_strong["ID"]))], min_distance=0)[1], peptide_lengths)
            filtered_weak_posarray_m2 = nan_post_length(filter_same(filtered_data_weak[list(map(lambda x: x.endswith("m2"), filtered_data_weak["ID"]))], min_distance=0)[1], peptide_lengths)
            filtered_garbage_posarray_m2 = nan_post_length(filter_same(filtered_data_garbage[list(map(lambda x: x.endswith("m2"), filtered_data_garbage["ID"]))], min_distance=0)[1], peptide_lengths)
            plot_complete_h(filtered_strong_posarray_m1, filtered_weak_posarray_m1, filtered_garbage_posarray_m1,
                            filtered_strong_posarray_m2, filtered_weak_posarray_m2, filtered_garbage_posarray_m2,
                            mutation_counts, width=60, height=20, renamings=renamings,
                            output_path=f"outputimages/minimum_distance_0_complete_result_{hla}_horizontal.svg")

        data = data_weak
        hlas = list(data["HLA"].sort_values().unique())
        candidates = ["TAF1B_m1", "TAF1B_m2", "ASTE1_m1", "ASTE1_m2", "AIM2_m1", "AIM2_m2", "TGFBR2_m1", "TGFBR2_m2", "RNF43_m1", "RNF43_m2"]
        for candidate in candidates:
            posarray = filter_same(data, min_distance=0, by_hla=candidate, length_total=200)
            posarrays.append(posarray)
        plot_many(posarrays, hlas, candidates[0:20], output_path="outputimages/multiple_hlas.png", width=40, height=80,
                  x_label="amino acid position [1]", cbar_label="# overlapping epitopes [1]",
                  plot_title="predicted binders by HLA-type, candidate peptide and position")

    plot_v2()
