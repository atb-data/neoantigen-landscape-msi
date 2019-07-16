import csv
from pathlib import Path
import unique_count
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import plottools
from mpl_toolkits.axes_grid1 import make_axes_locatable

def isNA(string):
  """Filter #N/A from incoming tables."""
  if string == "#N/A" or string == "NA":
    return True
  else:
    return False

def tryFloat(string):
  """Try to parse a string into a floating point number."""
  if isNA(string):
    return float("nan")
  else:
    return float(string)

def parse_colon_complete(csv_path, nancount_samples=20):
  """Parses the colon raw data, eliminating tumor samples with more than nancount_samples nans."""
  data = {}
  candidates = []
  with csv_path.open() as csv_file:
    header = next(csv_file)
    for line in csv_file:
      fields = line.split(";")
      tumor_id = fields[0]
      candidate = fields[2]
      m1 = float(fields[7]) + float(fields[10])
      m2 = float(fields[6]) + float(fields[9])
      wt = float(fields[8]) + float(fields[5]) + float(fields[11])
      if not (tumor_id in data):
        data[tumor_id] = {}
      data[tumor_id][candidate] = {
        "wt": wt,
        "m1": m1,
        "m2": m2
      }
      if not (candidate in candidates):
        candidates.append(candidate)
  cleaned_data = {}
  for key in data.keys():
    cleaned_data[key] = []
    for candidate in candidates:
      if candidate in data[key]:
        cleaned_data[key].append(1.0 - data[key][candidate]["wt"])
      else:
        cleaned_data[key].append(float("nan"))
    cleaned_data[key] = np.array(cleaned_data[key])
  denaned_data = {}
  for key in cleaned_data.keys():
    nancount = 0
    for idx, candidate in enumerate(candidates):
      if not (cleaned_data[key][idx] >= 0.0):
        nancount += 1
    if nancount < nancount_samples:
      denaned_data[key] = cleaned_data[key]
  return (denaned_data, candidates)

def parse_endometrium_complete(csv_path):
  """Parses the endometrium raw data."""
  data = {}
  candidates = ["MAPRE3"]
  with csv_path.open() as csv_file:
    header = next(csv_file)
    for line in csv_file:
      fields = line.split(";")
      tumor_id = fields[0]
      candidate = fields[1]
      m2 = float(fields[32])
      m1 = float(fields[33])
      wt = float(fields[34])
      if not (tumor_id in data):
        data[tumor_id] = {}
      data[tumor_id][candidate] = {
        "wt": wt,
        "m1": m1,
        "m2": m2
      }
      if not (candidate in candidates):
        candidates.append(candidate)
  cleaned_data = {}
  for key in data.keys():
    cleaned_data[key] = []
    for candidate in candidates:
      if candidate in data[key]:
        cleaned_data[key].append(1.0 - data[key][candidate]["wt"])
      elif candidate == "MAPRE3":
        cleaned_data[key].append(0.0)
      else:
        cleaned_data[key].append(float("nan"))
    cleaned_data[key] = np.array(cleaned_data[key])
  return (cleaned_data, candidates)

def parse_csv(csv_path, delimiter=",", force=None, nancount_rows=10, nancount_cols=20):
  """Parses a given raw-data csv file."""
  if force is None:
    force = ["TFAM", "CLOCK"]
  data = {}
  with csv_path.open() as csv_file:
    rd = csv.reader(csv_file, delimiter=delimiter)
    header = next(rd)
    rowlength = 0
    for line in rd:
      sample = line[0]
      vals = []
      currentnan = 0
      for idx in range(21, len(line)):
        if isNA(line[idx]):
          currentnan += 1
        floatopt = tryFloat(line[idx])
        val = 0.0
        if floatopt < 0.15:
          val = 0.0
        elif floatopt > 0.85:
          val = 1.0
        else:
          val = floatopt
        vals.append(val)
      if currentnan < nancount_rows:
        data[sample] = vals
        rowlength = len(data[sample])
    column_data = {}
    nancounts = np.zeros(rowlength)
    for sample in data.keys():
      for idx in range(len(data[sample])):
        if np.isnan(data[sample][idx]):
          nancounts[idx] += 1
    above_threshold = [idx for idx, count in enumerate(nancounts) if count > nancount_cols]
    candidates = [elem for idx, elem in enumerate(header[21:]) if (not idx in above_threshold) or elem in force]
    for sample in data.keys():
      row = []
      for idx in range(len(data[sample])):
        if (not idx in above_threshold) or header[21 + idx] in force:
          row.append(data[sample][idx])
      column_data[sample] = np.array(row)
  return (column_data, candidates)

def anonymize(tumor_list, prefix="CS"):
  """Anonymizes a list of tumor names."""
  return {tumor : tumor for tumor in tumor_list}

def rename_to_hugo(candidate, renamings):
  """Renames a candidate to correct HUGO nomenclature."""
  result = unique_count.rename_candidate(candidate, renamings=renamings)
  if result == "A2":
    result = "RNF43.C6"
  if result == "A3":
    result = "RNF43.G7"
  if result == "RNF43.A2":
    result = "RNF43.C6"
  if result == "RNF43.A3":
    result = "RNF43.G7"
  return result

def plot_data(data, endodata, cmap_code="seismic", x_label="candidate peptides",
              y_label="tumor samples",
              plot_title="mutation frequencies by peptide and tumor sample",
              cbar_label="percentage mutated [%]",
              vertical=False, mutsort=True, output="outputimages/out_heatmap.svg"):
  """Plots a given set of mutation data as a heatmap."""
  renamings = unique_count.parse_renamings("Renamings.csv", delimiter=",")
  block = []
  endoblock = []
  dictionary = data[0]
  endodictionary = endodata[0]
  candidates = data[1]
  endocandidates = endodata[1]
  keylist = list(dictionary.keys())
  endokeylist = list(endodictionary.keys())
  anonymized = anonymize(keylist)
  endoanonymized = anonymize(endokeylist, prefix="ES")
  endoyticks = [idx for idx in range(len(endokeylist))]

  if mutsort:
    keylist.sort(key=lambda x: -sum([elem for elem in dictionary[x] if elem > 0.0]))
    endokeylist.sort(key=lambda x: -sum([elem for elem in endodictionary[x] if elem > 0.0]))
    sumx = [0.0] * len(candidates)
    for sample in keylist:
      for idx in range(len(candidates)):
        if dictionary[sample][idx] > 0.0:
            sumx[idx] += dictionary[sample][idx]
    endosumx = [0.0] * len(endocandidates)
    for sample in endokeylist:
      for idx in range(len(endocandidates)):
        if endodictionary[sample][idx] > 0.0:
            endosumx[idx] += endodictionary[sample][idx]
    candidates_sorted = [elem[0] for elem in sorted(zip(candidates, sumx), key=lambda x: -x[1])]
  else:
    keylist.sort(key=lambda x: sum([1 for elem in dictionary[x] if not elem > -1.0]))
    endokeylist.sort(key=lambda x: sum([1 for elem in endodictionary[x] if not elem > -1.0]))
    sumx = [0] * len(candidates)
    for sample in keylist:
      for idx in range(len(candidates)):
        if not dictionary[sample][idx] > -1.0:
          sumx[idx] += 1
    candidates_sorted = sorted(candidates, key=lambda x: rename_to_hugo(x, renamings))

  merge_candidates = candidates_sorted

  for sample in keylist:
    sortedsample = []
    for cand in merge_candidates:
      sampleval = dictionary[sample][candidates.index(cand)]
      sortedsample.append(sampleval)
    block.append(sortedsample)
  for sample in endokeylist:
    sortedsample = []
    for cand in merge_candidates:
      if cand in endocandidates:
        sampleval = endodictionary[sample][endocandidates.index(cand)]
      else:
        sampleval = float("nan")
      sortedsample.append(sampleval)
    endoblock.append(sortedsample)

  xticks = [idx for idx in range(len(merge_candidates))]
  yticks = [idx for idx in range(len(keylist))]

  # Plot block:
  if vertical:
    array = np.array(block)
    endoarray = np.array(endoblock)
    cmap = matplotlib.cm.get_cmap(cmap_code)
    cmap = plottools.truncate_colormap(cmap, minval=0.0, maxval=0.8, n=100)
    cmap.set_bad(color="black")
    fig, ax = plt.subplots(nrows=3, ncols=1, gridspec_kw={"height_ratios": [1, len(keylist), len(endokeylist)], "width_ratios": [1]}, sharex=False, figsize=(35, 80))
    plt.subplots_adjust(hspace=0.01)
    ax[1].xaxis.labelpad = 20
    ax[1].yaxis.labelpad = 20
    ax[1].set_xticks(xticks)
    ax[1].set_yticks(yticks)
    ax[1].set_xticklabels([])
    ax[1].set_yticklabels([anonymized[key] for key in keylist])
    ax[2].set_xticks(xticks)
    ax[2].set_yticks(endoyticks)
    ax[2].set_xticklabels(candidates)
    ax[2].set_yticklabels([endoanonymized[key] for key in endokeylist])
    for idx in range(1, 3):
      for label in ax[idx].get_yticklabels():
        label.set_fontsize(40)
        label.set_fontweight("bold")
      for label in ax[idx].get_xticklabels():
        label.set_rotation(90)
        label.set_horizontalalignment("center")
        label.set_fontsize(40)
        label.set_fontweight("bold")
      ax[idx].patch.set(hatch="xxxxxxx", edgecolor="black")
    im = ax[1].imshow(array, cmap=cmap, vmax=1.0, vmin=0.0, aspect="auto")
    endoim = ax[2].imshow(endoarray, cmap=cmap, vmax=1.0, vmin=0.0, aspect="auto")
    cbar = fig.colorbar(im, cax=ax[0], orientation = "horizontal")
    cbar.ax.tick_params(labelsize=32)
    cbar.set_ticks([idx * 0.2 for idx in range(6)])
    cbar.set_ticklabels([idx * 20 for idx in range(6)])
    for label in cbar.ax.get_yticklabels():
      label.set_fontsize(40)
    cbar.set_label(cbar_label, labelpad=20, fontsize=40)
    cbar.ax.xaxis.set_ticks_position("top")

    for axis in ax:
      axis.set_anchor("W")

    plt.tight_layout()

    plt.savefig(output)
    plt.close()
  else:
    array = np.array(block).T
    endoarray = np.array(endoblock).T
    cmap = matplotlib.cm.get_cmap(cmap_code)
    cmap = plottools.truncate_colormap(cmap, minval=0.0, maxval=0.8, n=100)
    cmap.set_bad(color="black")
    fig, ax = plt.subplots(nrows=1, ncols=3, gridspec_kw={"height_ratios": [1], "width_ratios": [len(keylist), len(endokeylist), 1]}, sharex=False, figsize=(80, 35))
    plt.subplots_adjust(wspace=0.001)
    ax[0].xaxis.labelpad = 20
    ax[0].yaxis.labelpad = 20
    ax[0].set_yticks(xticks)
    ax[0].set_xticks(yticks)
    ax[0].set_yticklabels(list(map(lambda x: rename_to_hugo(x, renamings), merge_candidates)))
    ax[0].set_xticklabels([anonymized[key] for key in keylist])
    ax[1].set_yticks(xticks)
    ax[1].set_xticks(endoyticks)
    ax[1].set_yticklabels([])
    ax[1].set_xticklabels([endoanonymized[key] for key in endokeylist])
    for idx in range(0, 2):
      for label in ax[idx].get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(40)
        label.set_horizontalalignment("center")
        label.set_fontweight("bold")
      for label in ax[idx].get_yticklabels():
        label.set_fontsize(40)
        label.set_fontweight("bold")
      ax[idx].patch.set(hatch="xxxxxxx", edgecolor="black")
    im = ax[0].imshow(array, cmap=cmap, vmax=1.0, vmin=0.0, aspect="auto")
    endoim = ax[1].imshow(endoarray, cmap=cmap, vmax=1.0, vmin=0.0, aspect="auto")
    cbar = fig.colorbar(im, cax=ax[2], orientation = "vertical")
    cbar.ax.tick_params(labelsize=32)
    cbar.set_ticks([idx * 0.2 for idx in range(6)])
    cbar.set_ticklabels([idx * 20 for idx in range(6)])
    for label in cbar.ax.get_yticklabels():
      label.set_fontsize(40)
    cbar.set_label(cbar_label, labelpad=20, fontsize=40)

    for axis in ax:
      axis.set_anchor("W")

    plt.tight_layout()

    plt.savefig(output)
    plt.close()

if __name__ == "__main__":
  import sys

  data = parse_csv(Path("testfiles/CRCQMR.csv"),
                   delimiter=";", nancount_rows=18,
                   nancount_cols=20)
  endometrium = parse_endometrium_complete(Path("testfiles/ECQMR.csv"))
  for idx, key in enumerate(endometrium[1]):
    if key.endswith(".C6"):
      endometrium[1][idx] = "A2"
    if key.endswith(".G7"):
      endometrium[1][idx] = "A3"

  if int(sys.argv[1]):
    data_full = parse_csv(Path("testfiles/CRCQMR.csv"), nancount_rows=1000, nancount_cols=1000, delimiter=";")
    plot_data(data_full, endometrium, cmap_code="Blues", output="outputimages/SupFig2.svg", mutsort=False)
  else:
    plot_data(data, endometrium, cmap_code="Blues", output="outputimages/Fig1.svg", mutsort=True)

