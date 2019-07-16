# investigate bias
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from unique_count import read_data
from pathlib import Path
from scipy.stats import pearsonr
from excel_processing import read_mutation_files

lengths = {
  "RNF43_C6" : 6,
  "RNF43_G7" : 7,
  "ACVR2A" : 8,
  "GLYR1" : 8,
  "HPS1" : 8,
  "MAPRE3" : 8,
  "MYH11" : 8,
  "OR51E2" : 8,
  "SRRT" : 8,
  "TCF1" : 8,
  "TFE3" : 8,
  "WASF3" : 8,
  "ELAVL3" : 9,
  "EPHB2" : 9,
  "FLT3LG" : 9,
  "LMAN1" : 9,
  "NDUFC2" : 9,
  "RGS12" : 9,
  "RUFY2" : 9,
  "TCF7L2" : 9,
  "TTK" : 9,
  "ABCF1" : 10,
  "AIM2" : 10,
  "C4orf6" : 10,
  "CASP5" : 10,
  "RFC3" : 10,
  "SLC35F5" : 10,
  "SMAP1" : 10,
  "SPINK5" : 10,
  "TFAM" : 10,
  "TGFBR2" : 10,
  "TMEM97" : 10,
  "ASTE1" : 11,
  "CEP164" : 11,
  "MARCKS" : 11,
  "PTHLH" : 11,
  "SLC22A9" : 11,
  "TAF1B" : 11,
  "LTN1" : 11,
  "BANP" : 12,
}

top10 = [
  "TGFBR2", "SLC35F5", "TCF7L2",
  "LTN1", "CASP5", "MARCKS", "SLC22A9",
  "TTK", "MYH11"
]

potential_bias = [
  "TMEM97_m1",
  "FLT3LG_m1",
  "TFE3_m1",
  "SMAP1_m1",
  "MSH3_m1",
  "C1orf34_m1",
  "SPINK5_m1",
  "EPHB2_m1",
  "CEP164_m1",
  "ABCF1_m1",
  "TCF1_m1",
  "HPS1_m1",
  "MAPRE3_m1",
  "OR51E2_m1",
  "PRDM2_m1",
  "RGS12_m1",
  "RUFY2_m1",
  "SRRT_m1",
  "WASF3_m1"
]

mutations = read_mutation_files("testfiles/QMR1.xlsx",
                                "testfiles/QMR2.xlsx")
mutations["RNF43_C6_m1"] = mutations["RNF43.A2_m1"]
mutations["RNF43_G7_m1"] = mutations["RNF43.A3_m1"]
mutations["RNF43_C6_m2"] = mutations["RNF43.A2_m2"]
mutations["RNF43_G7_m2"] = mutations["RNF43.A3_m2"]
                            
data = read_data(Path("table_dump_colon.csv"), delimiter=";")
data = data[data["mm"] == 1]
data = data[data["strength"] == "weak"]
data = data[data["cohort"] == "European Caucasian"]
data = data[data["hla"] == "HLA-A02:01"]
data_names = list(data["candidate"])
data_gels = list(data["GEDS: 50 %"])
data_irs = list(data["IRS: 50 %"])
data_pmut = [
  data_irs[idx] / data_gels[idx]
  if data_gels[idx] != 0 else 0
  for idx in range(len(data_gels))
]
data_dict = {
  data_names[idx] : data_gels[idx]
  for idx in range(len(data_names))
}
data_dict_pmut = mutations

relevant_lengths = [
  lengths[name]
  for name in lengths
]
relevant_gels = [
  data_dict[name + "_m1"]
  if (name + "_m1") in data_dict
  else 0.0
  for name in lengths
]
relevant_pmut = [
  data_dict_pmut[name + "_m1"]
  if (name + "_m1") in data_dict_pmut
  else 0.0
  for name in lengths
]
color = [
  255 if (name + "_m1") in potential_bias else 0
  for name in lengths
]

print(pearsonr(relevant_lengths, relevant_gels))
print(pearsonr(relevant_lengths, relevant_pmut))

fig = plt.figure(figsize=(10, 5))
plt.scatter(relevant_lengths, relevant_gels)
plt.show()
fig = plt.figure(figsize=(10, 5))
plt.scatter(relevant_lengths, relevant_pmut)
plt.show()
