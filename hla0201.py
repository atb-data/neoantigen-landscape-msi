import os
import math
from pathlib import Path
import excel_processing
from unique_count import read_data

from scipy import stats
from matplotlib import pyplot as plt

def get_positive(path):
  with open(path) as f:
    result = []
    for line in f:
      name = line.strip()
      result.append(name)
  return result

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

def epitope_counts(m, strength=1):
  result = {}
  strengths = ["strong", "weak", "garbage"]
  for path in os.listdir("./outfiles"):
    for s in strengths[:strength + 1]:
      if path.startswith(f"epitope_counts_m{m}_{s}"):
        name = path.split("_")[4]
        with open(f"./outfiles/{path}") as f:
          counts = {}
          next(f)
          for line in f:
            hla, count = line.strip().split(";")
            count = int(count)
            counts[hla] = count
        if name in result:
          for hla in result[name]:
            result[name][hla] += counts[hla]
        else:
          result[name] = counts
  return result

def hla_percentage(hla, m=1, strength="weak", p=0.5):
  result = {}
  counts = epitope_counts(m, strength=strength)
  for name in counts:
    count = counts[name][hla]
    total = sum(counts[name][h] for h in counts[name])
    result[name] = count / total if total != 0 else 0
  return result

def hla_values(hla, m=1, strength="weak", p=0.5):
  result = {}
  counts = epitope_counts(m, strength=strength)
  for name in counts:
    count = counts[name][hla]
    if count > 0:
      result[name] = 1 - (1 - p) ** counts[name][hla]
    else:
      result[name] = 0
  return result

def positive_names(ec=False):
  name = "hlaA0201" if not ec else "hlaA0201-ec"
  positive = []
  negative = []
  with open(f"testfiles/{name}.dat") as f:
    for line in f:
      positive.append(line.strip())
  with open(f"testfiles/not-{name}.dat") as f:
    for line in f:
      negative.append(line.strip())
  return positive, negative

characteristics = excel_processing.tumor_characteristics("testfiles/CRCQMR.csv")
positive, negative = positive_names()
ec_positive, ec_negative = positive_names(ec=True)

def count_positive(filter):
  return len([tumor for tumor in characteristics if filter(tumor)])

pos_filter = lambda x: x in positive
neg_filter = lambda x: x in negative
ec_pos_filter = lambda x: x in ec_positive
ec_neg_filter = lambda x: x in ec_negative

positive_mutation = excel_processing.read_mutation_files(
  "testfiles/QMR1.xlsx", "testfiles/QMR2.xlsx",
  filter=pos_filter, method=excel_processing.read_mutated_mutations
)
negative_mutation = excel_processing.read_mutation_files(
  "testfiles/QMR1.xlsx", "testfiles/QMR2.xlsx",
  filter=neg_filter, method=excel_processing.read_mutated_mutations
)
positive_tumor_mutation = excel_processing.read_mutation_files(
  "testfiles/QMR1.xlsx", "testfiles/QMR2.xlsx",
  filter=pos_filter, method=excel_processing.read_tumor_mutations
)
negative_tumor_mutation = excel_processing.read_mutation_files(
  "testfiles/QMR1.xlsx", "testfiles/QMR2.xlsx",
  filter=neg_filter, method=excel_processing.read_tumor_mutations
)
positive_tumor_ec = excel_processing.read_mutation_files(
  "testfiles/ECQMR.xlsx",
  filter=ec_pos_filter, method=excel_processing.read_tumor_mutations
)
negative_tumor_ec = excel_processing.read_mutation_files(
  "testfiles/ECQMR.xlsx",
  filter=ec_neg_filter, method=excel_processing.read_tumor_mutations
)

print("N positive samples:", count_positive(pos_filter))
print("N negative samples:", count_positive(neg_filter))

strength = 1
data = read_data(Path("table_dump_colon.csv"), delimiter=";")
data = data[data["mm"] == 1]
data = data[data["strength"] == "weak"]
data = data[data["cohort"] == "European Caucasian"]
data = data[data["hla"] == "HLA-A02:01"]
data_names = list(data["candidate"])
data_gels = list(data["GEDS: 50 %"])
values = hla_values("HLA-A02:01", m=1, strength=strength, p=0.5)
percentage = hla_percentage("HLA-A02:01", m=1, strength=strength, p=0.5)
counts = epitope_counts(m=1, strength=strength)
counts = {
  name : counts[name]["HLA-A02:01"]
  for name in counts
}
candidates = [name[:-3] for name in positive_mutation]
candidates_ec = [name[:-3] for name in positive_tumor_ec]
print(candidates_ec)

valid = [
  "SLC22A9",
  "CASP5",
  "SLC35F5",
  "TGFBR2",
  "TTK",
  "MYH11",
  "TCF7L2",
  "MARCKS"
]

gels = {
    candidate : data_gels[data_names.index(f"{candidate}_m1")]
    for candidate in candidates
    if f"{candidate}_m1" in data_names
}

x = [values[candidate] if candidate in values else 0 for candidate in candidates]
x_p = [counts[candidate] if candidate in counts else 0 for candidate in candidates]
x_pp = [gels[candidate] if candidate in gels else 0 for candidate in candidates]
x_t = [
  values[candidate] if candidate in values else 0
  for candidate in candidates
  for tumor in positive
  if tumor in positive_tumor_mutation[candidate + "_m1"] and positive_tumor_mutation[candidate + "_m1"][tumor] > 0
]
x_t_p = [
  values[candidate] if candidate in values else 0
  for candidate in candidates
  for tumor in negative
  if tumor in negative_tumor_mutation[candidate + "_m1"] and negative_tumor_mutation[candidate + "_m1"][tumor] > 0
]

# EC
x_t_e = [
  values[candidate] if candidate in values else 0
  for candidate in candidates_ec
  for tumor in ec_positive
  if tumor in positive_tumor_ec[candidate + "_m1"] and positive_tumor_ec[candidate + "_m1"][tumor] > 0
]
x_t_e_p = [
  values[candidate] if candidate in values else 0
  for candidate in candidates_ec
  for tumor in ec_negative
  if tumor in negative_tumor_ec[candidate + "_m1"] and negative_tumor_ec[candidate + "_m1"][tumor] > 0
]

y = [positive_mutation[candidate + "_m1"] for candidate in candidates]
y_p = [negative_mutation[candidate + "_m1"] for candidate in candidates]

y_t = [
  positive_tumor_mutation[candidate + "_m1"][tumor]
  for candidate in candidates
  for tumor in positive
  if tumor in positive_tumor_mutation[candidate + "_m1"] and positive_tumor_mutation[candidate + "_m1"][tumor] > 0
]
y_t_p = [
  negative_tumor_mutation[candidate + "_m1"][tumor]
  for candidate in candidates
  for tumor in negative
  if tumor in negative_tumor_mutation[candidate + "_m1"] and negative_tumor_mutation[candidate + "_m1"][tumor] > 0
]

# EC
y_t_e = [
  positive_tumor_ec[candidate + "_m1"][tumor]
  for candidate in candidates_ec
  for tumor in ec_positive
  if tumor in positive_tumor_ec[candidate + "_m1"] and positive_tumor_ec[candidate + "_m1"][tumor] > 0
]
y_t_e_p = [
  negative_tumor_ec[candidate + "_m1"][tumor]
  for candidate in candidates_ec
  for tumor in ec_negative
  if tumor in negative_tumor_ec[candidate + "_m1"] and negative_tumor_ec[candidate + "_m1"][tumor] > 0
]

r_good, p_good = stats.pearsonr(x_t, y_t)
l_good, h_good = r_confidence_interval(r_good, 0.05, len(positive) * len(candidates))
r_bad, p_bad = stats.pearsonr(x_t_p, y_t_p)
l_bad, h_bad = r_confidence_interval(r_bad, 0.05, len(negative) * len(candidates))

y_data = [r_good, r_bad]
y_err = [(r_good - l_good, h_good - r_good), (r_bad - l_bad, h_bad - r_bad)]


fig, ax = plt.subplots(figsize=(5, 5))
ax.axhline(0.0, color="grey", ls="--")
ax.errorbar([1, 2],
            y_data,
            yerr=0,
            fmt="none",
            capsize=20,
            color="black")
ax.errorbar([1, 2],
            y_data,
            yerr=y_err,
            fmt="none",
            capsize=10,
            color="black")

ax.set_xticks([1, 2])
ax.set_xticklabels([
  "HLA-A*02:01",
  "other"
])
ax.set_xlim(0, 3)

plt.tight_layout()
plt.savefig(f"outputimages/hlaA0201-correlation.svg")

print("positive r, C_r, C_r, p:", r_good, l_good, h_good, p_good)
print("negative r, C_r, C_r, p:", r_bad, l_bad, h_bad, p_bad)
print("hla/gels r, p:", stats.pearsonr(x, x_pp))
