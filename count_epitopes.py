import sys

import unique_count

from math import sqrt
from math import atanh
from math import erf
from matplotlib import pyplot as plt
from scipy.stats import spearmanr
from pathlib import Path
from excel_processing import read_mutation_files
from excel_processing import parse_endometrium_mutations
import pandas as pd
from decimal import *

getcontext().prec = 100

def recursive_table_dump_header(table):  
  result = []
  if isinstance(table, dict):
    for key in table.keys():
      if not isinstance(table[key], dict):
        result += [str(key)]
      result += recursive_table_dump_header(table[key])
  return result
     
def recursive_table_dump_lines(table, lines):
  line = []
  if isinstance(table, dict):
    for key in table.keys():
      line += recursive_table_dump_lines(table[key], lines)
  else:
    line += [str(table)]
    lines.append(line)
  return line

def recursive_table_dump(table, output):
  with open(output, "w") as csv_file:
    lines = []
    header = ";".join(recursive_table_dump_header(table))
    recursive_table_dump_lines(table, lines)
    csv_file.write(header + "\n")
    for line in lines:
      csv_file.write(";".join(line) + "\n")

def count_by_hla_by_mutation(data, hlas, mutations, min_distance=0):
  """
  Reads a pandas dataframe of epitopes, counting distinct epitopes by HLA.
  """
  result = {}
  for mutation in mutations:
    result[mutation] = {}
    dictionaries, _ = unique_count.filter_same(data, min_distance=min_distance,
                                               by_hla=mutation, length_total=200)
    for hla in hlas:
      result[mutation][hla] = 0

    for candidate_dict in dictionaries:
      result[mutation][candidate_dict["HLA"]] += 1
  return result

def m1m2count(data):
  """
  Reads the counts of epitopes for -1 and -2 mutations by HLA in a dataset.
  """
  hlas = list(data["HLA"].sort_values().unique())
  mutations = list(data["ID"].sort_values().unique())
  m_1 = [mut for mut in mutations if mut.endswith("_m1")]
  m_2 = [mut for mut in mutations if mut.endswith("_m2")]
  m_1_counts = count_by_hla_by_mutation(data, hlas, m_1)
  m_2_counts = count_by_hla_by_mutation(data, hlas, m_2)
  return m_1_counts, m_2_counts, hlas, m_1, m_2

def m1m2count_merge(data):
  """
  Reads the total counts of epitopes for -1 and -2 mutations a dataset.
  """
  m_1_counts, m_2_counts, hlas, m_1, m_2 = m1m2count(data)
  m_1_total = {}
  m_2_total = {}
  for mutation in m_1:
    m_1_total[mutation] = 0
  for mutation in m_2:
    m_2_total[mutation] = 0
  for mutation in m_1:
    for hla in hlas:
      m_1_total[mutation] += m_1_counts[mutation][hla]
  for mutation in m_2:
    for hla in hlas:
      m_2_total[mutation] += m_2_counts[mutation][hla]
  return m_1_total, m_2_total, hlas, m_1, m_2

def std_normal_cdf(x):
  return 0.5 + 0.5 * erf(x / sqrt(2))

def fisher_transform(x):
  return atanh(x)

def fisher_Z(x, x0, N):
  return (fisher_transform(x) - fisher_transform(x0)) * sqrt(N - 3)

def pcc_test(pcc_val, N, two_sided=False):
  Z = fisher_Z(pcc_val, 0, N)
  cdf = std_normal_cdf(-abs(Z))
  if two_sided:
    return 2 * cdf
  else:
    return cdf

def pcc(x, y, indexes):
  """Pearson's correlation coefficient."""
  mean_x = sum([x[idx] for idx in indexes if x[idx] > 0.0]) / len(indexes)
  mean_y = sum([y[idx] for idx in indexes if y[idx] > 0.0]) / len(indexes)
  sum_sq_x = sum(map(lambda idx: (x[idx] - mean_x) ** 2 if x[idx] > 0.0 else 0.0, indexes))
  sum_sq_y = sum(map(lambda idx: (y[idx] - mean_y) ** 2 if y[idx] > 0.0 else 0.0, indexes))
  prod_sum = sum(map(lambda idx: (x[idx] - mean_x) * (y[idx] - mean_y) if y[idx] > 0.0 and x[idx] > 0.0 else 0.0, indexes))
  return (prod_sum + 1e-6) / sqrt(sum_sq_x * sum_sq_y + 1e-6)

def m1m2p1p2_correlate(data, probabilities):
  """
  Correlates the mutation probability with the number of epitopes for that
  mutation.
  """
  m_1_total, m_2_total, hlas, m_1, m_2 = m1m2count_merge(data)
  p_m = probabilities
  indexes_m_1 = [hla for key in p_m.keys() if key.endswith("_m1") and key in m_1]
  indexes_m_2 = [key for key in p_m.keys() if key.endswith("_m2") and key in m_2]
  for index in indexes_m_1:
    renamed = index[:-3] + "_m2"
    if renamed in m_2_total:
      continue
    m_2_total[renamed] = 0
  for index in indexes_m_2:
    renamed = index[:-3] + "_m1"
    if renamed in m_1_total:
      continue
    m_1_total[renamed] = 0
  indexes = [key[:-3] for key in indexes_m_1]
  indexes += [key[:-3] for key in indexes_m_2]
  indexes = list(sorted(set(indexes)))
  mq = [m_2_total[idx + "_m2"] / m_1_total[idx + "_m1"] for idx in indexes if p_m[idx + "_m1"] != 0 and m_1_total[idx + "_m1"] != 0]
  pq = [p_m[idx + "_m2"] / p_m[idx + "_m1"] for idx in indexes if p_m[idx + "_m1"] != 0 and m_1_total[idx + "_m1"] != 0]
  pcc_1 = pcc(m_1_total, p_m, indexes_m_1)
  pcc_2 = pcc(m_2_total, p_m, indexes_m_2)
  pcc_q = pcc(mq, pq, list(range(len(pq))))
  pval_1 = pcc_test(pcc_1, len(indexes_m_1), True)
  pval_2 = pcc_test(pcc_2, len(indexes_m_2), True)
  pval_q = pcc_test(pcc_q, len(pq), True)

  return pcc_1, pcc_2, pcc_q, pval_1, pval_2, pval_q

def dump_full(table, output):
  cohorts = [
    "European Caucasian",
    "USA African American",
    "USA Hispanic",
    "Japan",
    "Germany"
  ]
  strengths = ["strong", "weak", "garbage"]
  mms = [1, 2]
  result = []
  for cohort in cohorts:
    object_cohort = table[cohort]
    for strength in strengths:
      object_strength = object_cohort[strength]
      for mm in mms:
        object_mm = object_strength[mm]
        for candidate in object_mm.keys():
          if candidate == "mm":
            continue
          object_candidate = object_mm[candidate]
          for hla in object_candidate.keys():
            if not (hla in ["GEDS", "IRS", "candidate"]):
              object_hla = object_candidate[hla]
              object_eds = object_hla["EDS"]
              object_epds = object_hla["EPDS"]
              data = [cohort, strength, mm, candidate, hla]
              for idx in range(0, 100, 10):
                data.append(object_eds[idx * 0.01])
              for idx in range(0, 100, 10):
                data.append(object_epds[idx * 0.01])
              for idx in range(0, 100, 10):
                data.append(object_candidate["GEDS"][idx * 0.01])
              for idx in range(0, 100, 10):
                data.append(object_candidate["IRS"][idx * 0.01])
              result.append(data)
  header = ["cohort", "strength", "mm", "candidate", "hla"]
  for idx in range(0, 100, 10):
    header.append(f"EDS: {idx} %")
  for idx in range(0, 100, 10):
    header.append(f"EPDS: {idx} %")
  for idx in range(0, 100, 10):
    header.append(f"GEDS: {idx} %")
  for idx in range(0, 100, 10):
    header.append(f"IRS: {idx} %")
  with open(output, "w") as csv_file:
    csv_file.write(";".join(header) + "\n")
    for line in result:
      csv_file.write(";".join(list(map(lambda x: str(x), line))) + "\n")

def doGEDS(p, hlas, freqs, nepitopes):
  prod_A = Decimal(1.0)
  prod_B = Decimal(1.0)
  prob = Decimal(p)
  for raw_hla, freq in zip(list(freqs["Allele"]),
                           list(freqs["Allele Frequency"])):
    hla = f"HLA-{raw_hla}".replace("*", "")
    frequency = Decimal(freq)
    if hla in hlas:
      nh = nepitopes[hla] if hla in nepitopes.keys() else 0
      nh = Decimal(nh)
      pih = Decimal(1) - (Decimal(1) - prob) ** (nh)
      Fh = Decimal(1) - (Decimal(1) - frequency) ** 2
      term = Decimal(1) - Fh * pih
      if hla.startswith("HLA-A"):
        prod_A *= term
      else:
        prod_B *= term
  result = (Decimal(1) - prod_A) + (Decimal(1) - prod_B) - (Decimal(1) - prod_A) * (Decimal(1) - prod_B)
  return float(result)


if __name__ == "__main__":
  assert len(sys.argv) == 4
  DO_SORT = int(sys.argv[1])
  DO_POSTPROC = int(sys.argv[2])
  DO_ANALYSIS = int(sys.argv[3])

  # Eval:
  if DO_SORT:
    data_full = unique_count.read_data(Path("table_dump.csv"), delimiter=";")
    data_full = data_full.sort_values(["GEDS: 50 %"], ascending=False)
    data_wanted = data_full[data_full["strength"] == "weak"]
    data_wanted.to_csv("table_sorted_strong_GEDS.csv")

  # Get some numbers:
  if DO_POSTPROC:
    # Load datasets:
    data_strong = unique_count.read_data(Path("./NetMHCPanOutput/cMNR_peptides.csv.strong.csv"), delimiter=",")
    data_weak = unique_count.read_data(Path("./NetMHCPanOutput/cMNR_peptides.csv.weak.csv"), delimiter=",")
    data_garbage = unique_count.read_data(Path("./NetMHCPanOutput/cMNR_peptides.csv.garbage.csv"), delimiter=",")
    data_weak = pd.concat([data_strong, data_weak])
    data_garbage = pd.concat([data_weak, data_garbage])

    # Mutation data:
    mutation_counts_colon = read_mutation_files("testfiles/QMR1.xlsx",
                                                "testfiles/QMR2.xlsx")
    mutation_counts_endometrium = parse_endometrium_mutations(Path("testfiles/ECQMR.csv"))
    for mutation_name, mutation_counts in zip(
      ["colon", "endometrium"], [mutation_counts_colon, mutation_counts_endometrium]
    ):
      hlas = list(data_garbage["HLA"].sort_values().unique())

      # Epitope counts for all datasets:
      epitope_counts = {}
      epitope_counts["strong"] = m1m2count(data_strong)
      epitope_counts["weak"] = m1m2count(data_weak)
      epitope_counts["garbage"] = m1m2count(data_garbage)

      # MAX-EC & per-HLA binders
      for strength in ["strong", "weak", "garbage"]:
        max_ec = 0
        max_name = ""
        has_hla_epitope = {}
        hla_has_epitope = {}
        for hla in epitope_counts["garbage"][2]:
          hla_has_epitope[hla] = 0
        count = 0
        for mm in [1, 2]:
          mm_index = mm - 1
          ec = epitope_counts[strength][mm_index]
          for key in epitope_counts["garbage"][mm_index].keys():
            total_epitopes = 0
            has_hla_epitope[key] = {}
            count += 1
            for hla in epitope_counts["garbage"][2]:
              if key in ec:
                epitope_count_hla = ec[key][hla] if hla in ec[key] else 0
              else:
                epitope_count_hla = 0
              if epitope_count_hla > 0:
                has_hla_epitope[key][hla] = True
                hla_has_epitope[hla] += 1
              else:
                has_hla_epitope[key][hla] = False
              total_epitopes += epitope_count_hla
            if total_epitopes > max_ec:
              max_ec = total_epitopes
              max_name = key
        print(f"{strength}:")
        print(f"maximum number of epitopes for: {max_name} with: {max_ec}")
        for hla in hla_has_epitope.keys():
          val = (hla_has_epitope[hla] / count) * 100
          print(f"{val} % of candidates have at least one epitope for: {hla}")
        print("-----------------------------------------")
            

  if DO_ANALYSIS:
    # verbosity:
    concise = False

    # Load datasets:
    data_strong = unique_count.read_data(Path("./NetMHCPanOutput/cMNR_peptides.csv.strong.csv"), delimiter=",")
    data_weak = unique_count.read_data(Path("./NetMHCPanOutput/cMNR_peptides.csv.weak.csv"), delimiter=",")
    data_garbage = unique_count.read_data(Path("./NetMHCPanOutput/cMNR_peptides.csv.garbage.csv"), delimiter=",")
    data_weak = pd.concat([data_strong, data_weak])
    data_garbage = pd.concat([data_weak, data_garbage])

    # Mutation data:
    mutation_counts_colon = read_mutation_files("testfiles/QMR1.xlsx",
                                                "testfiles/QMR2.xlsx")
    mutation_counts_endometrium = parse_endometrium_mutations(Path("testfiles/ECQMR.csv"))
    for mutation_name, mutation_counts in zip(
      ["colon", "endometrium"], [mutation_counts_colon, mutation_counts_endometrium]
    ):

      # Hla-frequencies in the population:
      hlas = list(data_garbage["HLA"].sort_values().unique())
      hlas_A = [hla for hla in hlas if hla.startswith("HLA-A")]
      hlas_B = [hla for hla in hlas if hla.startswith("HLA-B")]

      cohort_names = [
        "European Caucasian",
        "USA African American",
        "USA Hispanic",
        "Japan",
        "Germany"
      ]
      hla_freq_paths = list(
        map(
          lambda x: "./testfiles/hla_population/HLA_Population - " + x + ".html.csv",
          cohort_names
        )
      )
      hla_frequencies = {}
      hla_any = {}
      for name, path in zip(cohort_names, hla_freq_paths):
        data_hla = unique_count.read_data(Path(path), delimiter=";")
        hla_frequencies[name] = data_hla

        # Cumulative hla frequencies.
        hla_any[name] = {}
        hla_any[name]["A"] = 0.0
        hla_any[name]["B"] = 0.0
        for raw_hla, frequency in zip(
          list(data_hla["Allele"]), list(data_hla["Allele Frequency"])
        ):
          hla = f"HLA-{raw_hla}".replace("*", "")
          if hla in hlas_A:
            hla_any[name]["A"] += frequency
          if hla in hlas_B:
            hla_any[name]["B"] += frequency

      # Epitope counts for all datasets:
      epitope_counts = {}
      epitope_counts["strong"] = m1m2count(data_strong)
      epitope_counts["weak"] = m1m2count(data_weak)
      epitope_counts["garbage"] = m1m2count(data_garbage)

      # Dump epitope counts per HLA:
      for strength in ["strong", "weak", "garbage"]:
        for mm in [1, 2]:
          mm_index = mm - 1
          ec = epitope_counts[strength][mm_index]
          for key in epitope_counts["garbage"][mm_index].keys():
            with open(f"./outfiles/epitope_counts_m{mm}_{strength}_{key}.csv", "w") as csv_file:
              csv_file.write("HLA;#Epitopes\n")
              for hla in epitope_counts["garbage"][2]:
                if key in ec:
                  epitope_count_hla = ec[key][hla] if hla in ec[key] else 0
                else:
                  epitope_count_hla = 0
                csv_file.write(f"{hla};{epitope_count_hla}\n")

      # Dump scores per HLA:
      scores = {}
      for name in cohort_names:
        score_object = {"cohort": name}
        hla_freq = hla_frequencies[name]
        valid_name = name.replace(" ", "_")
        for strength in ["strong", "weak", "garbage"]:
          score_object_strength = {"strength": strength}
          for mm in [1, 2]:
            mm_index = mm - 1
            ec = epitope_counts[strength][mm_index]
            score_object_mm = {"mm": mm}
            for key in epitope_counts["garbage"][mm_index].keys():
              if not (key in ec):
                continue
              ## Contains all scores for a mutation:
              score_object_candidate = {"candidate": key}
              with open(f"./outfiles/single_hla_scores_{mutation_name}_{valid_name}_m{mm}_{strength}_{key}.csv", "w") as csv_file:
                specs = ";".join([f"Specificity {x} %" for x in range(0, 100, 10)])
                csv_file.write(f"HLA;{specs}\n")
                all_epitopes = 0
                for raw_hla, frequency in zip(list(hla_freq["Allele"]),
                                              list(hla_freq["Allele Frequency"])):
                  hla = f"HLA-{raw_hla}".replace("*", "")
                  score_object_hla = {"hla": hla}
                  if (not (hla in ec[key])):
                    continue
                  if ec[key][hla] == 0.0 and concise:
                    continue
                  line = [hla]
                  line_eds = []
                  score_object_eds = {}
                  line_epds = []
                  score_object_epds = {}
                  for specificity in [x * 0.01 for x in range(0, 100, 10)]:
                    if key in ec:
                      # EDS
                      all_epitopes += ec[key][hla] if hla in ec[key] else 0
                      freq = (1.0 - (1.0 - frequency) ** 2)
                      prob = 1.0 - (1.0 - specificity) ** (ec[key][hla]) if hla in ec[key] else 0.0
                      EDS = freq * prob
                      score_object_eds[specificity] = EDS
                      EPDS = EDS * mutation_counts[key] if key in mutation_counts else float("nan")
                      score_object_epds[specificity] = EPDS
                      line.append(str(EDS))
                    else:
                      line.append("0.0")
                  score_object_hla["EDS"] = score_object_eds
                  score_object_hla["EPDS"] = score_object_epds
                  score_object_candidate[hla] = score_object_hla
                  csv_file.write(";".join(line) + "\n")
                # GEDS
                score_object_geds = {}
                score_object_irs = {}
                for specificity in [x * 0.01 for x in range(0, 100, 10)]:
                  merged_frequency = hla_any[name]["A"] + hla_any[name]["B"] - hla_any[name]["A"] * hla_any[name]["B"]
                  merged_freq_term = 1.0 - (1.0 - merged_frequency) ** 2
                  merged_prob_term = 1.0 - (1.0 - specificity) ** (all_epitopes)
                  GEDS = doGEDS(specificity, hlas, hla_freq, ec[key])
                  score_object_geds[specificity] = GEDS
                  IRS = GEDS * mutation_counts[key] if key in mutation_counts else float("nan")
                  score_object_irs[specificity] = IRS
                score_object_candidate["GEDS"] = score_object_geds
                score_object_candidate["IRS"] = score_object_irs
              score_object_mm[key] = score_object_candidate
            score_object_strength[mm] = score_object_mm
            score_object[strength] = score_object_strength
        scores[name] = score_object

      dump_full(scores, f"table_dump_{mutation_name}.csv")
