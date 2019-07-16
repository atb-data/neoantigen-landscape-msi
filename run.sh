
# Predict putative HLA binders for a set of HLAs and a set of peptide
# sequences given in "NetMHCPanInput/hlas" and "NetMHCPanHLA/", respectively:
mkdir NetMHCPanOutput
python csv_to_fasta.py NetMHCPanInput/hlas NetMHCPanHLA/ NetMHCPanOutput/

# Parse allele-frequency net html files:
for file in `ls testfiles/hla_population`; do
  python parse_allele_frequencies.py $file
done

# Compute GELS scores:
mkdir outfiles
python count_epitopes.py

# Cluster peptides and tumor samples by mutation pattern:
python cluster_genes.py

# Plot GELS against mutation frequencies for each peptide,
# and compute correlations.
mkdir outputimages
python scatter_plots.py

# Generate whisker plot for Fig. 4 D:
python correlation_box.py

# Generate HLA binder plots for main text and supplements:
python unique_count.py

# Generate heatmap for Fig. 1:
python heatmap_huge_csv.py 0
# Generate heatmap for Sup.Fig. 2:
python heatmap_huge_csv.py 1

# Generate scatter plot for Sup.Fig. 5:
python bias_test.py
