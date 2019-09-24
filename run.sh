# Set environment variables for NetMHCpan before running this script:
# export NETMHCpan=path/to/netmhcpan/binary
# and make sure that netMHCpan is in your PATH.

# Predict putative HLA binders for a set of HLAs and a set of peptide
# sequences given in "NetMHCPanInput/hlas" and "NetMHCPanHLA/", respectively:
mkdir -p NetMHCPanOutput
python csv_to_fasta.py NetMHCPanInput/ NetMHCPanHLA/hlas NetMHCPanOutput/

# Parse allele-frequency net html files:
for file in testfiles/hla_population/*; do
  python parse_allele_frequencies.py "$file"
done

# Compute GELS scores:
mkdir -p outfiles
python count_epitopes.py 0 0 1
python count_epitopes.py 0 1 0
python count_epitopes.py 1 0 0

# Cluster peptides and tumor samples by mutation pattern:
python cluster_genes.py

# Plot GELS against mutation frequencies for each peptide,
# and compute correlations.
mkdir -p outputimages
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
