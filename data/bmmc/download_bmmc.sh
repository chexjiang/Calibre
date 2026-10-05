# Luecken BMMC (bone marrow mononuclear cells)
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122

mkdir -p GSE194122
cd GSE194122

# Download the data
wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz
gunzip GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz