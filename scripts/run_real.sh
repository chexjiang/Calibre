# for celltype in Bcell Myeloid T; do
#   for method in epilink_poisson epilink_nb scmultimap signac; do
#     Rscript scripts/run_real.R \
#       Dogma_CTRL "$celltype" "$method" \
#       > "results/real/Dogma_CTRL/${celltype}_${method}.log" 2>&1
#   done
# done
# Simulation DA
Rscript --vanilla scripts/run_da_simulation.R --workers=16 
Rscript --vanilla scripts/plot_da_simulation.R
  
# BMMC: (Bcell Dendritic Erythroid Myeloid T)
Rscript scripts/run_real.R bmmc Erythroid epilink_nb
Rscript scripts/run_real.R bmmc Erythroid epilink_poisson
Rscript scripts/run_real.R bmmc Myeloid signac

# PBMC: (Bcell Myeloid T)
Rscript scripts/run_real.R pbmc T epilink_nb
Rscript scripts/run_real.R pbmc T epilink_poisson
Rscript scripts/run_real.R pbmc Myeloid scmultimap 
Rscript scripts/run_real.R pbmc Bcell signac

# Dogma
Rscript scripts/run_real.R Dogma_STIM Bcell epilink_nb
Rscript scripts/run_real.R Dogma_STIM Bcell epilink_poisson

Rscript scripts/run_real.R Dogma_CTRL Bcell epilink_nb
Rscript scripts/run_real.R Dogma_CTRL Bcell epilink_poisson

Rscript scripts/run_real.R Dogma_STIM T scmultimap
Rscript scripts/run_real.R Dogma_CTRL T scmultimap

Rscript scripts/run_real.R Dogma_STIM T signac
Rscript scripts/run_real.R Dogma_CTRL T signac


Rscript scripts/run_real.R NEAT CD4T epilink_nb
Rscript scripts/run_real.R NEAT CD4T epilink_poisson
Rscript scripts/run_real.R NEAT CD4T scmultimap
Rscript scripts/run_real.R NEAT CD4T signac

Rscript scripts/run_real.R SHARE LCL epilink_nb
Rscript scripts/run_real.R SHARE LCL epilink_poisson
Rscript scripts/run_real.R SHARE LCL scmultimap
Rscript scripts/run_real.R SHARE LCL signac

# Brain ("Excitatory neurons", Astrocytes, Oligodendroglia)
Rscript scripts/run_real.R Brain "Excitatory neurons" epilink_nb
Rscript scripts/run_real.R Brain "Excitatory neurons" epilink_poisson
Rscript scripts/run_real.R Brain "Excitatory neurons" scmultimap
Rscript scripts/run_real.R Brain "Excitatory neurons" signac

# SatpathyK562 K562
Rscript scripts/run_real.R SatpathyK562 K562 epilink_nb
Rscript scripts/run_real.R SatpathyK562 K562 epilink_poisson
Rscript scripts/run_real.R SatpathyK562 K562 scmultimap
Rscript scripts/run_real.R SatpathyK562 K562 signac

# Arthritis (T, B, Myeloid, Endothelial, Fibroblast)
Rscript scripts/run_real.R Arthritis Endothelial epilink_nb
Rscript scripts/run_real.R Arthritis Endothelial epilink_poisson

Rscript scripts/run_real.R Arthritis T scmultimap
Rscript scripts/run_real.R Arthritis B scmultimap
Rscript scripts/run_real.R Arthritis Endothelial scmultimap
Rscript scripts/run_real.R Arthritis Fibroblast scmultimap
Rscript scripts/run_real.R Arthritis Myeloid scmultimap

Rscript scripts/run_real.R Arthritis Endothelial signac

# Brain2 from scMultimap (Microglia, Astrocytes, "Inhibitory neurons", "Excitatory neurons")
Rscript scripts/run_real.R Brain2 "Inhibitory neurons" epilink_poisson \
  --condition-effect --overwrite
Rscript scripts/run_real.R Brain2 "Inhibitory neurons" epilink_poisson \
  --condition AD --overwrite
Rscript scripts/run_real.R Brain2 Microglia epilink_poisson \
  --condition CTRL --donor-effect --overwrite
  

Rscript scripts/run_real.R Brain2 "Inhibitory neurons" scmultimap \
  --donor-effect --overwrite
Rscript scripts/run_real.R Brain2 "Inhibitory neurons" scmultimap \
  --condition AD --overwrite
Rscript scripts/run_real.R Brain2 Microglia scmultimap \
  --no-donor-effect --overwrite
  
  
Rscript scripts/run_real.R Brain2 "Inhibitory neurons" signac  
Rscript scripts/run_real.R Brain2 Astrocytes signac

# Brain2 differential association (Ctrl vs AD; outputs under Brain2_diff)
# 0. No donor fixed effect + no donor robust
Rscript scripts/run_real_diff.R Brain2 Microglia epilink_poisson \
  --no-donor-robust --overwrite

# 1. No donor fixed effect + donor robust
Rscript scripts/run_real_diff.R Brain2 Microglia epilink_poisson --overwrite

# 2. donor fixed effect + donor robust
Rscript scripts/run_real_diff.R Brain2 Microglia epilink_poisson \
  --donor-effect --donor-robust --overwrite
  
Rscript scripts/run_real_diff.R Brain2 Microglia epilink_nb
Rscript scripts/run_real_diff.R Brain2 Astrocytes epilink_poisson
Rscript scripts/run_real_diff.R Brain2 Astrocytes epilink_nb
Rscript scripts/run_real_diff.R Brain2 Microglia scmultimap 
Rscript scripts/run_real_diff.R Brain2 Astrocytes scmultimap


# Pituitary
# "Corticotropes" "Gonadotropes"  "Lactotropes"   "Somatotropes"  "Stem cells"    "Thyrotropes"  
Rscript scripts/run_real.R Pituitary Corticotropes signac
Rscript scripts/run_real.R Pituitary Gonadotropes signac 
Rscript scripts/run_real.R Pituitary Lactotropes signac 
Rscript scripts/run_real.R Pituitary Somatotropes signac 
Rscript scripts/run_real.R Pituitary "Stem cells" signac 
Rscript scripts/run_real.R Pituitary Thyrotropes signac




Rscript scripts/run_crispr_validation.R Arthritis
Rscript scripts/run_pchic_validation.R Arthritis
Rscript scripts/run_hichip_validation.R pbmc


Rscript scripts/plot_crispr_validation.R pbmc
Rscript scripts/plot_pchic_validation.R Arthritis
Rscript scripts/plot_hichip_validation.R pbmc
