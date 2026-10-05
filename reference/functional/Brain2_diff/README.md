# Brain2 differential functional-validation references

- `anderson_2023_table_s3.xlsx`: Anderson et al., Cell Genomics (2023),
  Supplementary Table S3. It contains AD-versus-control DEGs within each cell
  type. The Microglia definition supplied by the authors is adjusted p-value
  below 0.01 and absolute log2 fold-change above 0.25.
- `mathys_2024_microglia_differential_expression.tsv.gz`: the
  `Mic_Immune_Mic` table extracted from Mathys et al., Nature (2024),
  Supplementary Table 9. `col_nm != 0` denotes concordant significant MAST and
  Nebula results; `col_nm` 1 and 2 indicate down- and upregulation.

The validation script uses the Mathys `nrad` (NIA-Reagan AD versus non-AD),
`allregions` results to match a disease-control comparison.
