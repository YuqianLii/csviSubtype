
<!-- README.md is generated from README.Rmd. Please edit that file -->

# csviSubtype

`csviSubtype` implements a GWAS-anchored and vascular-centered framework
for molecular subtype barcode discovery.

The current v0.1 package focuses on three core tasks:

1.  building scRNA-derived subtype gene features from CSS, CIS, LRS, and
    DEG inputs,
2.  fusing scRNA and GWAS gene-level evidence into subtype barcodes,
3.  summarizing robustness under ablation settings.

## Installation

Install from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("YuqianLii/csviSubtype")
```

For local development:

``` r
# install.packages(c("devtools", "data.table"))
devtools::load_all()
```

## Main functions

- `csvi_build_scrna_features()`
- `csvi_fuse_barcode()`
- `csvi_summarize_ablation()`

## Included example data

The package includes example input files under `extdata/`:

- `css_example.tsv`
- `cis_example.tsv`
- `lrs_example.tsv`
- `deg_example.tsv`
- `gwas_feature_example.tsv`

These files are intended to demonstrate the package workflow and output
structure.

## Minimal example

``` r
library(data.table)
library(csviSubtype)

css <- fread(system.file("extdata", "css_example.tsv", package = "csviSubtype"))
cis <- fread(system.file("extdata", "cis_example.tsv", package = "csviSubtype"))
lrs <- fread(system.file("extdata", "lrs_example.tsv", package = "csviSubtype"))
deg <- fread(system.file("extdata", "deg_example.tsv", package = "csviSubtype"))
gwas <- fread(system.file("extdata", "gwas_feature_example.tsv", package = "csviSubtype"))

scrna_res <- csvi_build_scrna_features(
  css = css,
  cis = cis,
  lrs = lrs,
  deg = deg,
  mode = "receptor"
)

fusion_res <- csvi_fuse_barcode(
  gwas = gwas,
  scrna_a1 = scrna_res$features_by_group$A1,
  scrna_b1 = scrna_res$features_by_group$B1
)

names(scrna_res)
#> [1] "css_group"         "cis_all"           "cis_top"          
#> [4] "lrs_use"           "features"          "features_by_group"
names(fusion_res)
#> [1] "gwas_gene_summary" "A1"                "B1"
```

## Inspect scRNA feature results

``` r
head(scrna_res$features)
#>             gene       logFC direction          M sum_contrib groupAB F_gene
#>           <char>       <num>    <char>      <num>       <num>  <char>  <num>
#> 1: 0610009B22RIK -0.59496479  up_in_A1 0.59496479           0      A1      0
#> 2: 0610009O20RIK -8.02181705  up_in_A1 8.02181705           0      A1      0
#> 3: 0610010F05RIK  2.83703947  up_in_B1 2.83703947           0      A1      0
#> 4: 0610010K14RIK  0.29458157  up_in_B1 0.29458157           0      A1      0
#> 5: 0610012G03RIK  0.04144461  up_in_B1 0.04144461           0      A1      0
#> 6: 0610030E20RIK -0.50770565  up_in_A1 0.50770565           0      A1      0
#>          q_M   q_F      q_sc scrna_tier scrna_percentile chainSupport
#>        <num> <num>     <num>     <char>            <num>       <lgcl>
#> 1: 0.3931844     0 0.3931844        low        0.3931844        FALSE
#> 2: 0.9702620     0 0.9702620       high        0.9702620        FALSE
#> 3: 0.8810481     0 0.8810481       high        0.8810481        FALSE
#> 4: 0.2084874     0 0.2084874        low        0.2084874        FALSE
#> 5: 0.0273268     0 0.0273268        low        0.0273268        FALSE
#> 6: 0.3464073     0 0.3464073        low        0.3464073        FALSE
table(scrna_res$features$scrna_tier)
#> 
#> high  low  mid 
#> 2492 6222 3730
```

## Inspect fused subtype ledgers

``` r
nrow(fusion_res$A1$full)
#> [1] 6236
nrow(fusion_res$A1$binary)
#> [1] 13
nrow(fusion_res$A1$strict)
#> [1] 1

nrow(fusion_res$B1$full)
#> [1] 6236
nrow(fusion_res$B1$binary)
#> [1] 13
nrow(fusion_res$B1$strict)
#> [1] 1
```

## View top binary entries

``` r
head(fusion_res$A1$binary[, .(
  gene_key, gwas_pp4_max, gwas_tier, q_sc, scrna_tier, zone_binary
)])
#> Key: <gene_key>
#>    gene_key gwas_pp4_max gwas_tier      q_sc scrna_tier zone_binary
#>      <char>        <num>    <char>     <num>     <char>      <char>
#> 1:     BANP    0.1603770       low 0.3499437        low         red
#> 2:    CELF1    0.7706769       mid 0.6042437        mid         red
#> 3:   GIGYF1    0.2002042       low 0.3968815        low         red
#> 4:     GNB2    0.4549226       low 0.3721267        low         red
#> 5:   MRPL38    0.2044441       low 0.4475165        low         red
#> 6:   NBEAL1    0.8002833      high 0.5675936        mid         red

head(fusion_res$B1$binary[, .(
  gene_key, gwas_pp4_max, gwas_tier, q_sc, scrna_tier, zone_binary
)])
#> Key: <gene_key>
#>    gene_key gwas_pp4_max gwas_tier      q_sc scrna_tier zone_binary
#>      <char>        <num>    <char>     <num>     <char>      <char>
#> 1:     BANP    0.1603770       low 0.3499437        low         red
#> 2:    CELF1    0.7706769       mid 0.6042437        mid         red
#> 3:   GIGYF1    0.2002042       low 0.3968815        low         red
#> 4:     GNB2    0.4549226       low 0.3721267        low         red
#> 5:   MRPL38    0.2044441       low 0.4475165        low         red
#> 6:   NBEAL1    0.8002833      high 0.5675936        mid         red
```

## Output structure

`csvi_build_scrna_features()` returns a list with:

- `css_group`
- `cis_all`
- `cis_top`
- `lrs_use`
- `features`
- `features_by_group`

`csvi_fuse_barcode()` returns a list with:

- `gwas_gene_summary`
- `A1`
- `B1`

Each subtype entry contains:

- `full`
- `binary`
- `strict`

## Notes

- `mode = "receptor"` is currently implemented in `v0.1`.
- `mode = "pathway"` is reserved for future expansion.
- Example files are provided for workflow demonstration and package
  testing.

## Citation

If you use `csviSubtype`, please cite the associated preprint/manuscript
when available.
