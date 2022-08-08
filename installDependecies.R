options(repos=structure(c(CRAN="https://cloud.r-project.org")))
install.packages(pkgs=c(
  "openxlsx",
  "rsq",
  "dplyr",
  "plyr",
  "matrixStats",
  "Rcpp",
  "tidyr",
  "readr"
),
repos="https://cloud.r-project.org",
Ncpus=4,
quiet = TRUE
)
