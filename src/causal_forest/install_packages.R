# R dependencies for the causal-forest scripts. Run once:
#   Rscript src/causal_forest/install_packages.R
pkgs <- c("grf", "tidyverse", "reticulate")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
# reticulate needs a Python with joblib + pandas to export the web-app tables:
#   reticulate::py_install(c("joblib", "pandas"))
