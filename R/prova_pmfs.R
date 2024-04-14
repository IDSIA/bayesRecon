dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)
source("reconc_TD.R")
source("PMF.R")

pmf1 = sample(1:30, size = 1e3, replace = TRUE)
pmf1 = pmf1 / sum(pmf1)
pmf2 = sample(1:30, size = 1e3, replace = TRUE)
pmf2 = pmf2 / sum(pmf2)

u = sample(0:1999, size = 1e4, replace = TRUE)

b = cond_biv_sampling(u, pmf1, pmf2)






outer(u, supp1, `-`)


