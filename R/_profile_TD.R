dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)

rm( list = ls()); gc(); #clean environment

source("PMF.R")
source("reconc_TD.R")
source("utils.R")
source("reconc.R")
source("hierarchy.R")
source("shrink_cov.R")

####

get_hier_M5 = function(save_ = FALSE, save_path = "../../../results/") {
  
  n_b = 3049
  n_u = 11
  n = n_u + n_b
  
  A = matrix(0, nrow = n_u, ncol = n_b)
  A[1,] = 1           # store
  A[2,1:565] = 1      # categories
  A[3,566:1612] = 1
  A[4,1613:3049] = 1
  A[5,1:416] = 1      # departments
  A[6,417:565] = 1
  A[7,566:1097] = 1
  A[8,1098:1612] = 1
  A[9,1613:1828] = 1
  A[10,1829:2226] = 1
  A[11,2227:3049] = 1
  
  S = rbind(A, diag(rep(1, n_b)))
  
  if (save_) {
    saveRDS(S, paste0(save_path, "S.rds"))
    saveRDS(A, paste0(save_path, "A.rds"))
  }
  
  return(list(A = A, S = S))
}

get_gauss_params_bott = function(fc_bottom) {
  
  mu_bott = c()
  sd_bott = c()
  for (fc_b in fc_bottom) {
    tab = fc_b$samples
    vals = as.integer(names(tab))
    mu_b = sum(vals * tab) / sum(tab)
    sd_b = (sum((vals-mu_b)**2 * tab) / sum(tab))**0.5
    mu_bott = c(mu_bott, mu_b)
    sd_bott = c(sd_bott, sd_b)
  }
  l = list(mu = mu_bott, sd = sd_bott)
  
  return(l)
}

get_gauss_params_upp = function(fc_upper) {
  
  mu = c()
  sd = c()
  for (fc in fc_upper) {
    mu = c(mu, fc$mu)
    sd = c(sd, fc$sigma)
  } 
  
  residuals.upper = lapply(fc_upper, "[[", "residuals")
  residuals.upper = t(do.call("rbind", residuals.upper))
  Sigma = schaferStrimmer_cov(residuals.upper)$shrink_cov
  
  l = list(mu = mu, sd = sd, Sigma = Sigma)
  return(l)  
}

PMF.from_tab = function(tab) {
  supp = names(tab)
  pmf = rep(0, max(as.integer(supp)))
  for (j in supp) {
    pmf[[as.integer(j)+1]] = tab[[j]] 
  }
  pmf = pmf / sum(pmf)
  return(pmf)
}

###

seed = 42
h = 1
STORE = "CA_1"
N_samples = 1e4
toll=1e-16
Rtoll=1e-7
smooth_bottom=TRUE
al_smooth=NULL
lap_smooth=FALSE
n_low_u_rprof = 7

CAT = c("HOBBIES", "HOUSEHOLD", "FOODS")
DEPT = c("HOBBIES_1", "HOBBIES_2", "HOUSEHOLD_1", "HOUSEHOLD_2", 
         "FOODS_1", "FOODS_2", "FOODS_3")
names_u = c(STORE, CAT, DEPT)

store_path = paste0("../../../Ricerca/Reconciliation/Codici/Mixed/results/", STORE)

A = get_hier_M5()$A
n_upp = nrow(A)
n_bott = ncol(A)
n = n_bott + n_upp
S = rbind(A, diag(n_bott))

str_bott = paste0(store_path,"/h=",h,"/base_fc_bottom.rds")
str_upp = paste0(store_path,"/h=",h,"/base_fc_upper.rds")
fc_bottom = readRDS(str_bott)
tabs_bottom = lapply(fc_bottom, "[[", "samples")
pmf_bottom = lapply(tabs_bottom, PMF.from_tab) 

fc_upper = readRDS(str_upp)
upper_params = get_gauss_params_upp(fc_upper)

###
# All the upper

Rprof()
out = reconc_TD(S, pmf_bottom, upper_params,
                bottom_in_type = "pmf", N_samples = 1e4,
                return_pmf = TRUE, return_samples = TRUE)
Rprof(NULL)
summaryRprof()

###
# Only one upper

S_ = S[c(1,12:3060),]
upper_params_ = list()
upper_params_$mu = upper_params$mu[1]
upper_params_$Sigma = upper_params$Sigma[1,1]

Rprof()
out = reconc_TD(S_, pmf_bottom, upper_params,
                bottom_in_type = "pmf", N_samples = 1e4,
                return_pmf = FALSE, return_samples = TRUE)
Rprof(NULL)
summaryRprof()

###
# Many upper, but all "lowest-level"

S_ = S[c(5:11,12:3060),]
upper_params_ = list()
upper_params_$mu = upper_params$mu[5:11]
upper_params_$Sigma = upper_params$Sigma[5:11,5:11]

Rprof()
out = reconc_TD(S, pmf_bottom, upper_params,
                bottom_in_type = "pmf", N_samples = 1e4,
                return_pmf = FALSE, return_samples = TRUE)
Rprof(NULL)
summaryRprof()


# Possible speed-up:
# -vectorize .cond_biv_sampling
# Rcpp?



