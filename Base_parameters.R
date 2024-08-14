


# This is an R file containing base settings, functions and packages to load.
# Both Sim_Cont_Many.R and Sim_Cont_Results.Rmd loads and runs this code. 



# packages needed
library("tidyverse")
library("dplyr")
library("readxl")
library("metafor")
library("meta")
library("pimeta")
library("plotrix")
library("abind")
library("sn")             # for skewed normal distribution
library("moments")        # for skewed normal distribution
library("stringr")        # for str_sub() function
library("EnvStats")       # for the bimodal distribution which is a mixture of 2 normals



# Setting the simulation. Only relevant for Sim_Cont_Many.R. 

distribution <- "Normal"
#distribution <- "Skewed_Normal_50"
#distribution <- "Skewed_Normal_75"
#distribution <- "Skewed_Normal_99"
#distribution <- "Bimodal"
#distribution <- "Uniform"

# Setting the PI expected coverage level
PI_level <- 0.95


# Setting the values for the simualtion parameters
# 
mu<- 0                                    # mean of random effects

r <- 1000                                   #number of meta-analysis to simulate

k <- c(3, 4, 5, 7, 10, 15, 20, 30, 100)                # number of studies (K)

n_numbers <- c(30, 50, 100, 200, 500, 1000, 2000)     # number of patients in the MA (N)

n_pattern <- c(paste0(n_numbers, "_all"), "_mixed")   # number of patients is either constant (_all) or mixed (vector of 50, 100, 500 repeated until K) 

t2 <-  c( 0.1, 0.2, 0.3, 0.5 ,1, 2, 5)             # tau2 values

within_study_sigma2 <- 10     # variance of the underlying variable (V) 




# Set random effects distribution based on distribution setting,
# The functions generate n random numbers from the distribution with given expected value and standard deviation
if (distribution == "Normal") {
  re_distribution <- function(n, mean, sd) rnorm(n=n, mean = mean, sd=sd)
  
} else if (distribution == "Uniform") {
  re_distribution <- function(n, mean, sd) runif(n=n, min = mean - sqrt(3*sd^2), max = mean + sqrt(3*sd^2))
  
# Skewed normal distributions
# last 2 number defines skewness (multiplied by 100), the cp2dp function transforms the given mean, sd, skewness to the parameters of rsn function
} else if (substr(distribution, 1, 13) == "Skewed_Normal") {
  re_distribution <- function(n, mean, sd)  rsn(n=n, dp =cp2dp( c(mean, sd, as.numeric(str_sub(distribution, - 2, - 1))/100 ), "SN") ) 
  
# Bimodal distribution
}  else if (distribution == "Bimodal") {
  re_distribution <- function(n, mean, sd) rnormMix(n=n, 
                                                    p.mix = 0.5,                # mixture of 2 normals with equal weights
                                                    mean1 =  mean - 0.8 * sd ,    # 0.8 is arbitrary
                                                    mean2 =  mean + 0.8 * sd ,    # expected value will be mu
                                                    sd1 = sqrt( (sd^2) / 5),    # define sd1^2 as tau2/5 
                                                    sd2 = sqrt(  (sd^2 -  0.5 * ((sd^2) / 5) - 0.5 * ( ( (mean - 0.8 * sd) -mean)^2 ) - 0.5*( ((mean + 0.8 * sd)-mean)^2 ))/0.5  )  ) # sd2 is bound to be this to have var as predefined tau2
}





# Set quantile function based on random effects distribution to compute true PI bounds
# The functions generate the quantile of the distribution with given probability (p), expected value and standard deviation
if (distribution == "Normal") {
  re_quantile <- function(p, mean, sd) qnorm(p, mean = mean, sd=sd)
  
} else if (distribution == "Uniform") {
  re_quantile <- function(p, mean, sd) qunif(p, min = mean - sqrt(3*sd^2), max = mean + sqrt(3*sd^2))
  
  # last 2 number defines skewness (multiplied by 100), the cp2dp function transforms the given mean, sd, skewness to the parameters of rsn function  
} else if (substr(distribution, 1, 13) == "Skewed_Normal") {
  re_quantile <- function(p, mean, sd)  qsn(p, dp =cp2dp( c(mean, sd, as.numeric(str_sub(distribution, - 2, - 1))/100 ), "SN") ) 
  
}  else if (distribution == "Bimodal") { 
  re_quantile <- function(p, mean, sd) qnormMix(p, 
                                                p.mix = 0.5,                # mixture of 2 normals with equal weights
                                                mean1 =  mean - 0.8 * sd ,    # 0.8 is arbitrary
                                                mean2 =  mean + 0.8 * sd ,    # expected value will be mu
                                                sd1 = sqrt( (sd^2) / 5),    # define sd1^2 as tau2/5 
                                                sd2 = sqrt(  (sd^2 -  0.5 * ((sd^2) / 5) - 0.5 * ( ( (mean - 0.8 * sd) -mean)^2 ) - 0.5*( ((mean + 0.8 * sd)-mean)^2 ))/0.5  )  ) # sd2 is bound to be this to have var as predefined tau2
}





# Define cumulative distribution function (CDF) for the computation of coverages from PI lower and upper estimates
# q is the given quantile
CDF_normal <- function(q, mean, sd, lower.tail = T) pnorm(q=q, mean = mu, sd = sd, lower.tail = lower.tail )
# CDF for the uniform distribution
CDF_uniform <- function(q, mean, sd, lower.tail = T) punif(q=q, min = mean - sqrt(3*sd^2), max = mean + sqrt(3*sd^2), lower.tail = lower.tail )

# CDF for the skewed normal distributions with skewness parameter 0.5, 0.75, 0.99, only for lower tail
CDF_skewed_normal_50 <- function(q, mean, sd) psn(x=q, dp =cp2dp( c(mean, sd, 0.5 ), "SN") ) 

CDF_skewed_normal_75 <- function(q, mean, sd) psn(x=q, dp =cp2dp( c(mean, sd, 0.75 ), "SN") ) 

CDF_skewed_normal_99 <- function(q, mean, sd) psn(x=q, dp =cp2dp( c(mean, sd, 0.99 ), "SN") ) 


# CDF for the bimodal distribution, only for lower tail
CDF_bimodal <- function(q, mean, sd) pnormMix(q=q, 
                                              p.mix = 0.5,                # mixture of 2 normals with equal weights
                                              mean1 =  mean - 0.8 * sd ,    # 0.8 is arbitrary
                                              mean2 =  mean + 0.8 * sd ,    # expected value will be mu
                                              sd1 = sqrt( (sd^2) / 5),    # define sd1^2 as tau2/5 
                                              sd2 = sqrt(  (sd^2 -  0.5 * ((sd^2) / 5) - 0.5 * ( ( (mean - 0.8 * sd) -mean)^2 ) - 0.5*( ((mean + 0.8 * sd)-mean)^2 ))/0.5  )  ) # sd2 is bound to be this to have var as predefined tau2



# dimensions of array for data generation
dim_generate <- c(r,max(k),length(k) ,length(n_pattern), length(t2) )

# dimensions of array for results
dim_results <- c(r ,length(k) ,length(n_pattern), length(t2)  )

# predefine dimnames for arrays
dimnames_gen <- list(
  "sim_number"=1:r, 
  "study_number" = 1:max(k), 
  "all_study"=k, 
  "n_pattern"=n_pattern, 
  "tau2" = t2
  )

dimnames_results <- list(
  "sim_number"=1:r, 
  "all_study"=k, 
  "n_pattern"=n_pattern, 
  "tau2" = t2
  )


dimnames_results_2<- list(
  "all_study"=k, 
  "n_pattern"=n_pattern, 
  "tau2" = t2
  )


# varnames for results
varlist_results <- c(  "true_PI_lower", "true_PI_upper", "teta_new",
                       "I2", "tau2_REML", "tau2_DL",
                       
                       "HTS_DL_PI_lower", "HTS_DL_PI_upper", 
                       "HTS_REML_PI_lower", "HTS_REML_PI_upper",
                       "HTS_HK_REML_PI_lower", "HTS_HK_REML_PI_upper",
                       "HTS_DL_k_1_PI_lower", "HTS_DL_k_1_PI_upper", 
                       "HTS_DL_z_PI_lower", "HTS_DL_z_PI_upper", 
                       
                       "Wang_PI_lower", "Wang_PI_upper"
)

# function that computes PI bounds for HTS methods and Ensemble method and other necessary statistics
PI_function <- function(y_i, v_i, level = PI_level) {
  
  res <- list() # empty list
  
  res$k <- sum( !is.na(y_i)) # number of studies
  
  temp_DL <- rma.uni(yi = y_i,    # DL estimate from metafor package
                     vi = v_i,
                     method = "DL")
  
  res$tau2_DL <- temp_DL$tau2
  res$SE_DL <-temp_DL$se
  res$mean_DL <-  as.vector(temp_DL$beta)
  res$I2<- temp_DL$I2
    
  # quantiles of t(df=k-2), t(df=k-1) and standard normal dist. 
  res$t_k_2 <- qt(p=1-( (1-level)/2 ), df=res$k-2) 
  res$t_k_1 <- qt(p=1-( (1-level)/2 ), df=res$k-1)
  res$z <-  qnorm(p=1-( (1-level)/2 ))
  
  
  # testing if REML estimate gives an error (Fisher scoring algorithm did not converge)
  temp_REML <- try( rma.uni(yi = y_i,    # REML estimate from metafor package
                            vi = v_i,
                            method = "REML") )
  
  
  res$tau2_REML <- ifelse(inherits(temp_REML, "try-error"),NA, temp_REML$tau2 )   # if convergence failed, assign NA, if convergence succeeded, assign the estimate
  res$SE_REML <-ifelse(inherits(temp_REML, "try-error"),NA, temp_REML$se )   # if convergence failed, assign NA, if convergence succeeded, assign the estimate
  res$mean_REML <-ifelse(inherits(temp_REML, "try-error"),NA, temp_REML$beta )   # if convergence failed, assign NA, if convergence succeeded, assign the estimate
  
  res$SE_HK_REML <- sqrt( (1/ (res$k-1)) *  sum( (( (res$tau2_REML + v_i)^(-1) )* ((y_i-res$mean_REML)^2))/ sum( (res$tau2_REML + v_i)^(-1) , na.rm = T) , na.rm = T) )
  

  # HTS PI with DL tau2 estimate with t distribution df=k-2
  res$HTS_DL_PI <- c( res$mean_DL - res$t_k_2*sqrt(res$tau2_DL + res$SE_DL^2) , res$mean_DL + res$t_k_2*sqrt(res$tau2_DL + res$SE_DL^2 ) )
  
  # HTS PI with REML tau2 estimate
  res$HTS_REML_PI <- c( res$mean_REML - res$t_k_2*sqrt(res$tau2_REML + res$SE_REML^2) , res$mean_REML + res$t_k_2*sqrt(res$tau2_REML + res$SE_REML^2 ) )
  
  # HTS PI with HK SE and REML tau2 estimate
  res$HTS_HK_REML_PI <- c( res$mean_REML - res$t_k_2*sqrt(res$tau2_REML + res$SE_HK_REML^2) , res$mean_REML + res$t_k_2*sqrt(res$tau2_REML + res$SE_HK_REML^2 ) )
  
  # HTS PI with DL tau2 estimate with t distribution df=k-1
  res$HTS_DL_k_1_PI <- c( res$mean_DL - res$t_k_1*sqrt(res$tau2_DL + res$SE_DL^2) , res$mean_DL + res$t_k_1*sqrt(res$tau2_DL + res$SE_DL^2 ) )
  
  # HTS PI with DL tau2 estimate with stnorm distribution
  res$HTS_DL_z_PI <- c( res$mean_DL - res$z*sqrt(res$tau2_DL + res$SE_DL^2) , res$mean_DL + res$z*sqrt(res$tau2_DL + res$SE_DL^2 ) )
  

  
  # compute PI proposed by Wang (original, based on their code and xlsx)

  PI_limit_Wang_original <- function(quantile) {
    
    post_mean_sorted <- sort(res$mean_DL + sqrt( res$tau2_DL/ (res$tau2_DL + v_i) ) * (y_i-res$mean_DL))
    
    AA <- floor(quantile*(res$k+1))
    BB <- AA+1
    
    if (AA == 0) {
      mu_wang <- res$mean_DL
      sigma_Wang <- (post_mean_sorted[1]-mu_wang)/qnorm(1/(res$k+1))
    } else if (AA == res$k) {
      mu_wang <- res$mean_DL
      sigma_Wang <- (post_mean_sorted[res$k]-mu_wang)/qnorm(res$k/(res$k+1))
    }
    else {
      sigma_Wang <- (post_mean_sorted[BB]-post_mean_sorted[AA])/(qnorm(BB/(res$k+1))-qnorm(AA/(res$k+1)))
      mu_wang <- post_mean_sorted[AA]- qnorm(AA/(res$k+1))*sigma_Wang
    }
    
    if (res$tau2_DL > 0) {
      return(qnorm(quantile,mean=mu_wang, sd=sigma_Wang))
    } else {
      return(res$mean_DL)
    }
  }
  
  res$Wang_PI <- c(PI_limit_Wang_original((1-level)/2), PI_limit_Wang_original(1-( (1-level)/2 )) )
  
  
  return(res)
  
}
