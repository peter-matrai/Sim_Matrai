

# Parameters for simulation settings for code run on HPC to compute PI bounderies for the parametric bootstrap method 

distribution <- "Normal"
#distribution <- "Uniform"
#distribution <- "Skewed_Normal_50"
#distribution <- "Skewed_Normal_75"
#distribution <- "Skewed_Normal_99"
#distribution <- "Bimodal"

PI_level <- 0.95


#set.seed(run)                        # set seed for random number generation
mu<- 0                                    # mean of random effects

r <- 1000                                   #number of meta-analysis 

k <- c(3, 4, 5, 7, 10, 15, 20, 30, 100)                # number of studies

n_numbers <- c(30, 50, 100, 200, 500, 1000, 2000)

n_pattern <- c(paste0(n_numbers, "_all"), "_mixed")

t2 <-  c( 0.1, 0.2, 0.3, 0.5 ,1, 2, 5)             # t2 values

within_study_sigma2 <- 10

# set random effects distribution based on distribution setting
if (distribution == "Normal") {
  re_distribution <- function(n, mean, sd) rnorm(n=n, mean = mean, sd=sd)
  
} else if (distribution == "Uniform") {
  re_distribution <- function(n, mean, sd) runif(n=n, min = mean - sqrt(3*sd^2), max = mean + sqrt(3*sd^2))
  
  # last 2 number defines skewness (multiplied by 100), the cp2dp function transforms the given mean, sd, skewness to the parameters of rsn function
} else if (substr(distribution, 1, 13) == "Skewed_Normal") {
  re_distribution <- function(n, mean, sd)  rsn(n=n, dp =cp2dp( c(mean, sd, as.numeric(str_sub(distribution, - 2, - 1))/100 ), "SN") ) 
  
}  else if (distribution == "Bimodal") {
  re_distribution <- function(n, mean, sd) rnormMix(n=n, 
                                                    p.mix = 0.5,                # mixture of 2 normals with equal weights
                                                    mean1 =  mean - 0.8 * sd ,    # 0.8 is arbitrary
                                                    mean2 =  mean + 0.8 * sd ,    # expected value will be mu
                                                    sd1 = sqrt( (sd^2) / 5),    # define sd1^2 as tau2/5 
                                                    sd2 = sqrt(  (sd^2 -  0.5 * ((sd^2) / 5) - 0.5 * ( ( (mean - 0.8 * sd) -mean)^2 ) - 0.5*( ((mean + 0.8 * sd)-mean)^2 ))/0.5  )  ) # sd2 is bound to be this to have var as predefined tau2
}


# set quantile function based on random effects distribution to compute true PI bounds
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


# Define CDF-s for the computation of quantile coverages from PI lower and upper estimates

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



