

# This an R file with which r repetations can be generated and the PI limits can be calculated for them with the HTS methods and the Ensemble method
# m, the number of "runs" can be specified, this defines how many times this code will be run, each with r repetations 
# r (the number of independent repetations) and all other parameters (RE distribution, K, N, etc.) can be set in the Base_parameters.R file
# A basedir folder should be set, where the output files will be saved in separate folders
# In the basedir folder, an "RE_distributions" folder should be created beforehand, in which folders with the names of the distributions should be created 
# These are the folders where the output files will be saved by this code


# starting time of computation 
start_time <- Sys.time()

# This is where you can set how many times to run the code 
for (m in 1:5) {

# Set basedir folder, in which an "RE_distributions" folder should be created beforehand, in which folders with the names of the distributions should be created beforehand  

basedir <- "D:/Matrai_Peter/Sim_To_Article/"
  



# Nothing to set from this point on

  
run <- m

source(file=paste0(basedir, "Base_parameters.R") )

set.seed(run)


# create empty list to store generated data
gen_list <- list()

# create varnames
varlist_generate <- c("N", "n_T", "n_C", "sigma2_T", "sigma2_C", "s2_T", "s2_C", "sample_SE2", "true_SE2", "teta", "Y")


# assign empty arrays to varnames with predefined dimensions and dimnames
for (o in 1:length(varlist_generate)) {
  assign(varlist_generate[o],array(dim=dim_generate, dimnames = dimnames_gen) )
  gen_list[[o]] <- get(varlist_generate[o])
  
}

names(gen_list) <- varlist_generate


# generation of the input data (r repetations)
    for (t in 1:length(t2)) {                 # which t2
      for (g in 1:length(n_pattern)) {        # n of total patients
        for (h in 1:length(k)) {              # n of studies
          for (i in 1:r) {                    # which simulation
            
        #generation of total study size based on n_pattern
        if (sub(".*_", "", n_pattern[g]) == "all") {   # if n pattern is eq, patients number is the same in for all studies
          gen_list[["N"]][i, 1:k[h], h, g, t] <- as.numeric(sub("_.*", "", n_pattern[g]))

        } else if  (sub(".*_", "", n_pattern[g]) == "mixed") { # if n pattern is mixed, a random mean N is chosen for each study (temp_n) and then N is sampled from uniform distribution with mean N
          temp_n <- sample(n_numbers, k[h], replace = T)
          gen_list[["N"]][i, 1:k[h], h, g, t] <- rep_len( c(50, 100, 500) , length.out = k[h])
        }
          
            
        # generation of sample sizes in the 2 groups 
          gen_list[["n_T"]][i, 1:k[h], h, g, t] <- gen_list[["N"]][i, 1:k[h], h, g, t] * 0.5
          gen_list[["n_C"]][i, 1:k[h], h, g, t] <- gen_list[["N"]][i, 1:k[h], h, g, t] - gen_list[["n_T"]][i, 1:k[h], h, g, t]
          

        # generation of true within study variances in the 2 groups 
          gen_list[["sigma2_T"]][i, 1:k[h], h, g, t] <- within_study_sigma2
          gen_list[["sigma2_C"]][i, 1:k[h], h, g, t] <- within_study_sigma2
          
        # generate estimated within study variances (s2)  in the 2 groups based on scaled chi2 dist df=n-1
          gen_list[["s2_T"]][i, 1:k[h], h, g, t] <- (rchisq(n= k[h], df=gen_list[["n_T"]][i, 1:k[h], h, g, t][1:k[h]] -1)/ (gen_list[["n_T"]][i, 1:k[h], h, g, t][1:k[h]] -1) )* gen_list[["sigma2_T"]][i, 1:k[h], h, g, t][1:k[h]]
          gen_list[["s2_C"]][i, 1:k[h], h, g, t] <- (rchisq(n= k[h], df=gen_list[["n_C"]][i, 1:k[h], h, g, t][1:k[h]] -1)/ (gen_list[["n_C"]][i, 1:k[h], h, g, t][1:k[h]] -1) )* gen_list[["sigma2_C"]][i, 1:k[h], h, g, t][1:k[h]]                                  
                                            
        # compute sample and true study SE2
          gen_list[["sample_SE2"]][i, 1:k[h], h, g, t] <- (gen_list[["s2_T"]][i, 1:k[h], h, g, t] / gen_list[["n_T"]][i, 1:k[h], h, g, t]) + (gen_list[["s2_C"]][i, 1:k[h], h, g, t] / gen_list[["n_C"]][i, 1:k[h], h, g, t])
          gen_list[["true_SE2"]][i, 1:k[h], h, g, t]   <- (gen_list[["sigma2_T"]][i, 1:k[h], h, g, t] / gen_list[["n_T"]][i, 1:k[h], h, g, t]) + (gen_list[["sigma2_C"]][i, 1:k[h], h, g, t] / gen_list[["n_C"]][i, 1:k[h], h, g, t])
          
        # generate teta for each study, the expected value of the study
          gen_list[["teta"]][i, 1:k[h], h, g, t] <- re_distribution(n= k[h], mean = mu, sd=sqrt(t2[t]) )
          
        # generate Y, the study estimates for teta
          gen_list[["Y"]][i, 1:k[h], h, g, t] <- rnorm(n= k[h], mean = gen_list[["teta"]][i, 1:k[h], h, g, t], sd= sqrt(gen_list[["true_SE2"]][i, 1:k[h], h, g, t] ) )
          
        
      }
    }
  }
} 


# This is where each gen_list files will be saved in the basedir folder in the "RE_distributions" folder 
save(gen_list,
     file = paste0(basedir, "RE_distributions/", distribution, "/Gen_list_", distribution, "_Run_", run, ".RData") )



# varlist_results is defined in Base_parameters

res_list <- list()

# assign empty arrays to varnames with predefined dimensions and dimnames
for (o in 1:length(varlist_results)) {
  assign(varlist_results[o],array(dim=dim_results, dimnames = dimnames_results) )
  res_list[[o]] <- get(varlist_results[o])
}

names(res_list) <- varlist_results


# This is where the PI limits and other necessary values will be computed based on the previously generated and saved gen_list data
# The computed numbers will be saved in a list named "res_list"

    for (t in 1:length(t2)) {                 # which t2
      for (g in 1:length(n_pattern)) {        # n of total patients
        for (h in 1:length(k)) {              # n of studies
          for (i in 1:r) {                    # which simulation
            
            # generate true lower and upper PI bounds
            res_list[["true_PI_lower"]][i, h, g, t] <- re_quantile(p=(1-PI_level)/2 , mean = mu, sd=sqrt(t2[t]))
            res_list[["true_PI_upper"]][i, h, g, t] <- re_quantile(p=1-((1-PI_level)/2), mean = mu, sd=sqrt(t2[t]))
            
            # generate a new teta for each study to see if the PI contains it
            res_list[["teta_new"]][i, h, g, t] <- re_distribution(n= 1, mean = mu, sd=sqrt(t2[t]) )
            
            # Temporary list to save PI_function output 
            Temp_list<- PI_function(  gen_list[["Y"]][i, 1:k[h], h, g, t], gen_list[["sample_SE2"]][i, 1:k[h], h, g, t] )
            
            # Extract meta analysis estimates and results from temporary list
            res_list[["I2"]][i , h, g, t]  <- Temp_list$I2
            
            res_list[["tau2_REML"]][i , h, g, t]  <- Temp_list$tau2_REML
            res_list[["tau2_DL"]][i , h, g, t]  <- Temp_list$tau2_DL
            
            
            
            res_list[["HTS_DL_PI_lower"]][i , h, g, t]  <- Temp_list$HTS_DL_PI[1]
            res_list[["HTS_DL_PI_upper"]][i , h, g, t]  <- Temp_list$HTS_DL_PI[2]
            
            res_list[["HTS_REML_PI_lower"]][i , h, g, t]  <- Temp_list$HTS_REML_PI[1]
            res_list[["HTS_REML_PI_upper"]][i , h, g, t]  <- Temp_list$HTS_REML_PI[2]

            res_list[["HTS_HK_REML_PI_lower"]][i , h, g, t]  <- Temp_list$HTS_HK_REML_PI[1]
            res_list[["HTS_HK_REML_PI_upper"]][i , h, g, t]  <- Temp_list$HTS_HK_REML_PI[2]
            
            res_list[["HTS_DL_k_1_PI_lower"]][i , h, g, t]  <- Temp_list$HTS_DL_k_1_PI[1]
            res_list[["HTS_DL_k_1_PI_upper"]][i , h, g, t]  <- Temp_list$HTS_DL_k_1_PI[2]
            
            res_list[["HTS_DL_z_PI_lower"]][i , h, g, t]  <- Temp_list$HTS_DL_z_PI[1]
            res_list[["HTS_DL_z_PI_upper"]][i , h, g, t]  <- Temp_list$HTS_DL_z_PI[2]
            
            
            
            res_list[["Wang_PI_lower"]][i , h, g, t]  <- Temp_list$Wang_PI[1]
            res_list[["Wang_PI_upper"]][i , h, g, t]  <- Temp_list$Wang_PI[2]
            
        
      }
    }
  }
} 


# This is where the "res_list" file will be saved in the basedir folder in the "RE_distributions" folder 
save(res_list,
     file = paste0(basedir,"RE_distributions/",distribution, "/Res_list_Many_", distribution, "_Run_", run, ".RData") )


rm(list=setdiff(ls(), "start_time")) #will clear all objects but start time .
gc() #free up memory and report the memory usage.

}   




# end time of computation
end_time <- Sys.time()

# time needed for computation
end_time - start_time





