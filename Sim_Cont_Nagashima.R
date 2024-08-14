

# This is an R code to run only on HPC to compute PI boundaries with the parametric bootstrap method

start_time <- Sys.time()


# Set basedir in HPC
basedir <- "~/Normal/"
#basedir <- "~/Uniform/"
#basedir <- "~/Skewed_Normal_50/"
#basedir <- "~/Skewed_Normal_75/"
#basedir <- "~/Skewed_Normal_99/"
#basedir <- "~/Bimodal/"



# load necessary packages
library(pimeta, lib.loc = '~/local/R/4.2.1.2/library')
library(parallel)

# set number of cores for the HPC to use
Cores <- 1500

# Which run to evaluate (for this method only 1 run is simulated with r=1000 repetitions) 
 run <- "1"


# load base parameters and functions
source(file= paste0(basedir, "Base_parameters_HPC.R" ) )


# load previously generated data
load(file = paste0(basedir,"Gen_list_",distribution, "_Run_",run, ".RData",sep="") )

# assign(paste0("Gen_list_" ,distribution, "_Run_",run), gen_list)
# rm(gen_list)

# Number of simulations in this dist and run
dim_sim<- dim(gen_list[["Y"]])[1]
# set number of sim based on generated data in this dist and sim
dim_results[1]<- dim_sim
# set dimnames for first dim
dimnames_results$sim_number<- 1:dim_sim


varlist_2_Nag <- c( "Nagashima_PI_lower", 
                    "Nagashima_PI_upper", 
                    "Nagashima_true_PI_lower",
                    "Nagashima_true_PI_upper",
                    "Nagashima_teta_new")

res_list_Nag <- list()

# assign empty arrays to varnames with predefined dimensions and dimnames
for (o in 1:length(varlist_2_Nag)) {
  assign(varlist_2_Nag[o],array(dim=dim_results, dimnames = dimnames_results) )
  res_list_Nag[[o]] <- get(varlist_2_Nag[o])
}

names(res_list_Nag) <- varlist_2_Nag


    for (t in 1:length(t2)) {                 # which t2
      for (g in 1:length(n_pattern)) {        # n of total patients
        for (h in 1:length(k)) {              # n of studies
          for (i in 1:dim_sim) {              # which simulation
            
            #Compute Nagashima PI by calling pima function from pimeta package and save it to a temporary object
              Nagashima_object<- pima(
              y =gen_list[["Y"]][i, 1:k[h], h, g, t],
              #se = NULL,
              v = gen_list[["sample_SE2"]][i, 1:k[h], h, g, t],
              alpha = 0.05,
              method = "boot",
              #method = "HTS",
              #method = "HK",
              #method = "SJ",
              #method = "KR",
              #method = "CL",
              # method = c("boot", "HTS", "HK", "SJ", "KR", "CL", "APX"),
              B = 5000,
              parallel = Cores
            )
            
            res_list_Nag[["Nagashima_PI_lower"]][i , h, g, t]  <- Nagashima_object$lpi
            res_list_Nag[["Nagashima_PI_upper"]][i , h, g, t]  <- Nagashima_object$upi
            
            
            # generate true lower and upper PI bounds
            res_list_Nag[["Nagashima_true_PI_lower"]][i, h, g, t] <- re_quantile(p=(1-PI_level)/2 , mean = mu, sd=sqrt(t2[t]))
            res_list_Nag[["Nagashima_true_PI_upper"]][i, h, g, t] <- re_quantile(p=1-((1-PI_level)/2), mean = mu, sd=sqrt(t2[t]))
            
            # generate a new teta for each study to see if the PI contains it
            res_list_Nag[["Nagashima_teta_new"]][i, h, g, t] <- re_distribution(n= 1, mean = mu, sd=sqrt(t2[t]) )
        
      }
    }
  }
} 


save(res_list_Nag,
     file = paste0(basedir,"Res_list_Nag_",distribution, "_Run_",run, ".RData",sep=""))

# N of MA
N_MA <-  length(t2) * length(n_pattern) * length(k) * length(dim_sim)

# end time of computation
end_time <- Sys.time()

# time needed for computation
Time_needed <- end_time - start_time

writeLines(paste("N_MA=", N_MA, "Time=" , as.character(Time_needed)), paste0(basedir, "Time_needed_",distribution,"_" ,run, ".txt" )  )



