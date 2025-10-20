# This is a Readme file to the R codes available to the article 'Assessing the properties of the prediction interval in random-effects meta-analysis' published in Research Synthesis Methods  
  
- Code written by: Péter Mátrai, Zoltán Sipos  
- Code maintainer: Péter Mátrai (peter.matrai@pte.hu)  
- The R version the code was written in: version 4.1.3.
- The necessary R packages and their versions: see **Base_parameters.R** file
- We recommend using R studio to run the codes
- Last modified: 20/10/2025


## Content of the Github repository

- Base_parameters.R
- Sim_Cont_Many.R
- Sim_Cont_Results.Rmd
- PI_Sim_Matrai.Rproj
- README.md


## For what purpose you can use these codes
You can replicate our simulation with your own settings. You can freely modify the parameters, test other true effects distributions, test other PI methods or modify the ones we used. You can easily modify the code to define new performance measures and display the results. The details of our simulation methods are described in the article in section 4.  

We know that our codes are not very user friendly. This is because we developed it for ourselves, specifically for our research article, and designed it the way that we can test our ideas. However, we decided to make it public because we think it might be useful for other researchers.    


## How to run the codes

1. First you need to determine a **base directory** on your own computer and **download the 4 files** to this directory from this Github repository (Base_parameters.R, Sim_Cont_Many.R, Sim_Cont_Results.Rmd, PI_Sim_Matrai.Rproj).


2. You need to have a very specific folder system in your base directory, which you can create one of the two following ways:   
    * If you would like to generate 'input data' (simulated meta-analyses with predifined settings), first you need to create the following folder system in your base directory, **manually**:
      - create a folder in your base directory named **RE_distributions** with 6 empty folders in it with the following names:
        - **Normal**  
        - **Uniform**  
        - **Bimodal**  
        - **Skewed_Normal_50**  
        - **Skewed_Normal_75**  
        - **Skewed_Normal_99**  

    * Alternatively, **you can download some or all of the input files that we generated** and used in our simulation through the following link: [Our input files](https://mega.nz/folder/d7o0VI6I#IDIanHGoP21XiwcZEPaYhg).
      Download it as a zipped folder and then put the unzipped folder in your base directory, so that you have the same folder system in the base directory as descibed above. The whole content is 13.6 GB.
      You do not need to download the whole content, you can download any part of these and place them in the appropriate location in your base directory.    


3. Open the **PI_Sim_Matrai.Rproj** R project file in your R studio. Open the **Base_parameters.R** file in the program. At the top of this file, from line 11, you can see the **R packages needed** to run the codes. Make sure that these packages are installed.


4. Generating input data
    * If you want to generate your own input data, open the **Sim_Cont_Many.R** file in the R Studio. If you downloaded some or all our input data and you want to display the results computed only from those, you can skip this step and proceed to step 5.
    * In the **Base_parameters.R** file you can set from which distribution you want to generate input meta-analyses. You can set it between lines 28 and 33.  Define the distribution by giving the appropriate value to the variable named **distribution**, and **comment         out all the other lines**.
    * Also, in line 44, you can define the number of meta-analyses to simulate in one 'Run', by giving a value to the variable named **r**. The default is `r <- 1000`, but it may take hours to generate them depending on your system. You can set **r** to a lower value to           reduce the running time.  
    * By running the **Sim_Cont_Many.R** file, you can generate **r** meta-analyses with the chosen true effects distribution and compute the PI boundaries for these with the different PI methods investigated in the article. The only exception is the **parametric              bootstrap method**, which is an extremely computation intensive method as we describe it in section 4.4 in the article. If you want to display the results for this method, you need to download the input data for this method [from our online repository](https://mega.nz/folder/d7o0VI6I#IDIanHGoP21XiwcZEPaYhg), and put this file in the appropriate distribution folder in your **RE_distributions** folder. Each distribution folder contains an input file for this method with r=1000 repetitions. For example, in the folder named **Normal** this file is called **Res_list_Nag_Normal_Run_1.RData**. In the other distribution folders, it is named similarly, but the distribution name is different.  
    * In the **Sim_Cont_Many.R** file you can set how many times to run the generating process. This was necessary for us, because we had not enough memory to generate 5000 repetitions. Instead, we run the whole process 5 times with 1000 repetitions and merged them in a later step. We call these **Runs**. You can set the number of runs in line 26 as defining a for loop. The default is `for (m in 1:1) {`, which means that the process will be run one time. You only need to set more runs if you set **r** to a high number and you run out of memory.  
    * **The Sim_Cont_Many.R file has two .RData output.** The code saves them both in the distribution folder, that was set in the **Base_parameters.R** file, because it genererates true effects from this chosen distribution.
      
      - The first part of the code in the **Sim_Cont_Many.R** file generates the meta-analyses (e.g. number of studies in the 2 groups, standard errors, true effects, observed effects, etc.). The code generates these in a list containing 11 arrays, each array is a characteristic of the generated meta-analyses. This list is saved as an R.Data file, which has the general name: **Gen_list_distribution_Run_run.RData**. For example, if you set the normal distribution with 1 Run, the name of the file will be **Gen_list_Normal_Run_1.RData** and the program will place it in the folder called Normal.
     
      - The second part of the **Sim_Cont_Many.R** file computes the lower and upper PI boundaries using the investigated PI methods. This part of the code also computes some additional properties, e.g. the theoretical PI limits, the I square value, etc. These are also stored in arrays and put in a list. This list is also saved as an R.Data file, which has the general name: **Res_list_Many_distribution_Run_run.RData**. For example, if you set the skewed normal distribution with skewness parameter 75 (Skewed_Normal_75) with 1 Run, the name of the file will be **Res_list_Many_Skewed_Normal_75_Run_1.RData** and the program will place it in the folder called Skewed_Normal_75.  

    * If you have created the necessary folder system, and you have set everything in the **Base_parameters.R** file and the **Sim_Cont_Many.R** file, just run the whole code in the **Sim_Cont_Many.R** file. The 2 output files will be generated and placed in the appropriate folder. 


5. Creating a HTML file to display the results
    * Open the **Sim_Cont_Results.Rmd**, knitting this file will create a html file showing the results, based on your settings.
    * This file loads the necessary Res_list files from the appropriate folders, merges the Runs, computes the performance measures described in section 4.6 of the article, creates the performance plots and puts these together in a html file.   
    * You can set the name of the html in line 11 by replacing **Sim_PI** with any other name in the code: `paste0("Sim_PI", ".html")))`.
    * All the other necessary settings are in the **chunk named 'base_settings'**, between lines 51-155 .
    * You do not need to modify the basedir setting, it will be automatically the folder where your PI_Sim_Matrai.Rproj R project file is located.
    * Between lines 87-103 you can set for which distributions you want to display results and which Runs you want to use. You need to have the appropriate Res_list RData file(s) in the appropriate distribution folder(s), otherwise you will get an error. Set '1' if you want to show the given distribution and set '0' if you do not want. Also, set the Runs in the format 1:k, if you want to process Run_1, Run_2, ...Run_k. 
    * For example, if you want to show results for the normal distribution of 2 runs, you need to following setting:
      - `Normal_to_display <- 1`
      - `Normal_run_to_display <- 1:2`  
      You also need to have the **Res_list_Many_Normal_Run_1.RData** and **Res_list_Many_Normal_Run_2.RData** files in your folder for the Normal distribution.
    * Between lines 110-116, you can choose which methods you want to display. Set '1' to display the given method and set '0' otherwise.
    * In line 125 you can set whether you want to display the parametric bootstrap method. If you do, you need to download the Res_list file(s) for this method and place them in the appropriate distribution folder, as it is explained in the previous section.
    * You do not need to modify the settings between lines 129-139.
    * In line 145, you can specify whether you want to display the histograms of the coverage probabilities. Set TRUE of FALSE.
    * In line 149, you can specify whether you want to display the 5 performance plots. Set TRUE of FALSE. If you set TRUE, you need to have at least the **Gen_list_Normal_Run_1.RData** file in your 'Normal' folder, because the program computes the true mean standard errors from this file.
    * If you set many things to display, the creation of the html file may take long hours. You can reduce the running time by setting fewer things to display. 
    * If you have set everything you wanted, simply knit the document.
    * The html file containing all of our simulation results is available as an [Online Supplementary Material](https://mega.nz/file/UqIT1QwL#XUhSyea-oApoFYpPjqZGudkfu0kV9xCdJ0xhsKadYP0) for the article. (Reference number 35 in the article).

