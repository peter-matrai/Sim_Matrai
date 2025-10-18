# This a Readme file to the R codes available to the article 'Assessing the properties of the prediction interval in random-effects meta-analysis' published by Research Synthesis Methods  
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

## For what purpose you can use these codes
You can replicate our simulation with your own settings. You can freely modify the parameters, test other true effects distributions, test other PI methods or modify the ones we used. You can easily modify the code to define new performance measures and display the results. The details of our simulation are described in the article in section 4. 

## How to run the codes

1. First you need to determine a **base directory** on your own computer and **download the 4 files** to this directory from this Github repository (Base_parameters.R, Sim_Cont_Many.R, Sim_Cont_Results.Rmd, PI_Sim_Matrai.Rproj)

2. You need to have a very specific folder system in your base directory, which you can create one of the two following ways:   
    * If you would like to generate 'input data' (simulated meta-analyses with the predifined settings), compute the PI limits of the tested methods, and visualize the results, first you need to create the following folder system in your base directory, **manually**:
      - create a folder named **RE_distributions** with 6 empty folders in it with the following names:
        - **Normal**  
        - **Uniform**  
        - **Bimodal**  
        - **Skewed_Normal_50**  
        - **Skewed_Normal_75**  
        - **Skewed_Normal_99**  

    * Alternatively, **you can download all of the input files that we generated** and used in our simulation through the following link: [Our input files](https://mega.nz/folder/d7o0VI6I#IDIanHGoP21XiwcZEPaYhg).
      Download it as a zipped folder and then put the unzipped folder in your base directory, so that you have the same folder system in the base directory as descibed above. The whole content is 13.6 GB.
      You do not need to download the whole content, you can download any part of these and place them in the appropriate location in your base directory.    

3. Open the **PI_Sim_Matrai.Rproj** R project file in your R studio. Open the **Base_parameters.R** file in the program. At the top of this file, from line 11, you can see the **R packages needed** to run the codes. Make sure that these packages are installed.

4.
    * If you want to generate your own input data, open the **Sim_Cont_Many.R** file in the R Studio. If you downloaded our input data and want to display the results computed from those, you can skip this step and proceed to the next one.
    * In the **Base_parameters.R** file you can set from which distribution you want to generate input meta-analyses. You can set it between lines 28 and 33.  
      Define the distribution by giving the appropriate value to the variable named **distribution**, and **comment out all the other lines**. 

### nagyon sub
dgzdfgfdgf  
fdgdgdfgdfgdf  
dfgdfgdfgdfgdfgdf





