# Joint_Model_Gap_Time
**************************************************************************************************************************
*****************************************Programs for the Study*********************************************************

*******************Part I Simulation Setting I to Setting V*************************
1. Generate Datasets for setting I to V. R
    - Used to generate 200 datasets for each setting

2. Macro program setting I to V. sas
   - Estimate the model parameters under settings I to V

 
3. Summary the fitting results from setting I to IV.R
   -Summary the estimation of model parameters under settings I to IV
   - Results shown in Table S1-S4



*******************Part II Self Defined Function*************************
1. Lambdainv.R
    -Function used to find vector T such that Lambda(T)=s for a vector input s

2. myprod.R
    -Function used to product of two step function

3. mysimRec
   -Function using Cinlar's inversion Method to generate Non-homogeneous Poisson process

4. mysimRec_gap
   -Function using Cinlar's inversion Method to generate Non-homogeneous Poisson process for gap time and then convert to study time

5. datagen_M_gap.R
   -Function used to generate recurrent events M0 and M1

6. datagen_X_gap.R
   -Function used to generate datasets with simulated X and Z

7. datagenQ13_M_gap.R
   -Function used to generate recurrent events M_Q1 and M_Q3
   -Q1 and Q3 are two comparing categories for NDE/NIE

8. mymed.R
   -Function used to compute the NDE/NIE for simulation settings

9. mymed_Q13.R
   -Function used to compute the NDE/NIE with Z containing two categories Q1 and Q3

10. S.R
    - Used to compute the survival functions

11. S00.R
     - Compute the survival functions S00

12. S01.R
     - Compute the survival functions S01

13. S10.R
     - Compute the survival functions S10
 
14. S11.R
     - Compute the survival functions S11
 
   

*******************Part III Estimate NDE/NIE Under Simulation Setting I to IV*************************

1. Folder "NDE NIE from Simulation I"
    *3-1. SimulationI_parallel_gap.R
              - Function used to estimate NDE/NIE under simulation setting I
    *3-2. SimulationI_parallel_boot_gap.R
             -Function of the bootstrap for NDE/NIE under simulation setting I
     *3-3. SimulationI_summary_gap.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting I


2. Folder "NDE NIE from Simulation II"
    *4-1. SimulationII_parallel_gap.R
              - Function used to estimate NDE/NIE under simulation setting II
    *4-2. SimulationII_parallel_boot_gap.R
             -Function of the bootstrap for NDE/NIE under simulation setting II
    *4-3. SimulationII_summary_gap.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting II


3. Folder "NDE NIE from Simulation III"
    *5-1. SimulationIII_parallel_gap.R
              - Function used to estimate NDE/NIE under simulation setting III
    *5-2. SimulationIII_parallel_boot_gap.R
             -Function of the bootstrap for NDE/NIE under simulation setting III
    *5-3. SimulationIII_summary_gap.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting III


4. Folder "NDE NIE from Simulation IV"
    *6-1. SimulationIV_parallel_gap.R
              - Function used to estimate NDE/NIE under simulation setting IV
    *6-2. SimulationIV_parallel_boot_gap.R
             -Function of the bootstrap for NDE/NIE under simulation setting IV
    *6-3. SimulationIV_summary_gap.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting IV


5. Folder "NDE NIE from Simulation V"
    *7-1. SimulationV_parallel_gap.R
              - Function used to estimate NDE/NIE under simulation setting V
    *7-2. SimulationV_parallel_boot_gap.R
             -Function of the bootstrap for NDE/NIE under simulation setting V
    *7-3. SimulationV_summary_gap.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting V



*******************Part IV Real Data Analysis*************************

1. CPCRA study analysis.sas
    -Estimates the parameters under simulation setting I and II for CPCRA study

2. Folder "CPCRA NDE NIE Plot"

   * NDE_NIE_Gap_Trt.R 
      -Estimate NDE/NIE for treatment

   * NDE_NIE_Gap_CI_Trt.R
     -Boostrap of NDE/NIE for treatment

   * Summary_Gap_Trt.R
    -Summary the results and plot the estimates of NDE/NIE/TE with bootstrapped 95% CI for treatment

  * NDE_NIE_Gap_CD4.R 
      -Estimate NDE/NIE for CD4

   * NDE_NIE_Gap_CI_CD4.R
     -Bootstrap of NDE/NIE for CD4

   * Summary_Gap_CD4.R
    -Summary the results and plot the estimates of NDE/NIE/TE with bootstrapped 95% CI for CD4



 
