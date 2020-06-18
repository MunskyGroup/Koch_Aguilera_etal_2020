Codes for "Quantifying the spatiotemporal dynamics of IRES versus Cap translation with single-molecule resolution in living cells"
=======

Amanda Koch<sup>1</sup>, Luis Aguilera<sup>2</sup>, Tatsuya Morisaki<sup>1</sup>, Brian Munsky<sup>2</sup>, and Timothy J. Stasevich<sup>1,3</sup>. <br/>

<sup>1</sup> Department of Biochemistry and Molecular Biology, Colorado State University, Fort Collins, CO 80523, USA. <br/>
<sup>2</sup> Keck Scholars, Department of Chemical and Biological Engineering and School of Biomedical Engineering, Colorado State University, Fort Collins, CO 80523, USA. <br/>
<sup>3</sup> World Research Hub Initiative, Institute of Innovative Research, Tokyo Institute of Technology, Yokohama, Kanagawa 226-8503, Japan. <br/>


<div align="center">
    <img ![Image](../blob/master/GA.jpg?raw=true) width="600px"</img>
</div>

For questions about the codes, please contact:  Luis.aguilera@colostate.edu and brian.munsky@colostate.edu <br/>

---
This repository contains the codes necessary to reproduce figures from the above manuscript. All codes are implemented in Matlab 2019b. <br/>

## Code organization <br/>

The codes are organized in the following categories. <br/>

* **Simulations for the selected model and the optimized parameter values** <br/>
 run__Desktop_Selected_Model.m <br/>

* **Comparison for all tested models** <br/>
run__Desktop_Compare_All_Models.m <br/>
run__Cluster_Compare_All_Models.m <br/>

* **Parameter estimation routines** <br/>
 run__Cluster_Optmization.m <br/>
 run__Desktop_Optmization.m <br/>

* **Parameter uncertainty** <br/>
 run__Cluster_Parameter_Uncertainty.m <br/>

---
## Tested Models. <br/>
* **Model abbreviation. Model number. Model description.** <br/>
3S  --  **Model 1.** -- 3 States .  <br/>
3S_C  --  **Model 2.** -- 3 States + Cross-over.  <br/>
4S_DD  --  **Model 3.** -- 4 States .  <br/>
4S_DDC --  **Model 4.** --  4 States + Cross-over <br/>
4S_II  --  **Model 5.** --  4 States - Independent. <br/>
4S_IIC  --  **Model 6.** --  4 States - Independent + Cross-over . <br/>
4S_DI  --  **Model 7.** -- 4 States - CAP is independent of IRES, but IRES depends on CAP. <br/>
4S_DIC  --  **Model 8.** -- 4 States - CAP is independent of IRES, but IRES depends on CAP. + Cross-over . <br/>
4S_ID  -- **Model 9.** -- 4 States - IRES is independent of CAP, but CAP depends on IRES. <br/>
4S_IDC  --  **Model 10.** -- 4 States - IRES is independent of CAP, but CAP depends on IRES + Cross-over . <br/>
4S_Im1 --  **Model 11.** -- 4 States - IRES activation rates dependent and inactivation rates dependent.  <br/>
4S_Im1C --   **Model 12.** -- 4 States - IRES activation rates dependent and inactivation rates dependent + Cross-over . <br/>
4S_Im2 --  **Model 13.** -- 4 States - IRES activation rates dependent and inactivation rates independent. <br/>
4S_Im2C --   **Model 14.** -- 4 States - IRES activation rates dependent and inactivation rates independent + Cross-over . <br/>
---

## Experimental data. <br/>

* **Experimental data is located in the following directory:** <br/>
Source_Code/Experimental_Data <br/>

* **Fraction of spots in Cap, IRES, CAP-IRES, and non-translating** <br/>
PercentPerCell_OriginalTag.xlsx <br/>
NumbersPerCell_OriginalTag.xlsx  <br/>

* **Intensity Distribution** <br/>
Cap_Only_mRNA_Threshold.xls <br/>
SwitchTag_IRESOnly_mRNA_Threshold.xls <br/>

* **Harringtonine data** <br/>
CapIntPerCell_HT.xlsx <br/>
IRESIntPerCell_HT.xlsx <br/>

* **NaAs data** <br/>
blueInPurpleTotalInt_NaAs.xlsx <br/>
blueInWhiteTotalInt_NaAs.xlsx <br/>
greenInWhiteTotalInt_NaAs.xlsx <br/>
greenInYellowTotalInt_NaAs.xlsx <br/>

* **DTT data** <br/>
blueInPurpleTotalInt_DTT.xlsx <br/>
blueInWhiteTotalInt_DTT.xlsx <br/>
greenInWhiteTotalInt_DTT.xlsx <br/>
greenInYellowTotalInt_DTT.xlsx <br/>

---

## Gene Sequences. <br/>
Sequence_CAP_KDM5B.txt <br/>
Sequence_IRES_Kif18b.txt <br/>

---  

## Code implementation.<br/>

To reproduce the results given in Figure 5 of the manuscript use: <br/>
**run__Desktop_Selected_Model.m** <br/>
All results are stored in: Source_Code/**Selected_Model_13** <br/>

To reproduce the results given in Supplementary Figure 6 of the manuscript use: <br/>
**run__Desktop_Compare_all_Models.m** <br/>
All results are stored in: Source_Code/**Model_Comparison** <br/>

---  

## Folders with simulation results. <br/>
Source_Code/**Best_Fits_For_Each_Model** .- contains 14 subfolders with the simulations for all the tested models. <br/>
Source_Code/**Selected_Model_13** .- contains the simulation results for the selected model (Model 13, 4S_Im2). <br/>
Source_Code/**Sup_Selected_Model_13** .- contains the simulation results for the selected model (Model 13, 4S_Im2) the figures have a special format for supplementary information. <br/>
Source_Code/**Complex_Model_4** .- contains the simulation results for the model with the largest number of free parameters (Model 4, 4S_DDC). <br/>
Source_Code/**Model_Comparison** .- contains subfolders with the Likelihood calculations for all 14 tested models.  <br/>
Source_Code/**ParameterUncertainty_13** .- contains the parameter uncertainty analysis for the selected model (Model 13, 4S_Im2). <br/>

---  

## Cluster implementation.<br/>
The codes performing the optimization, uncertainty, and model comparison (run__Cluster_Optmization.m and, run__Cluster_Parameter_Uncertainty.m, and run__Cluster_Compare_All_Models.m ) were implemented on the W.M. Keck Compute Cluster. <br/>

---  

## Files to submit jobs to the cluster.<br/>

cluster_submit_Model_Comparison.sh  <br/>
cluster_submit_Optimization.sh  <br/>
cluster_submit_Parameter_Uncertainty.sh   <br/>

To submit a job to the cluster, log to your account and use the following command:<br/>
**"qsub cluster_submit_Parameter_Uncertainty.sh"**

---  
