# 23_Geuking_thymus_SC
Code for downstream analysis of Geuking lab SC thymus data.
Preprocessing.R is a QC script on SC level which is run first for each sample, resulting objects are saved separately for future use.
IntegrationAndID.R relies heavily on the Seurat library. The script integrates the two samples and identify all cell populations at low resolution clustering by comparing DEGs to the ImmGen data base and to the Kernfeld et al. (10.1016/j.immuni.2018.04.015) publication. 
