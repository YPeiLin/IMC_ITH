## Fine-scale Spatial heterogeneity of tumor microenvironment revealed by Imaging Mass Cytometry data![image](https://github.com/user-attachments/assets/dc667e95-a9de-463d-b03d-5173e4995c58)
 
# Author: Yupei lin
# File ITH-main: Main functions for ITH calculation 

## This code can be used to generate the ITH scores for each of cell subpopulation 
# of interest. It contains 5 major parts with 1-4 for ITH calculation and 5 for 
# data visualization
# ' @ myinf1 - detailed cell coordinates and cell type
# ' @ myoutf2 - output RDS files to store ITH 
# ' @ immune_cells - List of immune cell types
# ' @ nn number of subsampling to be iterated
#' 
# ' Example usage: df=calculate_ITH(nn,myList) 
#' 
