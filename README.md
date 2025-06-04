## Fine-scale Spatial heterogeneity of tumor microenvironment revealed by Imaging Mass Cytometry data. 
![image](https://github.com/user-attachments/assets/83c428b8-9b04-42ac-88e6-7c4d64392597)
 
Author: Yupei lin
### File ITH-main: Main functions for ITH calculation 
This code can be used to generate the ITH scores for each of cell subpopulation of interest. It contains 5 major parts with 1-4 for ITH calculation and 5 for data visualization
Example usage: df=calculate_ITH(nn,myList) 

### File ITH-main: IMC_position.R
This code can be used to visualize cell distribution on raw IMC image with each color represents a unique cell subtype 

### File ITH-classifier
This code can be used do 5-CV on combining all ITH as one classfier

### File simulation.R
Simulate patients with no s-ITH to test pipeline robustness, number of cells tested (1600,3200,4800,6400,8000)
