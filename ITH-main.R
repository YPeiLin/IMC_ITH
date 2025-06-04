# Title: Main functions for ITH calculation 
## Author: Yupei lin
##This code can be used to generate the ITH scores for each of cell subpopulation 
# of interest. It contains 5 major parts with 1-4 for ITH calculation and 5 for 
# data visualization
#' @ myinf1 - detailed cell coordinates and cell type
#' @ myoutf2 - output RDS files to store ITH 
#' @ immune_cells - List of immune cell types
#' @ nn number of subsampling to be iterated
#' 
#' Example usage: df=calculate_ITH(nn,myList) 
#' 
#--Part1: Calculate cell fractions for spatial subsampling
#--Part2: Calculate cell fractions for random subsampling
#--Part3: Calculate Variation metrics
#--Part4: Summarize final ITH
#--Part5: Data visualization for boxplots, violin plots and heatmaps

myinf1="cell_coordinate.rda"
myoutf2="path/to/output.rda"
immune_cells <- c(
  "CD4 T cells",
  "CD8 T cells",
  "Double-positive T cells",
)


##Clear working environment
rm(list=ls())
load(myinf1)

#Iterating all unique labels, stored in @alllab
tmp=NULL
for(i in 1:length(myList)){
  tmp = append(tmp,unique(myList[[i]]$cel.typ))
}
alllab=unique(tmp)
alllab

# Find indices of immune cell types
imm <- which(alllab %in% immune_cells)


#####Define a function to iteratively calculate ITH 
#@@ input-
# nn:patch size
# myList: patient with positional info
#@@ output-
# data: ITH score dataframe


calculate_ITH <- function(nn, myList) {
  
  #--Part1: calculate cell fractions for spatial subsampling
  cf.List = list(NULL)
  cat("Calculating Spatial-subsampling cell fractions \n")
  for(k in 1:length(myList))
  {
    # cat(k,'\n')
    data = myList[[k]]
    xmax = as.numeric(max(data$posx))
    ymax = as.numeric(max(data$posy))
    xwin = round(xmax/3)
    ywin = round(ymax/3)
    xlim = xmax-xwin
    ylim = ymax-ywin
    tmp = table(data$cel.typ)
    if(k==1)
    {
      # cel.typ = names(tmp)
      cel.typ = alllab
    }
    xx = tmp[cel.typ]
    xx[is.na(xx)]=0
    xx = xx/sum(xx)
    obs = as.vector(xx)
    if(k==1)
    {
      # mymat = matrix(0, length(cel.typ), nn)
      mymat = matrix(0, length(alllab), nn)
      row.names(mymat) = cel.typ
      colnames(mymat) = paste("R", 1:nn, sep="")
      mymat = as.data.frame(mymat)
    }	
    for(i in 1:nn)
    {
      randx = sample(1:xlim, 1)
      randy = sample(1:ylim, 1)
      se = which((data$posx>=randx) & (data$posx<(randx+xwin)) & (data$posy>=randy) & (data$posy<(randy+ywin)))
      sdat = data[se,]
      tmp = table(sdat$cel.typ)
      xx = tmp[cel.typ]
      xx[is.na(xx)]=0
      xx = xx/sum(xx)
      mymat[,i] = xx	
    }
    res = cbind(obs, mymat)
    cf.List[[k]] = res
  }
  names(cf.List) = names(myList)
  cat("Finishing spatial-subsampling \n")
  #cf.list containting nn randomized patches with immune cell type fractions
  
  #--Part2: calculate cell fractions for random subsampling
  rd.List = list(NULL)
  for(k in 1:length(myList))
  {
    data = myList[[k]]
    if(k==1)
    {
      mymat = matrix(0, length(cel.typ), nn)
      row.names(mymat) = cel.typ
      colnames(mymat) = paste("R", 1:nn, sep="")
      mymat = as.data.frame(mymat)
    }
    sub.nn = nrow(data)/9
    
    for(i in 1:nn)
    {
      se = sample(1:nrow(data), sub.nn)
      sdat = data[se,]
      tmp = table(sdat$cel.typ)
      xx = tmp[cel.typ]
      xx[is.na(xx)]=0
      xx = xx/sum(xx)
      mymat[,i] = xx	
    }
    res = mymat
    rd.List[[k]] = res
  }
  names(rd.List) = names(myList)
  cat("Finish random-subsampling \n")
  
  
  #Part3- Calculate Variation metrics
  #For each patient and each cell type, find their obs, mean s-sub, sd1&2 s-sub, mean r-sub, sd1&2 r-sub
  res.List = list(NULL)
  for(k in 1:length(cf.List))
  {
    data = as.matrix(cf.List[[k]])
    
    IMM1 = apply(data[imm, ], 2, sum) #TODO
    data = rbind(data, IMM1)
    #spatial subsampling
    res = matrix(0, nrow(data), 7)		## why +2: (IMM1  & IMM2) 
    row.names(res) = row.names(data)
    colnames(res) = c("obs", "avg.sr","sd1.sr", "sd2.sr", "avg.rd","sd1.rd", "sd2.rd")
    res = as.data.frame(res)
    res[,1] = data[,1]
    data = data[, -1]
    res[,2] = apply(data, 1, mean, na.rm=T)
    res[,3] = apply(data, 1, sd, na.rm=T)
    
    xx = data-res[,1]#(cf - obs)
    xx = apply(xx^2, 1, sum, na.rm=T)/(nn-1)
    res[,4] = sqrt(xx)
    
    #random subsampling
    data = as.matrix(rd.List[[k]])
    IMM1 = apply(data[imm, ], 2, sum) #TODO
    data = rbind(data, IMM1)
    res[,5] = apply(data, 1, mean, na.rm=T)
    res[,6] = apply(data, 1, sd, na.rm=T)
    xx = data-res[,1]
    xx = apply(xx^2, 1, sum, na.rm=T)/(nn-1)
    res[,7] = sqrt(xx)
    
    myvec = c(0, ncol(res))
    myvec[1:2] = as.numeric(res["IMM1", 1:2]) # obs    avg.sr
    myvec[3] = sqrt(sum((res$sd1.sr[-nrow(data)])^2))
    myvec[4] = sqrt(sum((res$sd2.sr[-nrow(data)])^2))
    myvec[5] = res["IMM1", 5] #avg.rd
    myvec[6] = sqrt(sum((res$sd1.rd[-nrow(data)])^2))
    myvec[7] = sqrt(sum((res$sd2.rd[-nrow(data)])^2))
    IMM2 = myvec
    res = rbind(res, IMM2)
    row.names(res)[nrow(res)] = "IMM2"
    res.List[[k]] = res
  }
  names(res.List) = names(myList)
  
  
  #--Part4: Summarize final ITH
  cat("Calculating ITH \n")
  for(k in 1:length(res.List))
  {
    xx = res.List[[k]]
    if(k==1)
    {
      data = matrix(0, length(myList), nrow(xx)*2-1)
      row.names(data) = names(myList)
      tmp1 = row.names(xx)[1:(nrow(xx)-2)]
      tmp2 = paste(row.names(xx), "_ITH", sep="")
      colnames(data) = c(tmp1, "IMM", tmp2)
      data = as.data.frame(data)
    }
    myv1 = xx$obs
    myv1 = myv1[-length(myv1)]
    myv2 = xx$sd2.sr/xx$sd2.rd # myv2 = xx$sd1.sr/xx$sd1.rd
    data[k,] = c(myv1, myv2)
  }
  
  return(data)
}

#Call main function
df=calculate_ITH(nn,myList)
#Save output
save(df, file=myoutf2)


#####Function for Cell number Quality control
#@@ input-
# myList: patient with positional info
#@@ output-
# avg_spatial_cells: Average number of cells 

load(myoutf2)
data=ITH_df = df 
avg_spatial_cells <- mean(sapply(myList, function(d) {
  xmax = max(d$posx)
  ymax = max(d$posy)
  xwin = round(xmax / 3)
  ywin = round(ymax / 3)
  window_area = xwin * ywin
  total_area = xmax * ymax
  window_area / total_area * nrow(d)
})) 


#####Result visualization
#@@ input-
# data: ITH matrix from function 'calculate_ITH()'
#@@ output-
# Heatmaps and boxplots

#Libraries loading
library(pheatmap)
library(reshape2)
library(ggplot2)
raw.data=data
se = grep("_ITH", colnames(raw.data)) # Retrive only ITH columns
data = raw.data[, se]
tmp <- data[, imm]

#Color palette
pal <- c(
  "#D8811A", "#75A865", "darkorange","darkred", "#CAB08D", "#E7CA8E",
  "#A1A1C5", "#583201", "#ECE4B3", "#4D7696", "#D46135", 
  "#EBAFA4", "#5E8A89", "#ECD577", "#795F6A", "#5C6489","darkblue", "#578F12",
  "#E97F71", "#4290B7")

dat <- data.frame(group = rep(colnames(tmp), each = nrow(tmp)),
                  score = as.vector(as.matrix(tmp)))
g_order <- names(sort(apply(tmp, 2, median, na.rm = TRUE), decreasing = FALSE))
g_order= gsub("\\s*\\([^\\)]+\\)", "",g_order)
dat$group <- gsub("\\s*\\([^\\)]+\\)", "", dat$group)
dat$group <- factor(dat$group, levels = g_order)

#--Boxplot
g <- ggplot(dat, aes(x = group, y = score, fill = group)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  scale_fill_manual(values = pal) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(size = 0.5, fill = NA),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = 1, hjust = 1),
    title = element_text(face = "bold", size = 15),
    plot.title = element_text(hjust = 0.5),
    aspect.ratio = 3/4,
    legend.background = element_rect(color = NA, fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.position = "none"
  ) +
  labs(
    x = "Cell Subpopulations",
    y = "ITH Scores",
    title = "ITH Across Immune Cell Types"
  ) +
  ylim(0, 5)
print(g)

#--Violin plots
dat$group
dat$patient <- rownames(raw.data)  
ggplot(dat, aes(x = group, y = score, fill = patient)) +
  geom_violin(trim = FALSE, scale = "width", color = "black") +
  theme_minimal() +
  labs(x = "Cell Subpopulations", y = "ITH Score", title = "ITH per Cell Type by Patient") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
    legend.position = "right"
  )


#--Heatmap
myx <- tmp
cor_matrix <- cor(myx, use = "pairwise.complete.obs")
diag(cor_matrix) <- NA
cor_matrix[cor_matrix > 1] <- 1
cor_matrix[cor_matrix < -1] <- -1
colnames(cor_matrix) = gsub("\\s\\([^\\)]+\\)","",colnames(cor_matrix))
row.names(cor_matrix) = gsub("\\s\\([^\\)]+\\)","",row.names(cor_matrix))
pheatmap(
  mat = cor_matrix,
  color = colorRampPalette(c("lightyellow", "#D73027"))(50),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = FALSE,
  fontsize = 12,
  main = "Correlation of ITH Across Cell Types"
)

