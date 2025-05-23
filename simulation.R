rm(list=ls())
# myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/AAA_OtherData/Sorin2023_PMID36725934_IMC/done/IMC_Cell_position.rda"
# myDir1 = "/mount/amos1/home/ylin/IMC/output/IMC_Cell_position_updated.rda"
# myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/AAA_OtherData/Sorin2023_PMID36725934_IMC/LUAD_Clinical_416_Discovery.txt" 
# myinf3 = "/mount/amos1/home/ylin/IMC/output/IMC_patient-imagesize.txt"
myoutf1 = "/mount/amos1/home/ylin/IMC/output/simulated-5patient-noITH-1600.rds"

# 
# info = read.table(myinf2, sep="\t", header=T, row.names=1)
# row.names(info) = gsub("LUAD_", "", row.names(info))
# time = info$Survival.time
# event = info$Death
# info = cbind(time, event, info)



#Simulate 5 patients with no s-ITH, number of cells: 1600,3200,4800,6400,8000
#---------------------------------------------------------------------------------
set.seed(51)
iter <- 1
celtypelen =16
# x.all = y.all= 1000
# 
# #how to define patch
# patch.size.x = patch.size.y= x.all/10

cel.type.nam =c("Alt MAC" ,   "B cell"  ,     "Cancer" ,   "Cl MAC"  ,  
                "Cl Mo"   ,  "DCs cell" ,    "Endothelial cell" ,"Int Mo",   
                "Mast cell","Neutrophils","NK cell","Non-Cl Mo","T other",
                "Tc" ,"Th","Treg" )
myList=NULL
k=1600
x.cor.t= c(1:k)
y.cor.t= c(1:k)

for(k in 1:iter){

  myRes = NULL
  for(j in 1:length(x.cor.t)){
    for (i in 1:length(y.cor.t)){

    x.cor= i
    
    y.cor= j
    print(x.cor)
    print(y.cor)
    tmp = as.data.frame( cbind(x.cor,y.cor, "cancer"))
    myRes = as.data.frame(rbind(myRes,tmp))
    
    }
    cat(paste0(j,"-th row simulated \n\n"))
  }
colnames(myRes) = c("posx","posy","cel.type")
myRes$cel.type = rep(cel.type.nam,160000)
myList[[k]]=myRes
}
   
saveRDS(myList, myoutf1)
#-----
# original simulation method
# patch.size.x = patch.size.y= x.all/10
# 
# cel.type.nam =c("Alt MAC" ,   "B cell"  ,     "Cancer" ,   "Cl MAC"  ,  
#                 "Cl Mo"   ,  "DCs cell" ,    "Endothelial cell" ,"Int Mo",   
#                 "Mast cell","Neutrophils","NK cell","Non-Cl Mo","T other",
#                 "Tc" ,"Th","Treg" )
# myList=NULL
# k=10
# x.cor= sample((1:1000),celtypelen*k)
# y.cor= sample((1:1000),celtypelen*k)
# tmp = as.data.frame( cbind(x.cor,y.cor, cel.type.nam))
# myRes = as.data.frame(rbind(myRes,tmp))
# 
# for(k in 1:iter){
#   
#   myRes = NULL
#   for(j in 1:10){
#     for (i in 1:10){
#       start.row =  1+patch.size.x*(i-1)
#       start.col =  1+patch.size.y*(j-1)
#       end.row = start.row+patch.size.x
#       end.col = start.col+patch.size.y
#       print(paste0("start: [",start.row,",",start.col,"]"))
#       print(paste0("end: [",end.row,",",end.col,"]"))
#       cat("\n")
#       x.cor= sample((start.row:end.row),celtypelen*k)
#       y.cor= sample((start.col:end.col),celtypelen*k)
#       tmp = as.data.frame( cbind(x.cor,y.cor, cel.type.nam))
#       myRes = as.data.frame(rbind(myRes,tmp))
#       
#     }
#     cat(paste0(j,"-th row simulated \n\n"))
#   }
#   colnames(myRes) = c("posx","posy","cel.type")
#   myList[[k]]=myRes
# }


# myinf1= "/mount/amos1/home/ylin/IMC/output/100-IMC_sample_sITH_S9N199.rda"
#-------------------
# rm(list=ls())
nn = 100

myinf1 = "/mount/amos1/home/ylin/IMC/output/IMC_Cell_position_updated.rda"
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
    cel.typ = names(tmp)
  }
  xx = tmp[cel.typ]
  xx[is.na(xx)]=0
  xx = xx/sum(xx)
  obs = as.vector(xx)
  if(k==1)
  {
    mymat = matrix(0, length(cel.typ), nn)
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
#cf.list containting 100 randomized patches with 16 cell type fractions

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
  IMM1 = apply(data[-c(3,7), ], 2, sum) #Tumor infiltrating cells: Remove cancer and endothelial
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
  IMM1 = apply(data[-c(3,7), ], 2, sum)
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


#calculate final ITH
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
tmp = NULL
for(i in c(1:100)){
df=calculate_ITH(100,myList)
tmp= rbind(tmp, df[-c(1:17,34)])
}
tmp=round(tmp,digits = 3)

tmp =read.csv("/mount/amos1/home/ylin/IMC/output/simulation-ith_v3_3patient.txt")
mean1=apply(tmp[11:20,-1],2, mean)
mean1
tmp = myList[[1]]
df=tmp
# se= which(tmp$cel.typ=="Endothelial cell")
# df= tmp[se,]
names(df) = c("x","y","name")
df$name = as.factor(df$name)
p=ggplot(data=df, aes(x=x, y=y, label=name,col=name)) + 
  geom_point(size=2) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.background = element_rect(color=NA,fill=NA),axis.text.x=element_blank(),axis.text.y=element_blank())+labs(title = "Spatial image")

print(p)
# Compare ITH score result
# load(myinf1)
# df=calculate_ITH(100,myList)
# tmp1=df[1,-c(1:17)]
# tmp=rbind(tmp,tmp1)
# tmp1=round(tmp,digits=2)
# row.names(tmp1)=c(1:51)
# write.table(tmp1,"/mount/amos1/home/ylin/IMC/output/50-calculate_ith.csv",sep = ",",row.names = FALSE)

tmp=apply(df,2,sd)
tmp=as.data.frame(tmp)

######
rm(list=ls())
# Create a data frame with all possible positions
positions <- expand.grid(1:1000, 1:1000)
num_positions_to_select <- 100
  
cel.type.nam =c("Alt MAC" ,   "B cell"  ,     "Cancer" ,   "Cl MAC"  ,  
                  "Cl Mo"   ,  "DCs cell" ,    "Endothelial cell" ,"Int Mo",   
                  "Mast cell","Neutrophils","NK cell","Non-Cl Mo","T other",
                  "Tc" ,"Th","Treg" )
# Create a function to randomly select a position without repetition
selectRandomPosition <- function(positions) {
  
  if (nrow(positions) == 0) {
    return(NULL)  # No more positions available
  }
  
  # Randomly select one position
  selected_index <- sample(nrow(positions), 1)
  selected_position <- positions[selected_index, ]
  
  # Remove the selected position from the list
  positions <- positions[-selected_index, ]
  
  return(list(selected_position = selected_position, remaining_positions = positions))
}

# Select a specific number of positions

myRes = NULL
for( k in 1: length(cel.type.nam) ){
  cat("Randomizing ", cel.type.nam[k],"\n")

  for (i in 1:num_positions_to_select) {
    result <- selectRandomPosition(positions)
    
    if (is.null(result)) {
      cat("No more positions available.\n")
      break
    }
    
    selected_position <- result$selected_position
    positions <- result$remaining_positions
    
    cat("Selected Position:", selected_position[[1]], ", ", selected_position[[2]], "\n")
  
    x.cor= selected_position[[1]]
    y.cor= selected_position[[2]]
    tmp = as.data.frame( cbind(x.cor,y.cor, cel.type.nam[k]))
    myRes = as.data.frame(rbind(myRes,tmp))
    
  }
  
} 
which(duplicated(myRes[,c(1:2)])) #verify it's not duplicated location


colnames(myRes) = c("posx","posy","cel.type")
myList[[3]]=myRes

saveRDS(myList, "/mount/amos1/home/ylin/IMC/output/simulated-3patient-noITH-v3.rds")
