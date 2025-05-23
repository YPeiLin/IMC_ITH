# [2] generate the position information of all cells with original image size
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.1] 

rm(list=ls())
myDir1 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/AAA_OtherData/Sorin2023_PMID36725934_IMC/LUAD_IMC_CellType/"
myDir2 = "/mount/amos1/home/ylin/IMC/output/position"
myoutf1 = "/mount/amos1/home/ylin/IMC/output/IMC_Cell_position_updated.rda"

library(R.matlab)

#Reading each patient position file
myinf1 = dir(myDir1)
setwd(myDir1)
myList = list(NULL)
mynam = rep(0, length(myinf1))
xx = gsub("\\.mat", "", myinf1)
mytag = gsub("LUAD_", "", xx)

myinf2 = dir(myDir2)

for(ff in 1:length(myinf1))
{
  cat("\r", ff)
  setwd(myDir1)
  data <- readMat(myinf1[ff])
  #cell types 
  xx = data[[2]][,1]
  
  flag = rep("", length(xx))
  for(k in 1:length(xx))
  {
    flag[k] = ifelse(nrow(xx[[k]][[1]])>0, as.character(xx[[k]][[1]]), "")
  }
  setwd(myDir2)
  pos = read.table(myinf2[ff],sep = ",",header = TRUE)
  names(pos) = c("posx", "posy", "dist.x", "dist.y") 
  cel.typ = NULL
  index=which(flag=="")
  
  for(k in 1:nrow(pos))
  { 
    if(flag[k]=="") next
    cel.typ = c(cel.typ, flag[k])
  }
  pos = pos[-index,]
  tmp = data.frame(pos$posx, pos$posy, pos$dist.x, pos$dist.y, cel.typ)
  names(tmp) = c("posx","posy","dist.x","dist.y","cel.typ")
  myList[[ff]]=tmp 
}
names(myList) = mytag

save(myList, file=myoutf1)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.2] survival analysis
rm(list=ls())
# myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/AAA_OtherData/Sorin2023_PMID36725934_IMC/done/IMC_Cell_position.rda"
myDir1 = "/mount/amos1/home/ylin/IMC/output/IMC_Cell_position_updated.rda"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/AAA_OtherData/Sorin2023_PMID36725934_IMC/LUAD_Clinical_416_Discovery.txt" 

info = read.table(myinf2, sep="\t", header=T, row.names=1)
row.names(info) = gsub("LUAD_", "", row.names(info))
time = info$Survival.time
event = info$Death
info = cbind(time, event, info)

load(myDir1)
for(k in 1:length(myList))
{
  cat("\r", k)
  data = myList[[k]]
  tmp = table(data$cel.typ)
  if(k==1)
  {
    res = matrix(0, length(myList), length(tmp))
    row.names(res) = names(myList)
    colnames(res) = names(tmp)
    res = as.data.frame(res)
    res[k,] = tmp[colnames(res)]
    next
  }
  res[k,] = tmp[colnames(res)]
}

for(k in 1:ncol(res))
{
  tmp = res[,k]
  tmp[is.na(tmp)] =0
  res[,k] = tmp
}
xx = apply(res, 1, sum)
data = res/xx
## cell fractions
# data
#---------------------------------
comxx = intersect(row.names(info), row.names(data))
data = data[comxx,]
info = info[comxx,]
dim(data)

library(survival)

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  xx = xx[xx[, "time"]>0,]
  mycox = survreg(Surv(time, event)~mytf, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["mytf", "p"]
  mycox = coxph(Surv(time, event)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  lb1[k] = tmp[3]
  ub1[k] = tmp[4]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = colnames(data)
res = data.frame(name, survreg.pval1,  coxph.pval1,  HR= hr1)
xx = res[order(res[,3]), ]
xx


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.3] survival analysis  -- Ratios
rm(list=ls())
myinf1 = "/mount/amos1/home/ylin/IMC/output/IMC_Cell_position_updated.rda"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/AAA_OtherData/Sorin2023_PMID36725934_IMC/LUAD_Clinical_416_Discovery.txt" 


info = read.table(myinf2, sep="\t", header=T, row.names=1)
row.names(info) = gsub("LUAD_", "", row.names(info))
time = info$Survival.time
event = info$Death
info = cbind(time, event, info)

load(myinf1)

for(k in 1:length(myList))
{
  cat("\r", k)
  data = myList[[k]]
  tmp = table(data$cel.typ)
  if(k==1)
  {
    res = matrix(0, length(myList), length(tmp))
    row.names(res) = names(myList)
    colnames(res) = names(tmp)
    res = as.data.frame(res)
    res[k,] = tmp[colnames(res)]
    # next
  }
  res[k,] = tmp[colnames(res)]
}
xx = apply(res, 2, mean, na.rm=T)
xx
cel.typ = gsub(" ", "", names(xx))

#---------------------------------
nn = ncol(res)*(ncol(res)-1)
myres = matrix(0, nrow(res), nn)
row.names(myres) = row.names(res)
myres = as.data.frame(myres)
mynam = NULL
count = 0

for(i in 1:ncol(res))
  for(j in 1:ncol(res))
  {
    if(i==j) next
    count = count + 1
    mynam = c(mynam, paste(cel.typ[i], cel.typ[j], sep="__"))
    myres[, count] = res[,i]/res[,j]
  }
colnames(myres) = mynam

data = myres
comxx = intersect(row.names(info), row.names(data))
data = data[comxx,]
info = info[comxx,]
dim(data)

library(survival)

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  xx = xx[xx[, "time"]>0,]
  mycox = survreg(Surv(time, event)~mytf, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["mytf", "p"]
  mycox = coxph(Surv(time, event)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  lb1[k] = tmp[3]
  ub1[k] = tmp[4]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = colnames(data)
res = data.frame(name, survreg.pval1,  coxph.pval1,  HR= hr1)
xx = res[order(res[,3]), ]
xx
#-----------------------------------------------------------------------------

# [7.1] survival analysis
rm(list=ls())
myinf1= "/mount/amos1/home/ylin/IMC/output/100-IMC_sample_sITH_S9N199.rda"
# myinf1= "/mount/amos1/home/ylin/IMC/output/100-IMC_sample_sITH_adaptive.rda"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/AAA_OtherData/Sorin2023_PMID36725934_IMC/LUAD_Clinical_416_Discovery.txt" 

info = read.table(myinf2, sep="\t", header=T, row.names=1)
row.names(info) = gsub("LUAD_", "", row.names(info))
time = info$Survival.time
event = info$Death
#event = info$Progression
info = cbind(time, event, info)
raw.info = info

load(myinf1)

for(k in 1:length(myList))
{
  xx = myList[[k]]
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
  myv2 = xx$sd1.sr/xx$sd1.rd
  # myv2 = xx$sd2.sr/xx$sd2.rd
  data[k,] = c(myv1, myv2)
}

