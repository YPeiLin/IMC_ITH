#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [1] Sample-specific CCI -- KNN method
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.1] calcualte the z-score for cell-cell interactions for all pairs of cell types in each IMC image
rm(list=ls())
myinf1 = "/mount/amos1/home/ylin/IMC/output/IMC_Cell_position_updated.rda"

myoutf1 = "/mount/amos1/home/ylin/IMC/output/IMM-specific-CCI-KNN5_Zscore.rda"

cel.typ = c("Cancer", "Tc", "Endothelial cell", "Cl MAC", "B cell", "Th", "Alt MAC", "Treg", "T other", "Neutrophils", "Cl Mo", "Mast cell", "Non-Cl Mo", "Int Mo", "NK cell", "DCs cell")

load(myinf1)
CCI.List = list(NULL)
knn = 10

for(k in 1:length(myList))
{
  cat("\r", k)
  k=1
  data = myList[[k]]
  dis = dist(data[,1:2], method="euclidean", upper = TRUE)
  dis = as.matrix(dis)
  for(i in 1:ncol(dis))
  {
    vec = dis[,i]
    vec[i] = Inf
    tmp = rank(vec, ties="random")
    xx = rep(0, length(tmp))
    xx[tmp<=knn] = 1
    dis[,i] = xx
  }
  
  #--------	
  obs.knn = matrix(0, length(cel.typ), length(cel.typ))
  row.names(obs.knn) = colnames(obs.knn) = cel.typ
  for(i in 1:length(cel.typ))
  {
    # cancer
    se = which(data$cel.typ%in%cel.typ[i])
    if(length(se)==0)
    {
      obs.knn[i,] = 0
      next
    }
    tmp = dis[,se, drop=F]
    for(j in 1:length(cel.typ))
    {
      se = which(data$cel.typ%in%cel.typ[j])
      obs.knn[i,j] = sum(tmp[se,])
    }
  }
  
  #--------	
  pmu.knn1 = pmu.knn2 = matrix(0, length(cel.typ), length(cel.typ))
  pnn = 100
  for(p in 1:pnn)
  {
    pm.cel.typ = sample(data$cel.typ)
    for(i in 1:length(cel.typ))
    {
      se = which(pm.cel.typ%in%cel.typ[i])
      if(length(se)==0)
      {
        pmu.knn1[i,] = pmu.knn1[i,] + 0
        pmu.knn2[i,] = pmu.knn2[i,] + 0
        next
      }
      tmp = dis[,se, drop=F]
      for(j in 1:length(cel.typ))
      {
        se = which(pm.cel.typ%in%cel.typ[j])
        xx = sum(tmp[se,])
        pmu.knn1[i,j] = pmu.knn1[i,j] + xx
        pmu.knn2[i,j] = pmu.knn2[i,j] + xx^2
      }
    }
    pmu.avg = pmu.knn1/pnn
    pmu.var = (pmu.knn2 - pmu.knn1^2/pnn)/(pnn-1)
    pmu.std = sqrt(pmu.var)
  }
  zscore = (obs.knn-pmu.avg)/pmu.std
  CCI.List[[k]] = zscore
}
names(CCI.List) = names(myList)

save(CCI.List, file = myoutf1)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.2] correlated the distance with prognosis
rm(list=ls())
myinf1 = "/mount/amos1/home/ylin/IMC/output/IMM-specific-CCI-KNN5_Zscore.rda"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/AAA_OtherData/Sorin2023_PMID36725934_IMC/LUAD_Clinical_416_Discovery.txt" 


info = read.table(myinf2, sep="\t", header=T, row.names=1)
row.names(info) = gsub("LUAD_", "", row.names(info))
time = info$Survival.time
event = info$Death
info = cbind(time, event, info)
raw.info = info

cel.typ = c("Cancer", "Tc", "Endothelial cell", "Cl MAC", "B cell", "Th", "Alt MAC", "Treg", "T other", "Neutrophils", "Cl Mo", "Mast cell", "Non-Cl Mo", "Int Mo", "NK cell", "DCs cell")

load(myinf1)
myList = list(NULL)
for(k in 1:length(cel.typ))
{
  cat("\r",k)
  res = matrix(0, length(CCI.List), length(cel.typ))
  row.names(res) = names(CCI.List)
  colnames(res) = cel.typ
  res = as.data.frame(res)
  for(i in 1:length(CCI.List))
  {
    xx = CCI.List[[i]]
    res[i,] = as.numeric(xx[k,])
  }
  myList[[k]] = res
}
names(myList) = cel.typ

library(survival)

tmp = matrix(0, length(myList), length(cel.typ))
row.names(tmp) = names(myList)
colnames(tmp) = cel.typ
hr.mat = pv.mat = nn.mat = tmp

for(s in 1:length(myList))
{
  cat("\r", s)
  data = myList[[s]]
  info = raw.info
  comxx = intersect(row.names(info), row.names(data))
  data = data[comxx,]
  info = info[comxx,]
  
  for(k in 1:ncol(data))
  {
    mytf = as.numeric(data[,k])
    se = which(abs(mytf)=="Inf")
    mytf[se] = NA
    se = which(!is.na(mytf))
    if(length(se)<nrow(data)*0.2)
    {
      hr.mat[s,k] = pv.mat[s,k] = nn.mat[s,k] = NA
      next
    }
    nn.mat[s,k] = length(se)
    
    xx = cbind(mytf, info)
    xx = xx[xx[, "time"]>0,]
    mycox = coxph(Surv(time, event)~mytf, xx) 
    mycox = summary(mycox)
    pv.mat[s,k] = mycox$coefficients[5]
    tmp = mycox$conf.int
    hr.mat[s,k] = tmp[1]
  }
}
xx = round(-log10(pv.mat),1)
xx = xx * (I(hr.mat>1)-0.5)*2

Cancer   Tc Endothelial cell Cl MAC B cell   Th Alt MAC Treg
Cancer              0.3  0.4             -0.7   -0.6    2.1  0.3    -1.1 -1.2
Tc                  0.7 -0.6              0.5    0.0    0.2 -0.7     2.5  0.9
Endothelial cell   -0.3  0.4              0.8   -0.4    1.1  0.7     2.0  0.4
Cl MAC             -0.4  0.0             -0.3    0.2    0.2  0.4    -0.7  0.9
B cell              2.1  0.2              1.3    0.1   -1.3 -0.3     0.4  1.2
Th                  0.7 -0.7              1.1    0.5   -0.3 -0.6     1.5  1.1
Alt MAC            -1.2  2.6              2.1   -1.0    0.2  1.6    -0.3  1.0
Treg               -0.7  0.8              0.6    0.6    1.2  0.7     1.2  0.6
T other            -0.9  0.8              1.0    0.0    0.9  0.6     1.5  1.9
Neutrophils        -0.8 -1.3             -2.1   -0.5    0.9 -1.3    -1.3 -3.5
Cl Mo              -0.6  0.1             -0.2   -0.2    1.1  1.2     0.1  0.7
Mast cell           0.0 -0.2             -0.1   -0.7    2.1  0.5    -0.1 -0.1
Non-Cl Mo          -2.3  2.3              1.1   -0.7    0.4  1.2     0.3  1.3
Int Mo             -1.5  0.5              0.2   -1.0    0.1  0.6     1.1  0.2
NK cell             0.2 -0.1              0.6    0.2    1.5  1.0    -1.6  0.7
DCs cell            0.0  0.0              1.1    0.1    1.1  0.7    -1.0  0.1
T other Neutrophils Cl Mo Mast cell Non-Cl Mo Int Mo NK cell
Cancer              -1.3        -1.0  -1.3      -0.1      -2.8   -1.9     0.0
Tc                   0.4        -1.0   0.1      -0.1       1.6    0.3    -0.1
Endothelial cell     0.7        -2.3  -0.2       0.0       1.1    0.0     0.9
Cl MAC               0.0        -0.2  -0.1      -0.7      -0.8   -0.7     0.6
B cell               0.9         0.7   0.9       2.2       0.6    0.0     0.7
Th                   0.7        -1.1   1.0       0.5       0.9    0.4     1.0
Alt MAC              2.1        -1.6   0.1      -0.2       0.4    1.1    -1.0
Treg                 2.5        -3.5   0.3       0.0       0.8    0.4     0.9
T other              0.1        -1.5   0.3       0.1       0.1    0.2    -0.3
Neutrophils         -1.6         1.4  -0.6      -0.8      -0.5   -1.0    -0.5
Cl Mo                0.3        -0.5   0.2       0.0       0.7    1.6     0.8
Mast cell            0.3        -0.9  -0.1      -0.1       0.2   -0.6     0.1
Non-Cl Mo            0.1        -0.4   0.4       0.1       1.3    0.7     1.2
Int Mo               0.0        -1.3   1.3      -0.4       0.7    0.8     1.2
NK cell              0.0        -0.8   0.5       0.0       1.4    0.7    -1.6
DCs cell            -0.6         0.0  -0.1      -0.4       1.1    0.8    -1.0
DCs cell
Cancer               -0.1
Tc                    0.0
Endothelial cell      1.7
Cl MAC                0.1
B cell                2.1
Th                    0.7
Alt MAC              -0.9
Treg                  0.1
T other              -0.3
Neutrophils          -0.3
Cl Mo                 0.0
Mast cell            -0.5
Non-Cl Mo             1.3
Int Mo                0.6
NK cell              -0.6
DCs cell             -0.3

xx["Alt MAC", ]

Cancer               Tc Endothelial cell           Cl MAC 
-1.2              2.6              2.1             -1.0 
B cell               Th          Alt MAC             Treg 
0.2              1.6             -0.3              1.0 
T other      Neutrophils            Cl Mo        Mast cell 
2.1             -1.6              0.1             -0.2 
Non-Cl Mo           Int Mo          NK cell         DCs cell 
0.4              1.1             -1.0             -0.9 

