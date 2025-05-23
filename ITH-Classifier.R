# Title: Classifier pipeline; 
## Author: Yupei lin
rm(list=ls())
##Trainingtestingsplit
library(caret)
library(glmnet)
library(xgboost)
library(cvTools)
library(pROC)
library(gplots)
library(ggpubr)
library(patchwork)
#############################[1-4]Training#############################

#####[3]5-fold cross validation####
myinf5 = "/mount/amos1/home/ylin/IMC/output/CAF_IMC_sample_sITH.rda"
myinf6="/mount/amos1/home/ylin/IMC/output/lung_meta.rda"
load(myinf5)
load(myinf6)
info <- unique(data.frame(
  patientID = cell_info$patientID,
  event=cell_info$event,
  OS=cell_info$OS,
  DFS=cell_info$DFS,
  # eversmoke =   as.factor(cell_info$Smok),
  grade =   as.factor(cell_info$grade),
  stage =   as.factor(cell_info$stage),
  gender =   as.factor(cell_info$gender),
  age =   cell_info$age,
  type = cell_info$Tumortype,
  smok = cell_info$smok
  ))
se = grep("_ITH", colnames(df))
selected <- df[,se]
selected<- selected[,-c(32,33)]

tmp=t(selected)
se=which(colnames(tmp)%in%info$patientID)
tmp=tmp[,se]
df = as.data.frame(t(tmp))

#----5-CV model initialization
Clinicaldata <- data.frame(matrix(ncol = 0, nrow = nrow(df)))
Clinicaldata$finalgenematrix=df[, colnames(selected)]
Clinicaldata$finalgenematrix <- apply(Clinicaldata$finalgenematrix, 2, as.numeric)
Clinicaldata$eventclass=info$type
Clinicaldata$smok=info$smok
# Clinicaldata$eventclass=info$event
se=which(is.na(Clinicaldata$eventclass))
Clinicaldata=Clinicaldata[-se,] #remove no vital status
se=which((Clinicaldata$eventclass%in%c(0,1)))
Clinicaldata=Clinicaldata[se,]

#smoke
Clinicaldata$smok = as.numeric(Clinicaldata$smok)
se= which(Clinicaldata$smok %in% c(0 ,1, 2))

Clinicaldata$Model_glmnet <- 0
Clinicaldata$Model_xgboostdart <- 0
Clinicaldata$Model_xgboostlinear <- 0
# Outer loop for repeated cross-validation
# for (repeat_idx in 1:n_repeats) {
#   set.seed(15 + repeat_idx)
Clinicaldata$folds <- createFolds(Clinicaldata$eventclass, k = 5, list = FALSE)
#   

# INACCURATE: Replace missing values in finalgenematrix with 0
# se= which(is.na(Clinicaldata$finalgenematrix)) #7309
# Clinicaldata$finalgenematrix[se] <- 0
# se= which(is.infinite(Clinicaldata$finalgenematrix)) #5
# Clinicaldata$finalgenematrix[se] <- 0
#Missing value sets to mean
for(i in 1:ncol(Clinicaldata$finalgenematrix)){
  se=which(is.na(Clinicaldata$finalgenematrix[,i]))
  se2= which(is.infinite(Clinicaldata$finalgenematrix[,i]))
  cat(colnames(Clinicaldata$finalgenematrix)[i])
  cat("\n NA has", length(se),"\n\n")
  Clinicaldata$finalgenematrix[c(se,se2),i]<-colMeans(na.omit(Clinicaldata$finalgenematrix))[i]
}  

#Loop for 5CV
myList.1 <- list()
myList.2 <- list()
myList.3 <- list()
for(FOLDS in 1:5){
      ################# model 1
      ## Train GLMNET Model
      set.seed(15)

      GLM <- with(Clinicaldata[Clinicaldata$folds != FOLDS,c("finalgenematrix","eventclass")],{
        data.glmnet.lambda <- cv.glmnet(x = as.matrix(finalgenematrix), y = as.numeric(eventclass), type.measure="auc", family = "binomial", nfolds = 5)$lambda.min
        data.glmnet <- glmnet(x = as.matrix(finalgenematrix), y = as.numeric(eventclass), family = "binomial", lambda = data.glmnet.lambda,alpha = 0.5)
      })
      ## Make GLMNET prediction on samples not included in the training
      Clinicaldata$Model_glmnet[Clinicaldata$folds == FOLDS] <- predict(GLM,as.matrix(Clinicaldata$finalgenematrix), type="response",s=GLM$lambda)[Clinicaldata$folds == FOLDS,1] 
      myList.1[[FOLDS]]=Clinicaldata$Model_glmnet
      # 
      # ################# model 2
      # ##Train xgboost Model
      set.seed(16)
      xgboost1 <- with(Clinicaldata[Clinicaldata$folds != FOLDS,c("finalgenematrix","eventclass")],{

        gboostModeldart <- xgboost(data = as.matrix(Clinicaldata$finalgenematrix), label = Clinicaldata$eventclass,booster = "dart", max.depth = 10, eta = 0.05, nthread = 2, nrounds = 1000, objective = "binary:logistic",verbose = 0)

      })
      
      ## Make xgboost prediction on samples not included in the training
      Clinicaldata$Model_xgboostdart[Clinicaldata$folds == FOLDS] <- predict(xgboost1,as.matrix(Clinicaldata$finalgenematrix))[Clinicaldata$folds == FOLDS]
      myList.2[[FOLDS]]=Clinicaldata$Model_xgboostdart
      
      #################model 3
      ##Train xgboost Model
      set.seed(17)
      xgboost2 <- with(Clinicaldata[Clinicaldata$folds != FOLDS,c("finalgenematrix","eventclass")],{
        
        gboostModellinearboost <- xgboost(data = as.matrix(Clinicaldata$finalgenematrix), label = Clinicaldata$eventclass,booster = "gblinear", max.depth = 10, eta = 0.05, nthread = 2, nrounds = 1000, objective = "binary:logistic",verbose = 0)
        
      })
      ## Make xgboost prediction on samples not included in the training
      Clinicaldata$Model_xgboostlinear[Clinicaldata$folds == FOLDS] <- predict(xgboost2,as.matrix(Clinicaldata$finalgenematrix))[Clinicaldata$folds == FOLDS] 
      
      myList.3[[FOLDS]]=Clinicaldata$Model_xgboostlinear
    }

# Average predictions over all repeats
Clinicaldata$Model_glmnet <- rowMeans(sapply(myList.1, c))
Clinicaldata$Model_xgboostdart <- rowMeans(sapply(myList.2, c))
Clinicaldata$Model_xgboostlinear <- rowMeans(sapply(myList.3, c))

myoutf1="/mount/amos1/home/ylin/IMC/output/model-result.rds"
save.image(file = myoutf1)#save RDS file

#####[4]Training Model Evaluation####
# Evaluate model performance using the glmnet,xgboostdart,xgboostlinear
# se=which(label$braf_result=="Positive")
# # se=which(label$braf_result=="Negative")
# label=label[se,]
label = as.data.frame(Clinicaldata$eventclass)
rownames(label)=rownames(Clinicaldata$finalgenematrix)
thr = 0.5
model = "Model_xgboostlinear"
# model = "Model_xgboostdart"
# model = "Model_xgboostlinear"
# 1."Model_glmnet_binary"
# 2."Model_xgboostdart_binary"
# 3."Model_xgboostlinear_binary"
Clinicaldata$Model_glmnet_binary=NULL
Clinicaldata$Model_xgboostdart_binary=NULL
Clinicaldata$Model_xgboostlinear_binary=NULL

#------Class label and AUCcurve 
#Model1. glm
se= which(Clinicaldata$Model_glmnet>thr)
Clinicaldata[se,'Model_glmnet_binary']=1
Clinicaldata[-se,'Model_glmnet_binary']=0
roc_curve <- roc(Clinicaldata$eventclass, Clinicaldata$Model_glmnet)
plot(roc_curve, main = paste0("ROC Model_glmnet AUC: ", round(auc(roc_curve),3))) #0.716
round(auc(roc_curve),3) 

#Model2.xgboostdart
se= which(Clinicaldata$Model_xgboostdart>thr)
Clinicaldata[se,'Model_xgboostdart_binary']=1
Clinicaldata[-se,'Model_xgboostdart_binary']=0
roc_curve <- roc(Clinicaldata$eventclass, Clinicaldata$Model_xgboostdart)
plot(roc_curve, main = paste0("ROC Model_xgboostdart AUC: ", round(auc(roc_curve),3)))
round(auc(roc_curve),3) 

#Model3.xgboostlinear
se= which(Clinicaldata$Model_xgboostlinear>thr)
Clinicaldata[se,'Model_xgboostlinear_binary']=1
Clinicaldata[-se,'Model_xgboostlinear_binary']=0
roc_curve <- roc(Clinicaldata$eventclass, Clinicaldata$Model_xgboostlinear)
plot(roc_curve, main = paste0("ROC Model_xgboostlinear", round(auc(roc_curve),3))) #0.89

#------dotplot
test.df<-NULL
test.df=data.frame(predicted_class_binary=Clinicaldata[,paste0(model,"_binary")],
                   predicted_binary = as.numeric(Clinicaldata[,model]),binary_event=as.factor(Clinicaldata$eventclass) )

# rownames(test.df)=rownames(df)[se]
test.df$binary_event <- factor(test.df$binary_event, levels = c(0, 1), labels = c("ADC", "SQC"))
test.df$predicted_class_binary <- factor(test.df$predicted_class_binary , levels = c(0, 1), labels =  c("ADC", "SQC"))


pal=c("grey","#00BA38")#Order: Non-event, event, event+BRAF+, event-BRAF+
g=ggplot(test.df, aes(y=predicted_binary,x = binary_event))
g=g+geom_boxplot(outlier.shape = NA)
g=g+geom_jitter()
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g=g+theme(title=element_text(face="bold",size=20))
g=g+ggtitle(paste0("Training ", model,  " ", length(selected), " AUC: ",round(auc(roc(Clinicaldata$eventclass, Clinicaldata[,model])),3)))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="Black"))
g=g+theme(axis.text.x = element_text(size=10,color="Black"))
g=g+theme(aspect.ratio = 1/1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title =  element_blank())
# g=g+theme(legend.position = "None") #remove color legend
# g=g+scale_color_manual(values = pal)
g = g + geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")
g = g + labs(x = NULL, y = "Event-Mut Model") 
g
#####
#Boxplots
# Set significance level
significance_level <- 0.05

for(i in 1:ncol(Clinicaldata$finalgenematrix)) {
  df = NULL
  df$score = Clinicaldata$finalgenematrix[,i]
  df$group <- factor(Clinicaldata$eventclass, levels = c(0, 1), labels = c("ADC", "SQC"))
  
  df = as.data.frame(df)
  
  # Conduct t-test
  t_test_result <- t.test(score ~ group, data = df)
  
  # Check if the p-value is significant
  if (t_test_result$p.value < significance_level) {
    pal = brewer.pal(3, "Pastel2")
    
    g_order = names(sort(tapply(df$score, df$group, median, na.rm = TRUE), decreasing = TRUE))
    df$group <- factor(df$group, levels = g_order)
    
    g = ggplot(df, aes(x = group, y = score, fill = group)) + 
      geom_boxplot(color = "black") +
      scale_fill_manual(values = pal) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(size = 0.5, fill = NA),
            axis.text.y = element_text(size = 20, color = "black"),
            axis.text.x = element_text(size = 10, color = "black", angle = 30, vjust = 0.5, hjust = 1),
            plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
            aspect.ratio = 16/9,
            legend.background = element_rect(color = NA, fill = NA),
            legend.title = element_blank(),
            legend.position = "none",
            legend.text = element_text(size = 15)) +
      ylim(0, 5) +
      xlab("Subtype") + 
      ylab("ITH score") +
      labs(title = colnames(Clinicaldata$finalgenematrix)[i]) +
      stat_compare_means(method = "t.test")
    
    print(g)
  }
}

#----
significance_level <- 0.05
plots <- list()  # List to store significant plots

for (i in 1:ncol(Clinicaldata$finalgenematrix)) {
  df <- data.frame(score = Clinicaldata$finalgenematrix[, i],
                   group = factor(Clinicaldata$eventclass, levels = c(0, 1), labels = c("ADC", "SQC")))
  
  # Conduct t-test
  t_test_result <- t.test(score ~ group, data = df)
  
  # Check if the p-value is significant
  if (t_test_result$p.value < significance_level) {
    pal = brewer.pal(3, "Pastel2")
    
    g_order = names(sort(tapply(df$score, df$group, median, na.rm = TRUE), decreasing = TRUE))
    df$group <- factor(df$group, levels = g_order)
    
    g <- ggplot(df, aes(x = group, y = score, fill = group)) + 
      geom_boxplot(color = "black") +
      scale_fill_manual(values = pal) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(size = 0.5, fill = NA),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.text.x = element_text(size = 10, color = "black"),
            plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
            aspect.ratio = 16/9,  # Adjust aspect ratio as needed
            legend.background = element_rect(color = NA, fill = NA),
            legend.title = element_blank(),
            legend.position = "none",
            legend.text = element_text(size = 10)) +
      ylim(0, 5) +
      xlab(paste0("p= ", sprintf("%.2e", t_test_result$p.value))) + 
      ylab("ITH score") +
      labs(title = colnames(Clinicaldata$finalgenematrix)[i]) 
    
    # Add the plot to the list
    plots[[i]] <- g
  }
}

# Remove NULL entries
plots <- Filter(Negate(is.null), plots)

# Combine and display all significant plots in one canvas
if (length(plots) > 0) {
  combined_plot <- wrap_plots(plots) + 
    plot_layout(ncol = 6, guides = "collect")
  print(combined_plot)
} else {
  message("No significant plots to display.")
}



#######
#####
#Smoking Boxplots
# Set significance level
significance_level <- 0.05

for(i in 1:ncol(Clinicaldata$finalgenematrix)) {
  df = NULL
  df$score = Clinicaldata$finalgenematrix[,i]
  df$group <- factor(Clinicaldata$smok, levels = c(0, 1,2), labels = c("None", "Current", "Former"))
  
  df = as.data.frame(df)
  
    
    # g_order = names(sort(tapply(df$score, df$group, median, na.rm = TRUE), decreasing = TRUE))
    # df$group <- factor(df$group, levels = g_order)
    
    g = ggplot(df, aes(x = group, y = score, fill = group)) + 
      geom_boxplot(color = "black") +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(size = 0.5, fill = NA),
            axis.text.y = element_text(size = 20, color = "black"),
            axis.text.x = element_text(size = 10, color = "black", angle = 30, vjust = 0.5, hjust = 1),
            plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
            aspect.ratio = 16/9,
            legend.background = element_rect(color = NA, fill = NA),
            legend.title = element_blank(),
            legend.position = "none",
            legend.text = element_text(size = 15)) +
      ylim(0, 5) +
      xlab("Subtype") + 
      ylab("ITH score") +
      labs(title = colnames(Clinicaldata$finalgenematrix)[i]) 
      # stat_compare_means(method = "t.test")
    
    print(g)
  
}

#----
significance_level <- 0.05
plots <- list()  # List to store significant plots

for (i in 1:ncol(Clinicaldata$finalgenematrix)) {
  df <- data.frame(score = Clinicaldata$finalgenematrix[, i],
                   group = factor(Clinicaldata$smok, levels = c(0, 1), labels = c("ADC", "SQC")))
  
  # Conduct t-test
  # t_test_result <- t.test(score ~ group, data = df)
  
  # Check if the p-value is significant
  # if (t_test_result$p.value < significance_level) {
    pal = brewer.pal(3, "Pastel2")
    
    g_order = names(sort(tapply(df$score, df$group, median, na.rm = TRUE), decreasing = TRUE))
    df$group <- factor(df$group, levels = g_order)
    
    g <- ggplot(df, aes(x = group, y = score, fill = group)) + 
      geom_boxplot(color = "black") +
      scale_fill_manual(values = pal) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(size = 0.5, fill = NA),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.text.x = element_text(size = 10, color = "black"),
            plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
            aspect.ratio = 16/9,  # Adjust aspect ratio as needed
            legend.background = element_rect(color = NA, fill = NA),
            legend.title = element_blank(),
            legend.position = "none",
            legend.text = element_text(size = 10)) +
      ylim(0, 5) +
      # xlab(paste0("p= ", sprintf("%.2e", t_test_result$p.value))) + 
      ylab("ITH score") +
      labs(title = colnames(Clinicaldata$finalgenematrix)[i]) 
    
    # Add the plot to the list
    plots[[i]] <- g
  # }
}

# Remove NULL entries
plots <- Filter(Negate(is.null), plots)

# Combine and display all significant plots in one canvas
if (length(plots) > 0) {
  combined_plot <- wrap_plots(plots) + 
    plot_layout(ncol = 6, guides = "collect")
  print(combined_plot)
} else {
  message("No significant plots to display.")
}
