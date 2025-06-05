# clustering function (retrieved from melanoma paper)
library(pheatmap)
library(FlowSOM)
library(Rphenograph)
library(neighbouRhood)
library(dplyr)
library(ggplot2)
library(tibble)
runFlowsomPheno <- function(exp_df, markers, xdim=100, ydim=100, phenok=15){
  "xdim/ydim: nodes for FlowSOM,
  phenok: k NN for phenograph"
  set.seed(123)
  Do_FlowSOM = TRUE
  if(Do_FlowSOM){
    timestart<-Sys.time()
    
    # run SOM clustering:
    FlowSOM_input_data <- exp_df[,markers]
    print(dim(FlowSOM_input_data))
    map <- SOM(data= as.matrix(FlowSOM_input_data),
               xdim=xdim,
               ydim=ydim,
               silent = F)
    
    SOM_result<-map$mapping[,1]
    print(length(unique(SOM_result)))
    
    # run metacluster:
    FlowSOM_combined<-data.frame(FlowSOM=SOM_result,
                                 FlowSOM_input_data)
    
    metacluster_result<-metaClustering2(FlowSOM_combined,
                                        clustername = "FlowSOM",
                                        metaClustering_method = "metaClustering_PhenoGraph", 
                                        k_value=phenok, #<-- k for phenograph
                                        elbow_test=F, #<-- elbow test to determine k
                                        seed=123)   
    
    timeend<-Sys.time()
    runningtime<-timeend-timestart
    print(runningtime)
  }
  return(metacluster_result)
}

metaClustering2 <- function(indataframe=NULL,
                            clustername=NULL,
                            metaClustering_method=NULL,
                            usecol=NULL,
                            elbow_test=T,
                            k_value=NULL,
                            #tsne parameters：
                            view_tsne=T,
                            seed=NULL,
                            perplexity=15,
                            max_iter=1500,
                            ...)
{
  
  if(0){
    
    indataframe=FlowSOM_combined
    usecol=NULL
    #tsne parameters：
    view_tsne=T
    perplexity=15
    max_iter=1500
    
    clustername = "FlowSOM"
    metaClustering_method = "metaClustering_PhenoGraph" 
    k_value=10 
    elbow_test=F 
    seed=123
  }
  cat("Get expression matrix of cluster centers...\n")
  if (is.null(usecol))
  {usecol<-c(1:ncol(indataframe))}
  clusterparaid<-which(colnames(indataframe)==clustername)
  usecol<-union(usecol, c(clusterparaid))
  cluster_center_expr<-data.frame(indataframe[,usecol]) %>%
    dplyr::group_by_at(clustername) %>%
    dplyr::summarise_if(is.numeric,median)
  
  cluster_abundance<-data.frame(indataframe[,usecol]) %>%
    dplyr::group_by_at(clustername) %>%
    dplyr::summarise(num=n())
  
  cat("MetaClustering using method:",metaClustering_method,"...\n")
  if(elbow_test==T){
    cat("Drawing elow curve...\n")
  }
  #Findvalue<-DetermineNumberOfClusters(data=cluster_center_expr,max=20,kstep = 2,method=metaClustering_method,plot = T)
  
  if (metaClustering_method == 'metaClustering_PhenoGraph'){
    cc_metacluster <- metaClustering_PhenoGraph(cluster_center_expr[,-clusterparaid],k = k_value)
  }else{
    cc_metacluster <- MetaClustering(data=cluster_center_expr[,-clusterparaid],
                                     method=metaClustering_method)
  }
  
  cat("\nMetaclustering cluster centers is finished.\n")
  cat("Start to mapping metaclusters to single cells...\n")
  
  Cluster_arrange<-data.frame(cluster=cluster_center_expr[,clustername],
                              metacluster=cc_metacluster)
  
  cluster_arrange_fun<-function( cluster_id ){
    cellcluster <- subset(Cluster_arrange,Cluster_arrange[,1]==cluster_id)$metacluster
    return(cellcluster)
  }
  metacluster_result <- apply(as.matrix(indataframe[,colnames(indataframe)==clustername]),1,cluster_arrange_fun)
  
  if(view_tsne==T)
  {
    
    cat("Summarise metacluster information...\n")
    cat("Start to visualise metaclusters with tSNE...\n")
    #tsne (two important parameters：perplexity and max_iter)
    if(is.null(seed)){
      seed<-ceiling(runif(1,min=0,max=1)*10000)
      cat("Seed is not specified, randomly set to: ",seed,"\n")
      set.seed(seed)
    }else{
      cat("Seed is set to: ",seed,".\n")
      set.seed(seed)
    }
    
    tsne_result <- Rtsne(cluster_center_expr[,-c(1)], initial_dims = ncol(cluster_center_expr[,-c(1)]),
                         dims = 2, check_duplicates = FALSE, pca = F, perplexity=15,max_iter=1500)$Y
    colnames(tsne_result)<-c("tsne_1","tsne_2")
    
    combine_data_plot<-data.frame(cluster_center_expr,
                                  metacluster=cc_metacluster,
                                  tsne_result,
                                  num=cluster_abundance$num)
    
    centers<-combine_data_plot %>%
      dplyr::group_by(metacluster)  %>%
      dplyr::summarise(tsne_1=median(tsne_1),tsne_2=median(tsne_2))
    
    
    #visualization using ggplot2
    
    combine_data_plot$metacluster<-as.factor(combine_data_plot$metacluster)
    
    mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), 
                     legend.key = element_rect(fill = "white", colour = "white"), 
                     legend.background = (element_rect(colour= "white", fill = "white")))
    
    klab="Cluster number(k_value)"
    if(metaClustering_method=="metaClustering_PhenoGraph") klab="PhenoGraph_k(k_value)"
    ptsnemap<-ggplot(combine_data_plot)+
      geom_point(aes(x=tsne_1,y=tsne_2,colour=metacluster,size=num),alpha=0.7)+
      guides(colour = guide_legend(ncol = 2, bycol = T))+
      scale_size_continuous(range = c(0.1, 5))+
      labs(title = paste0(klab,": ",k_value))+
      mytheme+
      geom_text(data=centers,aes(x=tsne_1,y=tsne_2),label=centers$metacluster,colour="black",size=5)
    
    print(ptsnemap)
  }
  
  cat("Metaclustering is finished successufully.\n")
  return(list(metacluster_result=metacluster_result, tsne = tsne_result, plotdf = combine_data_plot))
}


metaClustering_PhenoGraph<-function(data,k=30){
  #as.numeric(membership(Rphenograph(data,k=k)))
  as.numeric(Rphenograph(data,k=k)[[2]]$membership)
}

renameMarker <- function(oldname){
  oldname <- ifelse(oldname == 'HLA_DR','HLA-DR', oldname)
  oldname <- ifelse(oldname == 'Alpha_SMA','SMA', oldname)
  oldname <- ifelse(oldname == 'E_Cadherin','E-cadherin', oldname)
  oldname <- ifelse(oldname == 'CA_IX','CAIX', oldname)
  oldname <- ifelse(oldname == 'Collagen_I','Collagen-I', oldname)
  oldname <- ifelse(oldname == 'Ki67','Ki-67', oldname)
  oldname <- ifelse(oldname == 'CD274_PDL1','PDL1', oldname)
  oldname <- ifelse(oldname == 'CD134_OX40','OX40', oldname)
  oldname <- ifelse(oldname == 'CD279_PD1','PD1', oldname)
  oldname <- ifelse(oldname == 'CD223_LAG3','LAG3', oldname)
  oldname <- ifelse(oldname == 'CD366_TIM3','TIM3', oldname)
  oldname <- ifelse(oldname == 'CD278_ICOS','ICOS', oldname)
  oldname <- ifelse(oldname == 'CD278','ICOS', oldname)
  return(oldname)
}

panel = read.csv('/mount/amos1/home/ylin/IMC/IMC_melanoma/data/input/meta/panel_metadata.csv', stringsAsFactors = FALSE)
markers <- renameMarker(panel$markers[panel$clustering_type_20==1])
a <- read.csv('/mount/amos1/home/ylin/IMC/IMC_melanoma/data/input/sc_data.csv',row.names = 1, check.names = FALSE)

result_test <- runFlowsomPheno(a, markers, xdim=100, ydim = 100)
a['cluster1'] <- (result_test[['metacluster_result']])
#write.csv(a, file = '../data/output/sc_data.csv')

marker_total <- renameMarker(panel$markers[panel$clustering_total==1])
colnames(a) <- renameMarker(colnames(a))
plotHeatmap(a, marker_total, 'flowsom100pheno15')