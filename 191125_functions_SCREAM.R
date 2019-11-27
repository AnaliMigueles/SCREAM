####### Functions

##### Function 1: Make a summary of the number of genes of each scRNA-seq cell type 
summary_background<-function(single_cell, list_clouds, output) {
  #selected<-single_cell[single_cell$cluster %in% list_clouds,]
  #clouds<-as.vector(selected$cluster)
  #clouds<-sort(unique(clouds))
  #list_clouds<-sort(as.character(unique(selected_smedwi$cluster)))
  #single_cell<-selected_smedwi
  
  clouds<-list_clouds
  
  total<-single_cell[single_cell$cluster %in% list_clouds,]
  clusterx_single<-total
  #Ahora creamos la matrix donde vamos a ir guardando todo
  row_total<-sort(unique(as.character(total$ID)))
  col_total<-clouds
  
  total_x<-matrix(data=1, ncol = length(col_total), nrow = length(row_total))
  colnames(total_x)<-col_total
  rownames(total_x)<-row_total
  #Teniendo ahora la referencia, podemos empezar a llenar la matrix
  
  for (i in 1:length(row_total)) {
    x<-total[total$ID%in% (row_total[i]),]
    sub_x<-sort(unique(as.character(x$cluster)))
    
    values_total<-match(col_total,sub_x)
    values_total [!is.na(values_total)]<-1
    values_total [is.na(values_total)]<-0
    total_x [i,]<-values_total
  }
  
  summary<-t(as.matrix(colSums(total_x)))
  rownames(summary)<-"#genes in each single cell cluster"
  
  
  matriz<-matrix(data = NA, ncol = length(list_clouds), nrow = 1)
  colnames(matriz)<-list_clouds
  rownames(matriz)<-c("Number of genes in each single cell cluster")
  for(i in 1:length(list_clouds)) {
    ##Quitar duplicados de las tablas 
    clusterx_single1<- clusterx_single[ clusterx_single$cluster %in% list_clouds[i],] 
    clusterx_single1<-clusterx_single1[!duplicated(clusterx_single1$ID),]
    num<-length(which(clusterx_single1$cluster==list_clouds[i]))
    matriz[1,i]<-num
  }
  
  
  if(output==1){
    return(matriz)
  }
  if(output==2) {
    return(total_x)
  }  
}


##### Function 2 Enrichment

enrichment<- function (background, summary, summary_A, cells_FULL, clusterA) {
  
  mhyper_A<-matrix(data=NA, ncol = dim(summary)[2], nrow=1)
  colnames(mhyper_A)<-cells_FULL
  rownames(mhyper_A)<-c("A") #Here you can put the name of the cluster
  p_values<-vector()
  
  universe<-dim(background)[1] #background
  for (i in 1:dim(summary)[2]) {
    white<-summary[i]
    black<-universe-white
    white_inSample<-summary_A[i]
    sample_size<-dim(clusterA)[1]
    
    math<-phyper(white_inSample-1, white, black, sample_size, lower.tail = FALSE)
    p_values<-append(p_values, math)
  }
  mhyper_A[1,]<-p_values
  return(mhyper_A)
}



######## Function 3: Overlap 

overlap<-function (binary_m, summary_cluster, summary_general, background, clusterx, pvalue, enrichment_logical_adj, list_clouds, output) {
  A_list<-vector()
  binary_m_df<-as.data.frame(binary_m) 
  for (i in 1:dim(binary_m)[2]){
    x<-subset(binary_m_df, binary_m_df[,i]==1)
    x_list<-rownames(x)
    A_list<-append(A_list, x_list)
  }
  
  m_A<-matrix(data = NA, ncol = 2, nrow = length(A_list))
  colnames(m_A)<-c("ID", "cluster")
  m_A[,1]<-A_list
  m_A[,2]<-rep(c(colnames(summary_cluster)), c(summary_cluster[1,]))
  
  ####### CLUSTER A
  df<-as.data.frame(m_A)
  populations<-as.character(unique(df$cluster))
  genes_populations<-matrix(data = NA, ncol = length(populations), nrow = length(populations))
  colnames(genes_populations)<-populations
  rownames(genes_populations)<-populations
  vector<-vector()
  
  for (i in 1:length(populations)) {
    for (j in 1:length(populations)) {
      selected1<-df[df$cluster %in% populations[i],]
      n1<-unique(selected1$ID)
      set1<-n1
      selected2<-df[df$cluster %in% populations[j],]
      n2<-unique(selected2$ID)
      set2<-n2
      overlap <- sum(set1 %in% set2)/length(set1) #what fraction of set1 is in set2
      genes_populations[j,i] <- overlap
      
    }
  }
  
  A_genes_populations_percentage<-  genes_populations * 100
  
  
  
  df<-as.data.frame(m_A)
  populations<-as.character(unique(df$cluster))
  genes_populations<-matrix(data = NA, ncol = length(populations), nrow = length(populations))
  colnames(genes_populations)<-populations
  rownames(genes_populations)<-populations
  vector<-vector()
  for (i in 1:length(populations)) {
    for (j in 1:length(populations)) {
      selected1<-df[df$cluster %in% populations[i],]
      n1<-unique(selected1$ID)
      set1<-n1
      selected2<-df[df$cluster %in% populations[j],]
      n2<-unique(selected2$ID)
      set2<-n2
      overlap <- sum(set1 %in% set2) #what fraction of set1 is in set2
      genes_populations[j,i] <- overlap
      
    }
  }
  
  
  A_genes_populations_interseccion<- genes_populations
  
  
  ######## Function 4: Enrichment of the overlap of genes across scRNA-seq cell types
  
  summary_clusterhypo <-as.data.frame(summary_general)
  keeps <- colnames(A_genes_populations_percentage)
  summary_clusterhypo<-summary_clusterhypo[keeps]
  summary_clusterhypo<-as.matrix(summary_clusterhypo)
  mhyper_A<-matrix(data=NA, ncol = dim(summary_clusterhypo)[2], nrow=dim(summary_clusterhypo)[2])
  colnames(mhyper_A)<-colnames(A_genes_populations_percentage)
  rownames(mhyper_A)<-colnames(A_genes_populations_percentage)
  
  universe<-dim(background)[1]
  for (i in 1:dim(summary_clusterhypo)[2]) {
    p_values<-vector()
    for (j in 1:dim(summary_clusterhypo)[2]) {
      white<-summary_clusterhypo[j]
      black<-universe-white
      white_inSample<-A_genes_populations_interseccion[j,i]
      sample_size<-dim(clusterx)[1]
      
      math<-phyper(white_inSample-1, white, black, sample_size, lower.tail = FALSE)
      p_values<-append(p_values, math)
    }
    mhyper_A[,i]<-p_values
  }
  
  #### Heatmap ###### 
  all_names<-colnames(summary_general)
  your_name<-colnames(mhyper_A)
  diff<-setdiff(all_names, your_name)
  
  the_name<-append(your_name, diff)
  
  enrichment_logical<-as.matrix((mhyper_A))
  a<-dim(enrichment_logical)[1]
  k<-vector()
  enrichment_logical2<-matrix(data=NA, ncol = length(the_name), nrow = length(the_name))
  colnames(enrichment_logical2)<-the_name
  rownames(enrichment_logical2)<-the_name
  for (i in 1:a) {
    k<-enrichment_logical[,i]
    k<-append(k, rep(1,length(diff)))
    k<-p.adjust(k )
    enrichment_logical2[,i]<-k
  }
  
  
  if(a!=length(the_name)) {
    for(j in (i+1):dim(enrichment_logical2)[1]) {
      k<-vector()
      k<-rep(1,dim(enrichment_logical2)[1])
      enrichment_logical2[,j]<-k
    }
    
  }
  hypo_enrichment_A<-enrichment_logical2
  
  ####### Cluster A
  
  ###### 
  enrichment<-as.data.frame(hypo_enrichment_A)
  enrichment<-enrichment[order(rownames(enrichment)),]
  enrichment<-enrichment[ ,order(colnames(enrichment))]
  clusterA_summary<-matrix(data = NA, ncol = 4, nrow = length(list_clouds))
  colnames(clusterA_summary)<-c("Cluster", "Population", "Enrichment(adj pvalue)","Possible identities")
  
  clusterA_summary[,1]<-rep("Cluster", length(names))
  clusterA_summary[,2]<-list_clouds
  clusterA_summary[,3]<-enrichment_logical_adj[1,]
  
  for(i in 1:length(list_clouds)) {
    n<-subset(enrichment, enrichment[,i]<=pvalue)
    n<-rownames(n)
    if (length(n)!=0){
      if(length(n)==1){
        if (n==as.character(clusterA_summary[i,2])) {
          clusterA_summary[i,4]<-paste(n, "is/are enriched")
        } 
      }else {
        clusterA_summary[i,4]<-paste(n, collapse  = ",")
      }
      
    } else {
      clusterA_summary[i,4]<-paste("No other possible enrichments")
    }
  }
  
  if(output==1) {
    return(clusterA_summary)
  }
  if(output==2) {
    return(A_genes_populations_interseccion)
  }
  if(output==3) {
    return(mhyper_A)
  }
  if(output==4) {
    return(hypo_enrichment_A)
  }
  if(output==5) {
    return(A_genes_populations_percentage)
  }
} #end of function


