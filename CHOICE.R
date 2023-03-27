



CHOICE<- function(User_working_dir,BlastSynWorking,query_length, distance_filter, Min_CDS_size_filter, Max_CDS_size_filter, ref_g) {
 
 CHOICE_input_file<-BlastSynWorking

 CHOICE_output_list <- list() 
 CHOICE_cluster_list<-list() 

 index_genome<-1 

 for(g_index in unique(CHOICE_input_file$Genome)){

   Target_g<- CHOICE_input_file[which(CHOICE_input_file$Genome==g_index), ]

    if(nrow(Target_g)>=2){
     


     Total_clusters<-length( unique( cutree( hclust(dist(Target_g[,6]) ), h = distance_filter) ) ) 
     
     hclust_similarity<-list() 
     hclust_size<-list() 
     hclust_members<-matrix(nrow=Total_clusters,ncol=1)

      for(cluster_index in 1:Total_clusters) {
       hclusters <-Target_g[ which( cutree(hclust(dist(Target_g[,6])), k =Total_clusters, h = distance_filter ) == cluster_index), ]
       hclust_similarity[cluster_index] <- mean(hclusters[,9])  
       hclust_size[cluster_index] <- sum(hclusters[,8])   
       hclust_members[cluster_index,1]<-nrow(hclusters) 
       }

     MeanSimilarity_cdsSize<-matrix(nrow=Total_clusters,ncol=2) 

      for(cluster_index in 1:Total_clusters) {
       MeanSimilarity_cdsSize[cluster_index,1]<-hclust_similarity[[cluster_index]]
       MeanSimilarity_cdsSize[cluster_index,2]<-hclust_size[[cluster_index]]
       }

     MeanSimilarity_cdsSize[,2]<-MeanSimilarity_cdsSize[,2]/query_length
     Size_filter<-which(MeanSimilarity_cdsSize[,2]>=Min_CDS_size_filter & MeanSimilarity_cdsSize[,2]<=Max_CDS_size_filter)
     
     if(length(Size_filter)>1){                                                               
      Similarity_filter<-which( max( MeanSimilarity_cdsSize[Size_filter,][,1] ) == MeanSimilarity_cdsSize[,1] ) 
      } else {
        Similarity_filter<-which( max( MeanSimilarity_cdsSize[Size_filter,][1] ) == MeanSimilarity_cdsSize[,1]  )
        }

     Target_cluster<-intersect(Size_filter,Similarity_filter)

     if(length(Target_cluster)==0){                             
      Ideal_Size <- 1.0
      Target_cluster<- which( abs(MeanSimilarity_cdsSize[,2] - Ideal_Size) == min( abs(MeanSimilarity_cdsSize[,2] - Ideal_Size) ) )
      } else if(length(Target_cluster)>1){
        Target_cluster<-which( max(MeanSimilarity_cdsSize[Target_cluster,][,1]) == MeanSimilarity_cdsSize[,1] )
        }

     CHOICE_output_list[[index_genome]]<-Target_g[ which( cutree( hclust( dist(Target_g[,6])), k =Total_clusters, h = distance_filter ) == Target_cluster), ]

     CHOICE_genome_summary<-matrix(nrow=Total_clusters,ncol=9)
     CHOICE_genome_summary[ ,1]<-rep(g_index, Total_clusters)
     CHOICE_genome_summary[1,2]<-nrow(Target_g)
     CHOICE_genome_summary[1,3]<-distance_filter/1000
     CHOICE_genome_summary[1,4]<-Total_clusters
     CHOICE_genome_summary[ ,5]<-1:nrow(MeanSimilarity_cdsSize)
     CHOICE_genome_summary[ ,6]<-hclust_members
     CHOICE_genome_summary[ ,7]<-round(MeanSimilarity_cdsSize[,1], digits = 3)
     CHOICE_genome_summary[ ,8]<-round(MeanSimilarity_cdsSize[,2], digits = 3)
     CHOICE_genome_summary[Target_cluster,9]<-c("Selected")

     colnames(CHOICE_genome_summary)<-c("Genomes","TotalHSPs","HeightCut(kb)","TotalClusters","ClusterIndex","Members","MeanSimilarity",paste("TotalLength","/",ref_g,sep=''),"CandidateCluster")
     CHOICE_output_Table<-as.data.frame(CHOICE_genome_summary)
     CHOICE_cluster_list[[index_genome]]<-CHOICE_output_Table 
    
     } else {
       CHOICE_output_list[[index_genome]]<-CHOICE_input_file[ which(CHOICE_input_file$Genome==g_index), ] 
       }

      index_genome<-index_genome+1

                }

   Filtered_By_CHOICE <-as.data.frame(rbindlist(CHOICE_output_list))
   
   if(length(CHOICE_cluster_list)!=0){
    CHOICE_summary_table <-as.data.frame(rbindlist(CHOICE_cluster_list))
    }

   write.table(Filtered_By_CHOICE, file=paste(User_working_dir,"Filtered_By_CHOICE.txt",sep=''), sep= "\t",quote = FALSE,row.names = FALSE )
   write.table(CHOICE_summary_table, file=paste(User_working_dir,"CHOICE_summary.txt",sep=''), sep= "\t",quote = FALSE,row.names = FALSE )

 #  my_list <- list(Filtered_By_CHOICE,CHOICE_summary_table)
 #  return(my_list)

}
#######################

########## example's parameters
User_working_dir<-"/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/script/For_github_upload/"

Gene<-"TraesCS4A02G058900"

BlastSynWorking <-read.table(paste(User_working_dir,Gene,"_Haplotype_syn",sep=''),header=T)

query_length<- 513  # reference CDS size
distance_filter<-20*1000 # distance filter
Min_CDS_size_filter <- 0.75 # CDS/query_length min value
Max_CDS_size_filter <- 1.25 # CDS/query_length max value
ref_g<-"IWGSC"                # reference genome
###########

CHOICE(User_working_dir,BlastSynWorking,query_length, distance_filter, Min_CDS_size_filter, Max_CDS_size_filter, ref_g)
