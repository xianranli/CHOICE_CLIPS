



CLIPS <- function(gDNAs_blast_,User_working_dir) {
 
 genomes <- as.vector(unique(gDNAs_blast_[,1]));

 self_idx <- which(gDNAs_blast_[,1] == gDNAs_blast_[,2])
 g_selfs <- gDNAs_blast_[self_idx,]

  g_selfs <- g_selfs[(g_selfs[,7] - g_selfs[,9] + g_selfs[,8] - g_selfs[,10] != 0),]

 
 gDNAs_blast_0 <- gDNAs_blast_[-self_idx,]

 n_g <- length(genomes)  
 b_matrix <- diag(n_g)
 for (g1 in 1:(n_g - 1)) {
  g1_self <- subset(g_selfs, g_selfs[,1] == genomes[g1])
   for (g2 in (g1 + 1): n_g ) {
    g2_self <- subset(g_selfs, g_selfs[,1] == genomes[g2]) 
    b_matrix[g1, g2] <- Slope_by_HSPs(gDNAs_blast_0, genomes[g1], genomes[g2], g1_self, g2_self)
    b_matrix[g2, g1] <- Slope_by_HSPs(gDNAs_blast_0, genomes[g2], genomes[g1], g2_self, g1_self)
   }
 }
 colnames(b_matrix) <- genomes
 rownames(b_matrix) <- genomes

 write.table(round(b_matrix,2), file = paste(User_working_dir, 'slop_matrix_output_from_CLIPS.txt', sep = ''),row.names=TRUE,col.names=TRUE,quote = FALSE,sep="\t")

# h_c <- hclust( dist(round(b_matrix,2) ) )
# return (h_c)
}

Slope_by_HSPs <- function(gDNAs_, g1_, g2_, g1_self_, g2_self_) {
 gDNAs <- subset(gDNAs_, gDNAs_[,1] == g1_ & gDNAs_[,2] == g2_)
 exact_match <- unique(c(which(gDNAs[,7] %in% g1_self_[,7] & gDNAs[,8] %in% g1_self_[,8]), which(gDNAs[,9] %in% g2_self_[,7] & gDNAs[,10] %in% g2_self_[,8])))
 if (length(exact_match) > 0) { gDNAs <- gDNAs[-exact_match,]}
 if (nrow(g1_self_) > 0 & nrow(g2_self_) > 0 & nrow(gDNAs) > 1) { gDNAs <- Remove_CrossOverMatch(gDNAs, g1_self_, g2_self_) }
 q_hits <- c(gDNAs[,7], gDNAs[,8]); s_hits <- c(gDNAs[,9,], gDNAs[,10])
 b_ <- round(lm(q_hits ~ s_hits)$coeff[2], 2)
 return (b_)
}

Remove_CrossOverMatch <- function(gDNAs_, g1_self_, g2_self_) {
 rep_tag <- c()
 for (i in 1:nrow(gDNAs_)) {
  for (j in 1:nrow(g1_self_)) {
   g1_s_e <- sort(c(g1_self_[j,7], g1_self_[j,8], gDNAs_[i,7], gDNAs_[i,8]))
   d_s <- abs(g1_s_e[1] - g1_s_e[2]) + abs(g1_s_e[3] - g1_s_e[4]);
   if (d_s < 20) { rep_tag <- append(rep_tag, i); break} ## 20 can be changed
  }
   ## check the subject columns might be reduant, but just be safe
  for (k in 1:nrow(g2_self_)) {
   g2_s_e <- sort(c(g2_self_[k,7], g2_self_[k,8], gDNAs_[i,9], gDNAs_[i,10]))
   d_s <- abs(g2_s_e[1] - g2_s_e[2]) + abs(g2_s_e[3] - g2_s_e[4]);
   if (d_s < 20) { rep_tag <- append(rep_tag, i); break} ## 20 can be changed
  }
 }
 if (length(rep_tag) > 0) { gDNAs_ <- gDNAs_[-unique(rep_tag),] }
 return(gDNAs_);
}

###
User_working_dir<-"/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/script/For_github_upload/"

Gene<-"TraesCS4A02G058900"

gDNAs_blast <- read.table(paste(User_working_dir, Gene, '_Haplotype-Self_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);

CLIPS(gDNAs_blast,User_working_dir)

#library(dendextend)
#as.dendrogram(h_c) %>% set("labels_cex", 0.9) %>% highlight_branches %>% plot(main = "CLIPS",ylab = "Height",horiz = FALSE);