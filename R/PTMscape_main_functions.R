### do every thing without decoy




# pssm_generation ---------------------------------------------------------



pssm_generation = function(candidate_Rds_name, 
                           center_position,
                           output_label)
{
  
  candidate = readRDS(candidate_Rds_name)
  
  pos_window <- candidate %>% dplyr::filter(label == "positive") %>%
    dplyr::select(window)
  
  ### to check if size and center position aligns
  if(center_position == (nchar(pos_window$window[1])+1)/2) 
  {
    noc_pos_window = pos_window$window
    str_sub(noc_pos_window, center_position, center_position) <- ""
    
  }else{
    warning("window size error!")
  }

  #neg_windows = noc_candi_window
  pos_windows = noc_pos_window
  
  
  positive_base_matrix  = get_base_freq_matrix(pos_windows)
  
  pssm=positive_base_matrix
  pssm[is.infinite(pssm)]=0
  
  saveRDS(pssm, file = paste0(output_label,"_pssm.Rds"))
  
  pssm_df = data.frame(AA =c("A","R","N","D","C","Q","E","G","H","I",
                             "L","K","M","F","P","S","T","W","Y","V","X","U"), pssm)
  
  write.table(pssm_df, paste0(output_label,"_pssm.tsv"),
              quote = F, sep = "\t", row.names = F)
  
  
}



# window_formation --------------------------------------------------------

window_formation = function(mod_site, flanking_size,
                            prot_seqs, prot_ids, positive_info, 
                            output_label)
{
  candidate = get_all_candidate(mod_site, flanking_size,prot_seqs, prot_ids, positive_info)
  saveRDS(candidate, file = paste0(output_label,"_candidate.Rds"))
  write.table(candidate, paste0(output_label,"_candidate.tsv"),
              row.names = F, quote = F, sep = "\t")
  
}


window_formation_no_positive = function(mod_site, flanking_size,
                            prot_seqs, prot_ids, 
                            output_label)
{
  candidate = get_all_candidate_no_positive(mod_site, flanking_size,prot_seqs, prot_ids)
  saveRDS(candidate, file = paste0(output_label,"_candidate.Rds"))
  write.table(candidate, paste0(output_label,"_candidate.tsv"),
              row.names = F, quote = F, sep = "\t")
  
}


# pssm_feature_extraction -------------------------------------------------


## there is no way to get the protein mean of pssm features


pssm_feature_extraction = function(candidate_Rds_name,  
                                   pssm_Rds_name,
                                   center_position,
                                   output_label)
{
  candidate = readRDS(candidate_Rds_name)
  pssm = readRDS(pssm_Rds_name)
  
  
  
  candi_window <- candidate %>% dplyr::filter(label == "negative") %>%
    dplyr::select(window)
 # cat("candi size: ",nrow(candi_window),"\n")
  
  ### to check if size and center position aligns
  if(center_position == (nchar(candi_window$window[1])+1)/2) 
  {
    ### remove the center site
 
    noc_candi_window = candi_window$window
    str_sub(noc_candi_window, center_position, center_position) <- ""
    
  }else{
    warning("window size error!")
  }
  
  
  noc_candi_pssm = t(sapply(1:length(noc_candi_window), function(x) {
    get_pssm_feature(pssm,noc_candi_window[x])
  }))
  
  saveRDS(noc_candi_pssm, file = paste0(output_label, "_noc_candi_pssm_matrix.Rds"))
  
  
  
  
  
  
  pos_window <- candidate %>% dplyr::filter(label == "positive") %>%
    dplyr::select(window)
  #cat("pos size: ",nrow(pos_window),"\n")
  
  if(nrow(pos_window)>0)
  {
    
    if(center_position == (nchar(pos_window$window[1])+1)/2) 
    {
      ### remove the center site
      noc_pos_window = pos_window$window
      str_sub(noc_pos_window, center_position, center_position) <- ""
      
    }else{
      warning("window size error!")
    }
    
    noc_pos_pssm = t(sapply(1:length(noc_pos_window), function(x) {
      get_pssm_feature(pssm,noc_pos_window[x])
    }))
    
    saveRDS(noc_pos_pssm, file = paste0(output_label, "_noc_pos_pssm_matrix.Rds"))
  
  }

}



# aaindex_feature_extraction ----------------------------------------------

### input: Rdata for candidate and decoy, the file of aaindex properties

### output: Rdata of the aaindex matrix 

aaindex_feature_extraction = function( aaindex_property, candidate_Rds_name, 
                                       center_position,
                                       output_label)
{
  ### for each protein get the protein mean of the features
  
  
  candidate = readRDS(candidate_Rds_name) 

  
  
  candi = candidate%>%dplyr::filter(label == "negative")
  candi_window <- candidate %>% dplyr::filter(label == "negative") %>%
    dplyr::select(window)
  
  
  pos = candidate%>%dplyr::filter(label == "positive")
  pos_window <- candidate %>% dplyr::filter(label == "positive") %>%
    dplyr::select(window)
  rm(candidate)
  
  
  ### to check if size and center position aligns
  if(center_position == (nchar(candi_window$window[1])+1)/2) 
  {
    noc_candi_window = candi_window$window
    str_sub(noc_candi_window, center_position, center_position) <- ""
    
    noc_candi_cluster = get_matrix_all(noc_candi_window,
                                       aaindex_property)
    
    saveRDS(noc_candi_cluster,file = paste0(output_label,"_noc_candi_cluster_matrix.Rds"))

   # cat("get window aaindex","\n")
    
    
  }else{
    warning("window size error!")
  }
  
  
  if(nrow(pos_window)>0)
  {
    
    if(center_position == (nchar(pos_window$window[1])+1)/2) 
    {
      noc_pos_window = pos_window$window
      str_sub(noc_pos_window, center_position, center_position) <- ""
      
      
      
      noc_pos_cluster = get_matrix_all(noc_pos_window,
                                       aaindex_property)
      
      saveRDS(noc_pos_cluster,file = paste0(output_label,"_noc_pos_cluster_matrix.Rds"))

      #cat("get window aaindex","\n")
      
      
    }else{
      warning("window size error!")
    }
  }
  
  
}



# spider_feature_extraction -----------------------------------------------


spider_feature_joining_without_mean = function(candidate_Rds_name,
                                               extracted_spider_Rds_name,
                                               spider_which_retain,
                                               spider_which_logit,
                                               center_position,
                                               #spider_which_center,
                                               output_label)
  
{


  protID_pos_spider_site_specific = readRDS(extracted_spider_Rds_name)
  
  ###for each position there is 4 spider properties.
  
  ### so the four properties for the center site are in column 
  
  center_start = 4*(center_position-1)+1
  
  spider_which_center = c(center_start:(center_start+3))
  
  
  candidate = readRDS(candidate_Rds_name) 
  
  
  pos = candidate%>%dplyr::filter(label == "positive")
  candi = candidate%>%dplyr::filter(label == "negative")
  
  #cat("see mem change","\n")
  #cat(mem_change(rm(candidate)),"\n")
  
  #### pos
  
  pos_site_specific = dplyr::left_join(pos,protID_pos_spider_site_specific,by = c("protID", "pos"))
  
  rm(pos)
  
  #### the only correct way to get from the back 
  
  last_col = ncol(pos_site_specific)
  
  get_col = c((last_col-100+1):last_col)
  
  pos_site_specific_matrix = as.matrix(pos_site_specific[,get_col])

  rm(pos_site_specific)
  
  noc_pos_structure = pos_site_specific_matrix[, -spider_which_center]
  
  saveRDS(noc_pos_structure, file = paste0(output_label, "_noc_pos_structure_matrix.Rds"))
  
  rm(noc_pos_structure)
  
  
 
  # cat("pos_processed","\n")
  
  #### candi
  candi_site_specific = dplyr::left_join(candi,protID_pos_spider_site_specific,by = c("protID", "pos"))
  
  rm(candi)
  candi_site_specific_matrix = as.matrix(candi_site_specific[,get_col])
  
  rm(candi_site_specific)
  
  noc_candi_structure = candi_site_specific_matrix[, -spider_which_center]
  
  saveRDS(noc_candi_structure, file = paste0(output_label, "_noc_candi_structure_matrix.Rds"))
  rm(noc_candi_structure)
  
  
  #cat("candi_processed","\n")
  
  
}



# feature_combining -------------------------------------------------------



combine_all_features = function(pos_aaindex = NULL,
                                candi_aaindex,
                                pos_spider = NULL,
                                candi_spider,
                                pos_pssm = NULL,
                                candi_pssm,
                                output_label)
  
{
  
  noc_not_na_candi_cluster = readRDS(candi_aaindex)
  noc_not_na_candi_structure = readRDS(candi_spider)
  noc_not_na_candi_pssm = readRDS(candi_pssm)


  noc_not_na_candi_feature = cbind(noc_not_na_candi_cluster,
                                   noc_not_na_candi_structure,
                                   noc_not_na_candi_pssm)
  saveRDS(noc_not_na_candi_feature, file = paste0(output_label, "_noc_not_na_candi_feature.Rds"))
  rm(noc_not_na_candi_feature)
  
  
  
  if(!is.null(pos_aaindex))
  {
    noc_not_na_pos_cluster = readRDS(pos_aaindex)
    noc_not_na_pos_structure = readRDS(pos_spider)
    noc_not_na_pos_pssm = readRDS(pos_pssm)
    
    noc_not_na_pos_feature = cbind(noc_not_na_pos_cluster,
                                   noc_not_na_pos_structure,
                                   noc_not_na_pos_pssm)
    saveRDS(noc_not_na_pos_feature, file = paste0(output_label, "_noc_not_na_pos_feature.Rds"))
    rm(noc_not_na_pos_feature)
    
  }
  
}




combine_all_features_no_spider = function(pos_aaindex = NULL,
                                candi_aaindex,
                                pos_pssm = NULL,
                                candi_pssm,
                                output_label)
  
{
  
  noc_not_na_candi_cluster = readRDS(candi_aaindex)
  noc_not_na_candi_pssm = readRDS(candi_pssm)
  
  
  noc_not_na_candi_feature = cbind(noc_not_na_candi_cluster,
                                   noc_not_na_candi_pssm)
  saveRDS(noc_not_na_candi_feature, file = paste0(output_label, "_noc_not_na_candi_feature.Rds"))
  rm(noc_not_na_candi_feature)
  
  
  
  if(!is.null(pos_aaindex))
  {
    noc_not_na_pos_cluster = readRDS(pos_aaindex)
    noc_not_na_pos_pssm = readRDS(pos_pssm)
    
    noc_not_na_pos_feature = cbind(noc_not_na_pos_cluster,
                                   noc_not_na_pos_pssm)
    saveRDS(noc_not_na_pos_feature, file = paste0(output_label, "_noc_not_na_pos_feature.Rds"))
    rm(noc_not_na_pos_feature)
    
  }
  
}




### more process with the negative data set , actually the negative dataset is either decoy or candidate

### and the process can be two forms now, 1 knn cleaning, 2 just randomly sample the same number as the positive set

# negative_selection ------------------------------------------------------



balanced_size_sampling = function(pos_matrix_Rds_name, neg_matrix_Rds_name, 
                                  output_label)
{
  
  pos_matrix = readRDS(pos_matrix_Rds_name)
  neg_matrix = readRDS(neg_matrix_Rds_name)
  
  pos_size = nrow(pos_matrix)
  neg_size = nrow(neg_matrix)
  
  set.seed(123)
  
  chosen_neg = sample(neg_size, pos_size)
  
  
  balanced_neg = neg_matrix[chosen_neg, ]
  
  saveRDS(balanced_neg, file = paste0(output_label, "_balanced_neg_feature.Rds"))
  
  
  
}



balanced_size_sampling_after_combining_pos_neg = function(feature_matrix_Rds_name,
                                                          feature_label_Rds_name,
                                                          output_label)
  
{
  
  feature_matrix = readRDS(feature_matrix_Rds_name)
  feature_label = readRDS(feature_label_Rds_name)
  
  
  pos_ind = which(feature_label == 1)
  neg_ind = which(feature_label == 0)
  
  
  pos_feature_matrix = feature_matrix[pos_ind,]
  neg_feature_matrix = feature_matrix[neg_ind,]
  
  
  
  rm(feature_matrix)
  
  pos_size = nrow(pos_feature_matrix)
  neg_size = nrow(neg_feature_matrix)
  
  
  if (pos_size<neg_size)
  {
    
    
    set.seed(123)
    
    chosen_neg_ind = sample(neg_size, pos_size)
    chosen_neg = neg_feature_matrix[chosen_neg_ind, ]
    
    
  }else{
    chosen_neg = neg_feature_matrix
  }
  
  
  balanced_feature_matrix = rbind(pos_feature_matrix, chosen_neg)
  
  balanced_feature_label = c(rep(1,pos_size),rep(0, nrow(chosen_neg)))
  
  
  
  
  saveRDS(balanced_feature_matrix, file = paste0(output_label, "_balanced_feature_matrix.Rds"))
  saveRDS(balanced_feature_label, file = paste0(output_label, "_balanced_feature_label.Rds"))
  
  
  write.table(balanced_feature_label, paste0(output_label, "_balanced_feature_label.tsv"), quote = F,
              sep = "\t", row.names = F)
  
}



# pos_neg_combination -----------------------------------------------------




combine_pos_neg = function(pos_feature_matrix_Rds_name, neg_feature_matrix_Rds_name,
                           output_label)
{
  
  pos_feature_matrix = readRDS(pos_feature_matrix_Rds_name)
  neg_feature_matrix = readRDS(neg_feature_matrix_Rds_name)
  
  
  ### rbind data 
  
  pos_neg_feature_matrix = rbind(pos_feature_matrix, neg_feature_matrix)
  
  
  
  col_mean =  colMeans(pos_neg_feature_matrix, na.rm = T)
  
  for(i in 1:ncol(pos_neg_feature_matrix))
  {
   this_na =  which(is.na(pos_neg_feature_matrix[,i]))
   pos_neg_feature_matrix[this_na,i] = col_mean[i]
   
  }
  
  
  saveRDS(pos_neg_feature_matrix, file = paste0(output_label, "_feature_matrix.Rds"))
  
  write.table(pos_neg_feature_matrix, paste0(output_label, "_feature_matrix.tsv"),
              row.names = F, col.names = F, sep = "\t", quote = F)
  
  rm(pos_neg_feature_matrix)
  
  
  pos_neg_feature_label = c(rep(1, nrow(pos_feature_matrix)), rep(0, nrow(neg_feature_matrix)))
  saveRDS(pos_neg_feature_label, file = paste0(output_label, "_feature_label.Rds"))
  write.table(pos_neg_feature_label, paste0(output_label, "_feature_label.tsv"),
              row.names = F, col.names = F, sep = "\t", quote = F)
  
  rm(pos_neg_feature_label)
  
}





# n_fold_cv_without_decoy -------------------------------------------------



construct_n_fold_cv_without_decoy = function(pos_candi_feature_matrix_Rds_name,pos_candi_feature_label_Rds_name,
                                             n, output_label)
{
  
  # pos_candi_feature_matrix_Rds_name = "glcnac_s_wp_feature_matrix.Rds"
  # pos_candi_feature_label_Rds_name = "glcnac_s_wp_feature_label.Rds"
  # n = 2
  # output_label = "tg"
  # 
  pos_candi_feature_matrix = readRDS(pos_candi_feature_matrix_Rds_name)
  pos_candi_feature_label = readRDS(pos_candi_feature_label_Rds_name)
  
  
  
  ### so the output should be n pairs of training/test sets
  
  ### for simplicity I want to use Caret package, it is not always good to build your own wheel
  
  

  pos_size = length(which(pos_candi_feature_label==1))
  pos_feature_matrix = pos_candi_feature_matrix[1:pos_size,]
  pos_seq = c(1:pos_size)
  
  candi_size = length(which(pos_candi_feature_label==0))
  candi_feature_matrix = pos_candi_feature_matrix[((pos_size+1):(pos_size+candi_size)),]
  candi_seq = c(1:candi_size)
  
  
  
  set.seed(123)
  pos_folds = caret::createFolds(pos_seq, n, FALSE )
  
  set.seed(123)
  candi_folds = caret::createFolds(candi_seq, n, FALSE)
  
  
  
  for(i in 1:n)
  {

    test_pos_ind = which(pos_folds==i)
    test_candi_ind = which(candi_folds==i)
    
    saveRDS(test_pos_ind, file = paste0(output_label,"_test_pos_ind_",i,".Rds"))
    
    saveRDS(test_candi_ind, file = paste0(output_label,"_test_candi_ind_",i,".Rds"))
    
    #cat(test_pos_ind,"\n")
    #cat(test_candi_ind,"\n")
    
    train_pos_ind = pos_seq[-test_pos_ind]
    train_candi_ind = candi_seq[-test_candi_ind]
    
    saveRDS(train_pos_ind, file = paste0(output_label,"_train_pos_ind_",i,".Rds"))
    
    saveRDS(train_candi_ind, file = paste0(output_label,"_train_candi_ind_",i,".Rds"))
    
    
    
    #cat(train_pos_ind,"\n")
    #cat(train_candi_ind,"\n")
    
    test_feature = rbind(pos_feature_matrix[test_pos_ind,],candi_feature_matrix[test_candi_ind,])
    
    # write.table(test_feature, paste0(output_label,"_test_feature_" ,i,".tsv"),
    #            row.names = F, col.names = F, sep = "\t", quote = F)
    saveRDS(test_feature, file = paste0(output_label,"_test_feature_" ,i,".Rds"))
    rm(test_feature)
    
    
    
    
    train_feature = rbind(pos_feature_matrix[train_pos_ind,],candi_feature_matrix[train_candi_ind,])
    #write.table(train_feature, paste0(output_label,"_train_feature_" ,i,".tsv"),
    #    row.names = F, col.names = F, sep = "\t", quote = F)
    
    saveRDS(train_feature, file = paste0(output_label,"_train_feature_" ,i,".Rds"))
    rm(train_feature)
    
    
    
    test_label = c(rep(1,length(test_pos_ind)),
                   rep(0, length(test_candi_ind)))
    write.table(test_label, paste0(output_label,"_test_label_" ,i,".tsv"),
                row.names = F, col.names = F, sep = "\t", quote = F)
    saveRDS(test_label, file = paste0(output_label,"_test_label_" ,i,".Rds"))
    
    
    
    
    train_label = c(rep(1,length(train_pos_ind)),
                    rep(0, length(train_candi_ind)))
    
    write.table(train_label, paste0(output_label,"_train_label_" ,i,".tsv"),
                row.names = F, col.names = F, sep = "\t", quote = F)
    
    saveRDS(train_label, file = paste0(output_label,"_train_label_" ,i,".Rds"))
    
    
    
  }
  
  
  
}




# scale the data ----------------------------------------------------------

scale_train_test = function(feature_train_name, feature_test_name, 
                            n_fold,
                            upper_bound, lower_bound,
                            output_label_fortrain, output_label_fortest)
{
  
  for(i in 1:n_fold)
  {
    feature_train_matrix = readRDS(paste0(feature_train_name, "_",i,".Rds"))
    
    feature_test_matrix = readRDS(paste0(feature_test_name, "_",i,".Rds"))
    

    
    get_train_range = scale_train(feature_train_matrix, upper_bound, lower_bound,
                                  paste0(output_label_fortrain,"_",i))
    
    scale_test(feature_test_matrix, get_train_range, upper_bound, lower_bound,
               paste0(output_label_fortest,"_",i))

  }
  
  
  
}



scale_train_test_single = function(feature_train_name, feature_test_name,
                                   upper_bound, lower_bound,
                                   output_label_fortrain, output_label_fortest)
{
  
  
  feature_train_matrix = readRDS(feature_train_name)
  
  feature_test_matrix = readRDS(feature_test_name)
  

  
  get_train_range = scale_train(feature_train_matrix, upper_bound, lower_bound,
                                output_label_fortrain)
  
  scale_test(feature_test_matrix, get_train_range, upper_bound, lower_bound,
             output_label_fortest)

  
}




scale_train_single = function(feature_train_name,
                                   upper_bound, lower_bound,
                                   output_label_fortrain)
{
  
  
  feature_train_matrix = readRDS(feature_train_name)
  
 it = scale_train(feature_train_matrix, upper_bound, lower_bound,
                                output_label_fortrain)
  
  
}




# libsvm formating --------------------------------------------------------



libsvm_formating_single = function(train_feature_name,
                                   train_label_name,
                                   test_feature_name,
                                   test_label_name,
                                   output_label)
{
  
  
  
  feature_train_out_Rds_name = paste0(output_label, "_train_feature")
  feature_train_out_tsv_name = paste0(output_label, "_train_feature")
  
  feature_test_out_Rds_name = paste0(output_label, "_test_feature")
  feature_test_out_tsv_name = paste0(output_label, "_test_feature")
  
  
  
  libsvm_formating(train_feature_name, train_label_name,
                   feature_train_out_Rds_name, feature_train_out_tsv_name)
  
  
  libsvm_formating(test_feature_name, test_label_name,
                   feature_test_out_Rds_name, feature_test_out_tsv_name)
  
  
  
  
  
}








# PCA plots ---------------------------------------------------------------




plot_plsda_for_two_types = function(pos_feature_Rds_name,
                                    candi_feature_Rds_name,
                                    output_label)
  
{
  
  
  # pos_feature_Rds_name = "py/py_nr_size25_nms_noc_not_na_pos_feature.Rds"
  # candi_feature_Rds_name = "py/py_nr_size25_nms_noc_not_na_candi_feature.Rds"
  # output_label  = "zpy"
  # output_label = "zpy"
  # 
  pos_feature = readRDS(pos_feature_Rds_name)
  candi_feature = readRDS(candi_feature_Rds_name)
  
  
  pos_size = nrow(pos_feature)
  candi_size = nrow(candi_feature)
  
  
  ### limit pos size to 5000
  if(pos_size>5000)
  {
    set.seed(123)
    sel_pos = pos_feature[sample(pos_size,5000),]
    
  }else{
    sel_pos = pos_feature
  }
  
  pos_size = nrow(sel_pos)
  
  if(candi_size>=pos_size)
  {
    set.seed(123)
    sel_candi = candi_feature[sample(candi_size,pos_size),]
    
  }else{
    sel_candi = candi_feature
  }
  
  
  
  sel_feature = rbind(sel_pos, sel_candi)
  
  type_label = as.factor(c(rep("positive",nrow(sel_pos)), 
                           rep("negative", nrow(sel_candi))))
  
  
  get_plsda_pos_candi(sel_feature, type_label, paste0(output_label,"_plsda_plot_two.pdf"))
  
  
  
}






plot_plsda_for_two_types_with_score_selection= function(pos_feature_Rds_name,
                                    candi_feature_Rds_name,
                                    pos_score_Rds_name,
                                    candi_score_Rds_name,
                                    candi_score_threshold,
                                    output_label)
  
{
  
  
  # pos_feature_Rds_name = "py/py_nr_size25_nms_noc_not_na_pos_feature.Rds"
  # candi_feature_Rds_name = "py/py_nr_size25_nms_noc_not_na_candi_feature.Rds"
  # output_label  = "zpy"
  # output_label = "zpy"
  # 
  pos_feature = readRDS(pos_feature_Rds_name)
  candi_feature = readRDS(candi_feature_Rds_name)
  
  pos_score = readRDS(pos_score_Rds_name)
  candi_score = readRDS(candi_score_Rds_name)
  
  pos_size = nrow(pos_feature)
  candi_size = nrow(candi_feature)
  
  
  ### limit pos size to 5000
  if(pos_size>5000)
  {
    set.seed(123)
    
    sel_pos_ind = sample(pos_size,5000)
    sel_pos = pos_feature[sel_pos_ind,]
    
    sel_pos_score = pos_score[sel_pos_ind,]
    
  }else{
    sel_pos = pos_feature
    sel_pos_score = pos_score
  }
  
  pos_size = nrow(sel_pos)
  
  if(candi_size>=pos_size)
  {
    set.seed(123)
    
    sel_candi_ind = sample(candi_size,pos_size)
    sel_candi = candi_feature[sel_candi_ind,]
    sel_candi_score = candi_score[sel_candi_ind,]
    
    
  }else{
    sel_candi = candi_feature
    sel_candi_score = candi_score
  }
  
  
  
  sel_feature = rbind(sel_pos, sel_candi)
  
  candi_score_label = rep("not_selected",nrow(sel_candi))
  
  candi_score_label[which(candi_score>= candi_score_threshold)] = "selected"
  
  score_label = as.factor(c(rep("selected", nrow(sel_pos)), candi_score_label))
  
  type_label = as.factor(c(rep("positive",nrow(sel_pos)), 
                           rep("negative", nrow(sel_candi))))

  
  get_plsda_pos_candi_with_score_selection(sel_feature, type_label, score_label, paste0(output_label,"_plsda_plot_two_score.pdf"))
  
  
  
}











# AUC_calculation ---------------------------------------------------------


# MCC calculation for prediction score threshold --------------------------


### modifiy this so that it can hold k fold rather than 2-fold only

# get_score_threshold_for_whole_proteome_mode -----------------------------



get_score_threshold = function(prediction_score_file_names,
                               test_label_file_names,
                               specificity_level,
                               output_label)
  
{
  # prediction_score_path = "/data/ginny/liblinear-2.11/test_package/ps_cv_predict/"
  # test_label_path = "/data/ginny/test_package/pred/ps_cv_predict/"
  # 
  # prediction_score_file_names = "ps_0103_predict.txt"
  # test_label_file_names = "ps_0103_test_label.txt"
  # output_label = "ps_0103"
  # specificity_level = 0.99
  
  score_files = readRDS(prediction_score_file_names)
  label_files = readRDS(test_label_file_names)
  
  ### need to think about how I can combine them together
  all_pred_df = data.frame(rbindlist(lapply(1:length(score_files), function(i) {
    this_predict = data.table::fread(score_files[i], stringsAsFactors = F, header = T)
    this_label = data.table::fread(label_files[i], stringsAsFactors = F, header = F)
    pred_df = cbind(this_predict, this_label)
    
    old_cn = colnames(pred_df)
    new_cn = c("pred_label","pos_score","neg_score","true_label")
    
    new_cn[which(old_cn == "labels")] = "pred_label"
    new_cn[which(old_cn == "1")] = "pos_score"
    new_cn[which(old_cn == "0")] = "neg_score"
    
    colnames(pred_df) = new_cn
    
    return(pred_df %>% dplyr::select(pos_score, true_label))
  })), stringsAsFactors = F)
  
  
  candi_score = all_pred_df %>%
    dplyr::filter(true_label == 0) %>%
    dplyr::select(pos_score)
  
  positive_score = all_pred_df %>%
    dplyr::filter(true_label == 1) %>%
    dplyr::select(pos_score)
  
  
  saveRDS(candi_score$pos_score, file = paste0(output_label, "_candi_score.Rds"))
  saveRDS(positive_score$pos_score, file = paste0(output_label, "_positive_score.Rds"))

  
  #### then output AUC and specificity and MCC etc
  
  total_label = all_pred_df$true_label
  total_score = all_pred_df$pos_score
  
  
  #### combine all the positive together and all the candidate together
  # roc_test = roc(total_label, total_score)
  # cat("AUC", roc_test$auc,"\n")
  # 
  record_mcc = rep(0,999)
  record_spec = rep(0,999)
  
  for(i in 1:999)
  {
    cutoff = i/1000
    
    tp = sum(total_label==1 & total_score>cutoff)
    tn = sum(total_label==0 & total_score<=cutoff)
    fp = sum(total_label==0 & total_score>cutoff)
    fn = sum(total_label==1 & total_score<=cutoff)
    
    # cat(cutoff, tp, tn, fp, fn, "\n")
    
    this_spec = tn/(tn+fp)
    this_mcc = calculate_MCC(tp,fp, fn, tn)
    
    # cat(this_mcc, "\n")
    
    record_mcc[i] = this_mcc
    record_spec[i] = this_spec
    
  }
  
  
  
  max_mcc = max(record_mcc, na.rm = T)
  score_max_mcc = which.max(record_mcc)/1000
  
  get_spec = abs(record_spec-specificity_level)
  score_which_spec = which.min(get_spec)/1000
  
  cat("at specificity level wanted, score cutoff is: ",
      score_which_spec,"\n")
  #cat("at specificity level wanted, how many sites are predicted? ",
   #   length(which(candi_score$pos_score>score_which_spec)), "\n")
  
  cat("threshold at best MCC ", score_max_mcc, "\n" )
  cat("best MCC", max_mcc,"\n")
  
  
  
  btp = sum(total_label==1 & total_score>score_max_mcc)
  btn = sum(total_label==0 & total_score<=score_max_mcc)
  bfp = sum(total_label==0 & total_score>score_max_mcc)
  bfn = sum(total_label==1 & total_score<=score_max_mcc)
  sens = btp/(btp+bfn)
  spec = btn/(btn+bfp)
  
  #cat("sens and spec at best MCC ", sens,"\t", spec, "\n" )
  #cat("how many candidate predicted: ", length(which(candi_score>score_max_mcc)), "\n")
  
  # pdf(paste0(output_label,"_roc.pdf"), useDingbats = F)
  # plot(roc_test)
  # dev.off()
  # 
  pdf(paste0(output_label,"_candidate_score_hist.pdf"), useDingbats = F)
  hist(candi_score$pos_score,breaks = 50, main = "predicted_score_on_candidate_sites")
  dev.off()
  
  pdf(paste0(output_label,"_positive_score_hist.pdf"), useDingbats = F)
  hist(positive_score$pos_score,breaks = 50, main = "predicted_score_on_positive_sites")
  dev.off()
  
  pdf(paste0(output_label,"_both_score_dens.pdf"), useDingbats = F)
  
  pd = density(positive_score$pos_score)
  cd = density(candi_score$pos_score)
  
  ymax = max(c(pd$y, cd$y))
  
  plot(cd, ylim = c(0,ymax),
       main = "candidate_pos_score", col = "blue")
  lines(pd, col = "red")
  dev.off()
  
  
  return(score_which_spec)
  
}





assemble_window_score_cv = function(candidate_df_Rds,
                                 positive_index_file_names,
                                 candi_index_file_names,
                                 positive_score_Rds,
                                 candi_score_Rds,
                                 score_threshold,
                                 id_convert,
                                 output_label)
  
{
  # candidate_df_Rds = "ps_wp_52_candidate.Rds"
  # positive_index_file_names = "ps_wp_52_pos_index_names.Rds"
  # candi_index_file_names = "ps_wp_52_candi_index_names.Rds"
  # positive_score_Rds = "ps_wp_52_positive_score.Rds"
  # candi_score_Rds = "ps_wp_52_candi_score.Rds"
  # score_threshold = 0.683
  # output_label = "test_0125"

  candidate_df = readRDS(candidate_df_Rds)
  
  ### now get the order correct
  
  positive_score = readRDS(positive_score_Rds)
  candi_score = readRDS(candi_score_Rds)
  
  
  #### tidy up the order of index and score and combine them together
  ### only need to look at test indices
  
  positive_df = candidate_df%>%dplyr::filter(label == "positive")
  candi_df = candidate_df%>%dplyr::filter(label == "negative")
  
  ###start from here tomorrow
  
  ### get the order of the index and the order of score in the same way
  
  all_positive_ind = readRDS(positive_index_file_names)
  all_candi_ind = readRDS(candi_index_file_names)

  
  tidy_positive_ind = data.frame(rbindlist(
    lapply( 1:length(all_positive_ind), function(i) {
      this_positive_ind = readRDS(all_positive_ind[i])
      return(data.frame(ind = this_positive_ind, stringsAsFactors = F))
    })), stringsAsFactors = F)
  
  positive_df_order = positive_df[tidy_positive_ind$ind,]
  
  positive_df_order_score =  data.frame(positive_df_order, pred_score = positive_score)
  rm(positive_df_order)
  
  
  tidy_candi_ind = data.frame(rbindlist(
    lapply( 1:length(all_candi_ind), function(i) {
      this_candi_ind = readRDS(all_candi_ind[i])
      return(data.frame(ind = this_candi_ind, stringsAsFactors = F))
    })), stringsAsFactors = F)
  
  candi_df_order = candi_df[tidy_candi_ind$ind,]
  
  candi_df_order_score =  data.frame(candi_df_order, pred_score = candi_score)
  
  rm(candi_df_order)
  
  df_score = rbind(positive_df_order_score, candi_df_order_score)
  
  rm(positive_df_order_score, candi_df_order_score)
  
  colnames(id_convert) = c("protID","gene_name")
  
  
  df_score_label = df_score %>%
    dplyr::mutate(pred_label = "negative") %>%
    dplyr::mutate(pred_label = replace(pred_label, pred_score >= score_threshold, "positive")) %>%
    dplyr::mutate(prediction_label = pred_label) %>%
    dplyr::mutate(combined_label = replace(pred_label, label == "positive", "positive")) %>%
    dplyr::mutate(known_label = label)%>%
    dplyr::mutate(threshold = score_threshold)%>%
    dplyr::left_join(id_convert) %>%
    dplyr::arrange(protID, pos) %>%
    dplyr::select(protID, gene_name, pos, window, pred_score, threshold,prediction_label, known_label, combined_label)
  
  rm(df_score)
  
  
  write.table(df_score_label,paste0(output_label, "_window_score_df.tsv"), 
              sep = "\t", quote = F, row.names = F, na = "")
  
  saveRDS(df_score_label, file = paste0(output_label, "_window_score_df.Rds"))
  
  
  
}







assemble_window_score_target = function(prediction_score_file,
                                        predict_candidate_df_Rds,
                                        id_convert,
                                        score_threshold,
                                        output_label)
{
  #### simply combine
  
  # 
  # prediction_score_file = "/data/ginny/liblinear-2.11/test_package/ps_predict_test_predict.tsv"
  # 
  # predict_candidate_df_Rds = "ps_predict_candidate.Rds"
  # 
  # score_threshold = 0.683
  # 
  # output_label = "ps_predict"
  # 
  # 
  pred_score_df = data.table::fread(prediction_score_file, stringsAsFactors = F, header = T)
  
  
  ### I think it is safer to code it in a matching manner
  old_cn = colnames(pred_score_df)
  new_cn = c("pred_label","pos_score","neg_score")
  
  new_cn[which(old_cn == "labels")] = "pred_label"
  new_cn[which(old_cn == "1")] = "pos_score"
  new_cn[which(old_cn == "0")] = "neg_score"
  
  colnames(pred_score_df) = new_cn
  
  
  pred_df = readRDS(predict_candidate_df_Rds)
  
  
  
  
  colnames(id_convert) = c("protID","gene_name")
  
  
  pred_df_score_label = pred_df %>%
    dplyr::mutate(pred_score = pred_score_df$pos_score) %>%
    dplyr::mutate(pred_label = "negative") %>%
    dplyr::mutate(pred_label = replace(pred_label, pred_score >= score_threshold, "positive")) %>%
    dplyr::mutate(prediction_label = pred_label) %>%
    dplyr::mutate(combined_label = prediction_label)%>%
    dplyr::mutate(known_label = label)%>%
    dplyr::mutate(threshold = score_threshold) %>%
    dplyr::left_join(id_convert) %>%
    dplyr::arrange(protID, pos) %>%
    dplyr::select(protID, gene_name, pos, window, pred_score, threshold, prediction_label, known_label, combined_label)
  

  
  write.table(pred_df_score_label,paste0(output_label, "_window_score_df.tsv"), 
              sep = "\t", quote = F, row.names = F, na = "")
  
  saveRDS(pred_df_score_label, file = paste0(output_label, "_window_score_df.Rds"))
  
  
  
  
}










# select_before_feature_extraction ----------------------------------------

select_equal_size_candidate_decoy = function(candidate_Rds_name, output_label)
{
  candidate_df = readRDS(candidate_Rds_name)
  

  
  
  
  candi_df <- candidate_df%>%dplyr::filter(label == "negative")
  
  candi_nrow = nrow(candi_df)
  
  
  pos_df <- candidate_df%>%dplyr::filter(label == "positive")
  pos_nrow = nrow(pos_df)
  
  

  set.seed(123)
  
  choose_candi = sample(candi_nrow, pos_nrow)

  
  choose_candi_df = candi_df[choose_candi,]

  
  choose_candidate_df = rbind(pos_df, choose_candi_df)
  
  write.table(choose_candidate_df, paste0(output_label, "_candidate.tsv"),
              quote = F, row.names = F, sep = "\t")
  saveRDS(choose_candidate_df, paste0(output_label, "_candidate.Rds"))
  
  
  
  
  
}





### use R to control the terminal



extract_site_specific_features_new=function(feature_data, to_extract, to_logit, center_position)
{
  
  
  #feature_data = half_spider_matrix1
  
  # set.seed(123)
  # feature_data = matrix(abs(rnorm(10*250)),nrow = 10, ncol = 250)
  #center_position = 13
  
  #to_extract = c(1,6,8,9)
  #to_logit = c(8,9)
  
  
  logit_from_extract = which(to_extract %in%  to_logit)
  
  ##### extract first 
  total_feature_length = 2*(center_position-1)+1
  
  
  extractbs = matrix(rep(seq(0,  10*(total_feature_length-1),10),length(to_extract)),
                     nrow=total_feature_length,
                     ncol=length(to_extract))
  
  add_extractbs = sapply(1:length(to_extract), function(x) extractbs[,x]+to_extract[x])
  vec_add_extractbs = sort(as.vector(add_extractbs))
  
  
  extract_feature_data = feature_data[,vec_add_extractbs, with = F]
  
  rm(feature_data)
  
  ##### log second
  
  ef = length(to_extract)
  
  logbs = matrix(rep(seq(0,  ef*(total_feature_length-1),ef),length(logit_from_extract)),
                 nrow=total_feature_length,
                 ncol=length(logit_from_extract))
  
  add_logbs = sapply(1:length(logit_from_extract), function(x) logbs[,x]+logit_from_extract[x])
  vec_add_logbs = sort(as.vector(add_logbs))
  
  
  ### arrange a matrix to get all the columns need to be loggit
  
  extract_feature_data = as.matrix(extract_feature_data)
  tar = extract_feature_data[, vec_add_logbs]
  
  tar[which(tar<0.001)]=0.001
  tar[which(tar>0.999)]=0.999
  
  
  logitit = function(p1){return(log(p1/(1-p1)))} 
  
  logit_tar = logitit(tar)
  
  ### place back these columns to the orignal data 
  
  extract_feature_data[,vec_add_logbs] = logit_tar 
  
  return(extract_feature_data)
  
}




# window score assembly ----------------------------------------------


#### ok the following two functions are still in the test procedure


# prediction_annotation ---------------------------------------------------



domain_subcellular_mapping = function(pos_window_score_Rds,
                                      candi_window_score_Rds, 
                                      dm_df_Rds, 
                                      sc_df_Rds,
                                      output_label)
{


  # pos_window_score_Rds = "glcnac_s_pred/glcnac_s_pred_pos_window_score.Rds"
  # candi_window_score_Rds = "glcnac_s_pred/glcnac_s_pred_candi_window_score.Rds"
  # 
  # dm_df_Rds = "domain_df_pure.Rds"
  # sc_df_Rds = "subcellular_location_df_pure.Rds"
  # 
  # output_label = "glcnac_s_try"
  # output_label = "glcnac_s_try"
  # 



  
  dm_df = readRDS(dm_df_Rds)
  sc_df = readRDS(sc_df_Rds)
  
  
  pos_window_score = readRDS(pos_window_score_Rds)
  
  
  pos_each_domain = map_domain(dm_df, pos_window_score)
  pos_each_subcellular = map_subcellular_location(sc_df, pos_window_score)
  
  
  pos_info_df = data.frame(pos_window_score, 
                           domain = pos_each_domain,
                           subcellular = pos_each_subcellular,
                           stringsAsFactors = F)
  
  saveRDS(pos_info_df, file = paste0(output_label, "_pos_info_df.Rds"))
  write.table(pos_info_df, paste0(output_label,"_pos_info_df.tsv"),
              quote =  F, row.names = F, sep = "\t")
  
  
  rm(pos_info_df)
  rm(pos_window_score)
  rm(pos_each_domain)
  rm(pos_each_subcellular)
  
  candi_window_score = readRDS(candi_window_score_Rds)
  
  
  candi_each_domain = map_domain(dm_df, candi_window_score)
  candi_each_subcellular = map_subcellular_location(sc_df, candi_window_score)
  
  
  candi_info_df = data.frame(candi_window_score, 
                             domain = candi_each_domain,
                             subcellular = candi_each_subcellular,
                             stringsAsFactors = F)
  
  saveRDS(candi_info_df, file = paste0(output_label, "_candi_info_df.Rds"))
  write.table(candi_info_df, paste0(output_label,"_candi_info_df.tsv"),
              quote =  F, row.names = F, sep = "\t")
  
  
  
}




retrieve_domain_all_mod = function(mod_names, mod_pos_score, mod_candi_score, domain_name,
output_label)
{
    
    ### try to use rbindlist
    
    all_retrieve = data.frame(rbindlist(lapply(1:length(mod_names), function(i) {
        
        mod_name = mod_names[i]
        
        pos_info = readRDS(paste0(mod_name, "_pred_pos_info_df.Rds"))
        candi_info = readRDS(paste0(mod_name, "_pred_candi_info_df.Rds"))
        
        pos_score = mod_pos_score[i]
        candi_score = mod_candi_score[i]
        
        
        
        this_mod_retrieve = retrieve_for_each_mod(pos_info, candi_info, pos_score,candi_score,
        mod_name, domain_name)
        
        return(this_mod_retrieve)
    })))
    
    
    all_retrieve_df = all_retrieve %>%
    dplyr::arrange(protID, pos)
    
    
    saveRDS(all_retrieve_df, file = paste0(output_label, "_",domain_name,"_retrieve.Rds"))
    write.table(all_retrieve_df, paste0(output_label, "_",domain_name,"_retrieve.tsv"),
    quote = F, row.names= F, sep = "\t")
    
    
    
}







create_gglogo_plot = function(candidate_df_Rds, 
                              output_label)
{
  
  
  candidate_df = readRDS(candidate_df_Rds)
#  candidate_df = ps_candidate
#  output_label = "try_new"
  
  pos_windows = candidate_df %>% 
    dplyr::filter(label == "positive") %>%
    dplyr::select(window)
  delete_center = paste0(substr(pos_windows$window,1,12),substr(pos_windows$window,14,25))
  
  pos_windows$window = delete_center
  
  
  
 # pdf(paste0(output_label,"_logoPlot.pdf"), useDingbats = F)
  
  ggplot(data = ggfortify(pos_windows, "window", method = "frequency")) +      
    geom_logo(aes(x=position, y=info, group=element, 
                  label=element, fill=interaction(Water, Polarity)),
              alpha = 0.6)  +
    scale_fill_brewer(palette="Paired") +
    theme(legend.position = "bottom")
  
  #dev.off()
  
  ggsave(filename = paste0(output_label,"_logoPlot.pdf"),device = "pdf",
         width = 10, height = 8)
  
  
}





  
get_average_weights_of_two = function(first_train_model_file, second_train_model_file,
                                      full_feature, arrange_feature, output_label)
  
{
  weight_matrix = matrix(0, nrow = length(full_feature), ncol = 2)


    first_weight_file = readLines(first_train_model_file)
    
    first_ws = as.numeric(first_weight_file[7:length(first_weight_file)])
    
    weight_matrix[,1] = first_ws
    
    
    second_weight_file = readLines(second_train_model_file)
    
    second_ws = as.numeric(second_weight_file[7:length(second_weight_file)])
    
    weight_matrix[,2] = second_ws

    
    
  
  ave_weights = rowMeans(weight_matrix)
  
  ws_df = data.frame(name = full_feature,
                     weights = ave_weights,
                     abs_weights = abs(ave_weights), stringsAsFactors = F)
  
  name_arr_ws_df = arrange_feature%>%dplyr::left_join(ws_df, by = c("arrange_feature" ="name"))
  
  
  write.table(name_arr_ws_df,paste0(output_label,"_arrange_name_coeff_df.tsv"),
              sep = "\t", row.names = F, quote = F)

}







plot_weights = function(coef_file, plot_name)
{
  
  #coef_file = "/Users/ginny/PTMtopographer_2017_08/classifier_performance/model_weight/ps_arrange_name_coeff_df.tsv"
  
  coef_df = data.table::fread(coef_file, header = T, stringsAsFactors = F)
  
  col_hydrophobicity = rep("skyblue", 8)
  col_aaindex = rep("purple", 45)
  col_ASA = rep("green",24 )
  col_HSE = rep("yellow",24 )
  col_pC = rep("red",24 )
  col_pH = rep("orange",24 )
  col_pssm = rep("pink",24)
  
  cols = c(col_hydrophobicity, col_aaindex, col_ASA, col_HSE, col_pC, col_pH, col_pssm)
  
  
  #  max_y = max(abs(coef_df$weights))
  
  
  pdf(paste0(plot_name,"_weights_bar.pdf"), useDingbats = F)
  
  barplot(coef_df$weights, main = plot_name, col = cols, ylim = c(-0.55, 0.55))
  # abline (h = 0, col = "green", lty = 2)
  # abline(v= 8, col = "red", lty = 2)
  # abline(v= 53, col = "red", lty = 2)
  # abline(v= 101, col = "red", lty = 2)
  # abline(v= 149, col = "red", lty = 2)
  # 
  dev.off()
  
  
}








calculate_seq_pairs = function(anchor_mod_site, cross_mod_site, distance,
                               anchor_mod, cross_mod,
                               anchor_new_merge_Rds, cross_new_merge_Rds,
                               output_label)
{
  anchor_new_merge = readRDS(anchor_new_merge_Rds)
  cross_new_merge = readRDS(cross_new_merge_Rds)

  cn_anchor = colnames(anchor_new_merge)
  cn_anchor[which(grepl(paste0(anchor_mod, "_label"), cn_anchor))] = "anchor_label"
  colnames(anchor_new_merge) = cn_anchor
  
  cn_cross = colnames(cross_new_merge)
  cn_cross[which(grepl(paste0(cross_mod, "_label"), cn_cross))] = "cross_label"
  colnames(cross_new_merge) = cn_cross
  
  anchor_domain_list =unique(anchor_new_merge$domain)
  anchor_domain_list = unique( unlist(strsplit(anchor_domain_list, split = " ")))
  cross_domain_list = unique(cross_new_merge$domain)
  cross_domain_list = unique( unlist(strsplit(cross_domain_list, split = " ")))
  common_domain = intersect(anchor_domain_list, cross_domain_list)
  common_domain = common_domain[!is.na(common_domain)]
  
  output_df = data.frame(domain_name = common_domain, a = 0, b=0, c=0, d=0)
  
  
  domain_prot_df = data.frame(domain_name = common_domain, 
                              protIDs = character(length(common_domain)),
                              stringsAsFactors = F)
  
  for(i in 1:length(common_domain))
  {
    
    
    if(i%%100==0)
      cat(i, "\n")
    
    domain_prot_df$protIDs[i] = NA
    #which(common_domain == "CTP_synth_N")
    
    a =0;b=0;c=0;d=0
    ### because some domains are connected by blanks in a row
    
    this_domain = common_domain[i]
    this_anchor_domain = anchor_new_merge %>%
      dplyr:: filter(grepl(paste0("\\b",this_domain,"\\b") ,domain))
    this_cross_domain = cross_new_merge %>%
      dplyr:: filter(grepl(paste0("\\b",this_domain,"\\b"), domain))
    
    
    if(nrow(this_anchor_domain)>0 & nrow(this_cross_domain)>0)
    {
      ### crosstalk happen within the same protein and the same domain
      
      this_anchor_proteins = unique(this_anchor_domain$protID)
      this_cross_proteins = unique(this_cross_domain$protID)
      common_proteins = intersect(this_anchor_proteins, this_cross_proteins)
      
      if(length(common_proteins)>0)
      {
        domain_crosstalk_in_prot = c("")
        
        for(j in 1:length(common_proteins))
        {
          this_prot_anchor_domain = this_anchor_domain %>%
            filter(protID == common_proteins[j])
          this_prot_cross_domain = this_cross_domain %>%
            filter(protID == common_proteins[j])
          
          ### get the pairs of anchor cross positions for analysis
          anchor_cross_pair_pos =
            data.frame(anchor_position = numeric(), cross_position = numeric())
          
          anchor_pos = this_prot_anchor_domain %>%
            dplyr::select(pos)
          anchor_pos = anchor_pos$pos
          
          cross_pos = this_prot_cross_domain %>%
            dplyr::select(pos)
          cross_pos = cross_pos$pos
          
          for(p in 1:length(anchor_pos))
          {
            dis = abs(anchor_pos[p]-cross_pos)
            
            find_cross = which(dis>0 & dis<=distance)
            
            if(length(find_cross)>0)
            {
              anchor_position = rep(anchor_pos[p], length(find_cross))
              cross_position = cross_pos[find_cross]
              
              pos_pair = cbind(anchor_position, cross_position)
              anchor_cross_pair_pos = rbind(anchor_cross_pair_pos, pos_pair)
              
            }
            
          }
              
          if(nrow(anchor_cross_pair_pos)>0)
          {
            
            get_anchor_match = match(anchor_cross_pair_pos$anchor_position,
                                     this_prot_anchor_domain$pos)
            get_cross_match = match(anchor_cross_pair_pos$cross_position,
                                    this_prot_cross_domain$pos)
            
            anchor_match_label = this_prot_anchor_domain$anchor_label[get_anchor_match]
            cross_match_label = this_prot_cross_domain$cross_label[get_cross_match]
            
            add_to_a = which(anchor_match_label == T & cross_match_label == T)
            add_to_b = which(anchor_match_label == F & cross_match_label == T)
            add_to_c =  which(anchor_match_label == T & cross_match_label == F)
            add_to_d =  which(anchor_match_label == F & cross_match_label == F)
            a = a+length(add_to_a)
            b = b+length(add_to_b)
            c = c+length(add_to_c)
            d = d+length(add_to_d)
            
            if(length(add_to_a)>0)
            {
              domain_crosstalk_in_prot = c(domain_crosstalk_in_prot, common_proteins[j])
              domain_prot_df$protIDs[i] = paste(domain_crosstalk_in_prot, collapse = "    ")
            }
            
          }
        }
      }
      
    }
    
    output_df$a[i] = a
    output_df$b[i] = b
    output_df$c[i] = c
    output_df$d[i] = d
    
  }
  saveRDS(output_df, file = paste0(output_label,"_tbt_table.Rds"))
  write.table(output_df, paste0(output_label, "_tbt_table.tsv"), sep = "\t",
              quote = F, row.names = F)
  
  saveRDS(domain_prot_df, file = paste0(output_label, "_domain_prot_match.Rds"))
  write.table(domain_prot_df, paste0(output_label, "_domain_prot_match.tsv"), sep = "\t",
              quote = F, row.names = F)
  
  get_tbt = output_df %>%
    dplyr::group_by(domain_name) %>%
    dplyr::summarise(both_positive = sum(a), cross_positive = sum(b), anchor_positive = sum(c), both_negative = sum(d))

  colnames(get_tbt) = c("domain", "both_positive",paste0(cross_mod, "_positive"),
                        paste0(anchor_mod,"_positive"), "both_negative")


  ### calculate fisher's exact test on this
  fisher_p = rep(0, nrow(get_tbt))
  or = rep(0, nrow(get_tbt))

  for(i in 1:nrow(get_tbt))
  {

    fm = matrix(as.numeric(get_tbt[i,c(2:5)]), nrow = 2, ncol = 2, byrow = T)

    ft =  fisher.test(fm, alternative = "g")

    fisher_p[i] = ft$p.value
    or[i] = ft$estimate
  }

  tbt_p = get_tbt %>%
    dplyr::mutate(fisher_pvalue = fisher_p)%>%
    dplyr::mutate(odds_ratio = or)%>%
    dplyr::arrange(fisher_pvalue)

  write.table(tbt_p, paste0(output_label, "_test.tsv"),
              quote = F, row.names = F, sep = "\t")

}




calculate_seq_pairs_negative = function(compete_mod_site,
                                        anchor_mod, cross_mod,
                                        new_merge_Rds,
                                        output_label)
{
  
  new_merge = readRDS(new_merge_Rds)
 
  cn_merge = colnames(new_merge)
  cn_merge[which(grepl(paste0(anchor_mod, "_label"), cn_merge))] = "anchor_label"
  cn_merge[which(grepl(paste0(cross_mod, "_label"), cn_merge))] = "cross_label"
  
  colnames(new_merge) = cn_merge
  
  
  domain_list =unique(new_merge$domain)
  domain_list = unique( unlist(strsplit(domain_list, split = " ")))
  
  common_domain = domain_list[!is.na(domain_list)]
  
  output_df = data.frame(domain_name = common_domain, a = 0, b=0, c=0, d=0)
  
  
  domain_prot_df = data.frame(domain_name = common_domain, 
                              protIDs = character(length(common_domain)),
                              stringsAsFactors = F)
  
  for(i in 1:length(common_domain))
  {
    
    if(i%%100==0)
      cat(i, "\n")
    domain_prot_df$protIDs[i] = NA
    
    a =0;b=0;c=0;d=0
    
    this_domain = common_domain[i]
    this_domain_df = new_merge %>%
      dplyr:: filter(grepl(paste0("\\b",this_domain,"\\b") ,domain))
    
    if(nrow(this_domain_df)>0 )
    {
      
      anchor_match_label = this_domain_df$anchor_label
      cross_match_label = this_domain_df$cross_label
      
      a = length(which(anchor_match_label == T & cross_match_label == T))
      b = length(which(anchor_match_label == F & cross_match_label == T))
      c = length(which(anchor_match_label == T & cross_match_label == F))
      d = length(which(anchor_match_label == F & cross_match_label == F))
    
      have_cross = this_domain_df %>%
        dplyr::filter(anchor_label == T, cross_label == T) %>%
        dplyr::select(protID)
      
      if(nrow(have_cross)>0)
        domain_prot_df$protIDs[i] = paste(unique(have_cross$protID),
                                          collapse = "    ") 
      
    }
    output_df$a[i] = a
    output_df$b[i] = b
    output_df$c[i] = c
    output_df$d[i] = d
    
  }
  
  saveRDS(domain_prot_df, file = paste0(output_label, "_domain_prot_match.Rds"))
  write.table(domain_prot_df, paste0(output_label, "_domain_prot_match.tsv"), sep = "\t",
              quote = F, row.names = F)
  
  saveRDS(output_df, file = paste0(output_label,"_tbt_table.Rds"))
  write.table(output_df, paste0(output_label, "_tbt_table.tsv"), sep = "\t",
              quote = F, row.names = F)

  get_tbt = output_df %>%
    dplyr::group_by(domain_name) %>%
    dplyr::summarise(both_positive = sum(a), cross_positive = sum(b), anchor_positive = sum(c), both_negative = sum(d))
  
  colnames(get_tbt) = c("domain", "both_positive",paste0(cross_mod, "_positive"),
                        paste0(anchor_mod,"_positive"), "both_negative")
  
  ### calculate fisher's exact test on this
  fisher_p = rep(0, nrow(get_tbt))
  or = rep(0, nrow(get_tbt))
  
  for(i in 1:nrow(get_tbt))
  {
    
    fm = matrix(as.numeric(get_tbt[i,c(2:5)]), nrow = 2, ncol = 2, byrow = T)
    
    ft =  fisher.test(fm, alternative = "g")
    
    fisher_p[i] = ft$p.value
    or[i] = ft$estimate
  }
  
  tbt_p = get_tbt %>%
    dplyr::mutate(fisher_pvalue = fisher_p)%>%
    dplyr::mutate(odds_ratio = or)%>%
    dplyr::arrange(fisher_pvalue)
  
  write.table(tbt_p, paste0(output_label, "_test.tsv"),
              quote = F, row.names = F, sep = "\t")

}



calculate_individual_ptm = function(mod_type, 
                                    new_merge_Rds, 
                                   output_label)
{
  new_merge = readRDS(new_merge_Rds)
  
  cn_merge = colnames(new_merge)
  cn_merge[which(grepl(paste0(mod_type, "_label"), cn_merge))] = "anchor_label"
  colnames(new_merge) = cn_merge

  num_positive_ptm = sum(new_merge$anchor_label)
  
  domain_list = unique(new_merge$domain)
  domain_list = domain_list[!is.na(domain_list)]
  domain_list = unique( unlist(strsplit(domain_list, split = " ")))
  
  output_df = data.frame(domain = domain_list, a = 0, b=0, c=0, d=0)
  
  domain_prot_df = data.frame(domain_name = domain_list, 
                              protIDs = character(length(domain_list)),
                              stringsAsFactors = F)
  
  for(i in 1:length(domain_list))
  {
    
    if(i%%1000 == 0)
      cat(i, "\n")
    domain_prot_df$protIDs[i] = NA
  
    this_domain_rows = new_merge %>% 
      dplyr::filter(grepl(paste0("\\b",domain_list[i],"\\b"), domain))

    have_ptm_domain = this_domain_rows %>%
      dplyr::filter(anchor_label == T)
    
    if(nrow(have_ptm_domain)>0)
      domain_prot_df$protIDs[i] = paste(unique(have_ptm_domain$protID),
                                        collapse = "    ") 
    
    a= sum(this_domain_rows$anchor_label)
    
    b = num_positive_ptm - a 
    
    c = nrow(this_domain_rows) - a
    
    d = nrow(new_merge) - a - b - c
    
    output_df$a[i] = a
    output_df$b[i] = b
    output_df$c[i] = c
    output_df$d[i] = d
    
    
  }
  
  saveRDS(domain_prot_df, file = paste0(output_label, "_domain_prot_match.Rds"))
  write.table(domain_prot_df, paste0(output_label, "_domain_prot_match.tsv"), sep = "\t",
              quote = F, row.names = F)
  
  colnames(output_df) = c("domain", "InDomain_positive","OutDomain_positive",
                          "InDomain_negative", "OutDomain_negative")
  
  saveRDS(output_df, file = paste0(output_label,"_tbt_table.Rds"))
  write.table(output_df, paste0(output_label, "_tbt_table.tsv"), sep = "\t",
              quote = F, row.names = F)
  
  chisq_p = rep(0, nrow(output_df))
  or = rep(0, nrow(output_df))
  
  for(i in 1:nrow(output_df))
  {
    if(i%%100 == 0)
      cat(i, "\n")
    
    fm = matrix(as.numeric(output_df[i,c(2:5)]), nrow = 2, ncol = 2, byrow = T)
    
    ft =  chisq.test(fm, simulate.p.value = T)
    
    chisq_p[i] = ft$p.value
    or[i] = output_df[i,2]*output_df[i,5]/output_df[i,3]/output_df[i,4]
    
  }
  
  tbt_p = output_df %>%
    dplyr::mutate(chisq_pvalue = chisq_p)%>%
    dplyr::mutate(odds_ratio = or)%>%
    dplyr::arrange(chisq_pvalue)
  

  write.table(tbt_p, paste0(output_label, "_test.tsv"),
              quote = F, row.names = F, sep = "\t")
  
}





map_domain_subcellular = function(window_score_label_Rds,
                                  dm_df_Rds, 
                                  sc_df_Rds,
                                  output_label)
{
  
  
  # window_score_label_Rds = "ps_0103_window_score_df.Rds"
  # dm_df_Rds = "domain_df_pure.Rds"
  # sc_df_Rds = "subcellular_location_df_pure.Rds"
  # output_label = "ps_0103"
  
  
  dm_df = readRDS(dm_df_Rds)
  sc_df = readRDS(sc_df_Rds)
  
  
  window_score_label = readRDS(window_score_label_Rds)
  
  
  each_domain = map_domain(dm_df, window_score_label)
  each_subcellular = map_subcellular_location(sc_df, window_score_label)
  
  
  info_df = data.frame(window_score_label, 
                       domain = each_domain,
                       subcellular = each_subcellular,
                       stringsAsFactors = F)
  
  saveRDS(info_df, file = paste0(output_label, "_mapped_df.Rds"))
  write.table(info_df, paste0(output_label,"_mapped_df.tsv"),
              quote =  F, row.names = F, sep = "\t",na = "")
  
}




# chisq_test_for_single_ptm -----------------------------------------------



calculate_tbt_single_ptm = function(mapped_window_score_label_Rds, 
                                    output_label)
{  
  #mapped_window_score_label_Rds = "ps_0103_mapped_df.Rds"
  mapped_window_score_label = readRDS(mapped_window_score_label_Rds)
  
  
  
  colnames(mapped_window_score_label) = c("protID","gene_name","pos","window","pred_score","threshold",
                                 "prediction_label","known_label","combined_label",
                                 "domain","subcellular")
  
  cn_merge = colnames(mapped_window_score_label)
  cn_merge[which(grepl("combined_label", cn_merge))] = "anchor_label"
  colnames(mapped_window_score_label) = cn_merge
  
  ### recode positive negative to TURE and FALSE
  
  mapped_window_score_label = mapped_window_score_label %>%
    dplyr::mutate(anchor_label = replace(anchor_label, anchor_label == "positive", TRUE)) %>%
    dplyr::mutate(anchor_label = as.logical(replace(anchor_label, anchor_label == "negative", FALSE)))
  
  
  num_positive_ptm = sum(mapped_window_score_label$anchor_label)
  
  
  domain_list = unique(mapped_window_score_label$domain)
  domain_list = domain_list[!is.na(domain_list)]
  domain_list = unique(unlist(strsplit(domain_list, split = " ")))
  
  output_df = data.frame(domain = domain_list, a = 0, b=0, c=0, d=0)
  
  
  
  domain_prot_df = data.frame(domain_name = domain_list, 
                              protIDs = character(length(domain_list)),
                              stringsAsFactors = F)
  
  
  
  
  for(i in 1:length(domain_list))
  {
    
    # if(i%%1000 == 0)
    #   cat(i, "\n")
    # 
    
    domain_prot_df$protIDs[i] = NA
    
    
    this_domain_rows = mapped_window_score_label %>% 
      dplyr::filter(grepl(paste0("\\b",domain_list[i],"\\b"), domain))
    
    
    have_ptm_domain = this_domain_rows %>%
      dplyr::filter(anchor_label == T)
    
    if(nrow(have_ptm_domain)>0)
      domain_prot_df$protIDs[i] = paste(unique(have_ptm_domain$protID),
                                        collapse = "    ") 
    
    
    a= sum(this_domain_rows$anchor_label)
    
    b = num_positive_ptm - a 
    
    c = nrow(this_domain_rows) - a
    
    d = nrow(mapped_window_score_label) - a - b - c
    
    
    output_df$a[i] = a
    output_df$b[i] = b
    output_df$c[i] = c
    output_df$d[i] = d
    
    
  }
  
  
  
  saveRDS(domain_prot_df, file = paste0(output_label, "_domain_prot_match.Rds"))
  write.table(domain_prot_df, paste0(output_label, "_domain_prot_match.tsv"), sep = "\t",
              quote = F, row.names = F)
  
  
  
  colnames(output_df) = c("domain", "InDomain_positive","OutDomain_positive",
                          "InDomain_negative", "OutDomain_negative")
  
  saveRDS(output_df, file = paste0(output_label,"_tbt_table.Rds"))
  write.table(output_df, paste0(output_label, "_tbt_table.tsv"), sep = "\t",
              quote = F, row.names = F)
  
  chisq_p = rep(0, nrow(output_df))
  or = rep(0, nrow(output_df))
  adjusted_or = rep(0,nrow(output_df))
  for(i in 1:nrow(output_df))
  {
    # if(i%%100 == 0)
    #   cat(i, "\n")
    
    fm = matrix(as.numeric(output_df[i,c(2:5)]), nrow = 2, ncol = 2, byrow = T)
    
    set.seed(123)
    ft =  chisq.test(fm, simulate.p.value = T)
    
    chisq_p[i] = ft$p.value
    or[i] = output_df[i,2]*output_df[i,5]/output_df[i,3]/output_df[i,4]
    adjusted_or[i] = (output_df[i,2]+0.5)*(output_df[i,5]+0.5)/(output_df[i,3]+0.5)/(output_df[i,4]+0.5)
    
  }
  
 # qv = qvalue(p = chisq_p)
  
  
  tbt_p = output_df %>%
    dplyr::mutate(chisq_pvalue = chisq_p)%>%
    dplyr::mutate(odds_ratio = or)%>%
    #dplyr::mutate(qvalue = qv$qvalues)%>%
    dplyr::mutate(adjusted_odds_ratio = adjusted_or)%>%
    dplyr::arrange(chisq_pvalue)
  
  write.table(tbt_p, paste0(output_label, "_test.tsv"),
              quote = F, row.names = F, sep = "\t")
  
}


