
#' Generate feature data for whole proteome
#' 
#' This function generates feature data for whole proteome by mapping the positive PTM info to protein sequences, constructing windows and extracting 3 sets of features.
#' @param ptm_site The amino acid this PTM involves, in upper-case single letter representation.
#' @param flanking_size The number of residues surrounding each side of the center residue, the total window size will be 2*flanking_size+1, default to 12.
#' @param SPIDER A boolean variable indicating the usage of SPIDER3 features, default set to TRUE.
#' @param positive_info_file A text file containing the positive PTM sites info in required format.
#' @param protein_fasta_file A text file containing the proteins sequences of interest in Fasta format.
#' @param output_label The string to tag the output files.
#' @import stringr dplyr magrittr data.table
#' @export
#' @details This function outputs the features generated from input files.
#' @examples 
#' generate_feature_wp(ptm_site = "S",
#'              flanking_size = 12,
#'              SPIDER = T,
#'              positive_info_file = "ps_PSP.tsv",
#'              protein_fasta_file = "S_sp_fasta.tsv",
#'              output_label = "ps_0103")

generate_feature_wp = function(ptm_site, flanking_size=12, SPIDER = T,
                            positive_info_file, protein_fasta_file,
                            output_label)
{
  ### get aaindex feature data
  aaindex_cluster_order = data.table::fread("aaindex_cluster_order.tsv", stringsAsFactors = F)
  aais = aaindex_cluster_order$cluster_name
  
  
  ### get spider names
  structural_fc = rep(c("ASA","HSE_au","pC", "pH"),2*flanking_size)
       # flanking_size = 12
  structural_fn1 = as.vector(t(matrix(rep(seq(1:flanking_size),4),flanking_size,4)))
  structural_fn2 = as.vector(t(matrix(rep(flanking_size + 1 + seq(1:flanking_size),4),flanking_size,4)))
  structural_fn = c(structural_fn1, structural_fn2)
  
  fis = sapply(1:length(structural_fc), function(x){
    return(paste(structural_fc[x], structural_fn[x], sep = "_"))
  })
  
  ### get pssm names
  position_seq = c(seq(1:flanking_size), (flanking_size + 1 + seq(1:flanking_size)))
  pssm_name = paste0("pssm_", position_seq)
  
  full_feature = c(aais, fis, pssm_name)
  
  
  ### find the intersection of user provided protein sequences and the spider info 
  #### protein_fasta_file = "S_sp_fasta.tsv"
  all_sp = readLines(protein_fasta_file)
  sp_id = all_sp[c(T,F)]
  sp_seq = all_sp[c(F,T)]
  sp_uni_id = sapply(1:length(sp_id),function(i) strsplit(sp_id[i], split = "|", fixed = T)[[1]][2])
  
  ps_info = data.table::fread(positive_info_file, stringsAsFactors = F)
  
  if(SPIDER == TRUE)
  {
    
    spider_protID = data.table::fread("spider_protID.tsv", stringsAsFactors = F)
    
    not_na_id = intersect(sp_uni_id, spider_protID$x)
    which_sel = which(sp_uni_id%in%not_na_id)
    
    not_na_sp_uni_id = sp_uni_id[which_sel]
    not_na_sp_seq = sp_seq[which_sel]
    
    
    
    ### positive PTM sites from PSP
    #### positive_info_file = "ps_PSP.tsv"
    
    
    ##############################################################################################
    ############ FEATURE EXTRACTION
    ##############################################################################################
    window_formation(ptm_site,flanking_size,
                     not_na_sp_seq, not_na_sp_uni_id, ps_info,
                     output_label)
    
    
    aaindex_feature_extraction(aaindex_cluster_order,
                               paste0(output_label,"_candidate.Rds"),
                               flanking_size + 1,
                               output_label)
    
    
    ### ok I need to change the file name to flexible for different size scenerio
    
    
    
    if(flanking_size == 12)
    {
      spider_rds_name = "extracted_spider.Rds"
    }
    else if(flanking_size == 7)
    {
      spider_rds_name = "extracted_spider_7.Rds"
    }else if(flanking_size == 5)
    {
      spider_rds_name = "extracted_spider_5.Rds"
    } else
    {
      cat("SIZE ERROR!","\n")
    }
    
    spider_feature_joining_without_mean(paste0(output_label,"_candidate.Rds"),
                                        spider_rds_name,
                                        c(1,6,8,9),
                                        c(8,9),
                                        flanking_size + 1,  
                                        output_label)
    
    
    pssm_generation(paste0(output_label,"_candidate.Rds"),
                    flanking_size+1,
                    output_label)
    
    pssm_feature_extraction(paste0(output_label,"_candidate.Rds"),
                            paste0(output_label,"_pssm.Rds"),
                            flanking_size + 1,
                            output_label)
    
    
    
    combine_all_features(pos_aaindex = paste0(output_label, "_noc_pos_cluster_matrix.Rds"),
                         candi_aaindex = paste0(output_label, "_noc_candi_cluster_matrix.Rds"),
                         pos_spider = paste0(output_label, "_noc_pos_structure_matrix.Rds"),
                         candi_spider = paste0(output_label, "_noc_candi_structure_matrix.Rds"),
                         pos_pssm = paste0(output_label, "_noc_pos_pssm_matrix.Rds"),
                         candi_pssm = paste0(output_label, "_noc_candi_pssm_matrix.Rds"),
                         output_label)
    
    
    
  }else{
    

    
    not_na_sp_uni_id = sp_uni_id
    not_na_sp_seq = sp_seq
    
    ### positive PTM sites from PSP
    #### positive_info_file = "ps_PSP.tsv"
    
    
    ##############################################################################################
    ############ FEATURE EXTRACTION
    ##############################################################################################
    window_formation(ptm_site,flanking_size,
                     not_na_sp_seq, not_na_sp_uni_id, ps_info,
                     output_label)
    
    
    aaindex_feature_extraction(aaindex_cluster_order,
                               paste0(output_label,"_candidate.Rds"),
                               flanking_size + 1,
                               output_label)
    
    
    pssm_generation(paste0(output_label,"_candidate.Rds"),
                    flanking_size+1,
                    output_label)
    
    pssm_feature_extraction(paste0(output_label,"_candidate.Rds"),
                            paste0(output_label,"_pssm.Rds"),
                            flanking_size + 1,
                            output_label)
    
    
    combine_all_features_no_spider(pos_aaindex = paste0(output_label, "_noc_pos_cluster_matrix.Rds"),
                         candi_aaindex = paste0(output_label, "_noc_candi_cluster_matrix.Rds"),
                         pos_pssm = paste0(output_label, "_noc_pos_pssm_matrix.Rds"),
                         candi_pssm = paste0(output_label, "_noc_candi_pssm_matrix.Rds"),
                         output_label)
    
    
    
  }
  
}


#' Generate training and prediction datasets for whole proteome
#' 
#' This function generates training and prediction datasets for whole proteome to predict 1/n_fold of the proteome by the other n-1 fold data.
#' @param n_fold Number of folds used for training and prediction, default set to 2
#' @param lower_bound The lower bound of the scaled data range, default to -1.
#' @param upper_bound The upper bound of the scaled data range, default to 1.
#' @param positive_feature_file A Rds file containing the positive features(generated by g_feature_wp()).
#' @param negative_feature_file A Rds file containing the negative features(generated by g_feature_wp()).
#' @param output_label The string to tag the output files.
#' @import stringr dplyr magrittr caret data.table
#' @export
#' @details This function outputs formatted feature files ready for Liblinear training and prediction.
#' @examples 
#' generate_training_test_data_wp(n_fold = 2,
#'                         lower_bound = -1,
#'                         upper_bound = 1,
#'                         positive_feature_file = "ps_0103_noc_not_na_pos_feature.Rds",
#'                         negative_feature_file = "ps_0103_noc_not_na_candi_feature.Rds",
#'                         output_label = "ps_0103")

generate_training_test_data_wp = function(n_fold = 2, 
                                   lower_bound = -1, 
                                   upper_bound = 1,
                                   positive_feature_file,
                                   negative_feature_file,
                                   output_label)
{
  
    combine_pos_neg(positive_feature_file,
                  negative_feature_file,
                  output_label)
  
  ### it seems the seed is not correct
  
  construct_n_fold_cv_without_decoy(paste0(output_label,"_feature_matrix.Rds"),
                                    paste0(output_label,"_feature_label.Rds"),
                                    n_fold, 
                                    output_label)
  #paste0(output_label, "_",n_fold))
  
  test_label_file_names = rep("", n_fold)
  test_candi_index_file_names = rep("", n_fold)
  test_pos_index_file_names = rep("", n_fold)
  
  
  
  
  for(i in 1:n_fold)
  {
    test_label_file_names[i] = paste0(output_label, "_test_label_",i,".tsv")
    test_candi_index_file_names[i] = paste0(output_label, "_test_candi_ind_",i,".Rds")
    test_pos_index_file_names[i] = paste0(output_label, "_test_pos_ind_",i,".Rds")
    
  }
  
  saveRDS(test_label_file_names, file = paste0(output_label, "_test_label_names.Rds"))
  saveRDS(test_candi_index_file_names, file = paste0(output_label, "_candi_index_names.Rds"))
  saveRDS(test_pos_index_file_names, file = paste0(output_label, "_pos_index_names.Rds"))
  
  ### the following need to be changed to adjust to n folds
  
  training_feature_file_names  = rep("", n_fold)
  test_feature_file_names = rep("", n_fold)
  
  
  for(i in 1:n_fold)
  {
    
    balanced_size_sampling_after_combining_pos_neg(paste0(output_label, "_train_feature_",i,".Rds"),
                                                   paste0(output_label, "_train_label_",i,".Rds"),
                                                   paste0(output_label, "_train_",i))
    
    scale_train_test_single(paste0(output_label, "_train_",i,"_balanced_feature_matrix.Rds"),
                            paste0(output_label, "_test_feature_",i,".Rds"),
                            upper_bound,lower_bound,
                            paste0(output_label, "_train_",i,"_balanced_feature"),
                            paste0(output_label, "_test_",i, "_feature"))
    
    libsvm_formating_single(paste0(output_label, "_train_",i,"_balanced_feature_scale.Rds"),
                            paste0(output_label, "_train_",i,"_balanced_feature_label.Rds"),
                            paste0(output_label, "_test_",i,"_feature_scale.Rds"),
                            paste0(output_label, "_test_label_",i,".Rds"),
                            paste0(output_label,"_",i))
    
    training_feature_file_names[i] = paste0(output_label, "_",i,"_train_feature_svm.tsv")
    test_feature_file_names[i] = paste0(output_label, "_", i, "_test_feature_svm.tsv")
    
    
  }
  
  ### write the names of the training and test files into two Rds files
  
  saveRDS(training_feature_file_names, file = paste0(output_label, "_training_feature_names.Rds"))
  saveRDS(test_feature_file_names, file = paste0(output_label, "_test_feature_names.Rds"))
  
  
}


#' Process n_fold cv 
#'
#' Process the feature files with Liblinear training and prediction. Conduct n_fold cross validation with the files supplied.
#' @param liblinear_dir Absolute peth of Liblinear tool.
#' @param n_fold Number of folds used for training and prediction, default set to 2
#' @param feature_file_path Absolute path of the feature files.
#' @param training_file_names An Rds file containing the file names of the training feature data.
#' @param test_file_names An Rds file containing the file names of the test feature data.
#' @param output_label The string to tag the output files.
#' @param cvlog_path_name The path and name of the log files, which hold the details of Liblinear procedures.
#' @import stringr dplyr magrittr caret data.table
#' @export
#' @details This function call Liblinear library to perform n_fold training/prediction, the prediction score will be part of the output files.
#' @examples 
#' process_n_fold_cross_validation(liblinear_dir = "/data/ginny/liblinear-2.11/",
#'                                 n_fold = 2,
#'                                 feature_file_path = "/data/ginny/test_package/",
#'                                 training_file_names = "ps_0103_training_feature_names.Rds",
#'                                 test_file_names = "ps_0103_test_feature_names.Rds",
#'                                 output_label = "ps_0103",
#'                                 cvlog_path_name = "/data/ginny/test_package/cvlog.txt")
#' 
process_n_fold_cross_validation = function(liblinear_dir,
                                           n_fold,
                                           feature_file_path,
                                           training_file_names,
                                           test_file_names,
                                           output_label,
                                           cvlog_path_name)
{
   
  training_files = readRDS(training_file_names)
  test_files = readRDS(test_file_names)
  
  sink(cvlog_path_name, append = T)
  
  
  prediction_file_names = rep("", n_fold)
  
  for(i in 1:n_fold)
  {
    cat(i, "\n")
    
    
    train_feature_name = paste0(feature_file_path, training_files[i])
    test_feature_name = paste0(feature_file_path, test_files[i])
    model_name = paste0(feature_file_path, output_label, "_model_",i,".tsv")
    test_prediction_name = paste0(feature_file_path,output_label, "_predict_",i,".tsv" )
    
    prediction_file_names[i] = test_prediction_name
    
    
    cv_command = paste0(liblinear_dir,"train -s 2 -C ", train_feature_name," >> cv_log.txt")
    system(cv_command)
    
    tmp_5cv_file = readLines("cv_log.txt")
    
    it = tmp_5cv_file[length(tmp_5cv_file)]
    sit = strsplit(it, "=")[[1]][2]
    get_c =  as.numeric(strsplit(sit, split = " ")[[1]][2])
    
    
    training_command = paste0(liblinear_dir, "train -s 2 -c ", get_c," ",train_feature_name," ",model_name)
    system(training_command)
    
    
    predict_test_command = paste0(liblinear_dir, "predict -b 1 ",test_feature_name," ", model_name,
                                  " ",test_prediction_name," >> predict_test_log.txt")
    system(predict_test_command)
    
    system(paste0("rm ", "*log.txt") )
    
  }
  
  
  saveRDS(prediction_file_names, file = paste0(output_label, "_prediction_file_names.Rds"))
  
  
  sink()
  
  
}



#' Present Liblinear prediction results 
#'
#' Combine Liblinear prediction results and window, position informations of sites. Calculate score threshold at user specified sensitivity level.
#' @param positive_index_file_names An Rds file containing the names of indices for positive windows.
#' @param candi_index_file_names An Rds file containing the names of indices for candidate(negative) windows.
#' @param prediction_score_file_names An Rds file containing the file names of Liblinear predicted scores.
#' @param test_label_file_names An Rds file containing the file names of the label of the test data.
#' @param candidate_df_Rds An Rds file containing the data frame of all candiate sites information(generated by g_feature_wp)/
#' @param specificity_level A numerical number indicating the specificity user requires the classifier to achieve, default set to 0.99
#' @param output_label The string to tag the output files.
#' @import stringr dplyr magrittr data.table
#' @export
#' @examples 
#' present_prediction_wp(positive_index_file_names = "ps_0103_pos_index_names.Rds",
#'                                  candi_index_file_names = "ps_0103_candi_index_names.Rds",
#'                                  prediction_score_file_names = "ps_0103_prediction_file_names.Rds",
#'                                  test_label_file_names = "ps_0103_test_label_names.Rds",
#'                                  candidate_df_Rds = "ps_0103_candidate.Rds",
#'                                  specificity_level = 0.99,
#'                                  output_label = "ps_0103")

present_prediction_wp = function(positive_index_file_names,
                                 candi_index_file_names,
                                 prediction_score_file_names,
                                 test_label_file_names,
                                 candidate_df_Rds,
                                 specificity_level = 0.99,
                                 output_label)
{
  
  st =  get_score_threshold(prediction_score_file_names = prediction_score_file_names, 
                            test_label_file_names = test_label_file_names,
                            specificity_level = specificity_level,
                            output_label = output_label)
  
  cat("score threshold: ", st, "\n")
  
  
  uniprot_genename = fread("uniprotID_genename.tsv", header = T, stringsAsFactors = F)
  
  
  
  assemble_window_score_cv(candidate_df_Rds = candidate_df_Rds,
                           positive_index_file_names = positive_index_file_names,
                           candi_index_file_names = candi_index_file_names,
                           positive_score_Rds = paste0(output_label, "_positive_score.Rds"),
                           candi_score_Rds = paste0(output_label, "_candi_score.Rds"),
                           score_threshold = st,
                           id_convert = uniprot_genename,
                           output_label = output_label)
  
}


#'  A function to predict and annotate whole proteom PTM sites
#'  
#' @param ptm_site The target amino acid of the given PTM type, in upper-case single letter representation.
#' @param flanking_size The number of residues surrounding each side of the center residue, The total window size will be 2*flanking_size+1 (default to 12).
#' @param SPIDER A boolean variable indicating whether to use SPIDER3 features (default set to TRUE.)  
#' @param positive_info_file A text file containing the positive PTM sites in the required format. 
#' @param protein_fasta_file A text file containing the protein sequences of interest in fasta format.  
#' @param liblinear_dir The path for the Liblinear tool. 
#' @param n_fold The number of folds used for training and prediction in cross validation stage (default set to 2).
#' @param feature_file_path The path for the feature files.  
#' @param lower_bound The lower bound of the scaled data range (default to -1).
#' @param upper_bound The upper bound of the scaled data range (default to 1).
#' @param cvlog_path_name The path and name of the log files, which hold the details of Liblinear procedures.
#' @param specificity_level  A number ranges from 0 to 1 indicating the specificity user requires the classifier to achieve (default to 0.99).
#' @param output_label The string to tag the output files.
#' @import stringr dplyr magrittr data.table
#' @export
#' @details This function outputs the features generated from input files.
#' @examples 
#' predict_on_whole_proteome(ptm_site = "S",
#'                          flanking_size = 12,
#'                          SPIDER = T,
#'                          positive_info_file = "ps_PSP.tsv",
#'                          protein_fasta_file = "S_sp_fasta.tsv",
#'                          n_fold = 2,
#'                          lower_bound = -1,
#'                          upper_bound = 1,
#'                          liblinear_dir = "/data/ginny/liblinear-2.11/",
#'                          feature_file_path = "/data/ginny/test_package/",
#'                          cvlog_path_name = "/data/ginny/test_package/cvlog.txt",
#'                          specificity_level = 0.99,
#'                          output_label = "ps_0103")




predict_on_whole_proteome = function(ptm_site,
                                     flanking_size = 12,
                                     SPIDER = T,
                                     positive_info_file,
                                     protein_fasta_file,
                                     n_fold = 2,
                                     lower_bound = -1,
                                     upper_bound = 1,
                                     liblinear_dir,
                                     feature_file_path,
                                     cvlog_path_name,
                                     specificity_level = 0.99,
                                     output_label)
{
  
  generate_feature_wp(ptm_site = ptm_site,
                                  flanking_size = flanking_size,
                      SPIDER = SPIDER,
                                  positive_info_file = positive_info_file,
                                  protein_fasta_file = protein_fasta_file,
                                  output_label = output_label)
  cat("STEP1: Feature generated.", "\n")
  
  generate_training_test_data_wp(n_fold = n_fold, 
                                             lower_bound = lower_bound,
                                             upper_bound = upper_bound,
                                             positive_feature_file = paste0(output_label,
                                                                            "_noc_not_na_pos_feature.Rds"),
                                             negative_feature_file = paste0(output_label,
                                                                            "_noc_not_na_candi_feature.Rds"),
                                             output_label = output_label)
  
  cat("STEP2: Training and predict prepared.", "\n")
  
  process_n_fold_cross_validation(liblinear_dir = liblinear_dir,
                                  n_fold = n_fold,
                                  feature_file_path = feature_file_path,
                                  training_file_names = paste0(output_label,
                                                               "_training_feature_names.Rds"),
                                  test_file_names = paste0(output_label,
                                                           "_test_feature_names.Rds"),
                                  output_label = output_label,
                                  cvlog_path_name = cvlog_path_name)
  
  cat("STEP3: Liblinear processed.", "\n")
  
  
  present_prediction_wp(positive_index_file_names = paste0(output_label,"_pos_index_names.Rds"),
                                    candi_index_file_names = paste0(output_label,"_candi_index_names.Rds"),
                                    prediction_score_file_names = paste0(output_label,"_prediction_file_names.Rds"),
                                    test_label_file_names = paste0(output_label, "_test_label_names.Rds"),
                                    candidate_df_Rds = paste0(output_label,"_candidate.Rds"),
                                    specificity_level = specificity_level,
                                    output_label = output_label)
  
  cat("STEP4: Prediction results combined.", "\n")
  
  PTM_domain_mapping_enrichment(window_score_label_Rds = paste0(output_label,"_window_score_df.Rds"), 
                                output_label = output_label)

  cat("STEP5:Domain mapped and annotated.", "\n")
  cat("Finished!","\n")
}




