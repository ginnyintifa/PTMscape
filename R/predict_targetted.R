### the wrapper function for targetted mode


# t is targeted.


#' Generate feature data for targeted proteome.
#' 
#' This function generates feature data for targeted by mapping the positive PTM info to protein sequences, constructing windows and extracting 3 sets of features.
#' @param ptm_site The amino acid this PTM involves, in upper-case single letter representation.
#' @param flanking_size The number of residues surrounding each side of the center residue, the total window size will be 2*flanking_size+1, default to 12.
#' @param SPIDER A boolean variable indicating the usage of SPIDER3 features, default set to TRUE.
#' @param positive_info_file A text file containing the positive PTM sites info in required format.
#' @param known_protein_fasta_file A text file containing the proteins sequences of interest and known PTM sites in Fasta format.
#' @param predict_protein_fasta_file A text file containing the proteins sequences with PTM sites to be predicted in Fasta format.
#' @param output_label_training The string to tag the output files associated with training proteins.
#' @param output_label_predict The string to tag the output files associated with prediction proteins.
#' @import stringr dplyr magrittr data.table
#' @export
#' @details This function outputs the features generated from input files.
#' @examples 
#' generate_feature_t(ptm_site = "S",
#'             flanking_size = 12,
#'             SPIDER = T,
#'             positive_info_file = "known_ps.tsv",
#'             known_protein_fasta_file = "known_fasta.tsv",
#'             predict_protein_fasta_file = "predict_fasta.tsv",
#'             output_label_training = "ps_training",
#'             output_label_predict = "ps_predict")


generate_feature_t = function(ptm_site, flanking_size=12, 
                              SPIDER = T,
                                     positive_info_file, 
                                     known_protein_fasta_file, predict_protein_fasta_file,
                                     output_label_training, output_label_predict)
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
  
  
  
  #####
  known_all_sp = readLines(known_protein_fasta_file)
  known_sp_id = known_all_sp[c(T,F)]
  known_sp_seq = known_all_sp[c(F,T)]
  known_sp_uni_id = sapply(1:length(known_sp_id),function(i) 
    strsplit(known_sp_id[i], split = "|", fixed = T)[[1]][2])
  
  predict_all_sp = readLines(predict_protein_fasta_file)
  predict_sp_id = predict_all_sp[c(T,F)]
  predict_sp_seq = predict_all_sp[c(F,T)]
  predict_sp_uni_id = sapply(1:length(predict_sp_id),function(i) 
    strsplit(predict_sp_id[i], split = "|", fixed = T)[[1]][2])
  
  
  
  ps_info = data.table::fread(positive_info_file, stringsAsFactors = F)
  
  if(SPIDER == TRUE)
  {
    spider_protID = data.table::fread("spider_protID.tsv", stringsAsFactors = F)
    
    
    known_not_na_id = intersect(known_sp_uni_id, spider_protID$x)
    known_which_sel = which(known_sp_uni_id%in%known_not_na_id)
    
    known_not_na_sp_uni_id = known_sp_uni_id[known_which_sel]
    known_not_na_sp_seq = known_sp_seq[known_which_sel]
    
    #####
    
    predict_not_na_id = intersect(predict_sp_uni_id, spider_protID$x)
    predict_which_sel = which(predict_sp_uni_id%in%predict_not_na_id)
    
    predict_not_na_sp_uni_id = predict_sp_uni_id[predict_which_sel]
    predict_not_na_sp_seq = predict_sp_seq[predict_which_sel]
    
    
    ### positive PTM sites from PSP
    #### positive_info_file = "ps_PSP.tsv"
    
    
    ##############################################################################################
    ############ FEATURE EXTRACTION 
    ##############################################################################################
    window_formation(ptm_site,flanking_size,
                     known_not_na_sp_seq, known_not_na_sp_uni_id, ps_info,
                     output_label_training)
    
    
    aaindex_feature_extraction(aaindex_cluster_order,
                               paste0(output_label_training,"_candidate.Rds"),
                               flanking_size + 1,
                               output_label_training)
    
    
    
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
    
    
    spider_feature_joining_without_mean(paste0(output_label_training,"_candidate.Rds"),
                                        spider_rds_name,
                                        c(1,6,8,9),
                                        c(8,9),
                                        flanking_size + 1,  
                                        output_label_training)
    
    
    pssm_generation(paste0(output_label_training,"_candidate.Rds"),
                    flanking_size+1,
                    output_label_training)
    
    pssm_feature_extraction(paste0(output_label_training,"_candidate.Rds"),
                            paste0(output_label_training,"_pssm.Rds"),
                            flanking_size + 1,
                            output_label_training)
    
    
    
    combine_all_features(pos_aaindex = paste0(output_label_training, "_noc_pos_cluster_matrix.Rds"),
                         candi_aaindex = paste0(output_label_training, "_noc_candi_cluster_matrix.Rds"),
                         pos_spider = paste0(output_label_training, "_noc_pos_structure_matrix.Rds"),
                         candi_spider = paste0(output_label_training, "_noc_candi_structure_matrix.Rds"),
                         pos_pssm = paste0(output_label_training, "_noc_pos_pssm_matrix.Rds"),
                         candi_pssm = paste0(output_label_training, "_noc_candi_pssm_matrix.Rds"),
                         output_label = output_label_training)
    
    #### for uncharted proteins, only divide into windows, not need to record the label 
    
    
    ### think about this, only extract features 
    window_formation_no_positive(ptm_site,flanking_size,
                                 predict_not_na_sp_seq, predict_not_na_sp_uni_id,
                                 output_label_predict)
    ### the windows in the predition part are called "_predict.Rds"
    
    aaindex_feature_extraction(aaindex_cluster_order,
                               paste0(output_label_predict,"_candidate.Rds"),
                               flanking_size + 1,
                               output_label_predict)
    
    spider_feature_joining_without_mean(paste0(output_label_predict,"_candidate.Rds"),
                                        spider_rds_name,
                                        c(1,6,8,9),
                                        c(8,9),
                                        flanking_size + 1,  
                                        output_label_predict)
    
    
    pssm_feature_extraction(paste0(output_label_predict,"_candidate.Rds"),
                            paste0(output_label_training,"_pssm.Rds"),
                            flanking_size + 1,
                            output_label_predict)
    
    combine_all_features(candi_aaindex = paste0(output_label_predict, "_noc_candi_cluster_matrix.Rds"),
                         candi_spider = paste0(output_label_predict, "_noc_candi_structure_matrix.Rds"),
                         candi_pssm = paste0(output_label_predict, "_noc_candi_pssm_matrix.Rds"),
                         output_label = output_label_predict)
    
  }else{
    
    known_not_na_sp_uni_id = known_sp_uni_id
    known_not_na_sp_seq = known_sp_seq
    
    
    predict_not_na_sp_uni_id = predict_sp_uni_id
    predict_not_na_sp_seq = predict_sp_seq
    
    
    
    window_formation(ptm_site,flanking_size,
                     known_not_na_sp_seq, known_not_na_sp_uni_id, ps_info,
                     output_label_training)
    
    
    aaindex_feature_extraction(aaindex_cluster_order,
                               paste0(output_label_training,"_candidate.Rds"),
                               flanking_size + 1,
                               output_label_training)
    
    
    
    
    pssm_generation(paste0(output_label_training,"_candidate.Rds"),
                    flanking_size+1,
                    output_label_training)
    
    pssm_feature_extraction(paste0(output_label_training,"_candidate.Rds"),
                            paste0(output_label_training,"_pssm.Rds"),
                            flanking_size + 1,
                            output_label_training)
    
    
    
    combine_all_features_no_spider(pos_aaindex = paste0(output_label_training, "_noc_pos_cluster_matrix.Rds"),
                         candi_aaindex = paste0(output_label_training, "_noc_candi_cluster_matrix.Rds"),
                         pos_pssm = paste0(output_label_training, "_noc_pos_pssm_matrix.Rds"),
                         candi_pssm = paste0(output_label_training, "_noc_candi_pssm_matrix.Rds"),
                         output_label = output_label_training)
    
    #### for uncharted proteins, only divide into windows, not need to record the label 
    
    
    ### think about this, only extract features 
    window_formation_no_positive(ptm_site,flanking_size,
                                 predict_not_na_sp_seq, predict_not_na_sp_uni_id,
                                 output_label_predict)
    ### the windows in the predition part are called "_predict.Rds"
    
    aaindex_feature_extraction(aaindex_cluster_order,
                               paste0(output_label_predict,"_candidate.Rds"),
                               flanking_size + 1,
                               output_label_predict)
    
    
    pssm_feature_extraction(paste0(output_label_predict,"_candidate.Rds"),
                            paste0(output_label_training,"_pssm.Rds"),
                            flanking_size + 1,
                            output_label_predict)
    
    combine_all_features_no_spider(candi_aaindex = paste0(output_label_predict, "_noc_candi_cluster_matrix.Rds"),
                         candi_pssm = paste0(output_label_predict, "_noc_candi_pssm_matrix.Rds"),
                         output_label = output_label_predict)
    
    
    
    
  }
  
 
  }



### for this mode
### no need to do cross validation, train on known part and test on predict part





#' Generate training and prediction datasets for targeted proteome
#' 
#' This function generates training and prediction datasets for whole proteome to predict 1/n_fold of the proteome by the other n-1 fold data.
#' @param lower_bound The lower bound of the scaled data range, default to -1.
#' @param upper_bound The upper bound of the scaled data range, default to 1.
#' @param positive_training_feature_file A Rds file containing the positive training features(generated by g_feature_t()).
#' @param negative_training_feature_file A Rds file containing the negative training features(generated by g_feature_t()).
#' @param negative_predict_feature_file A Rds file containing the negative prediction features(generated by g_feature_t()).
#' @param output_label_training The string to tag the output files associated with training proteins.
#' @param output_label_predict The string to tag the output files associated with prediction proteins.
#' @import stringr dplyr magrittr caret data.table
#' @export
#' @details This function outputs formatted feature files ready for Liblinear training and prediction.
#' @examples 
#' 
#' generate_training_test_data_t(lower_bound = -1, upper_bound = 1,
#'positive_training_feature_file ="ps_training_noc_not_na_pos_feature.Rds",
#'negative_training_feature_file = "ps_training_noc_not_na_candi_feature.Rds",
#'negative_predict_feature_file = "ps_predict_noc_not_na_candi_feature.Rds",
#'output_label_training = "ps_training",
#'output_label_predict = "ps_predict")

generate_training_test_data_t = function(lower_bound = -1, upper_bound = 1,
                                                positive_training_feature_file, negative_training_feature_file,
                                                negative_predict_feature_file,
                                                output_label_training, output_label_predict)
{
  
  # output_label = "ubi_pred"
  combine_pos_neg(positive_training_feature_file,
                  negative_training_feature_file,
                  output_label_training)
  
  ### create label Rds file for the negative_predict_feature_file 
  npf = readRDS(negative_predict_feature_file)
  pred_label = rep(0, nrow(npf))
  saveRDS(pred_label, file = paste0(output_label_predict,"_test_label.Rds"))
  

  
  balanced_size_sampling_after_combining_pos_neg(paste0(output_label_training, "_feature_matrix.Rds"),
                                                 paste0(output_label_training, "_feature_label.Rds"),
                                                 paste0(output_label_training))
  
  scale_train_test_single(paste0(output_label_training, "_balanced_feature_matrix.Rds"),
                          negative_predict_feature_file,
                          upper_bound,lower_bound,
                          paste0(output_label_training, "_train_balanced_feature"),
                          paste0(output_label_predict, "_test_feature"))
  
  libsvm_formating_single(paste0(output_label_training, "_train_balanced_feature_scale.Rds"),
                          paste0(output_label_training, "_balanced_feature_label.Rds"),
                          paste0(output_label_predict, "_test_feature_scale.Rds"),
                          paste0(output_label_predict, "_test_label.Rds"),
                          output_label_predict)
  
  
 }




#' Process prediction with Liblinear 
#'
#' Process the feature files with Liblinear training and prediction. Train on known proteins and predict on unknown proteins
#' @param liblinear_dir Absolute path of Liblinear tool.
#' @param feature_file_path Absolute path of the feature files.
#' @param training_file_name A string indicating file name of the training feature data.
#' @param test_file_name A string indicating the file name of the test feature data.
#' @param output_label The string to tag the output files.
#' @param cvlog_path_name The path and name of the log files, which hold the details of Liblinear procedures.
#' @import stringr dplyr magrittr caret data.table
#' @export
#' @details This function call Liblinear library to perform n_fold training/prediction, the prediction score will be part of the output files.
#' @examples 
#' predict_with_liblinear(liblinear_dir = "/data/ginny/liblinear-2.11/",
#'                       feature_file_path = "/data/ginny/test_package/",
#'                       training_file_name = "ps_predict_train_feature_svm.tsv",
#'                       test_file_name = "ps_predict_test_feature_svm.tsv",
#'                       output_label = "ps_predict",
#'                       cvlog_path_name = "/data/ginny/test_package/cvlog.txt")



predict_with_liblinear = function(liblinear_dir,
                                  feature_file_path,
                                  training_file_name,
                                  test_file_name,
                                  output_label,
                                  cvlog_path_name)
{
  sink(cvlog_path_name, append = T)
  
  ### I will need to make sure the files are written in order 
  ### ok all the file of the file names are Rds format and the data structure is vector
  ###

    train_feature_name = paste0(feature_file_path, training_file_name)
    test_feature_name = paste0(feature_file_path, test_file_name)
    model_name = paste0(feature_file_path, output_label,"_model.tsv")
    test_prediction_name = paste0(feature_file_path,output_label, "_predict.tsv" )
    
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
    
  
  
  sink()
  
  
  
  
}








#' Present Liblinear prediction results 
#'
#' Combine Liblinear prediction results and window, position informations of sites. Calculate score threshold at user specified sensitivity level.
#' @param flag_for_score_threshold_chosen A string indicating whether use reference score threshold or get from the user supplied training data, defalt set to "reference".
#' @param score_threshold A numerical value between 0 to 1 indicating the reference score threshold (supply when in "reference")
#' @param ptm_site The amino acid this PTM involves, in upper-case single letter representation.
#' @param flanking_size The number of residues surrounding each side of the center residue, the total window size will be 2*flanking_size+1, default to 12.
#' @param SPIDER A boolean variable indicating the usage of SPIDER3 features, default set to TRUE.
#' @param positive_info_file A text file containing the positive PTM sites info in required format.
#' @param known_protein_fasta_file A text file containing the proteins sequences of interest and known PTM sites in Fasta format.
#' @param pred_candidate_df_Rds_name An Rds file containing the candidate data frame of the proteins to be predicted.
#' @param pred_score_file_name     A text file containing the predicted score for the proteins of interest.   
#' @param liblinear_dir Absolute path of Liblinear tool.
#' @param n_fold Number of folds used for training and prediction, default set to 2
#' @param feature_file_path Absolute path of the feature files.
#' @param lower_bound The lower bound of the scaled data range, default to -1.
#' @param upper_bound The upper bound of the scaled data range, default to 1.
#' @param cvlog_path_name The path and name of the log files, which hold the details of Liblinear procedures.
#' @param specificity_level A numerical number indicating the specificity user requires the classifier to achieve, default set to 0.99. Used only not in "reference" mode.
#' @param output_label The string to tag the output files in threshold getting purpose.
#' @import stringr dplyr magrittr data.table
#' @export
#' @examples 
#' present_prediction_t(flag_for_score_threshold_chosen = "cv",
#'                            score_threshold = NULL,
#'                            ptm_site = "S", flanking_size = 12,
#'                            pred_candidate_df_Rds_name = "ps_predict_candidate.Rds",
#'                            pred_score_file_name = "ps_predict_predict.tsv",
#'                            positive_info_file = "known_ps.tsv", 
#'                            known_protein_fasta_file = "known_fasta.tsv",
#'                            n_fold = 2,
#'                            lower_bound = -1,
#'                            upper_bound = 1,
#'                            liblinear_dir = "/data/ginny/liblinear-2.11/",
#'                            feature_file_path = "/data/ginny/test_package/",
#'                            output_path = "/data/ginny/test_package/",
#'                            cvlog_path_name = "/data/ginny/test_package/cvlog.txt",
#'                            specificity_level = 0.99,
#'                            output_label = "ps_target")

present_prediction_t = function(flag_for_score_threshold_chosen = "reference",
                                       score_threshold,
                                       ptm_site, flanking_size = 12,
                                       SPIDER = T,
                                       pred_candidate_df_Rds_name,
                                       pred_score_file_name,
                                       positive_info_file, 
                                       known_protein_fasta_file,
                                       n_fold = 2,
                                       lower_bound = -1,
                                       upper_bound = 1,
                                       liblinear_dir,
                                       feature_file_path,
                                       cvlog_path_name,
                                       specificity_level,
                                       output_label)
{
  
  uniprot_genename = fread("uniprotID_genename.tsv", header = T, stringsAsFactors = F)
  
  
  
  if(flag_for_score_threshold_chosen == "reference")
  {
    assemble_window_score_target(prediction_score_file = pred_score_file_name,
                                 predict_candidate_df_Rds = pred_candidate_df_Rds_name,
                                 score_threshold = score_threshold,
                                 id_convert = uniprot_genename,
                                 output_label = output_label)
    
  }else{
    generate_feature_wp(ptm_site = ptm_site,
                                    flanking_size = flanking_size,
                                    SPIDER = SPIDER,
                                    positive_info_file = positive_info_file,
                                    protein_fasta_file = known_protein_fasta_file,
                                    output_label = paste0(output_label, "_known_cv"))
    
    generate_training_test_data_wp(n_fold = n_fold, 
                                               lower_bound = lower_bound,
                                               upper_bound = upper_bound,
                                               positive_feature_file = paste0(output_label,"_known_cv",
                                                                              "_noc_not_na_pos_feature.Rds"),
                                               negative_feature_file = paste0(output_label, "_known_cv",
                                                                              "_noc_not_na_candi_feature.Rds"),
                                               output_label = paste0(output_label, "_known_cv"))
    process_n_fold_cross_validation(liblinear_dir = liblinear_dir,
                                    n_fold = n_fold,
                                    feature_file_path = feature_file_path,
                                    training_file_names = paste0(output_label,"_known_cv",
                                                                 "_training_feature_names.Rds"),
                                    test_file_names = paste0(output_label,"_known_cv",
                                                             "_test_feature_names.Rds"),
                                    output_label = paste0(output_label, "_known_cv"),
                                    cvlog_path_name = cvlog_path_name)
    
    st = get_score_threshold(prediction_score_file_names = paste0(output_label,"_known_cv",
                                                                  "_prediction_file_names.Rds"),
                             test_label_file_names = paste0(output_label,"_known_cv",
                                                            "_test_label_names.Rds"),
                             specificity_level = specificity_level,
                             output_label = paste0(output_label, "_known_cv"))
          

    assemble_window_score_target(prediction_score_file = pred_score_file_name,
                                 predict_candidate_df_Rds = pred_candidate_df_Rds_name,
                                 score_threshold = st,
                                 id_convert = uniprot_genename,
                                 output_label = output_label)
    
    
  }
  
  
  
}






#' Predict on targeted proteome.
#' 
#' @param ptm_site The target amino acid of the given PTM type, in upper-case single letter representation.
#' @param flanking_size The number of residues surrounding each side of the center residue, The total window size will be 2*flanking_size+1 (default to 12).
#' @param SPIDER A boolean variable indicating whether to use SPIDER3 features (default set to TRUE.)  
#' @param positive_info_file A text file containing the positive PTM sites in the required format. 
#' @param known_protein_fasta_file A text file containing the proteins sequences of interest and known PTM sites in fasta format.  
#' @param predict_protein_fasta_file A text file containing the proteins sequences with PTM sites to be predicted in fasta format. 
#' @param output_label_training The string to tag the output files associated with training proteins.
#' @param output_label_predict The string to tag the output files associated with prediction proteins.
#' @param liblinear_dir The path for the Liblinear tool. 
#' @param n_fold The number of folds used for training and prediction in cross validation stage (default set to 2).
#' @param feature_file_path The path for the feature files.  
#' @param lower_bound The lower bound of the scaled data range (default to -1).
#' @param upper_bound The upper bound of the scaled data range (default to 1).
#' @param cvlog_path_name The path and name of the log files, which hold the details of Liblinear procedures.
#' @param specificity_level  A number ranges from 0 to 1 indicating the specificity user requires the classifier to achieve (default to 0.99).
#' @param flag_for_score_threshold_chosen A string indicating whether use reference score threshold or get from the user supplied training data (default set to "reference").  
#' @param score_threshold A numerical value between 0 to 1 indicating the reference score threshold (required in "reference" mode).
#' @import stringr dplyr magrittr data.table
#' @export
#' @details This function outputs the features generated from input files.
#' @examples 
#' predict_on_targeted_proteome = function (ptm_site = "S", 
#'                                          flanking_size=12, 
#'                                          SPIDER = T,
#'                                          positive_info_file = "known_ps.tsv", 
#'                                          known_protein_fasta_file = "known_fasta.tsv",
#'                                          predict_protein_fasta_file = "predict_fasta.tsv",
#'                                          output_label_training = "ps_training",
#'                                          output_label_predict = "ps_predict",
#'                                          lower_bound = -1,
#'                                          upper_bound = 1,
#'                                          liblinear_dir = "/data/ginny/liblinear-2.11/",
#'                                          feature_file_path = "/data/ginny/test_package/",
#'                                          cvlog_path_name = "/data/ginny/test_package/cvlog.txt",
#'                                          specificity_level = 0.99,
#'                                          n_fold = 2,
#'                                          flag_for_score_threshold_chosen = "cv",
#'                                          score_threshold = NULL)



predict_on_targeted_proteome = function (ptm_site, flanking_size=12, 
                                         SPIDER = T,
                                         positive_info_file, 
                                         known_protein_fasta_file, predict_protein_fasta_file,
                                         output_label_training, output_label_predict,
                                         lower_bound = -1, upper_bound = 1,
                                         liblinear_dir,
                                         feature_file_path,
                                         cvlog_path_name,
                                         specificity_level,
                                         n_fold = 2,
                                         flag_for_score_threshold_chosen = "reference",
                                         score_threshold)
{
  
  generate_feature_t(ptm_site = ptm_site,
                     flanking_size = flanking_size, 
                     SPIDER = SPIDER,
                     positive_info_file = positive_info_file,
                     known_protein_fasta_file = known_protein_fasta_file,
                     predict_protein_fasta_file = predict_protein_fasta_file,
                     output_label_training = output_label_training,
                     output_label_predict = output_label_predict)
  
  
  cat("STEP1: Feature generated.", "\n")                   
  
  generate_training_test_data_t(lower_bound = lower_bound,
                         upper_bound = upper_bound,
                         positive_training_feature_file = paste0(output_label_training,"_noc_not_na_pos_feature.Rds"), 
                         negative_training_feature_file = paste0(output_label_training,"_noc_not_na_candi_feature.Rds"),
                         negative_predict_feature_file = paste0(output_label_predict, "_noc_not_na_candi_feature.Rds"),
                         output_label_training = output_label_training,
                         output_label_predict = output_label_predict)
  
  cat("STEP2: Training and predict prepared.", "\n")
  
  predict_with_liblinear(liblinear_dir = liblinear_dir,
                         feature_file_path = feature_file_path,
                         training_file_name = paste0(output_label_predict,"_train_feature_svm.tsv"),
                         test_file_name = paste0(output_label_predict,"_test_feature_svm.tsv"),
                         output_label = output_label_predict,
                         cvlog_path_name = cvlog_path_name)
  
  cat("STEP3: Liblinear processed.", "\n")
  
  
  present_prediction_t(flag_for_score_threshold_chosen = flag_for_score_threshold_chosen,
                              score_threshold = score_threshold,
                              ptm_site = ptm_site, 
                              flanking_size = flanking_size,
                              pred_candidate_df_Rds_name = paste0(output_label_predict,"_candidate.Rds"),
                              pred_score_file_name = paste0(output_label_predict,"_predict.tsv"),
                              positive_info_file = positive_info_file,
                              known_protein_fasta_file = known_protein_fasta_file,
                              n_fold = n_fold,
                              lower_bound = lower_bound,
                              upper_bound = upper_bound,
                              liblinear_dir = liblinear_dir,
                              feature_file_path = feature_file_path,
                              cvlog_path_name = cvlog_path_name,
                              specificity_level = specificity_level,
                              output_label = output_label_predict)
  
  
  cat("STEP4: Prediction results combined.", "\n")
  
  
  PTM_domain_mapping_enrichment(window_score_label_Rds = paste0(output_label_predict,"_window_score_df.Rds"), 
                                output_label = output_label_predict)
  cat("STEP5: Domain mapped and annotated.", "\n")
  
  
  ### insert code here to remove redundant files 
  
  ### keep only _test.tsv and _test.Rds
  
  
  
  to_delete_matches = c("feature", "[0-9]", "matrix","match","score","candi",
                        "pssm","tbt","names")
  
  for(i in 1:length(to_delete_matches))
  {
    
    match_string = paste0(output_label_training,"_*",to_delete_matches[i],"*")
    rm_cmd = paste0("find -type f -name '", match_string,  "' -delete")
    
    system(rm_cmd)
    
    match_string = paste0(output_label_predict,"_*",to_delete_matches[i],"*")
    rm_cmd = paste0("find -type f -name '", match_string,  "' -delete")
    
    system(rm_cmd)
    
  
    
  }
  
  
  
  
  
  
  cat("Finished!","\n")
}
          



