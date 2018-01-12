

get_phosphoSVM_site = function(s_ID,directory, mod_site,output_seq_name,output_site_name)
{
  
  both_id = intersect(s_ID$V1, sp_uni_id)
  which_in_sp = which(sp_uni_id%in%both_id)
  sel_seq = sp_seq[which_in_sp]
  sel_id = sp_uni_id[which_in_sp]
  
  sel_seq_id = rep(0,2*length(sel_id))
  
  for(i in 1:length(sel_id))
  {
    sel_seq_id[2*i-1] = sel_id[i]
    sel_seq_id[2*i] = sel_seq[i]
    
  }
  
  write.table(sel_seq_id, output_seq_name, quote = F, row.names = F, sep = "\t")
  
  
  
  ### so the sequence I'm uisng will be from both_id 
  ### extract the sites 
  
  
  # so I readin the filenames of all the seqs and select to read in the real files of the sel_id only
  
  
  sel_sites =lapply(1:length(sel_id),function(i) {
    make_name = paste0(sel_id[i], ".txt")
    get_site = data.table::fread(paste0(directory,make_name), 
                     stringsAsFactors = F, header = F)
    return(get_site)
    
  })
  
  
  ### change this to the format of files in Parse_PSP
  
  
  site_frame = data.frame(rbindlist(lapply(1:length(sel_id), function(i)
  {
    ptm_position = which(sel_sites[[i]]==1)
    nrow = length(ptm_position)
    ptm_ID=rep(sel_id[i], nrow)
    ptm_type = rep(mod_site, nrow)
    this_frame = data.frame(ptm_ID, ptm_position, ptm_type)
    return(this_frame)
  })))
  
  write.table(site_frame, output_site_name, quote = F, row.names = F, sep = "\t")
  
}




### a function to select the set of proteins for each different types of modification 


#whole_sp_seq_filename = "sp_0814_fasta.tsv"
#the_site = "S"

get_wanted_proteins = function(the_site, whole_sp_seq_filename)
{
  
  all_sp = readLines(whole_sp_seq_filename)
  sp_id = all_sp[c(T,F)]
  
  cat(length(sp_id), "\n")
  sp_seq = all_sp[c(F,T)]
  #sp_uni_id = sapply(1:length(sp_id),function(i) strsplit(sp_id[i], split = "|", fixed = T)[[1]][2])
  
 have_ind =  sapply(1:length(sp_seq), function(i) {
    grepl(the_site, sp_seq[i])
  })
  
  
  
  have_id = sp_id[have_ind]
  have_seq = sp_seq[have_ind]
  
  cat(length(have_id), "\n")
  
  
  out_fasta = rep("",2*length(have_id))
  
  for(i in 1:length(have_id))
  {
    
    out_fasta[2*i-1]= have_id[i]
    out_fasta[2*i] = have_seq[i]
    
    
  }
  
  write.table(out_fasta, paste0(the_site,"_sp_fasta.tsv"), row.names = F, quote = F, sep = "\t")
  
  
}




get_fasta_multi_lines_seq = function(all_sp)
{
  
  sp_id = rep("",0)
  sp_seq = rep("",0)
  seq = ""
  
  for(i in 1:(length(all_sp)+1))
  {
    if(substr(all_sp[i],1,1)==">"| i == (length(all_sp)+1))
    {
      sp_seq = c(sp_seq, seq)
      sp_id = c(sp_id,all_sp[i])
      seq = ""
    }else{
      seq = paste0(seq, all_sp[i])
    }
    
  }
  
  # remove the first line of seq and the last line of id
  sp_seq = sp_seq[-1]
  sp_id = sp_id[-length(sp_id)]
  
  return(list(sp_id,sp_seq))
  
}









get_z_transform = function(tmatrix)
{
  
  each_col_mean = apply(tmatrix, 2, mean)
  each_col_sd = apply(tmatrix, 2,sd)
  
  tmatrix_new=sapply(1:ncol(tmatrix), function(i)  {
    (tmatrix[,i]-each_col_mean[i])/each_col_sd[i]
  })
  
  
  return(tmatrix_new)
}


### build a function for each sequence 
build_candidate_df = function(cand_site, flank_size, one_seq, one_id, ps_info)
{
  
  one_seq_vec = unlist(strsplit(one_seq, split = ""))
  cand_pos = which(one_seq_vec == cand_site)
  
  ### get the candidate windows
  ### consider the head and tail senario 
  
  ### I want to change the seq so that it is easier
  
  modseq = paste0(paste0(rep("X",flank_size),collapse = ""),
                  one_seq, 
                  paste0(rep("X",flank_size),collapse = ""))
  ### since after add X of flank_size at the begining and end of the seq,
  ### the candidate position will need to be shift to the right of flank_size,
  ###so it should be (can_pos+flank_size-flank_size, can_pos+flank_size+flank_size)
  
  cand_windows = stringr::str_sub(modseq, cand_pos,cand_pos+2*flank_size)
  
  psp_label = rep("negative", length(cand_pos))
  
  sel_ps_info = ps_info %>% dplyr::filter(ptm_ID == one_id)
  
  psped = which(cand_pos %in% sel_ps_info$ptm_position)
  
  psp_label[psped] = "positive"
  
  cand_df = data.frame(protID = rep(one_id,length(cand_pos)),
                       cand_pos,
                       cand_windows,
                       psp_label, stringsAsFactors = F)
  
  return(cand_df)
  
}


#### the following is built for seqs without any positive site information

build_candidate_df_no_positive = function(cand_site, flank_size, one_seq, one_id)
{
  
  one_seq_vec = unlist(strsplit(one_seq, split = ""))
  cand_pos = which(one_seq_vec == cand_site)
  
  ### get the candidate windows
  ### consider the head and tail senario 
  
  ### I want to change the seq so that it is easier
  
  modseq = paste0(paste0(rep("X",flank_size),collapse = ""),
                  one_seq, 
                  paste0(rep("X",flank_size),collapse = ""))
  ### since after add X of flank_size at the begining and end of the seq,
  ### the candidate position will need to be shift to the right of flank_size,
  ###so it should be (can_pos+flank_size-flank_size, can_pos+flank_size+flank_size)
  
  cand_windows = stringr::str_sub(modseq, cand_pos,cand_pos+2*flank_size)
  
  psp_label = rep("negative", length(cand_pos))
  
  cand_df = data.frame(protID = rep(one_id,length(cand_pos)),
                       cand_pos,
                       cand_windows,
                       psp_label, stringsAsFactors = F)
  
  return(cand_df)
  
}





get_all_candidate = function(cand_site, flank_size, sp_seq, sp_uni_id, ps_info)
{
  all_cand = data.frame(rbindlist(lapply(1:length(sp_seq), function(i) {
    cdf = build_candidate_df(cand_site = cand_site, 
                             flank_size = flank_size, 
                             one_seq = sp_seq[i], 
                             one_id = sp_uni_id[i], 
                             ps_info = ps_info)
    ncdf = data.frame(prot = rep(i, nrow(cdf)),
                      protID = cdf$protID,
                      pos = cdf$cand_pos,
                      window = cdf$cand_windows,
                      label = cdf$psp_label, stringsAsFactors = F)
    # if (i %% 1000 ==0 )
    # {
    #   cat(i)
    #   cat("\n")
    # }
    
    return(ncdf)
  })), stringsAsFactors = F)
  
  return (all_cand)
  
}




get_all_candidate_no_positive = function(cand_site, flank_size, sp_seq, sp_uni_id)
{
  all_cand = data.frame(rbindlist(lapply(1:length(sp_seq), function(i) {
    cdf = build_candidate_df_no_positive(cand_site = cand_site, 
                             flank_size = flank_size, 
                             one_seq = sp_seq[i], 
                             one_id = sp_uni_id[i])
    ncdf = data.frame(prot = rep(i, nrow(cdf)),
                      protID = cdf$protID,
                      pos = cdf$cand_pos,
                      window = cdf$cand_windows,
                      label = cdf$psp_label, stringsAsFactors = F)
    # if (i %% 1000 ==0 )
    # {
    #   cat(i)
    #   cat("\n")
    # }
    # 
    return(ncdf)
  })), stringsAsFactors = F)
  
  return (all_cand)
  
}





# away_center_size < flank_size

build_decoy_df = function(cand_site, flank_size, awayc_size, one_seq, one_id)
{
  
  
  
  
  one_seq_vec = unlist(strsplit(one_seq, split = ""))
  other_pos = which(one_seq_vec != cand_site)
  
  modseq = paste0(paste0(rep("X",flank_size),collapse = ""),
                  one_seq, 
                  paste0(rep("X",flank_size),collapse = ""))
  
  
  other_windows = stringr::str_sub(modseq, other_pos,other_pos+2*flank_size)
  other_center_windows = stringr::str_sub(modseq, 
                                 other_pos + flank_size - awayc_size,
                                 other_pos + flank_size + awayc_size)
  
  
  in_window = which(grepl(cand_site, other_windows))
  outof_center = which(!(grepl(cand_site, other_center_windows)))
  
  in_decoy = intersect(in_window, outof_center)
  
  decoy_windows = other_windows[in_decoy]
  decoy_pos = other_pos[in_decoy]
  
  
  decoy_label = rep("decoy", length(decoy_pos))
  
  # sel_ps_info = ps_info %>% dplyr::filter(ptm_ID == one_id)
  # 
  # psped = which(cand_pos %in% sel_ps_info$ptm_position)
  # 
  # psp_label[psped] = "positive"
  
  
  decoy_df = data.frame(protID = rep(one_id,length(decoy_pos)),
                        decoy_pos,
                        decoy_windows,
                        decoy_label, stringsAsFactors = F)
  
  return(decoy_df)
  
}




get_all_decoy = function(cand_site, flank_size,awayc_size, sp_seq, sp_uni_id)
{
  all_decoy = data.frame(rbindlist(lapply(1:length(sp_seq), function(i) {
    ddf = build_decoy_df(cand_site = cand_site, 
                         flank_size = flank_size, 
                         awayc_size = awayc_size,
                         one_seq = sp_seq[i], 
                         one_id = sp_uni_id[i])
    nddf = data.frame(prot = rep(i, nrow(ddf)),
                      protID = ddf$protID,
                      pos = ddf$decoy_pos,
                      window = ddf$decoy_windows,
                      label = ddf$decoy_label, stringsAsFactors = F)
    # if (i %% 1000 ==0 | i == 20170)
    # {
    #   cat(i)
    #   cat("\n")
    # }
    # 
    return(nddf)
  })), stringsAsFactors = F)
  
  return (all_decoy)
  
}




### build the pure decoy to see
### without S

build_pure_decoy_df = function(cand_site, flank_size, one_seq, one_id)
{
  
  #one_seq = sp_seq[1]
  #flank_size= 5
  #cand_site = "S"
  
  
  
  one_seq_vec = unlist(strsplit(one_seq, split = ""))
  other_pos = which(one_seq_vec != cand_site)
  
  modseq = paste0(paste0(rep("X",flank_size),collapse = ""),
                  one_seq, 
                  paste0(rep("X",flank_size),collapse = ""))
  
  
  other_windows = stringr::str_sub(modseq, other_pos,other_pos+2*flank_size)
  
  
  out_window = which(grepl(cand_site, other_windows)==F)
  
  ### a tip maybe, when you take subset, take the ne included, do not do exclude ones not included.
  decoy_windows = other_windows[out_window]
  decoy_pos = other_pos[out_window]
  
  
  
  decoy_label = rep("pure_decoy", length(decoy_pos))
  
  # sel_ps_info = ps_info %>% dplyr::filter(ptm_ID == one_id)
  # 
  # psped = which(cand_pos %in% sel_ps_info$ptm_position)
  # 
  # psp_label[psped] = "positive"
  
  
  pure_decoy_df = data.frame(protID = rep(one_id,length(decoy_pos)),
                             decoy_pos,
                             decoy_windows,
                             decoy_label,
                             stringsAsFactors = F)
  
  return(pure_decoy_df)
  
}



get_all_pure_decoy = function(cand_site, flank_size, sp_seq, sp_uni_id)
{
  all_pure_decoy = data.frame(rbindlist(lapply(1:length(sp_seq), function(i) {
    ddf = build_pure_decoy_df(cand_site = cand_site, 
                              flank_size = flank_size, 
                              one_seq = sp_seq[i], 
                              one_id = sp_uni_id[i])
    nddf = data.frame(prot = rep(i, nrow(ddf)),
                      protID = ddf$protID,
                      pos = ddf$decoy_pos,
                      window = ddf$decoy_windows,
                      label = ddf$decoy_label,
                      stringsAsFactors = F)
    # if (i %% 1000 ==0 | i == 20170)
    # {
    #   cat(i)
    #   cat("\n")
    # }
    # 
    return(nddf)
  })), stringsAsFactors = F)
  
  return (all_pure_decoy)
  
}






#### I want to build the purest decoy , for each purest decoy window the window does not
#### contain any fragment of candidate and PSP windows


build_purest_decoy_df = function(cand_site, flank_size, one_seq, one_id)
{
  
  #one_seq = sp_seq[1]
  #flank_size= 5
  #cand_site = "S"
  
  clear_size = 2*flank_size
  
  
  one_seq_vec = unlist(strsplit(one_seq, split = ""))
  other_pos = which(one_seq_vec != cand_site)
  
  modseq = paste0(paste0(rep("X",clear_size),collapse = ""),
                  one_seq, 
                  paste0(rep("X",clear_size),collapse = ""))
  
  
  other_windows = stringr::str_sub(modseq, other_pos,other_pos+2*clear_size)
  
  
  out_window = which(grepl(cand_site, other_windows)==F)
  
  ### a tip maybe, when you take subset, take the ne included, do not do exclude ones not included.
  decoy_windows_clear = other_windows[out_window]
  decoy_pos = other_pos[out_window]
  
  
  
  decoy_windows = stringr::str_sub(decoy_windows_clear, flank_size+1, flank_size+1+2*flank_size)
  
  decoy_label = rep("purest_decoy", length(decoy_pos))
  
  # sel_ps_info = ps_info %>% dplyr::filter(ptm_ID == one_id)
  # 
  # psped = which(cand_pos %in% sel_ps_info$ptm_position)
  # 
  # psp_label[psped] = "positive"
  
  
  pure_decoy_df = data.frame(protID = rep(one_id,length(decoy_pos)),
                             decoy_pos,
                             decoy_windows,
                             decoy_label,
                             stringsAsFactors = F)
  
  return(pure_decoy_df)
  
}


get_all_purest_decoy = function(cand_site, flank_size, sp_seq, sp_uni_id)
{
  all_pure_decoy = data.frame(rbindlist(lapply(1:length(sp_seq), function(i) {
    ddf = build_purest_decoy_df(cand_site = cand_site, 
                              flank_size = flank_size, 
                              one_seq = sp_seq[i], 
                              one_id = sp_uni_id[i])
    nddf = data.frame(prot = rep(i, nrow(ddf)),
                      protID = ddf$protID,
                      pos = ddf$decoy_pos,
                      window = ddf$decoy_windows,
                      label = ddf$decoy_label,
                      stringsAsFactors = F)
    # if (i %% 1000 ==0 | i == 20170)
    # {
    #   cat(i)
    #   cat("\n")
    # }
    # 
    return(nddf)
  })), stringsAsFactors = F)
  
  return (all_pure_decoy)
  
}





### get the matrix of window-property
### physical-chemical features


get_matrix_all = function(all_windows, aaindex_com)
{
  ## in which order the AA is arranged?
  
  # I think it is actually better to provide AAorder myself 
  # to hard code it 
  
  AAorder = c("A","R","N","D","C","Q","E","G","H","I",
              "L","K","M","F","P","S","T","W","Y","V","X","U")
  
  #AAorder = c(colnames(aaindex_com)[4:ncol(aaindex_com)],"X","U")
  
  ## get the frequency matrix for all the windows
  aafreq_mat = sapply(1:length(all_windows), function(i) {
    one_window_vec = unlist(strsplit(all_windows[i], split = ""))
    freq = sapply(1:length(AAorder), function(x) length(which(one_window_vec == AAorder[x])))
    return(freq)
  })
  ## rearrange the property matrix
  aaprop_mat  = as.matrix(data.frame(aaindex_com, 
                                     X = rep(0,nrow(aaindex_com)),
                                     U = rep(0,nrow(aaindex_com)))%>%dplyr::select(AAorder))
  
  ## how many AA in each window?
  aasum_vec = apply(aafreq_mat[1:20,],2,sum)
  
  ## do calculation
  multi = aaprop_mat%*%aafreq_mat
  output =  sweep(multi, 2, aasum_vec, "/")
  
  output_matrix = as.matrix(t(output))
  #colnames(output_df) =  aaindex_com$header
  
  #write.table(output_df, outputname, quote =F, row.names = F, sep = "\t")
  
  
  
  return(output_matrix)
  
}




### get the frequency matrix for pssm calculation
get_base_freq_matrix = function(try_windows)
{
  
  AAorder = c("A","R","N","D","C","Q","E","G","H","I",
              "L","K","M","F","P","S","T","W","Y","V","X","U")
  
  freq_list = lapply(1:length(try_windows), function(i){
    
    one_window = unlist(strsplit(try_windows[i],
                                 split = ""))
    sapply(1:length(AAorder), function(x) {
      return(as.numeric(one_window ==AAorder[x]))
    })
  })
  
  count = t(Reduce('+',freq_list))
  
  len = length(try_windows)
  
  freq = apply(count, c(1,2), function(t) t/len)
  
  return(freq)
}


get_pssm_feature = function(pssm_matrix, this_window)
{
  
  
  AAorder = c("A","R","N","D","C","Q","E","G","H","I",
              "L","K","M","F","P","S","T","W","Y","V","X","U")
  
  this_window_vec = unlist(strsplit(this_window,split = ""))
  
  pssm_feature = sapply(1:length(this_window_vec), function(x) {
    
    which_row = which(AAorder == this_window_vec[x])
    
    return(pssm_matrix[which_row, x])
    
  })
  
  return(pssm_feature)
  
}


# a function to get the train and test

get_train_test = function(positive_feature, decoy_feature, 
                          outputname_train_data, outputname_train_label,
                          outputname_test_data, outputname_test_label)
{
  
  nr_pos = nrow(positive_feature)
  nr_decoy = nrow(decoy_feature)
  
  slice_pos = floor(nr_pos/2)
  slice_decoy = floor(nr_decoy/2)
  
  pos_ind= get_slice_rest_subset_index(nr_pos,
                                       slice_pos,
                                       (nr_pos-slice_pos))
  decoy_ind = get_slice_rest_subset_index(nr_decoy,
                                          slice_decoy,
                                          (nr_decoy-slice_decoy))
  
  
  ### I want to shuffle before getiing slice
  
  # positive_feature = all_positive_df
  # decoy_feature = all_decoy_df
  # 
  # nr_pos = nrow(positive_feature)
  # nr_decoy = nrow(decoy_feature)
  # 
  set.seed(123)
  shuf_positive_df = positive_feature[sample(nr_pos),]
  

  set.seed(123)
  shuf_decoy_df = decoy_feature[sample(nr_decoy),]
  
  
  
  ### getting slices and rests
  
  
  ps_train_data =as.matrix(rbind(shuf_positive_df[pos_ind[[1]],],
                                 shuf_decoy_df[decoy_ind[[1]],]))
  
  ps_train_label = c(rep(1, length(pos_ind[[1]])), 
                     rep(0, length(decoy_ind[[1]])))
  
  ps_test_data =as.matrix(rbind(shuf_positive_df[pos_ind[[2]],],
                                shuf_decoy_df[decoy_ind[[2]],]))
  
  ps_test_label = c(rep(1, length(pos_ind[[2]])), 
                    rep(0, length(decoy_ind[[2]])))
  
  
  
  write.table(ps_train_data, outputname_train_data, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  write.table(ps_train_label, outputname_train_label, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  
  
  write.table(ps_test_data, outputname_test_data, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  write.table(ps_test_label, outputname_test_label,
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  
  
}

## same size of data in training and testing 


get_subset_train_test = function(positive_feature, decoy_feature, 
                                 positive_size, decoy_size,
                          outputname_train_data, outputname_train_label,
                          outputname_test_data, outputname_test_label)
{
  
  nr_pos = nrow(positive_feature)
  nr_decoy = nrow(decoy_feature)
  
  #slice_pos = floor(nr_pos/2)
  #slice_decoy = floor(nr_decoy/2)
  
  pos_ind= get_slice_rest_subset_index(nr_pos,
                                       positive_size,
                                       positive_size)
  decoy_ind = get_slice_rest_subset_index(nr_decoy,
                                          decoy_size,
                                          decoy_size)
  
 
  set.seed(123)
  shuf_positive_df = positive_feature[sample(nr_pos),]
  
  
  set.seed(123)
  shuf_decoy_df = decoy_feature[sample(nr_decoy),]
  
  
  
  ### getting slices and rests
  
  
  ps_train_data =as.matrix(rbind(shuf_positive_df[pos_ind[[1]],],
                                 shuf_decoy_df[decoy_ind[[1]],]))
  
  ps_train_label = c(rep(1, length(pos_ind[[1]])), 
                     rep(0, length(decoy_ind[[1]])))
  
  ps_test_data =as.matrix(rbind(shuf_positive_df[pos_ind[[2]],],
                                shuf_decoy_df[decoy_ind[[2]],]))
  
  ps_test_label = c(rep(1, length(pos_ind[[2]])), 
                    rep(0, length(decoy_ind[[2]])))
  
  
  
  write.table(ps_train_data, outputname_train_data, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  write.table(ps_train_label, outputname_train_label, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  
  
  write.table(ps_test_data, outputname_test_data, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  write.table(ps_test_label, outputname_test_label,
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  
  
}




### get subset of dataset for crossvalidation

get_slice_rest_subset_index = function(seq_length, slice_length, rest_subset_length)
{
  # seq_length = 100
  # slice_length = 10
  # rest_subset_length = 10
  
  
  set.seed(123)
  slice_ind = sample(seq_length, slice_length)
  rest_ind = seq(1:seq_length)[-slice_ind]
  
  
  set.seed(123)
  rest_subset_ind = sample(rest_ind, rest_subset_length)
  
  cat(slice_ind[1:10])
  cat("\n")
  cat(rest_subset_ind[1:10])
  cat("\n")
  
  cat("see if training and test overlaps",length(intersect(slice_ind, rest_subset_ind)))
  cat("\n")
  
  
  return (list(slice_ind, rest_subset_ind))
}


###

get_subset_train_test_no_shuffle = function(positive_feature, decoy_feature, 
                                 positive_size, decoy_size,
                                 outputname_train_data, outputname_train_label,
                                 outputname_test_data, outputname_test_label)
{
  
  nr_pos = nrow(positive_feature)
  nr_decoy = nrow(decoy_feature)
  
  #slice_pos = floor(nr_pos/2)
  #slice_decoy = floor(nr_decoy/2)
  
  pos_ind= get_slice_rest_subset_index(nr_pos,
                                       positive_size,
                                       positive_size)
  decoy_ind = get_slice_rest_subset_index(nr_decoy,
                                          decoy_size,
                                          decoy_size)
  
  
 
  shuf_positive_df = positive_feature
  
  shuf_decoy_df = decoy_feature
  
  
  
  ### getting slices and rests
  
  
  ps_train_data =as.matrix(rbind(as.matrix(shuf_positive_df[pos_ind[[1]],]),
                                 as.matrix(shuf_decoy_df[decoy_ind[[1]],])))
  
  ps_train_label = c(rep(1, length(pos_ind[[1]])), 
                     rep(0, length(decoy_ind[[1]])))
  
  ps_test_data =as.matrix(rbind(as.matrix(shuf_positive_df[pos_ind[[2]],]),
                                as.matrix(shuf_decoy_df[decoy_ind[[2]],])))
  
  ps_test_label = c(rep(1, length(pos_ind[[2]])), 
                    rep(0, length(decoy_ind[[2]])))
  
  
  
  write.table(ps_train_data, outputname_train_data, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  write.table(ps_train_label, outputname_train_label, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  
  
  write.table(ps_test_data, outputname_test_data, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  write.table(ps_test_label, outputname_test_label,
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  
  
}



get_subset_train_test_no_shuffle_windows = function(positive_windows, decoy_windows, 
                                                    positive_size, decoy_size,
                                                    outputname_train_data, outputname_train_label,
                                                    outputname_test_data, outputname_test_label)
{
  
  # positive_windows = pos_notNA_windows
  # decoy_windows = decoy_notNA_windows
  # positive_size = half_pos_len
  # decoy_size = half_pos_len
  
  nr_pos = nrow(positive_windows)
  nr_decoy = nrow(decoy_windows)
  
  #slice_pos = floor(nr_pos/2)
  #slice_decoy = floor(nr_decoy/2)
  
  pos_ind= get_slice_rest_subset_index(nr_pos,
                                       positive_size,
                                       positive_size)
  decoy_ind = get_slice_rest_subset_index(nr_decoy,
                                          decoy_size,
                                          decoy_size)
  
  
  
  shuf_positive_df = positive_windows
  
  shuf_decoy_df = decoy_windows
  
  
  
  
  ### getting slices and rests
  
  
  ps_train_data =c(shuf_positive_df[pos_ind[[1]],],
                   shuf_decoy_df[decoy_ind[[1]],])
  
  ps_train_label = c(rep(1, length(pos_ind[[1]])), 
                     rep(0, length(decoy_ind[[1]])))
  
  ps_test_data =c(shuf_positive_df[pos_ind[[2]],],
                  shuf_decoy_df[decoy_ind[[2]],])
  
  ps_test_label = c(rep(1, length(pos_ind[[2]])), 
                    rep(0, length(decoy_ind[[2]])))
  
  
  
  write.table(ps_train_data, outputname_train_data, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  write.table(ps_train_label, outputname_train_label, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  
  
  write.table(ps_test_data, outputname_test_data, 
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  write.table(ps_test_label, outputname_test_label,
              quote = F, row.names = F,col.names = F,sep = "\t")
  
  
  
}






## remove redundancy (final implementation is in C++)
## I will revert this part to RCPP eventually 

## the input is the set of vectorized windows. 


get_non_redundant = function(al)
{
  seq = al[1]
  tocompare = al[-1]
  to_put = seq
  #last_discard = T
  while(length(tocompare)>0)
  {
    getsames = sapply(tocompare,function(x) sum(x == unlist(seq)))
    
    
    to_discard = which(getsames>=7)
    
    
    
    if(length(to_discard)>0)
    {
      beforetocompare = tocompare[-to_discard]
      
    }else{
      beforetocompare = tocompare
    }
    
    seq = beforetocompare[1]
    
    to_put = c(to_put, seq)
    
    
    tocompare = beforetocompare[-1]
  }
  
  return(to_put)
  
  
}


## already zoom into same protein
## for each window

# 
# 
# get_structural_each_window = function(position, protein_structure )
# {
#   # drop one of the three SS in the last step
#   structure_prop_matrix = as.matrix(protein_structure[,5:14])
#   
#   # get start and end position of each aa in the window
#   window_start = max(0, position-5)
#   window_end = min(nrow(protein_structure), position+5)
#   
#   this_sub = structure_prop_matrix[window_start:window_end,]
#   
#   prop_vect = apply(this_sub, 2, mean)
#   
#   return(prop_vect)
#   
# }
# 
# 
# # for the whole protein 
# # 
# # 
# # positions = ps_try$cand_pos
# # windows = ps_try$cand_window
# # 
# # 
# # get_structural_each_protein = function(protID, positions, windows, protein_structure)
# # {
# #   prot_length = length(positions)
# #   
# #   
# #   prop_vect_prot = sapply(1:prot_length, function(x)
# #     {
# #     return(get_structural_each_window(positions[x], windows[x], protein_structure))
# #     
# #     })
# #   
# #   return(t(prop_vect_prot))
# #   
# # }
# # 
# # 
# # 
# # 
# # get_structural_each_window(position, window, protein_structure)
# 
# 
# get_structure_all = function(ps_data, protein_structure_info)
# {
#   all_vec = sapply(1:nrow(ps_data), function(i){
#     
#     this_protein = protein_structure_info%>%dplyr::filter(prot_ID == ps_data$protID[i])
#     
#     if(i%%10000 == 0)
#     {
#       cat(i)
#       cat("\n")
#     }
#     
#     
#     if(nrow(this_protein)>0)
#     {
#       
#       whichcol = grep("pos", colnames(ps_data))
#       this_vec = get_structural_each_window(ps_data[[whichcol]][i],
#                                             this_protein)
#     }else{
#       this_vec = rep(NA, 10)
#     }
#     
#     
#     return(this_vec)                                      
#     
#   })
#   
#   return(t(all_vec))
#   
# }
# 
# 
# 
# 
# 
# 



### get all the combined spd33 files given a folder of files

get_spd33_combined = function(filenames)
{
  combined_data = data.frame(rbindlist(lapply ( 1:length(filenames), function(i) {
    tmp_data = data.table::fread(filenames[i], stringsAsFactors = F, header = T)
    tmp_prot_name = strsplit(filenames[i], split = "_")[[1]][2]
    
    # cat(i)
    # cat("\t")
    # cat(nrow(tmp_data))
    # cat("\n")
    prot_name = rep(tmp_prot_name, nrow(tmp_data))
    new_data = cbind(prot_name, tmp_data)
    
    return(new_data)
    
    
  })))
  
  return(combined_data)
}






### to change the original file name first

# 
# 
# colnames(AOPspd33) =c("prot_ID","position","aa",
#                       "secondary_structure",
#                       "ASA","Phi","Psi","Theta","Tau",
#                       "HSE_au","HSE_ad","PC","PH","PE")
# 

## the final format will be 
## prot_ID, position(of the center aa), average property of that window




### the following idea is to build the average property first and then do a simple matching 
### against the pos/candidate/decoy data 




get_structure_average = function(Aspd33)
{
  all_unique_protein = unique(Aspd33$prot_ID)
  
  tit =rbindlist(lapply(1:length(all_unique_protein), function(x)
  {
    
    # if(x%%1000 ==0)
    #   cat(x,"\n")
    # 
    this_protein_data = Aspd33%>%dplyr::filter(prot_ID==all_unique_protein[x])
    
    this_structure = sapply(1:nrow(this_protein_data), function(i) 
    {
      start = max(1,i-5)
      end = min(nrow(this_protein_data), i+5)
      
      sub_protein_data = as.matrix(this_protein_data[start:end, 5:14])
      
      structure_ave = apply(sub_protein_data, 2, mean)
      
      return(structure_ave)
      
    })
    
    return(as.data.frame(t(this_structure)))
    
  }))
  
  return(tit)
  
}



### 

### instead of averaging, get the 10*11 dimension feature for each site(corresponding to a window)
### if some sites are not surrounded by +/-5 I will just put x number of 0s 

get_structure_site_specific = function(Aspd33,flank_size, protID, pos)
{
  all_unique_protein = unique(Aspd33$prot_ID)
  
  tit =rbindlist(lapply(1:length(all_unique_protein),function(x)
  {
    # 
    # if(x%%1000 ==0)
    #   cat(x,"\n")
    # #x =1 
    this_protein_data = Aspd33%>%dplyr::filter(prot_ID==all_unique_protein[x])
    
    this_structure = sapply(1:nrow(this_protein_data), function(i) 
    {
      ### the assumption is that each prot seq starts from position 1 in this program
     
      #i = 1 
       pos = i
      start = max(1,pos-flank_size)
      end = min(max(this_protein_data$position), pos+flank_size)
    
      sets_to_fill_at_start = flank_size-abs(pos-start)
      sets_to_fill_at_end = flank_size-abs(end-pos)
      
      
      sub_protein_data = as.matrix(this_protein_data[start:end, 5:14])
      
      ### 5:14 is from the original columns of the table.
      
      property_vec = as.vector(t(sub_protein_data))
      

      vec_median_properties = as.vector(apply(sub_protein_data,2,median))
      
      
      full_property_vec = c(rep(vec_median_properties,sets_to_fill_at_start),
                            property_vec, 
                            rep(vec_median_properties,sets_to_fill_at_end))
      
      
      
      return(full_property_vec)
      
    })
    
    return(as.data.frame(t(this_structure)))
    
  }))
  
  
  tit = data.frame(protID, pos, tit, stringsAsFactors = F)
  

  return(tit)
  
}

  


get_prot_mean_structure = function(full_spider)
{
  
  
  all_prots = unique(full_spider$protID)
  
 get_all_mean = sapply(1:length(all_prots), function(x) {
    sub_spider <- full_spider%>%
      dplyr::filter(protID == all_prots[x])%>%
      dplyr::select(-c(protID, pos))
    
    this_prot_mean = colMeans(sub_spider)
    
    return(this_prot_mean)
    
  })
 
 get_all_mean_df = data.frame(protID = all_prots, t(get_all_mean))
 
 
 return(get_all_mean_df)
 
 
}



### modify this function so it can take the center positon as the argument


extract_site_specific_features=function(feature_data, to_extract, to_logit, center_position)
{
  
  
  
  #feature_data = not_na_pos_structure_site_specific[1:10,]
  
 # center_position = 8
  
#  to_extract = c(1,6,8,9)
 # to_logit = c(8,9)
  
  total_feature_length = 2*(center_position-1)+1
  
  
  logbs = matrix(rep(seq(0,  10*(total_feature_length-1),10),length(to_logit)),nrow=total_feature_length,ncol=length(to_logit))
  add_logbs = sapply(1:length(to_logit), function(x) logbs[,x]+to_logit[x])
  vec_add_logbs = as.vector(add_logbs)
  
  logit_feature_data = feature_data
  
  for(i in 1:length(vec_add_logbs))
  {
    p1 = sapply(1:nrow(feature_data), function(x) 
      min(max(feature_data[x,vec_add_logbs[i]],0.001),0.999) )
    new_p1 = log(p1/(1-p1))
    
    logit_feature_data[,vec_add_logbs[i]] = new_p1
    
  }
  

  bs = matrix(rep(seq(0,  10*(total_feature_length-1),10),length(to_extract)),nrow=total_feature_length,ncol=length(to_extract))
  add_bs = sapply(1:length(to_extract), function(x) bs[,x]+to_extract[x])
  vec_add_bs = sort(as.vector(add_bs))
  
  
  e_feature = logit_feature_data[ ,vec_add_bs]
  
  
  return(e_feature)
}







### draw histogram to show each property
get_hist = function(nn_pos_structure_sse,nn_candi_structure_sse,nn_decoy_structure_sse,
                    fig_name)
{
  structure_properties = colnames(nn_pos_structure_sse)
  
  
  pdf(fig_name,useDingbats = F)
  
  par(mfrow = c(2,2))
  
  
  
  for(i in 1:ncol(nn_pos_structure_sse))
  {
    
    dens_positive = density(nn_pos_structure_sse[[i]], adjust = 1)
    dens_candidate = density(nn_candi_structure_sse[[i]], adjust = 1)
    dens_decoy =  density(nn_decoy_structure_sse[[i]], adjust = 1)
    
    common_max_y = max(max(dens_positive$y), max(dens_candidate$y), max(dens_decoy$y))
    common_max_x = max(max(dens_positive$x), max(dens_candidate$x), max(dens_decoy$x)) 
    common_min_x = min(min(dens_positive$x), min(dens_candidate$x), min(dens_decoy$x))
    
    
    
    hist(nn_pos_structure_sse[[i]], breaks = 1000,
         main = structure_properties[i],
         xlab = c("Red: positve"),
         cex.main = 0.5,
         cex.lab = 0.5,
         border = rgb(1,0,0,alpha = 1),
         col = rgb(1,0,0, alpha = 1))
    hist(nn_candi_structure_sse[[i]], breaks = 1000,
         main = structure_properties[i],
         xlab = c("Blue: candidate"),
         cex.main = 0.5,
         cex.lab = 0.5,
         border = rgb(0,0,1, alpha = 1),
         col = rgb(0,0,1, alpha = 1))
    hist(nn_decoy_structure_sse[[i]], breaks = 1000,
         main = structure_properties[i],
         xlab = c("Green: decoy"),
         cex.main = 0.5,
         cex.lab = 0.5,
         border = rgb(0,1,0, alpha = 1),
         col = rgb(0,1,0, alpha = 1))
    
    
    plot(dens_positive, col = "red", 
         xlim = c(common_min_x,common_max_x),
         ylim = c(0,common_max_y),
         main = structure_properties[i],
         xlab = c("Red: positve, Green: decoy, Blue: candidate"),
         cex.main = 0.5,
         cex.lab = 0.5
    )
    lines(dens_candidate,col = "blue")
    lines(dens_decoy, col = "green")
    
    
  }
  
  
  dev.off()
  
  
  
  
}



ggbiplot2 = function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                     ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                     alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                     varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                     color = muted("red"), linetype = "solid", alpha_arrow =1,point_size = 1,
                     ...) 
{
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, xend = xvar, yend = yvar), 
                          arrow = arrow(length = unit(1/2, "picas")), 
                          color = color, linetype = linetype, alpha = alpha_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha, size = point_size,stroke = 0,
                          shape=16)
    }
    else {
      g <- g + geom_point(alpha = alpha, size = point_size,stroke = 0,
                          shape=16)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar, angle = angle, hjust = hjust), 
                       color = "darkred", size = varname.size)
  }
  return(g)
}


### do pca and draw plot


#library(devtools)
#install_github("ggbiplot", "vqv")
#library(ggbiplot)


get_pca_plot = function(get_wd, label, outputname)
{
  #label = as.factor(get_type_label$x)
  
  get_wd_pca = prcomp(as.matrix(get_wd), center = TRUE, scale. = TRUE)
  
  ### the following is to make some modification so that the variable labels wont be shown 
  
  dimnames(get_wd_pca$rotation)[[1]] <- 
    Reduce(function(x,y) paste0(x,y),    
           rep(" ",dim(get_wd_pca$rotation)[1]),   
           acc=TRUE)       
  
  
  pdf(outputname, useDingbats = F)
  
  cols <- c("negative" = "blue", "decoy" = "green", "positive" = "red")
  
  g <- ggbiplot2(get_wd_pca, obs.scale = 1, var.scale = 1, 
                 groups = label, ellipse = TRUE, 
                 circle = F,
                 alpha = 0.4,
                 alpha_arrow = 0,
                 point_size = 0.8)
  
  g <- g + scale_color_manual(values = cols)
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top') 
  
  g<- g + theme_bw()
  
  
  print(g)
  
  dev.off()
  
  
  
}


get_pca_plot_two = function(get_wd, label, outputname)
{
  #label = as.factor(get_type_label$x)
  
  get_wd_pca = prcomp(as.matrix(get_wd), center = TRUE, scale. = TRUE)
  
  ### the following is to make some modification so that the variable labels wont be shown 
  
  dimnames(get_wd_pca$rotation)[[1]] <- 
    Reduce(function(x,y) paste0(x,y),    
           rep(" ",dim(get_wd_pca$rotation)[1]),   
           acc=TRUE)       
  
  
  pdf(outputname, useDingbats = F)
  
  cols <- c("negative" = "blue", "positive" = "red")
  
  g <- ggbiplot2(get_wd_pca, obs.scale = 1, var.scale = 1, 
                 groups = label, ellipse = TRUE, 
                 circle = F,
                 alpha = 0.6,
                 alpha_arrow = 0,
                 point_size = 0.9)
  
  g <- g + geom_hline(yintercept = 0)
  g <- g + geom_vline(xintercept = 0)
  
  g <- g + scale_color_manual(values = cols)
  #g <- g + theme(legend.direction = 'horizontal', 
                 #legend.position = 'top') 
  g <- g + theme_bw()
  
  print(g)
  
  dev.off()
  
  
  
}





get_tsne_pos_candi = function(get_wd, label, outputname)
{
  
  cols = rep("blue", length(label))
  
  which_pos = which(label == "positive")
  
  cols[which_pos] = "red"
 
  pdf(outputname, useDingbats = F)
  
  set.seed(123)
  get_wd_tsne = Rtsne(as.matrix(get_wd),pca_center = T,pca_scale = F, perplexity = 40)
  
  plot(get_wd_tsne$Y, col = alpha(cols,0.5), pch = 16, cex = 0.5)
  

  dev.off()
  
  
  
}




#### I will try plsda now



get_plsda_pos_candi = function(get_wd, label, outputname)
{
  
  # get_wd = sel_feature
  # label = type_label 
  #outputname = paste0(output_label,"_plsda_plot_two.pdf")
  # 
  colnames(get_wd) = as.character(c(1:ncol(get_wd)))
  
  get_plsda = plsda(X = get_wd, Y = label, ncomp = 2)
  
  
  pdf(outputname, useDingbats = F)

  cols = c( "blue","red")
  
  plotIndiv(get_plsda, ind.names = FALSE, ellipse = TRUE,
            ellipse.level = 0.75,
            legend = TRUE,
            col = alpha(cols,0.5), pch = 16, cex = 0.5)
  
  
  dev.off()
  
  
  
}








get_plsda_pos_candi_with_score_selection = function(get_wd, label,score_label, outputname)
{
  
  # get_wd = sel_feature
  # label = type_label 
  #outputname = paste0(output_label,"_plsda_plot_two.pdf")
  # 
  colnames(get_wd) = as.character(c(1:ncol(get_wd)))
  
  get_plsda = plsda(X = get_wd, Y = label, ncomp = 2)
  
  
  pdf(outputname, useDingbats = F)
  
  cols = c( "blue","red")
  
  plotIndiv(get_plsda, ind.names = FALSE, ellipse = TRUE,
            ellipse.level = 0.75,
            legend = TRUE,
            col = alpha(cols,0.5), pch = 16, cex = 0.5)
  
  
  dev.off()
  
  
  
}




#it = grep("pos", colnames(ps_pos))




whichpart <- function(x, n) {
  
  xp <- sort(x, partial=n)[n]
  return(which(x <= xp))
}

get_k_keep= function(one_point, target_points, k, pos_ind,bigs)
{
  alldist = rowSums(abs(matrix(one_point-as.vector(t(target_points)),
                               nrow = bigs, byrow = T )))
  nearest_ind = whichpart(alldist,k)
  keep = sum(nearest_ind%in%pos_ind)==0
  return(keep)
}


#### scale the data 

scale_train = function(feature_data, upper_bound, lower_bound, Rds_label)
{
  
  mins = apply(as.matrix(feature_data), 2, min)
  maxs = apply(as.matrix(feature_data), 2, max)
  
  new_feature = sapply(1:ncol(feature_data), function(x)
    {
   ((feature_data[,x]-mins[x])/(maxs[x]-mins[x]))*(upper_bound-lower_bound)+lower_bound
    
  })
  
  saveRDS(new_feature, file = paste0(Rds_label, "_scale.Rds"))
  
  
 return(cbind(mins, maxs))
  
  
}


scale_test = function(feature_data, training_range, 
                      upper_bound, lower_bound,Rds_label)
{
  
  mins = training_range[,1]
  maxs = training_range[,2]
  
  new_feature = sapply(1:ncol(feature_data), function(x)
  {
    ((feature_data[,x]-mins[x])/(maxs[x]-mins[x]))*(upper_bound-lower_bound)+lower_bound
    
  })
  
  saveRDS((new_feature), file = paste0(Rds_label, "_scale.Rds"))
  
  
}





libsvm_formating = function(feature_matrix_Rds, feature_label_Rds, Rds_label, output_label)
{
  
  feature_matrix = readRDS(feature_matrix_Rds)
  feature_label = readRDS(feature_label_Rds)
  
  
  feature_col= ncol(feature_matrix)
  
  for(i in 1:feature_col)
  {
    feature_matrix[,i] = paste0(i, ":",feature_matrix[,i])
  }
  
  
  combined_data = cbind(feature_label, feature_matrix)
  
  saveRDS(combined_data, paste0(Rds_label, "_svm.Rds"))
  write.table(combined_data, paste0(output_label, "_svm.tsv"), quote = F,
              sep = "\t", row.names = F, col.names = F)
  
  
}

#r1+r2 = 10

divide_candidate_to_train_test = function(candidate_Rds, r1,r2,Rds_label, output_label)

{
  
  candidate_df = readRDS(candidate_Rds)
  
  pos_df <- candidate_df %>% dplyr::filter(label == "positive")
  candi_df <- candidate_df %>% dplyr::filter(label == "negative")
  
  pos_nrow = nrow(pos_df)
  candi_nrow = nrow(candi_df)
  
  train_pos_length = floor(pos_nrow/10)*r1
  
  train_candi_length = floor(candi_nrow/10)*r1
  
  
  pos_seq = c(1:pos_nrow)
  candi_seq = c(1:candi_nrow)
  
 set.seed(123)
 train_pos_sel =  sample(pos_nrow, train_pos_length)
      
 set.seed(123)
 train_candi_sel =  sample(candi_nrow, train_candi_length)
 
 
 test_pos_sel = pos_seq[-train_pos_sel]
 
 test_candi_sel = candi_seq[-train_candi_sel]
 
  
 train_df = rbind(pos_df[train_pos_sel,], candi_df[train_candi_sel,])
 
 test_df = rbind(pos_df[test_pos_sel,], candi_df[test_candi_sel,])
 
  
  saveRDS(train_df, paste0(Rds_label, "_train_part.Rds"))
  saveRDS(test_df, paste0(Rds_label, "_test_part.Rds"))
  
  write.table(train_df, paste0(output_label, "_train_part.tsv"), quote = F, row.names = F, sep = "\t")
  write.table(test_df, paste0(output_label, "_test_part.tsv"), quote = F, row.names = F, sep = "\t")
  
  
  
}






get_average_weights_of_ten = function(model_path,mod_name, n_fold, full_feature, arrange_feature)
  
{
  weight_matrix = matrix(0, nrow = length(full_feature), ncol = n_fold)
  
  
  for(i in 1:n_fold)
  {
    this_name = paste0(model_path, mod_name, "_nr_size25_nms_10_train_feature_",i,"_svm_model_c.tsv")
    
    weight_file = readLines(this_name)
    
    ws = as.numeric(weight_file[7:length(weight_file)])
    
    weight_matrix[,i] = ws
  }
  
  ave_weights = rowMeans(weight_matrix)
  
  ws_df = data.frame(name = full_feature,
                     weights = ave_weights,
                     abs_weights = abs(ave_weights), stringsAsFactors = F)
  
  name_arr_ws_df = arrange_feature%>%dplyr::left_join(ws_df, by = c("arrange_feature" ="name"))
  
  
  write.table(name_arr_ws_df,paste0(mod_name,"_arrange_name_coeff_df.tsv"),
              sep = "\t", row.names = F, quote = F)
  
  
  
  arr_ws_df <- ws_df %>% dpylr::arrange(desc(abs_weights))
  
  
  write.table(arr_ws_df,paste0(mod_name,"_coeff_df.tsv"),
              sep = "\t", row.names = F, quote = F)
  
  
  
  write.table(ws_df,paste0(mod_name,"_unsorted_coeff_df.tsv"),
              sep = "\t", row.names = F, quote = F)
  
}



#######

get_3_roc = function(prediction_path, label_path,
                         prediction_name1, prediction_name2, prediction_name3, 
                         label_name1, label_name2, label_name3, 
                         n_fold,output_fig_name)
{
  
# 
#   prediction_path = "/data/ginny/liblinear-2.11/cv_1103/"
#              label_path =    "cv_1103/"
#             prediction_name1 = "py_nr_size25_nms_10_test_feature"
#              label_name1 =    "py_nr_size25_nms_10_test_label"
#              prediction_name2 = "py_nr_size15_nms_10_test_feature"
#              label_name2 =    "py_nr_size15_nms_10_test_label"
#              prediction_name3 = "py_nr_size11_nms_10_test_feature"
#              label_name3 =    "py_nr_size11_nms_10_test_label"
#              
#                n_fold =  10
#                 output_fig_name = "py_10"



  
  total_score1 = c()
  total_label1 = c()
  
  total_score2 = c()
  total_label2 = c()
  
  total_score3 = c()
  total_label3 = c()
  
  for(i in 1:n_fold)
  {
    
    predict1_tsv = paste0(prediction_path, prediction_name1, "_",i,"_svm_predict.tsv")
    predict2_tsv = paste0(prediction_path, prediction_name2, "_",i,"_svm_predict.tsv")
    predict3_tsv = paste0(prediction_path, prediction_name3, "_",i,"_svm_predict.tsv")
    
    label1_tsv = paste0(label_path, label_name1, "_",i, ".tsv")
    label2_tsv = paste0(label_path, label_name2, "_",i, ".tsv")
    label3_tsv = paste0(label_path, label_name3, "_",i, ".tsv")
    
    
    
    output_name = paste0(output_fig_name, "_hist_roc_",i,".pdf")
    
    test_predict1 = data.table::fread(predict1_tsv,stringsAsFactors = F, header = T)
    test_label1 = data.table::fread(label1_tsv, stringsAsFactors = F)
    total_score1 = c(total_score1, test_predict1$'1')
    total_label1 = c(total_label1, test_label1$V1)
    
    
    test_predict2 = data.table::fread(predict2_tsv,stringsAsFactors = F, header = T)
    test_label2 = data.table::fread(label2_tsv, stringsAsFactors = F)
    total_score2 = c(total_score2, test_predict2$'1')
    total_label2 = c(total_label2, test_label2$V1)
    
    test_predict3 = data.table::fread(predict3_tsv,stringsAsFactors = F, header = T)
    test_label3 = data.table::fread(label3_tsv, stringsAsFactors = F)
    total_score3 = c(total_score3, test_predict3$'1')
    total_label3 = c(total_label3, test_label3$V1)
    
    
  }
 
  roc_total1 = roc(total_label1, total_score1)
  roc_total2 = roc(total_label2, total_score2)
  roc_total3 = roc(total_label3, total_score3)
  
  output_name_total = paste0(output_fig_name, "_3_roc.pdf")
  
  
 # pdf(output_name_total)#, useDingbats = F)
  plot(roc_total1, col = "red", main = output_fig_name)
  plot(roc_total2, add = T, col = "green")
  plot(roc_total3, add = T, col = "yellow")
  #dev.off()
  # 
  #cat("total auc: ", roc_total$auc)
  #cat("\n")
  
}




calculate_MCC = function(tp, fp, fn, tn)
{
  #fenzi nominator
  
  
   nominator = (tp/10000)*(tn/10000)-(fp/10000)*(fn/10000)

   #nominator =  100000000*nom
     
  real_p = tp+fp
  pred_p = tp+fn
  real_n = tn+fn
  pred_n = tn+fn
  
  c_pr = c(real_p, pred_p, real_n, pred_n)
  
  isZero <- c_pr == 0
  
  if(sum(isZero)==0)
  {
    p1 = 1/sqrt(tp+fn)
    p2 = 1/sqrt(tp+fp)
    p3 = 1/sqrt(tn+fp)
    p4 = 1/sqrt(tn+fn)
    
    
    MCC = nominator*p1*p2*p3*p4*100000000
    
  }else
  {
    MCC = 0
  }
  
  
  return(MCC)
  
  
}



map_domain = function(dm_df, pos_window_score)
{
  
  #pos_each_domain = rep(NA, nrow(pos_window_score))
  
  pos_window_score <- pos_window_score %>% arrange(protID, pos)
  
  ## it may be faster if I do by protein
  
  all_pos_protein = unique(pos_window_score$protID)
  
  all_prot_domains = rep("",0)
  
  for(i in 1:length(all_pos_protein))
  {
      # if(i %% 1000 ==0) cat(i,"\n")
  
    
    this_prot  = all_pos_protein[i]
    
    this_prot_positions = pos_window_score%>%
      dplyr::filter(protID == this_prot )
      this_pos = this_prot_positions$pos

    this_prot_domains = rep(NA, length(this_pos))
    
    if(this_prot %in% dm_df$protID)
    {
      this_dm_df = dm_df%>%dplyr::filter(protID == this_prot)
      this_start_positions = this_dm_df$start_position
      this_end_positions = this_dm_df$end_position
      
      tf_matrix = matrix(F, nrow = length(this_start_positions),
                         ncol = length(this_pos))
      
       for(x in 1:length(this_pos))
         {
            tf_matrix[,x] <- this_pos[x]>=this_start_positions & this_pos[x]<=this_end_positions
         }
  
  
      for( t in 1:ncol(tf_matrix))
        {
        if(sum(tf_matrix[,t])>0)
           {
             this_prot_domains[t] = paste(this_dm_df[tf_matrix[,t],]$domain_name, collapse  = " ")
           }
        }
      
    }
    
    all_prot_domains = c(all_prot_domains, this_prot_domains)

  }
  
  
  return(all_prot_domains)
  
}





map_domain_old = function(dm_df, pos_window_score)
{
  
  pos_each_domain = rep(NA, nrow(pos_window_score))
  
  for(i in 1:nrow(pos_window_score))
  {
    
    # 
    # if(i%%1000 == 0)
    #   cat(i, "\n")
    # 
    this_prot = pos_window_score[i,]$protID
    this_pos = pos_window_score[i,]$pos
    if(this_prot %in% dm_df$protID)
    {
      this_dm_df = dm_df%>%dplyr::filter(protID == this_prot)
      
      this_start = this_dm_df$start_position
      this_end = this_dm_df$end_position
      
      get_tf <- this_pos>this_start & this_pos<this_end
      
      if(sum(get_tf)>0)
      {
        pos_each_domain[i] = paste(this_dm_df[get_tf,]$domain_name, collapse  = " ")
        
      }
      
    }
    
    
  }
  
  
  return(pos_each_domain)
}



map_subcellular_location = function(sc_df, pos_window_score)
{
  sc_ps = left_join(pos_window_score, sc_df, by = "protID") %>%
    dplyr::select(location)
  
  return(sc_ps$location)
  
}




retrieve_for_each_mod = function(pos_info, candi_info, pos_score,candi_score,
                                 mod_name, domain_name)
{
  
  get_domain_pos = pos_info %>%
    dplyr::filter(prediction_score >= pos_score) %>%
    dplyr::filter(domain == domain_name) %>%
    dplyr::mutate(ptm = mod_name)
  
  get_domain_candi = candi_info %>%
    dplyr::filter(prediction_score >= candi_score) %>%
    dplyr::filter(domain == domain_name) %>%
    dplyr::mutate(ptm = mod_name)
  
  this_retrive = rbind(get_domain_pos, get_domain_candi)
  
  return(this_retrive)
  
}







