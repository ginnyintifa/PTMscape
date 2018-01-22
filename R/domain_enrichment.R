
#' Building up domain analysis.
#'
#' Map sites to protein domains and analyse the enrichment of PTM sites for each domain.
#' @param window_score_label_Rds An Rds file containing the window, positive and score of interested sites(genereated from present_predict_wp()).
#' @param output_label The string to tag the output files.
#' @import stringr dplyr magrittr qvalue data.table
#' @export
#' @examples
#' PTM_domain_mapping_enrichment(window_score_label_Rds = "ps_0103_window_score_df.Rds",
#'                               output_label = "ps_0103")



PTM_domain_mapping_enrichment = function(window_score_label_Rds,
output_label)
{
    ### domain mapping and significance test
    map_domain_subcellular(window_score_label_Rds = window_score_label_Rds,
    dm_df_Rds = "domain_df_pure.Rds",
    sc_df_Rds = "subcellular_location_df_pure.Rds",
    output_label = output_label)
    
    calculate_tbt_single_ptm(paste0(output_label, "_mapped_df.Rds"),
    output_label)
    
}



#' Analyse positive cross talk events in domains.
#'
#' Test for enrichment of two user provided PTM types inside of domains.
#' @param distance A numerical value indicating distance between two PTM sites (defaul set to 5).
#' @param anchor_mod A string indicating the anchor modification.
#' @param cross_mod A string indicating the other modification.
#' @param anchor_mapped_df_Rds An Rds file containing the window score file of anchor PTM mapped domain.  
#' @param cross_mapped_df_Rds An Rds file containing the window score file of the other PTM with mapped domain.
#' @param output_label The string to tag the output files.
#' @import stringr dplyr magrittr qvalue data.table
#' @export
#' @examples
#' calculate_tbt_positive_ptms(distance = 5,
#'                        anchor_mod = "ps",
#'                        cross_mod = "pt",
#'                        anchor_mapped_df_Rds = "ps_0103_mapped_df.Rds",
#'                        cross_mapped_df_Rds = "pt_0103_mapped_df.Rds",
#'                        output_label = "pt_ps_positive")

calculate_tbt_positive_ptms = function(distance = 5,
                                  anchor_mod, 
                                  cross_mod,
                                  anchor_mapped_df_Rds,
                                  cross_mapped_df_Rds,
                                  output_label)
{
  # distance = 5
  # anchor_mod = "acety"
  # cross_mod = "ubi"
  # anchor_mapped_df_Rds = "acety_52/acety_wp_52_mapped_df.Rds"
  # cross_mapped_df_Rds = "ubi_52/ubi_wp_52_mapped_df.Rds"
  # output_label = "new_acety_ubi_positive_52"
  # 
  # 
  anchor_mapped_df = readRDS(anchor_mapped_df_Rds)
  cross_mapped_df = readRDS(cross_mapped_df_Rds)
  # 
  ### change column names
  
  cn_anchor = colnames(anchor_mapped_df)
  cn_anchor[which(grepl("pred_label", cn_anchor))] = "anchor_label"
  colnames(anchor_mapped_df) = cn_anchor
  
  cn_cross = colnames(cross_mapped_df)
  cn_cross[which(grepl("pred_label", cn_cross))] = "cross_label"
  colnames(cross_mapped_df) = cn_cross
  
  
  
  
  
  anchor_mapped_df = anchor_mapped_df %>%
    dplyr::mutate(anchor_label = replace(anchor_label, anchor_label == "positive", TRUE)) %>%
    dplyr::mutate(anchor_label = as.logical(replace(anchor_label, anchor_label == "negative", FALSE)))
  
  
  cross_mapped_df = cross_mapped_df %>%
    dplyr::mutate(cross_label = replace(cross_label, cross_label == "positive", TRUE)) %>%
    dplyr::mutate(cross_label = as.logical(replace(cross_label, cross_label == "negative", FALSE)))
  
  
  anchor_domain_list =unique(anchor_mapped_df$domain)
  anchor_domain_list = unique( unlist(strsplit(anchor_domain_list, split = " ")))
  
  
  cross_domain_list = unique(cross_mapped_df$domain)
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
    
    
    this_domain = common_domain[i]
    this_anchor_domain = anchor_mapped_df %>%
      dplyr:: filter(grepl(paste0("\\b",this_domain,"\\b") ,domain))
    this_cross_domain = cross_mapped_df %>%
      dplyr:: filter(grepl(paste0("\\b",this_domain,"\\b"), domain))
    
    
    if(nrow(this_anchor_domain)>0 & nrow(this_cross_domain)>0)
    {
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
  
  
  
  
  get_tbt = output_df 
  colnames(get_tbt) = c("domain", "both_positive",paste0(cross_mod, "_positive"),
                        paste0(anchor_mod,"_positive"), "both_negative")
  
  
  
  ## then I need to combine accroding to domains groupby domains and summ a b c d
  
  ### calculate fisher's exact test on this
  fisher_p = rep(0, nrow(get_tbt))
  or = rep(0, nrow(get_tbt))
  adjusted_or = rep(0,nrow(get_tbt))
  
  for(i in 1:nrow(get_tbt))
  {
    
    fm = matrix(as.numeric(get_tbt[i,c(2:5)]), nrow = 2, ncol = 2, byrow = T)
    
    ft =  fisher.test(fm, alternative = "g")
    
    fisher_p[i] = ft$p.value
    or[i] = ft$estimate
    adjusted_or[i] = (get_tbt[i,2]+0.5)*(get_tbt[i,5]+0.5)/(get_tbt[i,3]+0.5)/(get_tbt[i,4]+0.5)
    
  }
  
  qv = qvalue(p = fisher_p)
  
  tbt_p = get_tbt %>%
    dplyr::mutate(fisher_pvalue = fisher_p)%>%
    dplyr::mutate(odds_ratio = or)%>%
    dplyr::mutate(qvalue = qv$qvalues)%>%
    dplyr::mutate(adjusted_odds_ratio = adjusted_or)%>%
    dplyr::arrange(fisher_pvalue)
  
  
  
  write.table(tbt_p, paste0(output_label, "_test.tsv"),
              quote = F, row.names = F, sep = "\t")
  
  
  
}




#' Analyse negative cross talk events in domains.
#'
#' Test for enrichment of two user provided PTM types inside of domains.
#' @param anchor_mod A string indicating the anchor modification.
#' @param cross_mod A string indicating the other modification.
#' @param anchor_mapped_df_Rds An Rds file containing the window score file of anchor PTM mapped domain.  
#' @param cross_mapped_df_Rds An Rds file containing the window score file of the other PTM with mapped domain.
#' @param output_label The string to tag the output files.
#' @import stringr dplyr magrittr qvalue data.table
#' @export
#' @examples
#' calculate_tbt_negative_ptms(anchor_mod = "ubi",
#'                             cross_mod = "acety",
#'                             anchor_mapped_df_Rds = "ubi_0103_mapped_df.Rds",
#'                             cross_mapped_df_Rds = "acety_0103_mapped_df.Rds",
#'                             output_label = "ubi_acety_negtive")




calculate_tbt_negative_ptms = function(anchor_mod,
                                       cross_mod,
                                       anchor_mapped_df_Rds,
                                       cross_mapped_df_Rds,
                                       output_label)
{
  
  
  # anchor_mod = "ubi"
  # cross_mod = "methy_k"
  # anchor_mapped_df_Rds = "ubi_wp_52_mapped_df.Rds"
  # cross_mapped_df_Rds = "methy_k_80_mapped_df.Rds"
  # output_label = "ubi_methy_k_negtive"
  # 
  anchor_mapped_df = readRDS(anchor_mapped_df_Rds)
  cross_mapped_df = readRDS(cross_mapped_df_Rds)

  
  cn_anchor = colnames(anchor_mapped_df)
  cn_anchor[which(grepl("pred_label", cn_anchor))] = "anchor_label"
  colnames(anchor_mapped_df) = cn_anchor
  
  cn_cross = colnames(cross_mapped_df)
  cn_cross[which(grepl("pred_label", cn_cross))] = "cross_label"
  colnames(cross_mapped_df) = cn_cross
  
  
  
  anchor_mapped_df = anchor_mapped_df %>%
    dplyr::mutate(anchor_label = replace(anchor_label, anchor_label == "positive", TRUE)) %>%
    dplyr::mutate(anchor_label = as.logical(replace(anchor_label, anchor_label == "negative", FALSE)))
  
  
  cross_mapped_df = cross_mapped_df %>%
    dplyr::mutate(cross_label = replace(cross_label, cross_label == "positive", TRUE)) %>%
    dplyr::mutate(cross_label = as.logical(replace(cross_label, cross_label == "negative", FALSE)))
  
  
  anchor_domain_list =unique(anchor_mapped_df$domain)
  anchor_domain_list = unique( unlist(strsplit(anchor_domain_list, split = " ")))
  
  
  cross_domain_list = unique(cross_mapped_df$domain)
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
    

    a =0;b=0;c=0;d=0
    

    this_domain = common_domain[i]
    this_domain_anchor_df = anchor_mapped_df %>%
      dplyr:: filter(grepl(paste0("\\b",this_domain,"\\b") ,domain))
    this_domain_cross_df = cross_mapped_df %>%
      dplyr:: filter(grepl(paste0("\\b",this_domain,"\\b") ,domain))
    
    ### each domain may be involved in different proteins
    
    if(nrow(this_domain_anchor_df)>0 )
    {
      
      anchor_match_label = this_domain_anchor_df$anchor_label
      cross_match_label = this_domain_cross_df$cross_label
      
      a = length(which(anchor_match_label == T & cross_match_label == T))
      b = length(which(anchor_match_label == F & cross_match_label == T))
      c = length(which(anchor_match_label == T & cross_match_label == F))
      d = length(which(anchor_match_label == F & cross_match_label == F))
      
      # cat(a,b,c,d,"\n")
      
      
      
      have_anchor = this_domain_anchor_df %>%
        dplyr::filter(anchor_label == T) %>%
        dplyr::select(protID)
      
    
      have_cross = this_domain_cross_df %>%
        dplyr::filter(cross_label == T) %>%
        dplyr::select(protID)
      
      anchor_cross_prot = intersect(have_anchor$protID, have_cross$protID)
      
      
      if(nrow(have_cross)>0)
        domain_prot_df$protIDs[i] = paste(unique(anchor_cross_prot),
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
  
  
  
  get_tbt = output_df 
  colnames(get_tbt) = c("domain", "both_positive",paste0(cross_mod, "_positive"),
                        paste0(anchor_mod,"_positive"), "both_negative")
  
  
  
  ## then I need to combine accroding to domains groupby domains and summ a b c d
  
  ### calculate fisher's exact test on this
  fisher_p = rep(0, nrow(get_tbt))
  or = rep(0, nrow(get_tbt))
  adjusted_or = rep(0,nrow(get_tbt))
  
  
  for(i in 1:nrow(get_tbt))
  {
    
    fm = matrix(as.numeric(get_tbt[i,c(2:5)]), nrow = 2, ncol = 2, byrow = T)
    
    ft =  fisher.test(fm, alternative = "g")
    
    fisher_p[i] = ft$p.value
    or[i] = ft$estimate
    adjusted_or[i] = (get_tbt[i,2]+0.5)*(get_tbt[i,5]+0.5)/(get_tbt[i,3]+0.5)/(get_tbt[i,4]+0.5)
    
  }
  
  qv = qvalue(p = fisher_p)
  
  tbt_p = get_tbt %>%
    dplyr::mutate(fisher_pvalue = fisher_p)%>%
    dplyr::mutate(odds_ratio = or)%>%
    dplyr::mutate(qvalue = qv$qvalues)%>%
    dplyr::mutate(adjusted_odds_ratio = adjusted_or)%>%
    dplyr::arrange(fisher_pvalue)
  

  
  write.table(tbt_p, paste0(output_label, "_test.tsv"),
              quote = F, row.names = F, sep = "\t")
  
  
}




