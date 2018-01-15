# PTMscape user guide

# Download and install

```
install.packages("devtools")
library("devtools")
devtools::install_github("ginnyintifa/PTMscape")
library(PTMscape)

```




# Package and library dependency

PTMscape requires installation of the following R packages
```{r, eval = F}
install.packages("data.table")
install.packages("dplyr")
install.packages("magrittr")
install.packages("caret")
install.packages("stringr")
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")

```

PTMscape requires installation of Liblinear library, it can be downloaded from [Liblinear](https://www.csie.ntu.edu.tw/~cjlin/liblinear/). After installation of Liblinear on your machine, a small [change](https://www.csie.ntu.edu.tw/~cjlin/liblinear/FAQ.html#training_and_prediction) of source code in the file linear.cpp(see below) is required so as to activate probability outputs for SVM.



```

int check_probability_model(const struct model *model_)

{

	return (model_->param.solver_type==L2R_LR ||

to

int check_probability_model(const struct model *model_)

{

	return 1;
	
```


# Prediction mode choices.

## Whole proteome prediction mode.

In this mode, it is assumed that user wants to discover all possible PTM sites in proteins of interest. Known PTM sites will be mapped to all the proteins. Then the whole data will be divided into k folds as specified by the user. Each time (k-1) folds of the data will be used as training data to fit a linear SVM model, subsequently the rest one fold of the data will be predicted by the model trained. This process is conducted for k times so that each fold of the data will be predicted once by the rest data. Finally known PTM sites plus predicted PTM sites (at a specificity level set by the user) will be the positive PTM sites.

## Targeted prediction mode.

In this mode, it is assumed that user provides a set of reliable PTM sites and a Fasta file contaning associated protein seqeunces(i.e. User is confident about the positive/negative designation of PTM sites in these protein sequences). The aim is to predict PTM events in uncharted proteins sequences with model trained from the reliable set. After prediction, the score threshold of positive/negative site will be selected by either refering to the threshold derived from large dataset(provided by PTMscape tool) or conducting cross validation within the known data and select the cutoff corresponding to user specified specificity level.  






# Source file preparation 
## User provided input files 

### 1. Known PTM sites

In this file user has to provide the Uniprot accession number of proteins and position of the PTM sites in required format. See **sample_PTM_info.tsv**.

### 2. Fasta file for proteins of interest

In whole proteome prediction mode, only one Fasta file is required. In targeted prediction mode, two Fasta files are needed. One consists proteins containing the reliable PTM sites information, the other consists protein sequences from which PTM sites are to be predicted. The format can be learned from **sample_fasta.tsv**.

## PTMscape provided input files

Several files need to be downloaded from this [Google drive link]() before running PTMscape.

### 1. Feature extraction 

Clustered AAindex parameters in **aaindex_cluster_order.tsv**.  
Extracted SPIDER3 position specific structural features for the whole eligible Swiss-prot human proteome in **extract_spider.Rds** (this file is applicable when flanking size is set to 12, we also provide **extracted_spider_7.Rds** and **extracted_spider_5.Rds** for flanking size 7 and 5 respectively).   
File **spider_protID.tsv** is needed in feature extraction step, it contains all the Uniprot accession IDs of proteins for which SPIDER3 features are extractable.

### 2. Domain analysis


File **domain_df_pure.Rds** contains domain specifications compiled from Pfam tool.  
Subcellular locations of each protein is in file **subcellusr_location_df_pure.Rds**, the information is retrieved from the database called COMPARTMENTS.




# Input parameter specification 

PTMscape requires several user specified parameters.


## Whole proteome prediction mode

```ptm_site```  The amino acid this PTM involves, in upper-case single letter representation.  
```flanking_size``` The number of residues surround each side of the center residue, the total window size will be 2*flanking_size+1, default to 12.  
```positive_info_file```  A text file containing the positive PTM sites info in required format.  
```protein_fasta_file```  A fext file containing the proteins sequences of interest in Fasta format.  
```liblinear_dir``` Absolute peth of Liblinear tool.  
```n_fold```  Number of folds used for training and prediction, default set to 2.  
```feature_file_path``` Absolute path of the feature files.  
```lower_bound``` The lower bound of the scaled data range, default to -1.  
```upper_bound``` The upper bound of the scaled data range, default to 1.  
```cvlog_path_name``` The path and name of the log files, which hold the details of Liblinear procedures.  
```specificity_level``` A number indicating the specificity user requires the classifier to achieve, default set to 0.99.  
```output_label```  The string to tag the output files.  


## Targeted prediction mode



```ptm_site```  The amino acid this PTM involves, in upper-case single letter representation.  
```flanking_size``` The number of residues surround each side of the center residue, the total window size will be 2*flanking_size+1, default to 12.  
```positive_info_file```  A text file containing the positive PTM sites info in required format.  
```known_protein_fasta_file```  A fext file containing the proteins sequences of interest and known PTM sites in Fasta format.  
```predict_protein_fasta_file```  A fext file containing the proteins sequences with PTM sites to be predicted in Fasta format.  
```output_label_training``` The string to tag the output files associated with training proteins.  
```output_label_predict```  The string to tag the output files associated with prediction proteins.  
```liblinear_dir``` Absolute peth of Liblinear tool.  
```n_fold```  Number of folds used for training and prediction, default set to 2.  
```feature_file_path``` Absolute path of the feature files.  
```lower_bound``` The lower bound of the scaled data range, default to -1.  
```upper_bound``` The upper bound of the scaled data range, default to 1.  
```cvlog_path_name``` The path and name of the log files, which hold the details of Liblinear procedures.  ```specificity_level``` A number between 0 and 1 indicating the specificity user requires the classifier to achieve, default set to 0.99. Used only not in "reference" mode.  
```flag_for_score_threshold_chosen``` A string indicating whether use reference score threshold or get from the user supplied training data, defalt set to "reference".  
```score_threshold``` A numerical value between 0 to 1 indicating the reference score threshold (supply when in "reference").  



## Crosstalk domain enrichment analysis

```distance``` A numerical value indicating distance between two PTM sites, defaul set to 5.  
```anchor_mod``` A string indicating the anchor modification.  
```cross_mod``` A string indicating the crosstalk modification.  
```anchor_mapped_df_Rds``` An Rds file containing the achor window score file with domain mapped.  
```cross_mapped_df_Rds``` An Rds file containing the cross window score file with domain mapped.  
```output_label``` The string to tag the output files.  

# Example script


```{r, eval=F}

predict_on_whole_proteome(ptm_site = "K",
                          flanking_size = 12,
                          positive_info_file = "acety_PSP.tsv",
                          protein_fasta_file = "K_sp_fasta.tsv",
                          n_fold = 2,
                          lower_bound = -1,
                          upper_bound = 1,
                          liblinear_dir = "/data/ginny/Liblinear_prep/",
                          feature_file_path = "/data/ginny/PTMscape_test/",
                          cvlog_path_name = "/data/ginny/PTMscape_test/cvlog.txt",
                          specificity_level = 0.99,
                          output_label = "acety_wp")


predict_on_targeted_proteome (ptm_site = "S", 
                              flanking_size=12, 
                              positive_info_file = "sub_known_ps.tsv", 
                              known_protein_fasta_file = "sub_known_fasta.tsv",
                              predict_protein_fasta_file = "sub_predict_fasta.tsv",
                              output_label_training = "new_target_ps_training",
                              output_label_predict = "new_target_ps_predict",
                              lower_bound = -1,
                              upper_bound = 1,
                              liblinear_dir = "/data/ginny/liblinear-2.11/",
                              feature_file_path = "/data/ginny/test_package/",
                              cvlog_path_name = "/data/ginny/test_package/cvlog.txt",
                              specificity_level = 0.99,
                              flag_for_score_threshold_chosen = "cv",
                              score_threshold = NULL)


calculate_tbt_two_ptms(distance = 5,
                       anchor_mod = "ps",
                       cross_mod = "pt",
                       anchor_mapped_df_Rds = "ps_0103_mapped_df.Rds",
                       cross_mapped_df_Rds = "pt_0109_mapped_df.Rds",
                       output_label = "pt_ps_positive")
  
  



```







