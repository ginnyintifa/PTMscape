
#Table of contents
- [Installation] (#ids1)
- [Functions] (#ids2)
- [Input files] (#ids3)
- [Input parameters] (#ids4)
- [Reference score threshold] (#ids5)
- [Example script] (#ids6)
- [Output description] (#ids7)



<div id = 'ids1'/>
# 1.Installation

`PTMscape` can be downloaded and installed in R with following code:


Installation of `qvalue` from `bioconductor` is required:

```{r, eval = F}

source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")

```

If `devtools` is not installed:

```{r, eval = F}
install.packages("devtools")
```

Install `PTMscape`:

```{r, eval = F}

library("devtools")
devtools::install_github("ginnyintifa/PTMscape")
library(PTMscape)

```


`PTMscape` requires installation of `Liblinear` library, it can be downloaded from [Liblinear website](https://www.csie.ntu.edu.tw/~cjlin/liblinear/). After installation of `Liblinear` on your machine, a small [change](https://www.csie.ntu.edu.tw/~cjlin/liblinear/FAQ.html#training_and_prediction) of source code in the file linear.cpp(see below) is required so as to activate probability outputs for SVM.



```

int check_probability_model(const struct model *model_)

{

	return (model_->param.solver_type==L2R_LR ||

to

int check_probability_model(const struct model *model_)

{

	return 1;
	
```

<div id = 'ids2'/>
# 2. Functions

### Predict PTM events
#### Whole proteome prediction

In this mode, it is assumed that user wants to discover all possible PTM sites in proteins of interest. Known PTM sites will be mapped to all the proteins. Then the whole data will be divided into k folds as specified by the user. Each time (k-1) folds of the data will be used as training data to fit a linear SVM model, subsequently the rest one fold of the data will be predicted by the model trained. This process is conducted for k times so that each fold of the data will be predicted once by the rest data. Finally known PTM sites plus predicted PTM sites (at a specificity level set by the user) will be the positive PTM sites. 

Function `predict_on_whole_proteome()` should be called. Input files and parameters will be described in the following sections.

#### Targeted prediction

In this mode, it is assumed that user provides a set of reliable PTM sites and a Fasta file contaning associated protein seqeunces(i.e. User is confident about the positive/negative designation of PTM sites in these protein sequences). The aim is to predict PTM events in uncharted proteins sequences with model trained from the reliable set. After prediction, the score threshold of positive/negative site will be selected by either refering to the threshold derived from large dataset(provided by `PTMscape` tool, see the chart below) or conducting cross validation within the reliable data and select the cutoff corresponding to user specified specificity level.  

Function `predict_on_targeted_proteome()` should be called. Input files and parameters will be described in the following sections.

### PTM crosstalk analysis

#### Positive crosstalk

Positive crosstalk in a protein domain is defined as two types of PTM happening on two residues which are near to each other (on the same protein). The distance between the two residues can be specified by the user. `PTMscape` takes the mapped files(these files are the output files of the predicting functions) of the two PTM types as input and perform fisher exact test to see if co-occurrence of PTMs are more frequent than expected. 

Function `calculate_tbt_positive_ptms()` should be called. Input files and parameters will be described in the following sections.

#### Negative crosstalk

Negative crosstalk in a protein domain is defined as two types of PTM happening on the same residue. The two PTM types may compete with each other for the chance of modifying the target site. `PTMscape` takes the mapped files(these files are the output files of the predicting functions) of the two PTM types as input and perform fisher exact test to see if competing  of two PTM events are more frequent than expected. 

Function `calculate_tbt_negative_ptms()` should be called. Input files and parameters will be described in the following sections.


<div id = 'ids3'/>
# 3.Input files
### User provided input files

* 1. Known PTM sites

User has to provide the Uniprot accession ID of proteins and position of the PTM sites in required format. See **sample_known_ps.tsv**.

* 2. Fasta file for proteins of interest

In whole proteome prediction mode, only one Fasta file is required. In targeted prediction mode, two Fasta files are needed. One consists proteins containing the reliable PTM sites information, the other consists protein sequences from which PTM sites are to be predicted. The format can be learned from **sample_known_fasta.tsv** and **sample_predict_fasta.tsv**.

### PTMscape provided input files

Several files need to be downloaded from this [webpage](http://137.132.97.109:59739/CSSB_LAB/) before running `PTMscape`.



* 1.  Clustered AAindex parameters in **aaindex_cluster_order.tsv**.  
* 2.  Extracted SPIDER3 position specific structural features for the whole eligible Swiss-prot human proteome in **extract_spider.Rds** (this file is applicable when flanking size is set to 12, we also provide **extracted_spider_7.Rds** and **extracted_spider_5.Rds** for flanking size 7 and 5 respectively).   
* 3.  File **spider_protID.tsv** is needed in feature extraction step, it contains all the Uniprot accession IDs of proteins for which SPIDER3 features are extractable.
* 4.  File **domain_df_pure.Rds** contains domain specifications compiled from Pfam tool.  
* 5.  Subcellular locations of each protein is in file **subcellusr_location_df_pure.Rds**, the information is retrieved from the database called COMPARTMENTS.


Please download these files and put them in the same working directory where you installed `PTMscape`.

<div id = 'ids4'/>
# 4.Input parameters

`PTMscape` requires several user specified parameters.

### Predict PTM events.

#### Whole proteome prediction

```ptm_site```  The amino acid this PTM involves, in upper-case single letter representation.  
```flanking_size``` The number of residues surrounding each side of the center residue, the total window size will be 2*flanking_size+1, default to 12.  
```SPIDER``` A boolean variable indicating the usage of SPIDER3 features, default set to TRUE.  
```positive_info_file```  A text file containing the positive PTM sites info in required format.  
```protein_fasta_file```  A text file containing the proteins sequences of interest in Fasta format.  
```liblinear_dir``` Absolute path of Liblinear tool.  
```n_fold```  Number of folds used for training and prediction, default set to 2.  
```feature_file_path``` Absolute path of the feature files.  
```lower_bound``` The lower bound of the scaled data range, default to -1.  
```upper_bound``` The upper bound of the scaled data range, default to 1.  
```cvlog_path_name``` The path and name of the log files, which hold the details of Liblinear procedures.  
```specificity_level``` A number ranges from 0 to 1 indicating the specificity user requires the classifier to achieve, default set to 0.99.  
```output_label```  The string to tag the output files.  


#### Targeted prediction



```ptm_site```  The amino acid this PTM involves, in upper-case single letter representation.  
```flanking_size``` The number of residues surround each side of the center residue, the total window size will be 2*flanking_size+1, default to 12.  
```SPIDER``` A boolean variable indicating the usage of SPIDER3 features, default set to TRUE.  
```positive_info_file```  A text file containing the positive PTM sites info in required format.  
```known_protein_fasta_file```  A text file containing the proteins sequences of interest and known PTM sites in Fasta format.  
```predict_protein_fasta_file```  A text file containing the proteins sequences with PTM sites to be predicted in Fasta format.  
```output_label_training``` The string to tag the output files associated with training proteins.  
```output_label_predict```  The string to tag the output files associated with prediction proteins.  
```liblinear_dir``` Absolute path of Liblinear tool.  
```n_fold```  Number of folds used for training and prediction, default set to 2.  
```feature_file_path``` Absolute path of the feature files.  
```lower_bound``` The lower bound of the scaled data range, default to -1.  
```upper_bound``` The upper bound of the scaled data range, default to 1.  
```cvlog_path_name``` The path and name of the log files, which hold the details of Liblinear procedures.  ```specificity_level``` A number between 0 and 1 indicating the specificity user requires the classifier to achieve, default set to 0.99. Used only not in "reference" mode.  
```flag_for_score_threshold_chosen``` A string indicating whether use reference score threshold or get from the user supplied training data, default set to "reference".  
```score_threshold``` A numerical value between 0 to 1 indicating the reference score threshold (required in "reference" mode).  



### PTM crosstalk analysis

#### Positive crosstalk

```distance``` A numerical value indicating distance between two PTM sites, defaul set to 5.  
```anchor_mod``` A string indicating the anchor modification.  
```cross_mod``` A string indicating the crosstalk modification.  
```anchor_mapped_df_Rds``` An Rds file containing the achor window score file with domain mapped.  
```cross_mapped_df_Rds``` An Rds file containing the cross window score file with domain mapped.  
```output_label``` The string to tag the output files.  

#### Negative crosstalk

```anchor_mod``` A string indicating the anchor modification.  
```cross_mod``` A string indicating the crosstalk modification.  
```anchor_mapped_df_Rds``` An Rds file containing the achor window score file with domain mapped.  
```cross_mapped_df_Rds``` An Rds file containing the cross window score file with domain mapped.  
```output_label``` The string to tag the output files.  

<div id = 'ids5'/>
# 5.Reference score threshold derived from PhosphoSitePlus PTM data



| PTM_type  | Score_at_spec0.9 | Score_at_spec0.95 | Score_at_spec0.99|
| ------------- | ------------- | ------------- | ------------- |
| phosphoS  | 0.577 |0.613|0.683|
| phosphoT  | 0.572 |0.607|0.675|
| phosphoY  | 0.575 |0.604|0.663|
| ubiquitinK  | 0.566 |0.586|0.621|
| acetylK | 0.563 |0.59|0.646|
| SUMOyK | 0.554|0.604|0.688|
| methylK| 0.584|0.621|0.688|
| methylR| 0.562|0.608|0.704|



<div id = 'ids6'/>
# 6.Example script

### Whole proteom prediction
Predict acetylation in all the proteins provided. A 2-fold cross validation will be conducted. Score threshold will be determined at specificity level 0.99.Two text output files will be produced. **ps_samle_wp_mapped_df.tsv** and **ps_sample_wp_test.tsv**.

```{r, eval=F}

predict_on_whole_proteome(ptm_site = "S",
                          flanking_size = 12,
                          positive_info_file = "sample_known_ps.tsv",
                          protein_fasta_file = "sample_known_fasta.tsv",
                          n_fold = 2,
                          lower_bound = -1,
                          upper_bound = 1,
                          liblinear_dir = "/data/ginny/Liblinear_prep/",
                          feature_file_path = "/data/ginny/PTMscape_test/",
                          cvlog_path_name = "/data/ginny/PTMscape_test/cvlog.txt",
                          specificity_level = 0.99,
                          output_label = "ps_sample_wp")

```




### Targeted prediction

Know phosphoS sites and proteins are used to train a linear SVM model. Select score threshold by cross validating within known PTM sites. A score corresponding to specificity 0.99 will be chosen as the threshold. Two text output files will be produced. **ps_sample_predict_mapped_df.tsv** and **ps_sample_predict_test.tsv**.

```{r, eval=F}
predict_on_targeted_proteome (ptm_site = "S", 
                              flanking_size=12, 
                              SPIDER = T,
                              positive_info_file = "sample_known_ps.tsv", 
                              known_protein_fasta_file = "sample_known_fasta.tsv",
                              predict_protein_fasta_file = "sample_predict_fasta.tsv",
                              output_label_training = "ps_sample_training",
                              output_label_predict = "ps_sample_predict",
                              lower_bound = -1,
                              upper_bound = 1,
                              liblinear_dir = "/data/ginny/liblinear-2.11/",
                              feature_file_path = "/data/ginny/test_package/",
                              cvlog_path_name = "/data/ginny/test_package/cvlog.txt",
                              specificity_level = 0.99,
                              n_fold = 2,
                              flag_for_score_threshold_chosen = "cv",
                              score_threshold = NULL)
```




Select score threshold by cross validating within known PTM sites. Score threshold is chosen by looking up to the reference table

```{r, eval=F}
predict_on_targeted_proteome (ptm_site = "S", 
                              flanking_size=12, 
                              SPIDER = T,
                              positive_info_file = "sample_known_ps.tsv", 
                              known_protein_fasta_file = "sample_known_fasta.tsv",
                              predict_protein_fasta_file = "sample_predict_fasta.tsv",
                              output_label_training = "ps_sample_training",
                              output_label_predict = "ps_sample_predict",
                              lower_bound = -1,
                              upper_bound = 1,
                              liblinear_dir = "/data/ginny/liblinear-2.11/",
                              feature_file_path = "/data/ginny/test_package/",
                              cvlog_path_name = "/data/ginny/test_package/cvlog.txt",
                              specificity_level = NULL,
                              n_fold = 2,
                              flag_for_score_threshold_chosen = "reference",
                              score_threshold = 0.683)
```



### Crosstalk analysis

Analyze positive crosstalk events between methylationK and phosphorylationS. The output file will be **sample_ps_ubi_positive_test.sv**. 


```{r, eval = F}

calculate_tbt_positive_ptms(distance = 5,
                        anchor_mod = "ps",
                        cross_mod = "ubi",
                        anchor_mapped_df_Rds = "sample_ps_mapped_df.Rds",
                        cross_mapped_df_Rds = "sample_ubi_mapped_df.Rds",
                        output_label = "sample_ps_ubi_positive")
```

Analyze negative crosstalk events between ubiquitination and methylationK. The output will be **sample_methy_k_ubi_negative.tsv**.

```{r, eval = F}
calculate_tbt_negative_ptms(anchor_mod = "methy_k",
                             cross_mod = "ubi",
                             anchor_mapped_df_Rds = "sample_methy_k_mapped_df.Rds",
                             cross_mapped_df_Rds = "sample_ubi_mapped_df.Rds",
                             output_label = "sample_methy_k_ubi_negtive")

```
<div id = 'ids7'/>
# 7.Output files description


Basically, user can ignore all the `.Rds` files produced by the functions as they are used internally by the functions in this package. After the function finishes the runing process, the following `.tsv` files should be examined.


**output_label_mapped_df.tsv** contains the predicted score for each modifiable site in proteins given.

**output_label_test.tsv** contains the significance score for each domain where the PTM takes place. 







