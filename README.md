
# Table of contents
- [Installation](https://github.com/ginnyintifa/PTMscape#1installation)
- [Functions](https://github.com/ginnyintifa/PTMscape#2-functions)
- [Input and feature files](https://github.com/ginnyintifa/PTMscape#3input-and-feature-files)
- [Input parameters](https://github.com/ginnyintifa/PTMscape#4input-parameters)
- [Reference score threshold](https://github.com/ginnyintifa/PTMscape#5reference-score-threshold-derived-from-phosphositeplus-ptm-data)
- [Example script](https://github.com/ginnyintifa/PTMscape#6example-script)
- [Description of output files](https://github.com/ginnyintifa/PTMscape#7output-files-description)


## 1.Installation

`PTMscape` can be downloaded and installed in R as follows. As a prerequisite, `devtools` must be installed:

```{r, eval = F}
install.packages("devtools")
```

Next, install `PTMscape`:

```{r, eval = F}

library("devtools")
devtools::install_github("ginnyintifa/PTMscape")
library(PTMscape)

```


Finally, `PTMscape` requires installation of `Liblinear` library, which can be downloaded from [Liblinear website](https://www.csie.ntu.edu.tw/~cjlin/liblinear/). After downloading and unzipping, a small [change](https://www.csie.ntu.edu.tw/~cjlin/liblinear/FAQ.html#training_and_prediction) of source code in the file `linear.cpp` (see below) is required so as to produce probability score for SVM prediction.


```

int check_probability_model(const struct model *model_)
{
        return (model_->param.solver_type==L2R_LR ||
                        model_->param.solver_type==L2R_LR_DUAL ||
                        model_->param.solver_type==L1R_LR);
}


to

int check_probability_model(const struct model *model_)

{

	return 1;
}
```

After the change is made, compile the files following the instructions below:

>On Unix systems, type `make` to build the `train` and `predict`.
programs. Run them without arguments to show the usages.

>On other systems, consult `Makefile` to build them (e.g., see
'Building Windows binaries' in this file) or use the pre-built
binaries (Windows binaries are in the directory `windows').

>This software uses some level-1 BLAS subroutines. The needed functions are
included in this package.  If a BLAS library is available on your
machine, you may use it by modifying the Makefile: Unmark the following line

        #LIBS ?= -lblas

>and mark

        LIBS ?= blas/blas.a



Note: we provide a modified version of `linear.cpp` in the repository.

## 2. Functions

### Prediction of  PTM events
#### Whole proteome scale prediction

In this mode of analysis, it is assumed that the user wants to discover all possible PTM sites in proteins of interest. Position information of known PTM sites will be mapped to all proteins. Then the whole data will be divided into k folds as specified by the user. Each time (k-1) folds of the data will be used as training data to fit a linear SVM model, and subsequently the remaining one fold of the data will be predicted by the model trained. This process is repeated for k times so that each fold of the data will be predicted once by the rest of the data. Finally known PTM sites and predicted PTM sites (at a specificity level set by the user) will be the regarded as positive PTM sites. 

For this analysis, the function `predict_on_whole_proteome()` should be called. Input files and parameters are described below.

#### Targeted prediction

In this mode of analysis, it is assumed that the user provides a set of reliable PTM sites and a fasta file containing associated protein sequences i.e. user is confident about the positive/negative designation of PTM sites in these protein sequences). The aim here is to predict PTM events in previously uncharted proteins sequences with model trained from the reliable set. After prediction, the score threshold of positive/negative site will be selected by either referring to the threshold derived from large dataset (provided by `PTMscape` tool, see the chart below) or by conducting cross validation within the reliable data and select the cutoff corresponding to user specified specificity level.  

For this mode of analysis, the function `predict_on_targeted_proteome()` should be called. Input files and parameters are described in the following sections.

### PTM crosstalk analysis

#### Positive crosstalk

Positive crosstalk in a protein domain is defined as two types of PTM occurring on two residues close to one another (e.g. 5AA apart). The distance between the two residues can be specified by the user. `PTMscape` takes the mapped files (these files are the output files of the prediction functions) of the two PTM types as input and performs Fisher exact test to see whether the frequency of co-occurrence of PTMs is higher than expected. 

The function `calculate_tbt_positive_ptms()` should be called for this analysis. Input files and parameters are described in the following sections.

#### Negative crosstalk

Negative crosstalk in a protein domain is defined as two types of PTM occurring on the same residue. The two PTM types may compete with each other for the chance of modifying the target site. `PTMscape` takes the mapped files (these files are the output files of the prediction functions) of the two PTM types as input and performs Fisher exact test to see whether the frequency of co-occurrence of PTMs is higher than expected. 

The function `calculate_tbt_negative_ptms()` should be called for this analysis. Input files and parameters are described in the following sections.


## 3.Input and feature files
### User provided input files

* 1. Known PTM sites

User has to provide the Uniprot accession ID of proteins and position of the PTM sites in the required format. See **sample_known_ps.tsv**.

* 2. Fasta file for proteins of interest

In the whole proteome scale prediction mode, only one fasta file is required. In targeted prediction mode, two fasta files are needed. One consists proteins containing the reliable PTM sites information, and the other consists protein sequences from which PTM sites are to be predicted. See the sample data: **sample_known_fasta.tsv** and **sample_predict_fasta.tsv**.

### PTMscape provided feature files

Several files need to be downloaded from the following [website](http://137.132.97.109:59739/CSSB_LAB/) before running `PTMscape`:


* 1.  Clustered AAindex parameters in **aaindex_cluster_order.tsv**.  
* 2.  Extracted SPIDER3 position specific structural features for the whole eligible Swiss-prot human proteome in **extract_spider.Rds** (this file is applicable when flanking size is set to 12 (window size 25). We also provide **extracted_spider_7.Rds** and **extracted_spider_5.Rds** for window size 15 and 11 respectively).   
* 3.  **spider_protID.tsv** is needed in the feature extraction step. It contains all the Uniprot accession IDs of proteins, for which SPIDER3 features can be extracted.
* 4. **uniprotID_genename.tsv** contains mapping between Uniprot accession numbers and gene names.
* 4.  **domain_df_pure.Rds** contains domain specifications compiled from the `Pfam` tool.  
* 5.  Subcellular locations of each protein is in the file named **subcellusr_location_df_pure.Rds**, the information is retrieved from the database called `COMPARTMENTS`.


Download these files and put them in the same working directory where you installed `PTMscape`.

## 4.Input parameters

`PTMscape` requires several user specified parameters.

### Prediction of  PTM events.

#### Whole proteome prediction

```ptm_site```  The target amino acid of the given PTM type, in upper-case single letter representation.  
```flanking_size``` The number of residues surrounding each side of the center residue, The total window size will be 2*flanking_size+1 (default to 12).  
```SPIDER``` A boolean variable indicating whether to use SPIDER3 features (default set to TRUE).  
```positive_info_file```  A text file containing the positive PTM sites in the required format.  
```protein_fasta_file```  A text file containing the protein sequences of interest in fasta format.  
```liblinear_dir``` The path for the Liblinear tool.  
```n_fold``` The number of folds used for training and prediction in cross validation stage (default set to 2).  
```feature_file_path``` The path for the feature files.  
```lower_bound``` The lower bound of the scaled data range (default to -1).  
```upper_bound``` The upper bound of the scaled data range (default to 1).  
```cvlog_path_name``` The path and name of the log files, which hold the details of Liblinear procedures.  
```specificity_level``` A number ranges from 0 to 1 indicating the specificity user requires the classifier to achieve (default to 0.99).  
```output_label```  The string to tag the output files.  


#### Targeted prediction


```ptm_site```  The target amino acid of the given PTM type, in upper-case single letter representation.  
```flanking_size``` The number of residues surrounding each side of the center residue, The total window size will be 2*flanking_size+1 (default to 12).  
```SPIDER``` A boolean variable indicating whether to use SPIDER3 features (default set to TRUE).  
```positive_info_file```  A text file containing the positive PTM sites in the required format.  
```protein_fasta_file```  A text file containing the protein sequences of interest in fasta format. 
```known_protein_fasta_file```  A text file containing the proteins sequences of interest and known PTM sites in fasta format.  
```predict_protein_fasta_file```  A text file containing the proteins sequences with PTM sites to be predicted in fasta format.  
```output_label_training``` The string to tag the output files associated with training proteins.  
```output_label_predict```  The string to tag the output files associated with prediction proteins.   
```liblinear_dir``` The path for the Liblinear tool.  
```n_fold``` The number of folds used for training and prediction in cross validation stage (default set to 2).  
```feature_file_path``` The path for the feature files.  
```lower_bound``` The lower bound of the scaled data range (default to -1).  
```upper_bound``` The upper bound of the scaled data range (default to 1).  
```cvlog_path_name``` The path and name of the log files, which hold the details of Liblinear procedures.  
```specificity_level``` A number ranges from 0 to 1 indicating the specificity user requires the classifier to achieve (default to 0.99).  
```flag_for_score_threshold_chosen``` A string indicating whether use reference score threshold or get from the user supplied training data (default set to "reference").  
```score_threshold``` A numerical value between 0 to 1 indicating the reference score threshold (required in "reference" mode).  




### PTM crosstalk analysis

#### Positive crosstalk

```distance``` A numerical value indicating distance between two PTM sites (default to 5).  
```anchor_mod``` A string indicating the anchor modification.  
```cross_mod``` A string indicating the other modification.  
```anchor_mapped_df_Rds``` An Rds file containing the window score file of anchor PTM mapped domain.  
```cross_mapped_df_Rds``` An Rds file containing the window score file of the other PTM with mapped domain.  
```output_label``` The string to tag the output files.  

#### Negative crosstalk

```anchor_mod``` A string indicating the anchor modification.  
```cross_mod``` A string indicating the other modification in competition with the anchor modification.  
```anchor_mapped_df_Rds``` An Rds file containing the window score file of anchor PTM mapped domain.  
```cross_mapped_df_Rds``` An Rds file containing the window score file of the other PTM with mapped domain.  
```output_label``` The string to tag the output files.  

## 5.Reference score threshold derived from PhosphoSitePlus PTM data



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



## 6.Example script

### Whole proteome prediction
Let us assume that we predict lysine acetylation on all proteins provided. A two-fold cross validation will be conducted. Score threshold will be determined at specificity level 0.99. Two text output files will be produced: **ps_samle_wp_mapped_df.tsv** and **ps_sample_wp_test.tsv**.

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


##### Running time needed for serine phosphorylation on whole proteome scale

The program was run on a Linux server with 32 Intel Xeon(R) E4-2630 v3 (2.40GHz) processors and 64Gb memory.




| Total time  | Feature generation | Training/prediction preparation | Liblinear processing| Prediction score annotation | Domain enrichment analysis |
| ------------- | ------------- | ------------- | ------------- |------------- | ------------- |
| 39min  | 9min | 17min |4min| 3min | 6min |





### Targeted prediction

Known phosphoS sites and proteins are used to train a linear SVM model. Select a score threshold by cross validating within known PTM sites. The score corresponding to specificity 0.99 will be chosen as the threshold. Two text output files will be produced: **ps_sample_predict_mapped_df.tsv** and **ps_sample_predict_test.tsv**.

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
                              liblinear_dir = "/data/ginny/Liblinear_prep/",
                              feature_file_path = "/data/ginny/PTMscape_test/",
                              cvlog_path_name = "/data/ginny/PTMscape_test/cvlog.txt",
                              specificity_level = 0.99,
                              n_fold = 2,
                              flag_for_score_threshold_chosen = "cv",
                              score_threshold = NULL)
```




Select a score threshold by looking up to the reference table.

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
                              liblinear_dir = "/data/ginny/Liblinear_prep/",
                              feature_file_path = "/data/ginny/PTMscape_test/",
                              cvlog_path_name = "/data/ginny/PTMscape_test/cvlog.txt",
                              specificity_level = NULL,
                              n_fold = 2,
                              flag_for_score_threshold_chosen = "reference",
                              score_threshold = 0.683)
```



### Crosstalk analysis

Let us assume that we analyze positive crosstalk events between methylationK and phosphorylationS. The output file will be **sample_ps_ubi_positive_test.sv**. 


```{r, eval = F}

calculate_tbt_positive_ptms(distance = 5,
                        anchor_mod = "ps",
                        cross_mod = "ubi",
                        anchor_mapped_df_Rds = "sample_ps_mapped_df.Rds",
                        cross_mapped_df_Rds = "sample_ubi_mapped_df.Rds",
                        output_label = "sample_ps_ubi_positive")
```

Meanwhile, this time we analyze negative crosstalk events between ubiquitination and methylationK. The output will be **sample_methy_k_ubi_negative.tsv**.

```{r, eval = F}
calculate_tbt_negative_ptms(anchor_mod = "methy_k",
                             cross_mod = "ubi",
                             anchor_mapped_df_Rds = "sample_methy_k_mapped_df.Rds",
                             cross_mapped_df_Rds = "sample_ubi_mapped_df.Rds",
                             output_label = "sample_methy_k_ubi_negtive")

```

## 7.Description of output files.


The user can ignore all the `.Rds` files produced by the functions as they are used internally by the functions in this package. After the function finishes the running process, the following `.tsv` files should be examined:


**output_label_mapped_df.tsv** contains the predicted score for each modifiable site in proteins given.

**output_label_test.tsv** contains the significance score of enrichment for each domain where the PTM takes place. 







