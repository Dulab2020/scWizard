# scWizard
## What is scWizard?
scWizard is a tool with graphical user interface (GUI) for cancer single cell RNA-seq data analysis.
The aim of scWizard is:
1)	providing comprehensive analysis for integration strategies of cancer scRNA-seq data.
2)	providing comprehensive classification and annotation for 47 cell subtypes within tumor microenvironment and specific for cancer research.
3)	providing accurate and efficient classification and annotation for scRNA-seq cell subtypes based on hierarchical model by deep neural network.
4)	integrating automated methods into easy-to-use workflows and point-and-click tool helping for researchers without programming skills.
## Starting to use scWizard
To start using scWizard, you can install a R package, it can provides a interactive UI.
### Installing scWizard
NOTE: Works with R v4.0.2 or greater.<br>
To install this package from Github, please use the code below.
Windows:
First install jags. Download address: https://nchc.dl.sourceforge.net/project/mcmc-jags/JAGS/4.x/Source/JAGS-4.3.0.tar.gz. And then run the following code in R
Linux:
First run conda install -c conda-forge jags. And then run the following code in R
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("Dulab2020/scWizard")
```
### Running scWizard
Windows:
The following commands are used to start the graphical user interface (GUI).
```R
scWizard::run_app()
```
Linux:
The following commands are used to start the graphical user interface (GUI).
```R
options(shiny.host='0.0.0.0')
options(shiny.launch.browser=F)
scWizard::run_app()
```
Then enter the “server IP address: port number” in the browser.

1.“Input Data” function provides a default sample dataset including 300 cells and users can upload their local dataset. The counts matrix for the default dataset is in right column. scWizard supports multiple formats of input data. 'Data-rds': The input data is a Seurat object; 'Data-counts': The input data is the counts matrix, the rows are genes, and the columns are samples, including csv, rds, and h5 format; 'Data-10X': The input is a folder containing 10x data. To facilitate the first usage of scWizard, the example_data.rds file of testing set including 10,000 single cells of various cell types is randomly chosen and uploaded in Figshare (https://figshare.com/s/8f568e156d943754915e), users can selectively download and quickly get started with scWizard.<br>
![Fig1](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/inputdata.png)

2.“Quality Control” function provides necessary parameters for users to choose. When using, users can click the "view quality" button to view the original data quality; Then select parameters for quality control according to the situation.  he results are automatically saved to the folder where the data is entered.<br>
![Fig2](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/qc.png)

3.“Batch processing” function provides two ways to remove batch effects, enabling uses to select batch information from the metadata, and providing the necessary parameters for each method. As for the selectivity of the batch correction, ScWizard displays the graphical result before and after the batch effect. Users can perform batch effect removal analysis before cell annotation to determine whether there are significant batches between datasets including multiple scRNA-seq datasets generated using different technologies, obtained from different patients or samples, or from different experiments. Both Harmony and Seurat are the methods to remove batch effects, in which Harmony performs well on datasets with common cell types and different techniques and the comparatively shorter runtime of Harmony also makes it suitable for initial data exploration of large datasets. The 'View Batch' button allows user to view the batch of raw data. <br>
![Fig2](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/batch_processing.png)

4.“Cell annotation” function includes three modules: cell annotation is for major cell annotation, cell classification is for subcluster for cell subtypes, cell subcluster annotation is for cell subtypes annotation. The parameter including the number of network nodes could be adjusted. Users can adjust the hidden layer nodes according to three criteria: the number of hidden neurons should be between the size of the input layer nodes and the size of the output layer nodes; The number of hidden neurons should be 2/3 the size of the input layer nodes plus 2/3 the size of the output layer nodes; the number of hidden neurons should be less than twice the size of the input layer nodes. Users also can use default parameters. ScWizard provides training data for four major cell types and 47 cell subtypes, users could also use and input their own data set for model training in scWizard by clicking “Choose custom trainset” and “Choose custom label for trainset” button. The training data must be in H5 format and rows represent genes and columns represent cells. In addition, before using the annotation function, users need to click the "install Python" button to install the python environment; At the same time, you need to download the trainset and copy it to the "app/www/python" folder under the scwizard installation directory.Trainset download  address:Link:https://pan.baidu.com/s/1_J_17OKRLo4BQB33C_Ux3Q; Extraction code: 0ccq<br>
![Fig3](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/cell_anno.png)

5.“GSVA” function is for calculating GSVA score for specific gene sets or signaling pathways. GMT file could be uploaded by users and parameter can be adjusted. The results is displayed by below, and users can save the results.<br>
![Fig4](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/GSVA.png)

6.“inferCNV” function is for analyzing large-scale chromosome copy number alterations (CNA). Users can upload the prepared gene location file for copy number variation analysis. See this page for details of parameters. The results are saved in the "inforcnv_res" folder of the current working directory.<br>
![Fig5](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/inferCNV.png)

7.“TF-SCENIC” function is for screening out transcription factors with significant regulatory roles. The parameter can be adjusted and the results are presented in the form of heat map, which is saved in the int and output folders under the current path. Before running, users need to download the reference data set from the following address and copy it to the 'cistarget' folder under the current working directory.download address: https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/ <br>
![Fig5](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/scenic.png)

8.“CellphoneDB” function is for analyzing potential ligand receptor communication signals between cells and the parameter can be adjusted. The Seurat object is automatically processed as input data, and the results are saved in the cellphonedb_out folder under the current path.<br>
![Fig6](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/cellphonedb.png)

9.“Correlation Analysis” function is for calculating correlation among two or more variables to measure the degree of correlation between them. Pearson and Spearman correlation method, the cell type as well as the gene could be chosen by users’ needs. The results are saved in the folder under the current path.<br>
![Fig7](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/cor.png)

10.“Monocle” function is for constructing cell lineage development according to the changes of gene expression levels of different cell subpopulations. The parameter can be adjusted. The results are displayed in figure and are saved in the folder under the current path.<br>
![Fig8](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/Monocle.png)

11.“Plot” function displays the results of the above analysis modules and provides plotting function. The users can also upload their own data for plotting and all results can be downloaded by the users.<br>
![Fig9](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/plot.png)
