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

1.The single cell data should be saved in Seurat object and input in RDS file format. “Input Data” function provides a default sample dataset including 300 cells and users can upload local dataset. The counts matrix for the default dataset is in right column.<br>
![Fig1](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/inputdata.png)

2.“Quality Control” function provides necessary parameters for users to choose. When using, users can click the "view quality" button to view the original data quality; Then select parameters for quality control according to the situation.  he results are automatically saved to the folder where the data is entered.<br>
![Fig2](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/qc.png)

3.“Batch processing” function provides two ways to remove batch effects, enabling uses to select batch information from the metadata, and providing the necessary parameters for each method. At the same time, ScWizard displays the graphical result before and after the batch effect, which can be saved by users.<br>
![Fig2](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/batch_processing.png)

4.“Cell annotation” function includes three modules: cell annotation is for major cell annotation, cell classification is for subcluster for cell subtypes, cell subcluster annotation is for cell subtypes annotation. The parameter including number of network nodes could be adjusted. ScWizard provides training data for major cell types and subtypes for four cell types, users can also use their own training data for cell annotation. In addition, before using the annotation function, users need to click the "install Python" button to install the python environment; At the same time, you need to download the trainset and copy it to the "app/www/python" folder under the scwizard installation directory.<br>
![Fig3](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/cell_anno.png)

5.“GSVA” function is for calculating GSVA score for specific gene sets or signaling pathways. GMT file could be uploaded by users and parameter can be adjusted. The results is displayed by below, and users can save the results.<br>
![Fig4](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/GSVA.png)

6.“inferCNV” function is for analyzing large-scale chromosome copy number alterations (CNA). Users can upload the prepared gene location file for copy number variation analysis. See this page for details of parameters. The results are saved in the "inforcnv_res" folder of the current working directory.<br>
![Fig5](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/inferCNV.png)

7.“TF-SCENIC” function is for screening out transcription factors with significant regulatory roles. The parameter can be adjusted and the results are presented in the form of heat map, which is saved in the int and output folders under the current path. Before running, you need to download the reference data set from the following address and copy it to the 'cistarget' folder under the current working directory.<br>
![Fig5](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/scenic.png)

8.“CellphoneDB” function is for analyzing potential ligand receptor communication signals between cells and the parameter can be adjusted. The Seurat object is automatically processed as input data, and the results are saved in the cellphonedb_out folder under the current path.<br>
![Fig6](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/cellphonedb.png)

9.“Correlation Analysis” function is for calculating correlation among two or more variables to measure the degree of correlation between them. Pearson and Spearman correlation method, the cell type as well as the gene could be chosen by users’ needs. The results are saved in the folder under the current path.<br>
![Fig7](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/cor.png)

10.“Monocle” function is for constructing cell lineage development according to the changes of gene expression levels of different cell subpopulations. The parameter can be adjusted. The results are displayed in figure and are saved in the folder under the current path.<br>
![Fig8](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/Monocle.png)

11.“Plot” function displays the results of the above analysis modules and provides plotting function. The users can also upload their own data for plotting and all results can be downloaded by the users.<br>
![Fig9](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/plot.png)
