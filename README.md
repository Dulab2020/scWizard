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
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("Dulab2020/scWizard")
```
### Running scWizard
The following commands are used to start the graphical user interface (GUI).
```R
scWizard::run_app()
```

1.The single cell data should be saved in Seurat object and input in RDS file format. “Input Data” function provides a default sample dataset including 300 cells and users can upload local dataset. The counts matrix for the default dataset is in right column.<br>
![Fig1](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/inputdata.png)

2.“Batch processing” function provides two ways to remove batch effects, enabling uses to select batch information from the metadata, and providing the necessary parameters for each method. At the same time, ScWizard displays the graphical result before and after the batch effect, which can be saved by users.<br>
![Fig2](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/batch_processing.png)

3.“Cell annotation” function includes three modules: cell annotation is for major cell annotation, cell classification is for subcluster for cell subtypes, cell subcluster annotation is for cell subtypes annotation. The parameter including number of network nodes could be adjusted. ScWizard provides training data for major cell types and subtypes for four cell types, users can also use their own training data for cell annotation.Note：If you are using it for the first time, please click ‘Install Python’ first.<br>
![Fig3](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/cell_anno.png)

4.“GSVA” function is for calculating GSVA score for specific gene sets or signaling pathways. GMT file could be uploaded by users and parameter can be adjusted. The results is displayed by below, and users can save the results.<br>
![Fig4](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/GSVA.png)

5.“TF-SCENIC” function is for screening out transcription factors with significant regulatory roles. The parameter can be adjusted and the results are presented in the form of heat map, which is saved in the int and output folders under the current path.<br>
![Fig5](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/scenic.png)

6.“CellphoneDB” function is for analyzing potential ligand receptor communication signals between cells and the parameter can be adjusted. The Seurat object is automatically processed as input data, and the results are saved in the cellphonedb_out folder under the current path.<br>
![Fig6](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/cellphonedb.png)

7.“Correlation Analysis” function is for calculating correlation among two or more variables to measure the degree of correlation between them. Pearson and Spearman correlation method, the cell type as well as the gene could be chosen by users’ needs. The results are saved in the folder under the current path.<br>
![Fig7](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/cor.png)

8.“Monocle” function is for constructing cell lineage development according to the changes of gene expression levels of different cell subpopulations. The parameter can be adjusted. The results are displayed in figure and are saved in the folder under the current path.<br>
![Fig8](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/Monocle.png)

9.“Plot” function displays the results of the above analysis modules and provides plotting function. The users can also upload their own data for plotting and all results can be downloaded by the users.<br>
![Fig9](https://github.com/Dulab2020/scWizard/blob/master/inst/app/www/pageimage/plot.png)
