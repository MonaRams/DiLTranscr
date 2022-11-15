# DiLTranscr
Dictionary learning for transcriptomic data

## Data download and setup 
The idea in this step is to (i), if neccessary, download data and to then (ii) setup the data such that meta data and transcriptomic count data is aligned and the data is composed of the same number of sample per type (for example, tissue types). These are all R scripts.

### 0. Package download for R scripts
In order for the other scripts to run, you will require some R packages. To download the packages, run the file RScripts./pckgDownloadsGT.R. Note, that this might ask you to choose wheather you want to update *all/some/none [a/s/n]* packages. Choose the preferred option by entering the respective character.

### 1. Data download
If you do not want to use your own data, we prepared a script to download two datasets from the GEO database: *GSE120795* and *GSE112004*. For the data download, run the script RScripts./geoDownloadGT.R. 

If you want to use other GEO datasets: You case also use this script to download other datasets from GEO database. You will therefore have to change the variable geoIDs in line 8 of RScripts./geoDownloadGT.R. Also, in lines 53-65, the column names of the meta data file that (i) entails the sample identifers (for example, titleCol <- "description" for dataset *GSE120795*) and (ii) entails information about the type (for example, typeCol <- "source_name_ch1" for dataset *GSE120795*) are defined. This has to be adjusted if you use other datasets. Note also, that for some GEO datsets, none of the meta data columns agrees with all the sample indetifiers from the dataset. Often, you only need to replace certain characters, for examples replace "-" with ".". This is all considered for the two examplary datasets.

As a results, this script will create three files in the Data directory: 1. The data file (*geoID* + *DatRAW.txt*), 2. the meta data file (geoID + *SampleInfoRAW.txt*), and 3. a list of the gene names (*geoID* + *GenesRAW.txt*). The order of the samples in the meta data file (each sample in one row) is now identical to the oder of the samples in the data file (each sample in one row) and the order of the gene names is the same as the genes in the data file (each gene in one row).

### 2. Data setup
For this step, for each dataset, you need to have two files in your Data directory:  1. The data file (*dataset name* + *DatRAW.txt*) and 2. the meta data file (*dataset name* + *SampleInfoRAW.txt*). For the data setup, run the script RScripts./setupDataGT.R. 

If you want to use datasets, other than the two exemplary GEO datasets, you need have the data and the meta data file in the Data directory, named dataset *name* + *DatRAW.txt* and dataset *name* + *SampleInfoRAW.txt*). The dataset has to be of shape *genes* x *samples* and the meta data file has to be of shape *samples* x *features*. Further, (i) the identifiers in the meta data file, column *title*, have to be identical to the column names of the data file and (ii) the column of the meta dataset with the type separating feature hase to be called *type*. You will have to adjust line 8 of this script as: geoIDS <- *dataset name*.

As a results, this script will create three files in the Data directory: 1. The data file (*DatPREPPED.txt*), 2. the meta data file (SampleInfoPREPPED.txt*), 3. the gene names (*GeneNamesPREPPED.txt*), and 4. a file which specifies the column name of the meta data file that entails information on the sample types (*TypeColsPREPPED.txt*). The prefix of these datasets is *the number of samples per type* *x* *the number of types* + *dataset name*.

## DLT run, evaluation & method comparison
In this step, the actual DLT method will be run. In addition, results for ICA, PCA, NMF, t-SNE, and UMAP will be computed. An evaluation is performed via clustering. The performed clustering method is kmeans. It is run for as many clusters, as there are different groups in the dataset (identified via the *type* column in the meta data file). The clustering is evaluated with the Adjusted rand index and the Adjusted mutual information. DLT is implemented in PythonScripts./normaliseGT.py. 

### Script details
In lines 46-51 of PythonScripts./normaliseGT.py, the datasets are automatically detected based on those files which are called *dataset name* + *DatPREPPED.txt* or *dataset name* + *DatPREPPED.csv*. You can also change the code here if you want to upload different datasets. When the data is uploaded, further files are required: (1) metadata file, named *dataset name* + *SampleInfoPREPPED.txt*, (2) type column specifier file, named *dataset name* + *TypeColPREPPED.txt*. 

In line 59, the methods are defined ("mthds <- ['DLT', 'PCA', 'ICA', 'NMF', 'tSNE', 'UMAP']"). If you want to run the code online for a subset of the methods, you may change this line. In the following lines, it is defined whether the dictionary or the coefficient matrix will contain positive elements only ("posDict <- False", "posCode <- False" - hence, in this setting postive and negative values are permitted for both matrices). Also, the value of labmda from the DLT functional is set to one (alphaParam <- 1). In line 65, the random seeds, the evaluations are run for, are defined. You can set this to any array you like - note though, that this variable has to be an array. 

The computation starts in line 70 (until line 877). It starts with uploading the data in lines 71-103. Further, it sets the number of dictiounary columns, as 1,...,*nGroups* (where *nGroups* is the number of different groups in the dataset, identified via the *type* column in the meta data file).  In lines, 126-149, normalisation is conducted. 

The actual method run starts in line 200. For all evaluated parameter values (compare line 217: *for col in dictCols:*), the dimension reduction is run in lines 230-329. Starting in line 296, the coefficient matrix is computed for all values 1,... *number of atoms*. Starting in line 431, kmeans is run and the clustering is evaluated. 

Lastly, after the methods are run for all evaluated parameter values, the results are stored in a files named *ARIs_* + *dataset name* + *RS_* + *current random seed* + *.txt*. Hence, for each dataset, for each evaluated random seed, one such file is stored. 

### Script output
As a results, PythonScripts./normaliseGT.py will create for each method (DLT/PCA/ICA/NMF/tSNE/UMAP) a directory called *method name* + *Results*. In this directory, there will be directories that contain the results for each dataset. Each of these directories is called *dataset name* + *_* + *value of labmda from the DLT functional* + *_* + *positive dictionary as Pdict/ positive coefficient matrix as Pcode* + *_* *name of transform algorithm*. In each of these directories, the learned dictionaries are stored as *dataset name* + *_* + *D/MBD* (depending on whether the mini batch Dictionary learning has been run) + *_* + *number of atoms/components*  + *Comp* + *_* + *value of random seed*.

The ARIs_... files contain the evaluation values. It contains 6 columns: In column 1, *RS*, the value of the random seed in given; in column 2, *Natoms*, the number of atoms are given; in column 3, *Nsparse*, the sparsity is given; in column 4, *ARI*, the ARI for the clustering of the respective decomposition is given; in column 5, *MI*, the respective Asjusted mutual information is given; in column 6, *RecError*, the reconstruction error is given.
