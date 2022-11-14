# DiLTranscr
Dictionary learning for transcriptomic data

## Data download and setup
The idea in this step is to (i), if neccessary, download data and to then (ii) setup the data such that meta data and transcriptomic count data is aligned and the data is composed of the same number of sample per type (for example, tissue types).

### 0. Package download for R scripts
In order for the other scripts to run, you will require some R packages. To download the packages, run the file pckgDownloadsGT.R. Note, that this might ask you to choose wheather you want to update *all/some/none [a/s/n]* packages. Choose the preferred option by entering the respective character.

### 2. Data download
If you do not want to use your own data, we prepared a script to download two datasets from the GEO database: *GSE120795* and *GSE112004*. For the data download, run the script geoDownloadGT.R. 

If you want to use other GEO datasets: You case also use this script to download other datasets from GEO database. You will therefore have to change the variable geoIDs in line 8 of geoDownloadGT.R. Also, in lines 53-65, the column names of the meta data file that (i) entails the sample identifers (for example, titleCol <- "description" for dataset *GSE120795*) and (2.) entails information about the type (for example, typeCol <- "source_name_ch1" for dataset *GSE120795*) are defined. This has to be adjusted if you use other datasets. Note also, that for some GEO datsets, none of the meta data columns agrees with all the sample indetifiers from the dataset. Often, you only need to replace certain characters, for examples replace "-" with ".". This is all considered for the two examplary datasets.

As a results, this script will create three files in the Data directory: 1. The data file (*geoID* + *DatRAW.txt*), 2. the meta data file (geoID + *SampleInfoRAW.txt*), and 3. a list of the gene names (geoID + "GenesRAW.txt"). The order of the samples in the meta data file (each sample in one row) is now identical to the oder of the samples in the data file (each sample in one row) and the order of the gene names is the same as the genes in the data file (each gene in one row).

### 3. Data setup
For this step, for each dataset, you need to have two files in your Data directory:  1. The data file (*dataset name* + *DatRAW.txt*) and 2. the meta data file (*dataset name* + *SampleInfoRAW.txt*). 

If you want to use datasets, other than the two exemplary GEO datasets, you need have the data and the meta data file in the Data directory, named dataset *name* + *DatRAW.txt* and dataset *name* + *SampleInfoRAW.txt*). The dataset has to be of shape genes x samples and the meta data file has to be of shape samples x features. Further, (i) the identifiers in the meta data file, column *title*, have to be identical to the column names of the data file and (ii) the column of the meta dataset with the type separating feature hase to be called *type*.

As a results, this script will create three files in the Data directory: 1. The data file (*DatPREPPED.txt*), 2. the meta data file (SampleInfoPREPPED.txt*), 3. the gene names (*GeneNamesPREPPED.txt*), and 4. a file which specifies the column name of the meta data file that entails information on the sample types (*TypeColsPREPPED.txt*). The prefix of these datasets is *the number of samples per type* *x* *the number of types* + *dataset name*.

