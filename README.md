# DBSOMA
Read me in progress...

This project is a java based tool for training a Self Organizing Map based on single cell RNA-seq data which is then used to evaluate gene lists for the degree of cluster formation by the DBSCAN algorithm. 
The goal is to allow users to in silico screen transcriptional perturbation data to make more accurate predictions 

Code herein was a written by Timothy Kunz. The SOM training and viewing implementation is a modified from a previous implementation which can be found https://github.com/seqcode/chromosom, with an associated publication https://doi.org/10.1016/j.ymeth.2020.07.002. Other code, including DBSCAN implementation and scripts are original code.

Below is the description of the overall workflow, a guide for using the sample data provided in this repository, and description of the methods used at each step.

## Table of Contents
- [Workflow](#Workflow)
- [Executable_Jars](#Executable_Jars)
- [Sample_Files](#Sample_files)
- [Preprocessing](#preprocessing)
- [Training](#Training)
- [Neighborhood_Calling](#Neighborhood_Calling)
- [Scanning](#Scanning)
- [Viewing](#Viewing)
- [Use](#Use)
- [Supportive_scripts](#Supportive_scripts)
- [License](#license)

## Workflow
The general workflow for use is:
  
  1) Start with single cell RNA seq data
  
  2) Preprocess the single cell data into 2 files required for SOM training:
  
    a) a pairwise gene x gene correlation matrix based on single cell sequencing data
    
    b) a list of gene names associated with each row/column of the square correlation matrix
  
  3) Train a self organizing map with the pairwise correlation matrix
  
  4) Scan the self organizing map (or maps) for pertubation gene lists which overlap with a set of desired states
  
  5) Interact with the resulting scan outputs to generate a list of perturbations predicted to effect genes of the desired state

  The SOM can also be used to view various gene lists before scanning, which can help to decide on target states

## Executable_Jars
Provided are a set of executable jar files which allow for use of the various functionalities described below
Source code for each is contained in this repository, but interaction with the code is not necessary for use
Files:
1) SingleCellParserTSV.jar
2) SparseParse.jar
3) FilePreviewer.jar
4) transposeTSV.jar
5) transposeNFixCSV.jar
6) SOMTrainer.jar
7) SOMViewer.jar
8) DBSOMA.jar


## Sample_Files
This directory contains a few files which can be used to test the various functionalities of SOMSCAN repository

Included files:
1) A randomly generated mock of single cell data to test the preprocessing 
2) A randomly generated mock of a correlation matrix and sample gene names to test training
3) A sample bash script for implementing the training
4) The output SOM from actual training on the single cell sequencing data set aquired from GSE114412
5) Several transcriptional state targets for scanning
6) Several selected perturbation response data sets which overlap with various target states
   Datasets aquired from https://maayanlab.cloud/sigcom-lincs/#/Download
7) An empty directory which will be used to output the scan results
8) A sample bash script for running the in silico scan method

## Preprocessing
The SOM training requires two files

Code in the package "correlationMatrixMaker" can be used to produce files required for SOM Training:
1) SingleCellParserTSV takes a tab delimeted file containing read count data with cell identifiers on columns and genes in rows a pairwise correlation matrix and associated list of gene names
2) SparseParse takes read count data in sparse matrix notation to produce a pairwise correlation matrix and associated list of gene names

Code in the package "file_interaction" can be used to pre-organize single cell data to transpose matricies:
1) transposeNFixCSV = convert a CSV file to tsv file
2) transposeTSV = transpose a TSV file to have cell identifiers on columns and genes in rows
3) FilePreviewer = functionality to read a first 2 lines of a file to ensure the matrix looks as expected


Code in the "CorMatAnalyzer" package can be used to interact with the matrix before beginning training
1) provides functionality to look up a gene or set of genes in a produced pairwise correlation matrix to ensure QC has not removed desired genes

## Training
The source code used for SOM training can be found in "som" package in this repository, and runs from the SOMaestro main method.

Training requires 3 files:
1) SOMtrainer.jar
2) A pairwise correlation matrix file
3) A line split list of genes corresponding to the matrix files

The training is able to run with multiple threads. Due to the size of the matrix files, training on a local machine is quite time consuming. It is recommended to run training on a cluster with many threads.

Training runs with the following command: 
    java -jar /path/to/SOMtrainer.jar /path/to/correlation_matrix /path/to/genes threads [no. of threads] {desired output matrix size}

An example bash script for running training "trainingscript.sh" in this repository

Training outputs 2 files to the directory from which SOMTrainer.jar was called. 
1) The output map, a line split list of SOM coordinates followed by the names of genes assigned to each node
2) An info file describing the running and results of the SOM

Once a SOM has been trained, it can be used by the viewer (#Viewing), or to scan gene lists and categorize the degree of clustering and overlap with a desired other list of genes (#Scanning). 

## Neighborhood_Calling
The NeighborhoodCalling class in the som_analysis package can be use to find the genes which make up a cluster observed by eye. This is useful for defining a target neighborhood, specifically when a given list of genes shows multiple discrete clusters and you'd like to specifically target one of those.

Requires 3 files:
1) NeighborhoodCalling.jar
2) Trained SOM file
3) list of genes to identify neighborhoods among

Training requires 3 arguments:
1) Int number of clusters
2) Path to Trained SOM file (MySom.som)
3) Path to gene list
4) Optional: List of specific genes which can be the center of a neighborhood

Calling runs with the following command: 
    java -jar /path/to/NeighborhoodCalling.jar (number of clusters) /path/to/MySom.som /path/to/GenesList (optional)/path/to/DesiredGeneCenters

## Scanning
The source code used for the scanning process is in the "mapScanning" package in this repository.

Scanning a SOM requires 4 directories:
1) A directory containing any number of trained SOM files
2) A directory containing a number of gene lists defining target states
3) A directory containing a number of gene lists to be scanned against the target states (perturbagen responses)
4) An empty directory were results will be populated
Optionally, add DBSCAN arguments:
5) Radius argument for DBSCAN (1 by default)
6) MinPts argument for DBSCAN (5 by default)

The following process is performed for each SOM in the directory

Scanning proceeds by first loading the SOM file. Next, lists of genes a projected onto the SOM and subjected to the DBSCAN algorithm. 
The DBSCAN implementation iterates as follows:
1) For each gene in the list, find its assigned node
2) increment the count of that assigned node and nodes within a radius (1 by default).
3) After that has been done for each gene, identify nodes with count >= min (5 by default). These nodes are _called_
4) Calculate _structure_ and _quality_ matrics
5) Calculate overlap of _called_ nodes between current list and target state lists

Next, 100,000 lists of genes are randomly generated and the same metrics are calculated for each.

This is done for each gene list, then the output directory is populated as follows:
1) Analysis_metadata.txt contains information about the parameters and files used in the scanning.
2) The values calculated for each input gene list are tabulated in SOM_Analysis.txt
3) The values calculated for each list with _strucutre_ >0.075 and _quality_ >0.057 are tabulated in Structure.txt
4) The values for list passing the above thresholds and also have _overlap_ >0.40 values with the first target gene list are tabulated in Overlap_Analysis.txt
5) The values for lists in which the product of _strucutre_ * _quality_ * _overlap_ above the product of the thresholds are tabulated in Hits_Analysis.txt
6) The alphabetically first listed file in the target states directory will automatically populate two directories, one with images of the intial projection of the list (Images) and the other with DBSCAN resultant projection (Trimmed), for up to 250 lists that are scored as hits. 

While both the Hits_Analysis and Overlap_Analysis outputs will contain high scoring gene lists, we recommend interacting the the Structure_Analysis.txt file directly for subsequent analysis. SOM_Analysis.txt can also be used, but is often too big for easy use in excel.

Scanning runs with the following command: 
    java -jar /path/to/DBSOMA.jar /path/to/MySom.som /path/to/PerturbationGeneListsDirectory /path/to/outputDirectory /path/to/TargetStateGeneListsDirectory (optional){DBSCAN radius} {DBSCAN minPts}

## Viewing
The Viewer, contained in the exectuable jar file "SOMviewer.jar" will first open a file explorer which users can select a trained SOM.

After openning a selected SOM, the window will display a projection of all genes in their assigned node, and offers several functionalities.

1) Save: this button will save a screenshot of the current projection in the directory from which "SOMviewer.jar" was launched. The name of the file will be whatever has been typed into the dilague box (gene...) as default
2) Search: This button will open a file explorer dilague box with which users can select a line seperated list of genes. The selected list will then be projected onto the map. This can be used to determine by the eye the degree of clustering a list of genes forms on the map
3) Search 2: This button will open a file explorer dilague box twice. Users select two line separated lists of genes, both of which will be projected on to the SOM. The first list selected will map in red, the second will map in blue
4) View swap: This button cycles through a few different appearances with which one can examine the current SOM projection
5) gene: This button will search for whatever gene name is currently written in the dialog box, and if found that gene's assigned node will be colored red.
6) Dialog box: Users can type here. Used to select the name of a saved file or select a gene to search for.

## Use
A guide to using the provided sample files:


## Supportive_scripts

This repository also contains several python scripts that were used for various data processing and endpoint analyses. However, their general usability is not optimized. I would recommend users personalize analyses based on their exact needs, using these scripts only as guides.
Scripts:


This Repository also contains R scripts used to generate GoTerm analysis plots 
Script:

This Repository contains Bash scripts which serve as examples for how to run the various jar files.

## License
MIT License

Copyright (c) 2024 Timothy Kunz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
