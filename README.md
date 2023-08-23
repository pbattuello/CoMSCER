# CoMSCER

COmparative Mutational Signature analysis on Coding and Extragenic Regions. The tool preforms a comparison of mutational signature analysis on two given conditions using different tools, reference datases and genomic coding types.

**_Aim:_**

Given two conditions and a set of SBS signatures, the tool allows the identification of the best workflow for mutational signature analysis based on technical and biological readouts.

_**Requirement:**_

Conda is required for the correct usage of the tool. All dependencies are available in the CoMSCER-env.yml file attached. For conda environment creation use:

_conda env create -f CoMSCER-env.yml_

_OR_

_mamba env create -f CoMSCER-env.yml_


**_Installation:_**

After having downloaded the folder "CoMSCER" do the folloging commands:

	cd CoMSCER
 	bash install.sh

  
**_Input:_**

The tool takes as input two tab separated matrices representing the two conditions to compare. 


Here is an example of the input matrix: 


<img width="914" alt="Screenshot 2023-08-02 at 11 15 20" src="https://github.com/pbattuello/CoMSCER/assets/108470251/d9ca3a6d-be34-422b-8de5-bcd04ef56452">



It contains the sample name, chromosome, 1-based position, reference genome (in the current version only GRCh38 is available, further reference genomes will be available in the next update), reference base, altered base, coding or non_coding.

Current version of the tool is able to evaluate Cosmic v2 and Cosmic v3.2 mutational signatures. Further datasets and custom signatures analysis will be available in the next release. 

Be careful that CoMSCER does not check if the mutations are coding or not coding but its evaluation is only based on the input file annotation.

Example of command:

_Syntax:_

./mutational_signature_benchmark_tool.sh matrix1 matrix2 output_dir SBS_cosmic_v2 SBS_cosmic_v3.2

_Example:_

./mutational_signature_benchmark_tool.sh matrix1 matrix2 /storage/test/signatures_result SBS5,SBS6,SBS7 SBS5,SBS6,SBS7a,SBS7b

_**Output:**_

In the output_dir you will find three directories:
   1. **SBS_signature_contributions**: in this this directory you will find all the relative signature contributions for each matrix given in input. Specific results for each tool and reference datasets will be present.
   2. **Cosine_Similarity**: in this  directory you will find all the cosine similarity values for matrix given in input. Specific results for each tool and reference datasets will be present.
   3. **Figures**: in this directory four plots will be present showing signature contributions and cosine similarity for each reference datasets.
   4. **data**: in this directory you will find the initial mutational matrixes, and a folder "_logs_" containing all log and error files of the pipeline.


  _Citations:_

  We thanks all the authors of the following works. All the following tools are used by CoMSCER:

  - Blokzijl F, Janssen R, van Boxtel R, Cuppen E. MutationalPatterns: comprehensive genome-wide analysis of mutational processes. Genome Med. 2018 Apr 25;10(1):33. doi: 10.1186/s13073-018-0539-0. PMID: 29695279; PMCID: PMC5922316.
  - Rosenthal R, McGranahan N, Herrero J, Taylor BS, Swanton C. DeconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution. Genome Biol. 2016 Feb 22;17:31. doi: 10.1186/s13059-016-0893-4. PMID: 26899170; PMCID: PMC4762164.
  - Díaz-Gay M, Vangara R, Barnes M, Wang X, Islam SMA, Vermes I, Narasimman NB, Yang T, Jiang Z, Moody S, Senkin S, Brennan P, Stratton MR, Alexandrov LB. Assigning mutational signatures to individual samples and individual somatic mutations with SigProfilerAssignment. bioRxiv [Preprint]. 2023 Jul 11:2023.07.10.548264. doi: 10.1101/2023.07.10.548264. PMID: 37502962; PMCID: PMC10369904.
  - A. Degasperi et al. Substitution mutational signatures in whole-genome-sequenced cancers in the UK population. Science, doi:10.1126/science.abl9283, 2022.
  - Bergstrom EN, Huang MN, Mahto U, Barnes M, Stratton MR, Rozen SG, Alexandrov LB. SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. BMC Genomics. 2019 Aug 30;20(1):685. doi: 10.1186/s12864-019-6041-2. PMID: 31470794; PMCID: PMC6717374.

  
	 


