# CoMSCER

COmparative Mutational Signature analysis on Coding and Extragenic Regions. The tool preforms a comparison of mutational signature analysis on two given conditions.

Conda is required for the correct usage of the tool. All dependencies are available in the CoMSCER.yml file attached. For conda environment creation use:

_conda env create -f CoMSCER.yml_

**_Input:_**

The tool takes as input two tab separated matrices of the two conditions to compare. 


Here is an example of the imput matrix: 



<img width="914" alt="Screenshot 2023-08-02 at 11 15 20" src="https://github.com/pbattuello/CoMSCER/assets/108470251/b9402db9-43ac-4d9e-aaa4-ce3175e696cf">



It contains the sample name, chromosome, 1 based position, reference genome (in the current version only GRCh38 is available, further reference genomes will be available in the next update), reference base, altered base, coding or non_coding.

Current version of the tool is able to evaluate Cosmic v2 and Cosmic v3.2 mutational signatures. Further datasets and custom signatures analysis will be available in the next release. 

Example of command:

_Syntax:_

./mutational_signature_benchmark_tool.sh matrix1 matrix2 output_dir SBS_cosmic_v2 SBS_cosmic_v3.2

_Example:_

./mutational_signature_benchmark_tool.sh matrix1 matrix2 /storage/test/signatures_result SBS5,SBS6,SBS7 SBS5,SBS6,SBS7a,SBS7b



