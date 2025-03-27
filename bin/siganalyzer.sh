#!/usr/bin/bash

#import sys
#import signatureanalyzer as sa

sample=$1
path_script=$2
output_dir=$3


#def main():      
signatureanalyzer $sample -t spectra -n 10 --reference cosmic2 --objective poisson -o $output_dir
#if __name__ == "__main__":


#main()

