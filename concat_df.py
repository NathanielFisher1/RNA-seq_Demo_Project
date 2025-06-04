#!/usr/bin/env python

# script for concatenating verse output files into a single dataframe

import argparse
import pandas
import os

parser = argparse.ArgumentParser()
 
parser.add_argument("-i", "--input", help='A list of the VERSE output filenames provided by snakemake', dest="input", required=True, nargs='+')
parser.add_argument("-o", "--output", help='The output file name and path provided by snakemake', dest="output", required=True)

args = parser.parse_args()

concat = pandas.concat([pandas.read_csv(df, sep='\t', header=0, names = ['gene', '{}'.format(os.path.basename(df.split('.')[0]))], index_col='gene') for df in args.input], axis=1)

concat.to_csv(args.output)
