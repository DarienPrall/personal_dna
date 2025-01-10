#----------------------------------------------------------------------
#
# PROJECT INFORMATION
# Author: Darien Prall
# Start Date: 01-03-2025 (mmddyyyy)
# End Date: 
# Project Description: A personal project to further understand Python and Genomic Analysis
#
# PROJECT STEPS
# Task 0: Import and Clean Data
# Task 1: Count Allele Frequencies
# Task 2: Group RSIDs by Chromosome
# Task 3: Filter RSIDs by Chromosome and Position Range
# Task 4: Visualize Chromosome SNP Distribution
# Task 5: Lookup RSID Information from a Database
# Task 6: Identify Homozygous and Heterozygous SNPs
# Task 7: Identify Disease
# Task 8: Final Report
#
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# 
# DATASET INFORMATION
# File Type: .csv
# Size
# - rows(677,436) 
# - columns(5)
# - total(3,387,180)
# Features
# - rsid        object
# - chromosome  int64
# - position    int64
# - allele1     object
# - allele2     object
#
#----------------------------------------------------------------------

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

file_path = ('/Users/darienprall/Documents/GitHub/datasets/My_DNA.csv')

#----------------------------------------------------------------------
# IMPORT AND CLEAN DATA
# Step 0.1: Import Data
# Step 0.2: Change all headers to lowercase and remove whitespace
# Step 0.3: Get names of all columns and their types
#----------------------------------------------------------------------

# Step 0.1: Import Data
data = pd.read_csv(file_path)
#print(data.columns)

# Step 0.2: Change all headers to lowercase and remove whitespace
data.columns = data.columns.str.lower().str.replace(' ', '_')

# Step 0.3: Get names of all columns and their types
#print(data.dtypes)

data_columns = data.columns
#print(f"All columns in the dataset: {data_columns}")


#----------------------------------------------------------------------
# TASK 1: COUNT ALLELE FREQUENCIES
#----------------------------------------------------------------------
# Step 1.1: Store allele_1 and allele_2 to seperate variables
allele_1_list = data['allele_1'].to_numpy()
allele_2_list = data['allele_2'].to_numpy()
#print(f"Allele_1 List: {allele_1_list}")
#print(f"Allele_2 List: {allele_2_list}")

# Step 1.2: Create a dictionary to store all allele1 and allele2 counts
allele_1_counts = {
    'A' : 0,
    'C' : 0,
    'G' : 0,
    'T' : 0
 }

allele_2_counts = {
    'A' : 0,
    'C' : 0,
    'G' : 0,
    'T' : 0
 }

# Step 1.3: Create a function that loops through the allele column and updates the dictionary count
def allele_1_counter(allele_list):
    for allele in allele_list:
        if allele == 'A':
         allele_1_counts['A'] += 1
        elif allele == 'C':
            allele_1_counts['C'] += 1
        elif allele == 'G':
            allele_1_counts['G'] += 1
        else:
            allele_1_counts['T'] += 1

def allele_2_counter(allele_list):
    for allele in allele_list:
        if allele == 'A':
         allele_2_counts['A'] += 1
        elif allele == 'C':
            allele_2_counts['C'] += 1
        elif allele == 'G':
            allele_2_counts['G'] += 1
        else:
            allele_2_counts['T'] += 1

def count_all_alleles():
    allele_1_counter(allele_list=allele_1_list)
    allele_2_counter(allele_list=allele_2_list)

count_all_alleles()

# Step 1.4: Calculate Allele Frequency
# ALLELE 1
allele_1_total = sum(allele_1_counts.values())
allele_1_a_frequency = allele_1_counts['A'] / allele_1_total
allele_1_c_frequency = allele_1_counts['C'] / allele_1_total
allele_1_g_frequency = allele_1_counts['G'] / allele_1_total
allele_1_t_frequency = allele_1_counts['T'] / allele_1_total

#A ALLELE 2
allele_2_total = sum(allele_2_counts.values())
allele_2_a_frequency = allele_2_counts['A'] / allele_2_total
allele_2_c_frequency = allele_2_counts['C'] / allele_2_total
allele_2_g_frequency = allele_2_counts['G'] / allele_2_total
allele_2_t_frequency = allele_2_counts['T'] / allele_2_total

allele_1_frequencies = {
    'A': allele_1_a_frequency,
    'C': allele_1_c_frequency,
    'G': allele_1_g_frequency,
    'T': allele_1_t_frequency
}

allele_2_frequencies = {
    'A': allele_2_a_frequency,
    'C': allele_2_c_frequency,
    'G': allele_2_g_frequency,
    'T': allele_2_t_frequency
}


# Step 1.5: Print Results
#print(f"Allele_1 Counts: {allele_1_counts}")
#print(f"Allele_1 Frequencies: {allele_1_frequencies}")
#print(f"Allele_2 Counts: {allele_2_counts}")
#print(f"Allele_2 Frequencies: {allele_2_frequencies}")

# Step 1.6: Vizualize Results
allele_1_chart_file_path = ('/Users/darienprall/Documents/GitHub/outputs')
allele_1_key = list(allele_1_counts.keys())
allele_1_values = list(allele_1_counts.values())
plt.bar(allele_1_key, allele_1_values)
plt.grid(True, which = 'both', axis = 'y', linestyle = '--', alpha = 0.7)
plt.title('Allele 1 Counts')
plt.xlabel('Allele')
plt.ylabel('Count')
#plt.savefig(os.path.join(allele_1_chart_file_path, 'allele_1_counts.png'))
#plt.show()


allele_2_chart_file_path = ('/Users/darienprall/Documents/GitHub/outputs')
allele_2_key = list(allele_2_counts.keys())
allele_2_values = list(allele_2_counts.values())
plt.bar(allele_2_key, allele_2_values)
plt.grid(True, which ='major', axis = 'y',linestyle = '--', alpha = 0.7)
plt.title('Allele 2 Counts')
plt.xlabel("Allele")
plt.ylabel('Count')
#plt.savefig(os.path.join(allele_1_chart_file_path, 'allele_2_counts.png'))
#plt.show()

#----------------------------------------------------------------------
# TASK 2: GROUP RSIDs BY CHROMOSOME
#----------------------------------------------------------------------

# Step 2.1: Extract the Chromosome and RSID Columns
chromosomes = data['chromosome'].to_numpy()
rsids = data['rsid'].to_numpy()
# Step 2.2: Use dictionary to group SNPs by Chromosome

chromosome_rsid_dictionary = {}

for rsid, chrom in zip(rsids, chromosomes):
    if chrom not in chromosome_rsid_dictionary:
        chromosome_rsid_dictionary[chrom] = []
    chromosome_rsid_dictionary[chrom].append(rsid)

# Step 2.3: Display results in a table and save to a new file

chromosomes_rsid_data = []

for chrom, rsids in chromosome_rsid_dictionary.items():
    for rsid in rsids:
        chromosomes_rsid_data.append([chrom, rsid])

dataframe = pd.DataFrame(chromosomes_rsid_data, columns=['Chromosome', 'RSID'])
chromosome_rsid_folder_path = ('/Users/darienprall/Documents/GitHub/outputs')
chromosome_rsid_filename = ('chromosome_rsid_data.csv')
chromosome_rsid_file_path = os.path.join(chromosome_rsid_folder_path, chromosome_rsid_filename)
#dataframe.to_csv(chromosome_rsid_file_path, index = False)
#print(f"File saved to: {chromosome_rsid_file_path}")

#for i in chromosome_rsid_dictionary:
#    print(f"Chromosome {i}: {chromosome_rsid_dictionary[i][:10]}")

#for i in chromosome_rsid_dictionary:
#    print(f"Chromosome {i} has {len(chromosome_rsid_dictionary[i])} RSIDs")


#----------------------------------------------------------------------
# TASK 3: FILTER RSIDs BY CHROMOSE AND POSITION RANGE
#----------------------------------------------------------------------

# Step 3.1: Prompt User for a Chromosome and Position

TOTAL_POSITIONS = 249208772
MAX_CHROMOSOME_COUNT = 23

def get_chromosome():
    while True:
        lookup_chromosome = input(f"Enter a position number (1 to {MAX_CHROMOSOME_COUNT}): ")
        if lookup_chromosome.isdigit():
            lookup_chromosome = int(lookup_chromosome)
            if lookup_chromosome >= 1 and lookup_chromosome <= MAX_CHROMOSOME_COUNT:
                return lookup_chromosome
            else:
                print(f"Please enter a chromosome number between 1 and {MAX_CHROMOSOME_COUNT}")
        else: 
            print(f"Invalid Input. Please enter a number.")

def get_position():
    while True:
        lookup_position = input(f"Enter a position (1 to {TOTAL_POSITIONS}): ")
        if lookup_position.isdigit():
            lookup_position = int(lookup_position)
            if lookup_position >= 1 and lookup_position <= TOTAL_POSITIONS:
                return lookup_position
            else:
                print(f"Invalid Entry. Please enter a number between 1 and {TOTAL_POSITIONS}")
        else:
            print(f"Invalid Input. Please enter a number.")

def collect_user_inputs():
    chromosome = get_chromosome()
    position = get_position()
    return chromosome, position

chromosome, position = collect_user_inputs()

#print(f"You chose {chromosome} for chromosome and {position} for position")


