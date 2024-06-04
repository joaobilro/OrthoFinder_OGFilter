# !/usr/bin/env python3
#
# -*- coding: utf-8 -*-
#
# OGF_Table2Folder.py
#
# Version 1.0
#
# This Python 3 script allows the user to filter the ortholog files, resulting from the output created 
# with Orthofinder 2.5.5 by Emms, D.M. et al. (2019), from Orthogroups.GeneCount.tsv to a new folder 
# based on the number of individual copies and on the percentage of missing taxa per file. This script 
# will make use of the file structure created by OrthoFinder 2.5.5, it might not work for other versions. 
#           _
#         ><_> 
#     
#        
# MIT License
# 
# Copyright (c) 2024 João Bilro (joaobilro)
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
import numpy as np
import os
import shutil

ogf = argparse.ArgumentParser(description="This Python 3 script allows the user to filter the orthologs resulting from an Orthofinder run"
                                            "based on the number of individual copies and the percentage of missing taxa, defined as thresholds."
                                            "Do not alter the files or structure created by OrthoFinder 2.5.5, because it will not work as intended.")

ogf.add_argument("--input", "-i", dest="input_orthofinder_folder", required=True, type=str, help="The full path to the main directory which has the OrthoFinder run(s) and fastas.")

ogf.add_argument("--output", "-o", dest="output_orthologs_folder", required=True, type=str, help="The full path to new folder that will contain the desired orthologs.")

ogf.add_argument("--missing", "-t", dest="missing_taxa", required=True, choices=range(0.1, 1.0), type=float, help="The proportion of missing taxa.")

ogf.add_argument("--copies", "-c", dest="maximum_copies", required=True, type=int, help="The maximum number of copies permitted.")

args = ogf.parse_args()

class Table2Folder:
    """Contains the functions necessary to parse the information stored in Orthogroups.GeneCount.tsv, and copies the desired orthologs to a new folder."""

    def __init__(self, orthofinder_directory, output_directory, missing_taxa, maximum_copies):
        self.of_dir = orthofinder_directory
        self.out_dir = output_directory
        self.missing = missing_taxa
        self.copies = maximum_copies
    
    def MainDir2OFolder(self):
        """This function will get to the OrthoFinder results folder starting from the main directory that was provided as input."""
        
        ### Check if there is a OrthoFinder folder
        if "OrthoFinder" in os.listdir(self.of_dir):
            os.chdir(os.path.join(self.of_dir, "OrthoFinder"))
        else:
            raise FileNotFoundError("Could not find OrthoFinder folder in {}.".format(self.of_dir))
        
        of_folder = os.chdir()

        ### Check how many runs there are, and choose the desired one     
        runs = os.listdir(of_folder)

        ### For only 1 run    
        if len(runs) == 1:
            os.chdir(os.path.join(of_folder, runs))
            print(f"Found OrthoFinder run folder. Retrieving orthologs...")


        ### For no runs
        elif len(runs) == 0:
            raise FileNotFoundError("Could not find any OrthoFinder runs in {}.".format(of_folder))
        
        ### For more than 1 run
        elif len(runs) > 1:
            print("Please choose the desired OrthoFinder run:")
            for index, run in enumerate(runs):
                print(f"{index + 1}: {run}")

            while True:
                try:
                    choice = int(input("Select the desired run by typing the corresponding number:"))  
                    if 1 <= choice <= len(runs):
                        selected_run = runs[choice - 1]
                        print(f"Found OrthoFinder run folder. Retrieving orthologs...")
                        break
                    else:
                        print(f"Out of bounds. Please enter a number between 1 and {len(runs)}.")
                
                except ValueError:
                    print("Invalid input. Please enter a number.")

            os.chdir(os.path.join(of_folder, selected_run))               

    def OFolder2List(self):
        """This function will get a list containing the desired orthologs starting from the OrthoFinder results folder."""

        ### Change to Orthogroups folder, where Orthogroups.GeneCount.tsv is located
        run_folder = os.getcwd()
        os.chdir(os.path.join(run_folder, "Orthogroups"))

        ### Define Orthogroups.GeneCount.tsv and check if it exists
        tsv_file_name = "Orthogroups.GeneCount.tsv"

        if not os.path.exists(tsv_file_name):
            raise FileNotFoundError(f"Could not find {tsv_file_name} in 'Orthogroups' folder.")
        
        with open(tsv_file_name, "r") as file:
            ortologs = file.readlines()
        
        ### Create a list with the desired orthologs and call the missing taxa and maximum copies arguments
        line_list = []
        missing = self.missing   ### Float
        missing_list = []
        copies = self.copies   ### Integer
        final_list = []

        for line in ortologs[1:]:
            ### Split each line by tabs, tab separated values file
            columns = line.strip().split('\t')

            ### Exclude the copy count column
            columns = columns[:-1]

            ### Append the remaining data to a list
            line_list.append(columns)

            ### Apply the thresholds defined by the user
            for entries in line_list:
                proportion = np.count_nonzero(entries == 0) / (len(entries) - 1)
                if proportion <= missing:
                    missing_list.append(entries)

            for entries in missing_list:
                ### Convert the numeric entries to integers so they can be matched to the copy number threshold
                integers = list(map(int, missing_list[1:]))
                if all(x <= copies for x in integers):
                    final_list.append(missing_list)

        ### Convert the final list into a .txt file containing solely the ortholog file names
        first_elements = [entry[0] for entry in final_list]

        ### Move one directory level up
        os.chdir(os.path.dirname(os.getcwd()))

        ### Write the previous variable to a new text file
        with open("FilteredOrthologs.txt", "w") as file:
            for element in first_elements:
                file.write(element + "\n")
        
        print(f"Orthologs succesfully retrieved. Copying to a new folder...")
                    

    def List2Folder(self):    
        """This function will copy and paste the desired orthologs to a new folder contained in the 'Orthogroups' folder."""

        ### Define current working directory
        current_directory = os.getcwd
        filtered_ogs = os.path.join(current_directory, "Filtered_Orthologs")

        ### Create the folder if there is none yet
        if not os.path.exists(filtered_ogs):
            os.makedirs(filtered_ogs)

        ### Read the file names from the .txt list
        with open("FilteredOrthologs.txt", "r") as file:
            og_names = file.read().splitlines()

        ### Initiate a count to keep track of how many files were copied
        og_count = 0
        
        ### Iterate through the orthologs
        for label in og_names:
            ### Check if any .fa file corresponds to the label
            file_path = os.path.join(current_directory, "Orthogroup_Sequences", label, ".fa")

            if os.path.exists(file_path):
                ### Copy the file to the Filtered_Orthologs directory
                destination = os.path.join(filtered_ogs, label, ".fa")
                shutil.copy2(file_path, destination)
                og_count += 1 
            else:
                print(f"File {label}.fa not found.")

        print(f"{og_count} orthologs were identified and successfully copied to the 'Filtered_Orthologs' folder.")

def main():
    ### Matching arguments with their intended variables
    orthofinder_directory = args.input_orthofinder_folder
    output_directory = args.output_orthologs_folder
    missing_taxa = args.missing_taxa
    maximum_copies = args.maximum_copies

    table2folder = Table2Folder(orthofinder_directory, output_directory, missing_taxa, maximum_copies)

    ### Start running the functions
    table2folder.MainDir2OFolder()
    table2folder.OFolder2List()
    table2folder.List2Folder()

if __name__ == "__main__":

    main()