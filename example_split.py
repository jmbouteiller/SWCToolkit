
#********* MST BE TESTED !!!!  *********
# NOT FUNCTIONAL YET. 
# ERROR: 
# AttributeError: module 'SWCToolkit.swcToolkit' has no attribute 'read_swc_file'

#import modules
#from swcToolkit import swcToolkit
#import os

import os
import sys

# Get the absolute path of the current file
current_file_path = os.path.abspath(__file__)
# Get the directory name from the file path
current_directory = os.path.dirname(current_file_path)

#############################################
# PARAMETERS
#############################################
# Define the input SWC file and output directory
# Note: Ensure the input file path is correct and the output directory exists
input_file_name = 'simpleSWCFile_2somas.swc'  # replace with your SWC file path
input_file = current_directory + '\\data\\input\\'+ input_file_name  # replace with your SWC file path
output_path = current_directory + '\\data\\output'  # replace with your desired output directory

print ("input_file: ", input_file)
print ("output_path: ", output_path)

# Change the current working directory to the SWC toolkit directory
path_to_swctoolkit = os.path.join(current_directory, 'swctoolkit')

# Optional: Verify the change
# print(f"Current working directory: {os.getcwd()}")

# Add the SWC toolkit directory to the system path (needed for interactive use)
# (not needed for terminal since we have a launch.json and settions.json)
sys.path.append(path_to_swctoolkit )  # Adjust the path to add SWC toolkit directory

print ("New path: ", sys.path)
print ("DONE ")

print ("Current directory: ", os.getcwd())

#import swcToolkit  # Import the SWC toolkit
# If the import fails, you can try the following alternative import statements:  

# from directory swctoolkit, import the class Swctoolkit
from swctoolkit import Swctoolkit    

# Create an instance of the SWC toolkit
swctk = Swctoolkit.Swctoolkit()  

swctk.split_swc(input_file=input_file, output_path=output_path)  # Call the split function
#swcToolkit.split_swc(swcToolkit, input_file=input_file, output_path=output_path)  # Call the split function


