# Mat-SimAFM

05.09.2022 - Update to v1_2 with speed improvements 

A Matlab app and code that generates Simulated AFM images based the great works by Romain Amyot and Holger Flechsig:

Amyot, Romain, and Holger Flechsig. "BioAFMviewer: An interactive interface for simulated AFM scanning of biomolecular structures and dynamics." PLoS computational biology 16.11 (2020): e1008444.

## Requirements: 
MATLAB ‘Bioinformatics Toolbox’ to read local and online .pdb files. 
Without the ‘Bioinformatics Toolbox’ 3D coordinates can be read using ‘load .txt, .csv…’ button

## App Installation:
1.	Download the entire Simulation AFM app folder
2.	Open Mat SimAFM.mlappinstall and install the app to MATLAB
3.	Open Mat SimAFM from the app list within MATLAB

## For code only version:
1.	Run ‘Mat_SimAFM_code_only.m’ in MATLAB editor. 

## How to Use App:
1.	Load coordinates of structure to be simulated using either the ‘Load pdb ID’, ‘Load local .pdb file’ or ‘Load .txt, .csv…’ buttons. See below for description of each option:

1.1	Load pdb ID: Enter the 4-character alphanumeric identifier of your protein of interest found at https://www.rcsb.org or from accession codes within publications and then, Press the ‘Load pdb ID’ button. Note that this button requires an internet connection and may timeout if the connection is poor.

1.2	Load local .pdb file: Press this button and then Select a .pdb file from your local files

1.3	Load .txt, .csv…: Press this button and then select a .txt, .csv or .xlsx file with the x, y and z coordinates arranged as 3 single columns without any headers.


## Exporting: 

Press ‘Export to Matlab’ to export image data and any line profiles to your Matlab workspace 

Press ‘Export Tiff’ to export a 32-bit greyscale image to a user defined location. Height information will be encoded in the image intensity in units of Angstrom. X, Y pixel dimension will match that of the AFM Surface image and will therefore depend on the user defined ‘Pixel Sampling’ value. 

19.08.2022 - Export Tiff error fixed 

To Export the current AFM Surface image as a JPEG, PNG or PDF go to the upper right-hand corner of the AFM Surface image until the save option and zoom icons appear. Select the arrow and then floppy disk icon for Save As… and save the image as required (these symbols may differ between MATLAB versions, this is for MATLAB 2022a).
