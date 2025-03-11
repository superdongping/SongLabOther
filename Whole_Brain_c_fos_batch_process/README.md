This repository contains MATLAB code designed for processing the registration and quantification of c-fos expression following the ABBA workflow (ABBA – Aligning Big Brains & Atlases; documentation available at ABBA Documentation). 
https://abba-documentation.readthedocs.io/en/latest/

The repository includes the following components:

Allen_Brain_Atlas_database.json: This file contains the Allen Brain Atlas names, IDs, and structural relationships.
Batch_Tif_Annotation_export.groovy: A Groovy script used with QuPath for the batch export of TIF images.
readImageJROI.m: A MATLAB script that facilitates the reading and parsing of JSON files.
tif_annotation_overlay_V2_2.m: A MATLAB script that processes the exported TIF images from QuPath along with the corresponding brain nucleus annotations obtained from the ABBA registration.
The output of the processing pipeline includes the identification of brain nuclei, along with quantification of c-fos positive cell counts and the corresponding areas for each region.
