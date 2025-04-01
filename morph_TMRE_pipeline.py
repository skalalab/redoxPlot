# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 17:36:07 2024
Mitochondria Morphology Analysis
@author: Wenxuan Zhao
"""
#%%
import os 
import shutil
import re
import cv2
import numpy as np
from skimage import measure
import matplotlib.pyplot as plt
# rt: the path to the folder that contains a list of folders, each of which contains raw data of a certain well on certain day
rt = "Z:\Ty\TMRE intensity data - Redox\A549_TMRE Data"
rt =  "/mnt/z/Wenxuan/Ty/day6/"
rt = r"Z:\Ty\RD-Redox project\TMRE intensity data - Redox\Panc1_TMRE_extracted-data"
os.chdir(rt)
# analysis_dir is the path to the folder that stores all the files needed for analysis (the raw data file , mask, intermediate analysis file, to the final csv data)
analysis_dir = os.path.join(rt, "analysis")
# this script will call imageJ to execute a macro file, thus we need to specify the path to the imageJ and to the macro
imagej_path = r"C:\Users\tweaver\DOWNLO~1\FIJI-W~1\Fiji.app\ImageJ-win64"
#imagej_path = "~/Fiji.app/ImageJ-linux64"
# imagej_path
macro_name = 'Auto_Threshold_ImageJ.ijm'
macro_path = os.path.join(rt, macro_name)     # Path to the .ijm macro file
imgFileSuffix = ".tif"
maskFileSuffix = ".tiff"

#%%
## Step 0: make a special folder called "analysis"; categorize and copy raw data to their belonged well_day and ROI
# rt: a list of well_day folders, each of which contains a list of roi folders. the ROI folders are named starting with a number that denotes the ROI number
well_days = {}
if not os.path.exists(analysis_dir):
    os.mkdir(analysis_dir) 

# iterates through all the well_day folders 
for d in os.listdir(rt):
    # assumes the rt has four types of folders: the one starting with . (system hidden folders); the analysis folder we just made; the sdts folder; besides all of the above are all well_day folders. 
    # Therefore, make sure your rt folder contains only the well_day folders 
    if os.path.isdir(d) and not d.startswith('.') and not d.endswith('analysis') and not d.endswith('sdts'):
        # make an analysis folder for the current well_day 
        analysis_dir_d = os.path.join(analysis_dir, d)
        well_days[d] = []
        if not os.path.exists(analysis_dir_d):
            os.mkdir(analysis_dir_d)

        # given a well_day folder iterates through all the ROI folders  
        for roi_dir in os.listdir(d):
            data_dir = os.path.join(d, roi_dir)
            if os.path.isdir(data_dir) and not data_dir.startswith('.'): 
                # assume ROI folders are named with the ROI number as prefix
                match = re.match(r'(\d+)', roi_dir)
                roi_num = int(match.group(0))
                well_days[d].append(roi_num)
                # make an analysis folder for the current ROI
                roi_dir = os.path.join(analysis_dir_d, str(roi_num))
                if not os.path.exists(roi_dir):
                    os.mkdir(roi_dir)
                # iterate through all the data file inside the ROI folder, and copy them 
                for f in os.listdir(data_dir): 
                    # assumes the image file that contains "Ch3" as the NADH intensity image
                    # copy it to its belonged well_day and ROI of the analysis folder and rename it to {well_day}_N{ROI_num}.tif(f)
                    if "Ch3" in f and f.endswith(imgFileSuffix):
                        name = d + '_N' + str(roi_num) + imgFileSuffix
                        try:
                            shutil.copy2(os.path.join(data_dir, f), roi_dir) 
                            os.rename(os.path.join(roi_dir, f) , os.path.join(roi_dir,name))       
                        except FileNotFoundError as e:
                            print(f"{name} not found error: {e}")
                        except PermissionError as e:
                            print(f"Permission error: {e}")
                        except Exception as e:
                            print(f"An error occurred: {e}")
                    # assumes the image file that contains "Ch1" as the Mitochrondria intensity image
                    # copy it to its belonged well_day and ROI of the analysis folder and rename it to {well_day}_M{ROI_num}.tif(f)
                    elif "Ch1" in f and f.endswith(imgFileSuffix):
                        name = d + '_M' + str(roi_num) + imgFileSuffix
                        try:
                            shutil.copy2(os.path.join(data_dir, f), roi_dir) 
                            os.rename(os.path.join(roi_dir, f) , os.path.join(roi_dir,name))       
                        except FileNotFoundError as e:
                            print(f"{name} not found error: {e}")
                        except PermissionError as e:
                            print(f"Permission error: {e}")
                        except Exception as e:
                            print(f"An error occurred: {e}")
                    # assumes the sdt file that contains "LifetimeData" as the lifetime data
                    # copy it to its belonged well_day and ROI of the analysis folder and rename it to {well_day}_sdt{ROI_num}.sdt
                    elif "LifetimeData" in f and f.endswith('.sdt'):
                        name = d + '_sdt' + str(roi_num) + ".sdt"
                        try:
                            shutil.copy2(os.path.join(data_dir, f), roi_dir) 
                            os.rename(os.path.join(roi_dir, f) , os.path.join(roi_dir,name))       
                        except FileNotFoundError as e:
                            print(f"{name} not found error: {e}")
                        except PermissionError as e:
                            print(f"Permission error: {e}")
                        except Exception as e:
                            print(f"An error occurred: {e}")
                        
#%% Step 1: for each well_day, for each ROI, get cell mask from nadh intensity image ({well_day}_N{ROI_num}.tif(f)) using Napari (output: ROI-level cell mask)
### open the roi folders and uses Napari to segment mask and stores the output to the same folder. 
# ******The name of each mask should be {well_day}_N{ROI_num}_cellmask.tiff*********

## Step 2: separate the ROI-level cell mask to single-cell mask

def splitNapariMask(analysis_dir, well_day, roi, mask="cellmask"):
    roi_mask = well_day + "_N" + str(roi) + "_" + mask + maskFileSuffix
    image_path = os.path.join(analysis_dir, well_day, str(roi), roi_mask)
    image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
    print("Original image shape:", image.shape)
    print("Original image channels:", image.shape[-1])
    # If the image still has an unexpected number of channels, try to extract the first channel
    if image.ndim == 3 and image.shape[-1] != 1:
        image = image[:, :, 0]
    else:
        image = image  
    output_directory = os.path.join(analysis_dir, well_day , str(roi))
    unique_values = np.sort(np.unique(image)[np.unique(image) != 0])
    print("Unique values in the image:", unique_values)
    for cell_value in unique_values:
        binary_mask = np.where(image == cell_value, 1, 0)
        # binary_mask = binary_mask.astype(np.uint8)
        roi_mask = roi_mask.split(".")[0]
        cv2.imwrite(os.path.join(output_directory, f'{roi_mask}_{cell_value}{maskFileSuffix}'), binary_mask)   

for well_day in well_days.keys():
    for roi in well_days[well_day]:
       splitNapariMask(analysis_dir, well_day, roi)

## Output: a list of tiff files that are named as:  {well_day}_N{ROI_num}_cellmask_{cell_number}.tif(f).                                             
# The three things that are enclosed in {} represents the three identification levels of a certain single cell.



#%% Step 3: get the mito mask by multplying mito intensity with each single-cell mask

def getMitoMask(analysis_dir, roi, well_day):
    current_dir = os.path.join(analysis_dir, well_day, str(roi))
    mito = well_day + "_M" + str(roi) + ".tif"
    cellmask_pattern = re.compile(r'^(?!\.).*_cellmask_(\d+)\.tiff$')
    cell_masks = []
    for f in os.listdir(current_dir):
        match = cellmask_pattern.search(f)
        if match:
            cell_masks.append(f)
            
    for cell_mask in cell_masks:
        image_path = os.path.join(current_dir, cell_mask)
        print(image_path)
        cell_mask_image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
        image_path = os.path.join(current_dir, mito)
        mito_image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
        output_directory = os.path.join(analysis_dir, well_day , str(roi))
        if cell_mask_image.shape == mito_image.shape:
            # Multiply the images
            cell_mask_image_f = cell_mask_image.astype(np.uint16)
            mito_image_f = mito_image.astype(np.uint16)
            multiplied_image = cv2.multiply(cell_mask_image_f, mito_image_f)
            cv2.imwrite(os.path.join(output_directory,cell_mask.replace("cellmask", "MitoI")), multiplied_image)
        else: 
            print("size mismatch: cell_mask and mito intensity_image")

for well_day in well_days.keys():
    for roi in well_days[well_day]:
        getMitoMask(analysis_dir, roi, well_day)

# output: {well_day}_N{ROI_num}_MitoI_{cell_number}.tif(f)
        
#%% Step 4: get the sharpened mito mask
# for each single-cell level mito intensity image, we get a mask
import subprocess

def sharpenMitoMask(analysis_dir, roi, well_day):
    current_dir = os.path.join(analysis_dir, well_day, str(roi))
    png_dir = os.path.join(current_dir, "png")
    if not os.path.exists(png_dir):
        os.mkdir(png_dir)
        
    mitoI_pattern = re.compile(r'^(?!\.).*_MitoI_(\d+)\.tiff$')
    mitoIs = []
    for f in os.listdir(current_dir):
        match = mitoI_pattern.search(f)
        if match:
            mitoIs.append(f)
    
    for mitoI in mitoIs:
        input_path = os.path.join(current_dir, mitoI)
        output_path = os.path.join(current_dir,  mitoI.replace("MitoI", "MitoMask"))
        mask_path = os.path.basename(input_path)
        input_args = input_path + ',' +  output_path + ',' + mask_path  
        
        # Construct the full command
        command = [imagej_path, '-macro', macro_path, input_args]
        
        # Run the command
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred: {e}")
            
       # Step 4.2: remove the small objects
        image = cv2.imread(output_path)
        labels = measure.label(image)
        regions = measure.regionprops(labels)
        threshold = 10 # 5, 15, 20? need to be optimized
        region_areas = np.array([region.area for region in regions])
        cleaned_image = np.where(
            np.isin(labels, np.where(region_areas >= threshold)[0] + 1),
            image,
            0
        )
        cv2.imwrite(output_path, cv2.cvtColor(cleaned_image, cv2.COLOR_BGR2GRAY))
        print("Small objects removed")
for well_day in well_days.keys():
    for roi in well_days[well_day]:
        sharpenMitoMask(analysis_dir, roi, well_day)
        
#%% match nuclei label with cell label

def match_label(analysis_dir, roi, well_day):
    cur_dir = os.path.join(analysis_dir, well_day, str(roi))
    roi_nuclei_mask = os.path.join(cur_dir, well_day + "_N" + str(roi) + "_" + "nucleimask.tiff")
    roi_cell_mask = os.path.join(cur_dir, well_day + "_N" + str(roi) + "_" + "cellmask.tiff")
    nuclei_mask = cv2.imread(roi_nuclei_mask, cv2.IMREAD_UNCHANGED)
    cell_mask = cv2.imread(roi_cell_mask, cv2.IMREAD_UNCHANGED)
    if nuclei_mask is None or cell_mask is None: 
        return 
    cell_mask_labels = np.sort(np.unique(cell_mask)[np.unique(cell_mask) != 0])
    nuclei_mask_labels = np.unique(nuclei_mask)
    for cell_label in cell_mask_labels: 
        cell_binary_mask = (cell_mask == cell_label).astype(int)
        for nuclei_label in nuclei_mask_labels: 
            nuclei_binary_mask = (nuclei_mask == nuclei_label).astype(int)
            nuclei_area = np.sum(nuclei_binary_mask)
            mask_intersection = cell_binary_mask * nuclei_binary_mask
            intersection_area = np.sum(mask_intersection > 0)
        # Check if the intersection area is at least 70% of the cell are
            if intersection_area < 0.7 * nuclei_area:
                continue
            else: 
                # minus pixels with mitochrondia 
                print(os.path.join(cur_dir, well_day + "_N" + str(roi) + "_MitoMask_" + str(cell_label) +".tiff"))
                mitoMask_image = cv2.imread(os.path.join(cur_dir, well_day + "_N" + str(roi) + "_MitoMask_" + str(cell_label) +".tiff"), cv2.IMREAD_UNCHANGED)
                mito_binary_mask = (mitoMask_image > 0).astype(int)
                nuclei_binary_mask -=  mito_binary_mask
                nuclei_binary_mask [nuclei_binary_mask < 0] = 0
    
                cv2.imwrite(os.path.join(cur_dir, well_day + "_N" + str(roi) + "_nucleimask_" + str(cell_label) +".tiff"), nuclei_binary_mask)
            
            
for well_day in well_days.keys():
    for roi in well_days[well_day]:
        match_label(analysis_dir, roi, well_day)
    
#%% Step 5: TMRE analysis
def tmreAnalysis(analysis_dir, roi, well_day, csvPath):
    current_dir = os.path.join(analysis_dir, well_day, str(roi))
    mitoI_pattern = re.compile(r'^(?!\.).*_MitoI_(\d+)\.tiff$')
    mitoIs = []

    # step 1: register nucleus with their corresponing cell
    # for l in 
    for f in os.listdir(current_dir):
        match = mitoI_pattern.search(f)
        if match:
            mitoIs.append(f)
    rows = []
    for mitoI in mitoIs:
        mitoI_path = os.path.join(current_dir, mitoI)
        mitoI_image = cv2.imread(mitoI_path, cv2.IMREAD_UNCHANGED)
        nuclei_mask_path =  os.path.join(current_dir, mitoI.replace("MitoI", "nucleimask"))
        # try find the corresponding nuclei mask
        nuclei_mask_image = cv2.imread(nuclei_mask_path, cv2.IMREAD_UNCHANGED)
      
        if nuclei_mask_image is None:
            continue 
        nuclei_binary_mask = (nuclei_mask_image > 0).astype(int)
        mitoMask_path = os.path.join(current_dir, mitoI.replace("MitoI", "MitoMask"))
       
        mitoMask_image = cv2.imread(mitoMask_path, cv2.IMREAD_UNCHANGED)
        
        mito_binary_mask = (mitoMask_image > 0).astype(int)
        mito_result = np.multiply(mitoI_image, mito_binary_mask)
        cv2.imwrite(mitoI_path.replace("MitoI", "TMRE_Mito"), mito_result)
        mito_area = np.count_nonzero(mito_result > 0)
        mito_intensity_sum = np.sum(mito_result)
        if mito_result.dtype == 'uint8':
            mito_intensity_sum /= 255
        elif mito_result.dtype == 'uint16':
            mito_intensity_sum /= 65535
        elif mito_result.dtype == 'uint32':
            mito_intensity_sum /= (2**32-1)
        try:
            mito_intensity_avg = mito_intensity_sum / mito_area
        except ZeroDivisionError:
            print ("Image is all black")

        # nucleus tmre analysis
      
        nuclei_area = np.sum(nuclei_binary_mask)
        nuclei_result = np.multiply(mitoI_image, nuclei_binary_mask)
        cv2.imwrite(mitoI_path.replace("MitoI", "TMRE_Nuclei"), nuclei_result)
        nuclei_intensity_sum = np.sum(nuclei_result)

        if nuclei_result.dtype == 'uint8':
            nuclei_intensity_sum /= 255
        elif nuclei_result.dtype == 'uint16':
            nuclei_intensity_sum /= 65535
        elif nuclei_result.dtype == 'uint32':
            nuclei_intensity_sum/= (2**32-1)
        try:
            nuclei_intensity_avg = nuclei_intensity_sum / nuclei_area
        except ZeroDivisionError:
            print ("Image is all black")
            
        # # method 2: use nth percentile 
        # n = 40
        # subset_array = mitoI_image[(nuclei_mask == label)]
        # plt.hist(subset_array , bins=100, alpha=0.75, edgecolor='black')
        # plt.show()
        # nth_percentile = np.percentile(subset_array, n)
        # if nuclei_result.dtype == 'uint8':
        #     nth_percentile /= 255
        # elif nuclei_result.dtype == 'uint16':
        #     nth_percentile /= 65535
        
        # DOI: 10.2337/db06-0757
        RT = 61.5
        delta_psi_1 = -1 * RT * np.log(mito_intensity_avg / nuclei_intensity_avg)
     #   delta_psi_2 = -1 * RT * np.log(mito_intensity_avg / nth_percentile)
        nuclei_area_threshold = 160 # adjustable by user
        if ( nuclei_area >= nuclei_area_threshold): 
            row = [mitoI, mito_area, mito_intensity_sum, mito_intensity_avg, nuclei_area, \
                   nuclei_intensity_sum, nuclei_intensity_avg, delta_psi_1]
                
            rows.append(row)
    return rows

import pandas as pd
from datetime import date
current_date = date.today()
for well_day in well_days.keys():
    csvPath = os.path.join(analysis_dir, well_day, f"{current_date}_tmre.csv")
    headers = ['file_name', 'mito_area', 'mito_intensity_sum', 'mito_intensity_avg', "nuclei_area", \
    "nuclei_intensity_sum", "nuclei_intensity_avg", "Delta Psi 1"]
    roi_rows = []
    for roi in well_days[well_day]:
        cells= tmreAnalysis(analysis_dir, roi, well_day, csvPath)
        roi_rows += cells 
    if len(roi_rows) > 0: 
        df = pd.DataFrame(roi_rows, columns=headers)
        df.to_csv(csvPath, index=False)
        
#%% Plot

import pandas as pd
import matplotlib.pyplot as plt

tmre = pd.read_csv(os.path.join(rt, "analysis", "A1", "2024-11-22_tmre.csv"), header=None, names = ['file_name', 'mito_area', 'mito_intensity_sum', 'mito_intensity_avg', "nuclei_area", \
"nuclei_intensity_sum", "nuclei_intensity_avg", "Delta Psi"])

plt.figure(figsize=(10, 6))
plt.boxplot(tmre["Delta Psi"])
plt.title(f'Box Plot of Psi')
plt.ylabel("Delta Psi")
plt.grid(True)
plt.show()
