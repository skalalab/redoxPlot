
import os 
import shutil
import re
import cv2
import numpy as np
from skimage import measure
rt = "/home/useradmin/Downloads/tmre_redox"

def tmreAnalysis(cell_line, treatment, roi):
    current_dir = os.path.join(cell_line, treatment, str(roi))
    mitoI_pattern = re.compile(r'^(?!\.).*_MitoI_(\d+)\.tiff$')
    mitoIs = []
    treatment = "FCCP" if treatment == "A1" else "control" if treatment == "A2" else "Oligo"
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
        # try find the corresponding nuclei mask
        nuclei_mask_path =  os.path.join(current_dir, mitoI.replace("MitoI", "nucleimask"))
        nuclei_mask_image = cv2.imread(nuclei_mask_path, cv2.IMREAD_UNCHANGED)
        cell_mask_path = os.path.join(current_dir, mitoI.replace("MitoI", "cellmask"))
        cell_mask_image = cv2.imread(cell_mask_path, cv2.IMREAD_UNCHANGED)
        mitoMask_path = os.path.join(current_dir, mitoI.replace("MitoI", "MitoMask"))
        mitoMask_image = cv2.imread(mitoMask_path, cv2.IMREAD_UNCHANGED)
        if nuclei_mask_image is None or cell_mask_image is None or mitoMask_image is None:
            continue 
        nuclei_binary_mask = (nuclei_mask_image > 0).astype(int)
        cell_binary_mask = (cell_mask_image > 0).astype(int)
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
        cell_area = np.sum(cell_binary_mask)
        mito_content = mito_area / cell_area
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
                   nuclei_intensity_sum, nuclei_intensity_avg, delta_psi_1, cell_area, mito_content, cell_line, treatment]
                
            rows.append(row)
    return rows

import pandas as pd
from datetime import date
current_date = date.today()

cell_lines = ["A549", "HeLa", "HFF", "HPDE", "IMR90", "MCF10A", "MCF7", "Panc1"]
treatment_dir = ["A1", "A2", "A3"]
data_rows = []
headers = ['file_name', 'mito_area', 'mito_intensity_sum', 'mito_intensity_avg', "nuclei_area", \
        "nuclei_intensity_sum", "nuclei_intensity_avg", "Delta Psi 1", "cell_area", "mito_content", "cell_line", "treatment"]

for cell_line in cell_lines:
    for treatment in treatment_dir:
        treatment_path = os.path.join(rt, cell_line, treatment)
        if os.path.exists(treatment_path):
            roi_dirs = [d for d in os.listdir(treatment_path) if os.path.isdir(os.path.join(treatment_path, d)) and d.isdigit()]
            for roi in roi_dirs:
                cells = tmreAnalysis(cell_line, treatment, roi)
                data_rows += cells

csvPath = os.path.join(rt, f"tmre_analysis.csv")
if len(data_rows) > 0: 
    df = pd.DataFrame(data_rows, columns=headers)
    df.to_csv(csvPath, index=False)