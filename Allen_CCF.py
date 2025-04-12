
from pathlib import Path
import os
import numpy as np

from allensdk.core.reference_space_cache import ReferenceSpaceCache

# output_dir = '.'
# 
# reference_space_key = os.path.join('annotation', 'ccf_2017')
# resolution = 10
# rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest=Path(output_dir) / 'manifest.json')
# annotation, meta = rspc.get_annotation_volume()
# os.listdir(Path(output_dir) / reference_space_key)
# 
# rsp = rspc.get_reference_space()
# 
# # A complete mask for one structure
# # ... id goes all > nose > barrel > lower limb > mouth > upper limb > trunk > unassigned
# ROI_mask_S1_L1 = rsp.make_structure_mask([793, 558, 981, 1030, 878, 450, 1006, 182305693])
# ROI_mask_S1_L1 = np.argwhere(ROI_mask_S1_L1 == 1)
# ROI_mask_S1_L23 = rsp.make_structure_mask([346, 838, 201, 113, 657, 854, 670, 182305697])
# ROI_mask_S1_L23 = np.argwhere(ROI_mask_S1_L23 == 1)
# ROI_mask_S1_L4 = rsp.make_structure_mask([865, 654, 1047, 1094, 950, 577, 1086, 182305701])
# ROI_mask_S1_L4 = np.argwhere(ROI_mask_S1_L4 == 1)
# ROI_mask_S1_L5 = rsp.make_structure_mask([921, 702, 1070, 1128, 974, 625, 1111, 182305705])
# ROI_mask_S1_L5 = np.argwhere(ROI_mask_S1_L5 == 1)
# ROI_mask_S1_L6a = rsp.make_structure_mask([686, 889, 1038, 478, 1102, 945, 9, 182305709])
# ROI_mask_S1_L6a = np.argwhere(ROI_mask_S1_L6a == 1)
# ROI_mask_S1_L6b = rsp.make_structure_mask([719, 929, 1062, 510, 2, 1026, 461, 182305713])
# ROI_mask_S1_L6b = np.argwhere(ROI_mask_S1_L6b == 1)

def generate_roi_masks(output_dir = '.', reference_space_key = 'annotation/ccf_2017', resolution = 10):
    # Initialize the ReferenceSpaceCache
    rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest=Path(output_dir) / 'manifest.json')
    
    # Get the reference space
    rsp = rspc.get_reference_space()
    
    # Define the structure mask function
    def make_mask(structure_ids):
        mask = rsp.make_structure_mask(structure_ids)
        return np.argwhere(mask == 1)
    
    # Structure IDs and corresponding mask names
    structure_info = {
        "ROI_mask_S1_L1": [793, 558, 981, 1030, 878, 450, 1006, 182305693],
        "ROI_mask_S1_L23": [346, 838, 201, 113, 657, 854, 670, 182305697],
        "ROI_mask_S1_L4": [865, 654, 1047, 1094, 950, 577, 1086, 182305701],
        "ROI_mask_S1_L5": [921, 702, 1070, 1128, 974, 625, 1111, 182305705],
        "ROI_mask_S1_L6a": [686, 889, 1038, 478, 1102, 945, 9, 182305709],
        "ROI_mask_S1_L6b": [719, 929, 1062, 510, 2, 1026, 461, 182305713]
    }
    
    # Generate masks for each structure
    masks = {}
    for mask_name, structure_ids in structure_info.items():
        masks[mask_name] = make_mask(structure_ids)
    
    # Return the masks as a dictionary
    return masks


