#!/usr/bin/env python3
import sys
sys.dont_write_bytecode = True
from smb_disk import smb_disk, extract_files
from pathlib import Path as path

# File with info on how connect to our laboratory (remote) disk
disk_setups_fname = path.home()/"Documents"/"lab_disk_credentials.txt"
# Data folder on this disk
data_directory = path("Data")/"2023"/"231109"/"231109_caen_archive"
# Where extract the data to
destination = path.home()/"Documents"/"hdda"/"Data"/"231109"
# For renaming 'f1', 'f2', etc. to verbose folder names
zip_renaming = {"f41":"231109_1ph_LArN2_X-ray_14mm_coll_filt3_20kV_850V_46V",
                "f40":"231109_1ph_LArN2_X-ray_14mm_coll_filt3_18kV_850V_46V",
                "f39":"231109_1ph_LArN2_X-ray_14mm_coll_filt3_16kV_850V_46V",
                "f38":"231109_1ph_LArN2_X-ray_14mm_coll_filt3_14kV_850V_46V",
                "f37":"231109_1ph_LArN2_X-ray_14mm_coll_filt3_12kV_850V_46V",
                "f36":"231109_1ph_LArN2_X-ray_14mm_coll_filt3_10kV_850V_46V",
                "f35":"231109_1ph_LArN2_X-ray_14mm_coll_filt3_8kV_850V_46V",
                "f34":"231109_1ph_LArN2_X-ray_14mm_coll_filt3_0kV_850V_46V",
                "f33":"231109_1ph_LArN2_X-ray_14mm_coll_0kV_850V_46V",
                "f32":"231109_1ph_LArN2_X-ray_14mm_coll_8kV_850V_46V",
                "f31":"231109_1ph_LArN2_X-ray_14mm_coll_10kV_850V_46V",
                "f30":"231109_1ph_LArN2_X-ray_14mm_coll_12kV_850V_46V",
                "f29":"231109_1ph_LArN2_X-ray_14mm_coll_14kV_850V_46V",
                "f28":"231109_1ph_LArN2_X-ray_14mm_coll_16kV_850V_46V",
                "f27":"231109_1ph_LArN2_X-ray_14mm_coll_18kV_850V_46V",
                "f26":"231109_1ph_LArN2_X-ray_14mm_coll_20kV_850V_46V",
                "f25":"231109_1ph_LArN2_X-ray_6mm_coll_20kV_850V_46V_1",
                "f24":"231109_1ph_LArN2_X-ray_6mm_coll_20kV_850V_46V",
                "f23":"231109_1ph_LArN2_X-ray_6mm_coll_18kV_850V_46V",
                "f22":"231109_1ph_LArN2_X-ray_6mm_coll_16kV_850V_46V",
                "f21":"231109_1ph_LArN2_X-ray_6mm_coll_14kV_850V_46V",
                "f20":"231109_1ph_LArN2_X-ray_6mm_coll_12kV_850V_46V",
                "f19":"231109_1ph_LArN2_X-ray_6mm_coll_10kV_850V_46V",
                "f18":"231109_1ph_LArN2_X-ray_6mm_coll_8kV_850V_46V",
                }

def select_files(fname):
    return fname.endswith(".zip") and fname[:-4] in zip_renaming

# Get all .zip files in data_directory on remote disk
disk = smb_disk(disk_setups_fname)
if disk.ls(data_directory):
    local_files = disk.copy_all_files(destination, fname_selector=select_files)
    extract_files(local_files, destination, name_mapping=zip_renaming)

