#!/usr/bin/env python3
import sys
sys.dont_write_bytecode = True
from smb_disk import smb_disk, extract_files
from pathlib import Path as path

# File with info on how connect to our laboratory (remote) disk
disk_setups_fname = path.home()/"Documents"/"lab_disk_credentials.txt"
# Data folder on this disk
data_directory = path("Data")/"2022"/"221124"/"221124_caen_archive"
# Where extract the data to
destination = path.home()/"Documents"/"hdda"/"Data"/"221124"
# For renaming 'f1', 'f2', etc. to verbose folder names
zip_renaming = {}

def select_files(fname):
    return fname.endswith(".zip") and not "_Q_" in fname

# Get all .zip files in data_directory on remote disk
disk = smb_disk(disk_setups_fname)
if disk.ls(data_directory):
    local_files = disk.copy_all_files(destination, fname_selector=select_files)
    extract_files(local_files, destination, name_mapping=zip_renaming)

