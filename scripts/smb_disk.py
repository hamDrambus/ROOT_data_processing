#!/usr/bin/env python3
import sys
sys.dont_write_bytecode = True
import shutil
from pathlib import Path as path
from smb.SMBConnection import SMBConnection
from smb import smb_structs
import socket
import zipfile

class smb_disk:
    def __init__(self, credentials_fname):
        credentials = {}
        self.is_connected = False
        self.current_dir = path()
        with open(credentials_fname) as myfile:
            for line in myfile:
                try:
                    line = line.strip()
                    if not line:
                        continue
                    name, var = tuple(line.split("=", maxsplit=1))
                    credentials[name.strip()] = var.strip()[1:-1]
                except Exception as e:
                    continue 
        self.username = credentials["user"]
        self.password = credentials["password"]
        self.server_ip = credentials["ip"]
        self.shared_name = credentials["shared_name"]
        # Port is either 445 or 139
        self.port = int(credentials["port"]) if "port" in credentials else 445
        try:
            self.hostname = socket.gethostbyaddr(self.server_ip)
            self.hostname = self.hostname[0]
        except socket.herror as e:
            print(f"Uknknown host for {self.server_ip}")
            self.hostname = None
        if self.hostname:
            self.connection = SMBConnection(self.username, self.password, "python-smb-remote", self.hostname, use_ntlm_v2=True)
            self.is_connected = self.connection.connect(self.server_ip, self.port)
        if not self.is_connected:
            print(f"Error: connection to smb://{self.user}@{self.server_ip}:{self.port}/{self.shared_name} has failed.")

    def ls(self, directory):
        if not self.is_connected:
            return False
        try:
            attr = self.connection.getAttributes(self.shared_name, directory.as_posix())
        except Exception as e:
            print("ls:", e)
            return False
        if not attr.isDirectory:
            print(f"ls: cannot access smb://'{self.server_ip}/{self.shared_name}/{directory}': No such directory")
            return False
        self.current_dir = directory
        return True

    def pwd(self):
        return self.current_dir

    def copy_files(self, file_list, dest):
        """ Copy file_list files from remote disk to dest directory.

        Returns the list of successfuly copied files.
        """
        copied_files = []
        if not self.is_connected:
            return copied_files
        for fname in file_list:
            dest_fname = dest / fname.name
            dest_fname.parent.mkdir(parents=True, exist_ok=True)
            with open(dest_fname, "wb") as file:
                self.connection.retrieveFile(self.shared_name, fname.as_posix(), file)
                print(f"Retrieved \"{fname.name}\"")
                copied_files.append(dest_fname)
        return copied_files
    
    def copy_all_files(self, dest,
                       search_attr=smb_structs.SMB_FILE_ATTRIBUTE_ARCHIVE,
                       fname_selector=lambda x: x.endswith("*.zip")):
        """ Copy all files in the current directory on the remote disk to dest directory.
        Files can be filteted first by search_attr which defines which files are found on the smb disk.
        Additionally files can be filtered by fname_selector.

        Returns the list of successfuly copied files.
        """
        if not self.is_connected:
            return []
        files = self.connection.listPath(self.shared_name, self.current_dir.as_posix(), search=search_attr)
        files = [self.current_dir / f.filename for f in files if fname_selector(f.filename)]
        files = self.copy_files(files, dest)
        return files


def extract_files(file_list, dest, name_mapping={}):
    for fname in file_list:
        with zipfile.ZipFile(fname, 'r') as zip_ref:
            name = fname.name
            dest_name = name_mapping[name] if name in name_mapping else name
            if dest_name:
                dest_dir = (dest/dest_name).with_suffix('')
                shutil.rmtree(dest_dir, ignore_errors=True)
                dest_dir.mkdir(parents=True, exist_ok=True)
                zip_ref.extractall(dest_dir)
                print(f"Extracted \"{name}\" to \"{dest_dir}\"")
        fname.unlink()
