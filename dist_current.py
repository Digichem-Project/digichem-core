#!/usr/bin/env python3

# Script for archiving the current source tree for distribution to production environments.

from pathlib import Path
import subprocess
import distutils.dir_util
from datetime import datetime
import itertools
import shutil

import silico

# The name of the archive.
archive_name = silico.name + "-" + silico.version

# Where we're going to write to
dest = Path("../dist", archive_name)

# Check we're not overwriting an old archive.
if Path(str(dest) + ".tar.gz").exists():
    print("Archive with name {} already exists".format(dest))
    exit(-1)

# This is all just a quick hack for now...

# The date that the packaged version was last edited, which is right now.
edited_date = datetime.today()

# First copy our source to an appropriate folder name
distutils.dir_util.copy_tree("./", str(dest))

version_file_contents = ""
version_file_path = Path(dest, "silico/__init__.py")
with open(version_file_path, "r") as version_file:
    # Read everything.
    version_file_contents = version_file.readlines()
    
    # Modify those lines that we want to.
    for line_num, line in enumerate(version_file_contents):
        new_line = line
        
        # Replace our last modified string.
        if "_last_updated_string = " in line:
            new_line = '_last_updated_string = "{}"\n'.format(edited_date.strftime("%d/%m/%Y"))
        
        # Replace our line.
        version_file_contents[line_num] = new_line
        

with open(version_file_path, "w") as version_file:
    # Write everything.
    version_file.writelines(version_file_contents)

# Use git --archive to package
subprocess.run(
    ["git", "archive", "HEAD", "--prefix={}/".format(archive_name), "-o", "../{}.tar.gz".format(archive_name)],
    cwd = str(dest)
)

# Delete the folder.
shutil.rmtree(str(dest))
