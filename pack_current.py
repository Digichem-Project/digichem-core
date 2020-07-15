#!/usr/bin/env python3

# Script for archiving the current source tree for distribution to various virtual machines etc.

import silico
from pathlib import Path
import subprocess
import distutils.dir_util
import shutil
from datetime import datetime

# The name of the archive.
archive_name = silico.name + "-" + silico.version + "-source"

# Where we're going to write to
dest = Path("../archives", archive_name)

# Check we're not overwriting an old archive.
if Path(str(dest) + ".tar.gz").exists():
	print("Archive with name {} already exists".format(dest))
	exit(-1)

# This is all just a quick hack for now...

# The date that the packaged version was last edited, which is right now.
edited_date = datetime.today()

# Add a new header.
#changelog = ["[{} ({})]\n".format(silico.version_number, edited_date.strftime("%d/%m/%Y"))]

# First get our new changes
#with open("updates") as updates_file:
#	changelog.extend(updates_file.readlines())
	
# Add a blank line (?)
#changelog.append("\n")
	
# Then read in the old changelog.
#with open("changelog") as changelog_file:
#	changelog.extend(changelog_file.readlines())
	
# Now write a new changelog with both.
#with open("changelog", "wt")  as changelog_file:
#	changelog_file.writelines(changelog)
	
# Remove the updates file.
#os.remove("updates")

# First copy our source to an appropriate folder name
distutils.dir_util.copy_tree("./", str(dest))

# Modify our version info to remove the development label.
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
		else:
			pass
			# And our dev label.
			#new_line = line.replace("development = True", "development = False")
		
		# Replace our line.
		version_file_contents[line_num] = new_line
		
	#version_file_contents = [line.replace("development = True", "development = False") for line in version_file_contents]
	
with open(version_file_path, "w") as version_file:
	# Write everything.
	version_file.writelines(version_file_contents)

# Use the tar program.
#subprocess.run(
#	["tar", "-czf", str(dest) + ".tar.gz", str(dest)]
#)
subprocess.run(
    ["tar", "-czf", archive_name + ".tar.gz", archive_name],
    cwd = str(dest.parent)
)


# Delete the folder.
shutil.rmtree(str(dest))

# # And create a new blank updates file.
# with open("updates", "wt") as updates_file:
# 	pass
