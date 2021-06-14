# General imports.
from datetime import datetime
from pathlib import Path
import pkg_resources
import os
import sys
import PIL.Image

# Silico imports.
import silico.logging
from .base import init_obabel


PIL.Image.MAX_IMAGE_PIXELS = None

# The name of this package.
name = "silico"
# Brief description.
description = "Silico Computational Chemistry Package"
# Version information.
major_version = 0
minor_version = 19
revision = 5
version_number = "{}.{}.{}".format(major_version, minor_version, revision)
# The full version number of this package.
version = "{}.{}.{}{}".format(major_version, minor_version, revision, "-pre.{}".format(prerelease) if development else "")
# The bloke who wrote this.
author = "Oliver Lee"
# Program date (when we were last updated).
_last_updated_string = "11/02/2022"
last_updated = datetime.strptime(_last_updated_string, "%d/%m/%Y")

# The name of the command to launch new instances of silico.
silico_cmd = "silico" if not development else "silico-dev"


# Decide on whether we are frozen or not.
# The sys attribute 'frozen' is our flag, '_MEIPASS' is the dir location.
# https://pyinstaller.readthedocs.io/en/stable/runtime-information.html#run-time-information
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    frozen = True
else:
    frozen = False
    

########################
# Function Definitions #
########################

def default_template_directory():
    return Path(pkg_resources.resource_filename('silico', 'data/templates'))
