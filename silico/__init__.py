from datetime import datetime
from pathlib import Path
import pkg_resources
import os
import sys
import silico.base
import PIL.Image

from .base import init_obabel

PIL.Image.MAX_IMAGE_PIXELS = None

# The name of this package.
name = "silico"
# Brief description.
description = "Silico Computational Chemistry Package"
# Whether this is a development version.
development = True
# Version information.
major_version = 0
minor_version = 18
revision = 7
version_number = "{}.{}.{}".format(major_version, minor_version, revision)
# The full version number of this package.
version = "{}{}".format(version_number, "-dev" if development else "")
# The bloke who wrote this.
author = "Oliver Lee"
# Program date (when we were last updated).
_last_updated_string = "11/02/2021"
last_updated = datetime.strptime(_last_updated_string, "%d/%m/%Y")

# The name of the command to launch new instances of silico.
silico_cmd = "silico" if not development else "silico-dev"

# The name of the logger that will be used by silico.
logger_name = "silico"

# Init the logger.
silico.base.init_logger(logger_name)

# Decide on whether we are frozen or not.
# The sys attribute 'frozen' is our flag, '_MEIPASS' is the dir location.
# https://pyinstaller.readthedocs.io/en/stable/runtime-information.html#run-time-information
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    frozen = True
else:
    frozen = False

def default_template_directory():
    return Path(pkg_resources.resource_filename('silico', 'data/templates'))
