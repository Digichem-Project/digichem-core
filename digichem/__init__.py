"""Computational chemistry management"""

from datetime import datetime
from pathlib import Path
import pkg_resources
import os
import sys

import silico.logging


####################
# Package metadata.#
####################

# Version information.
major_version = 1
minor_version = 0
revision = 0
prerelease = 30
# Whether this is a development version.
development = False
# Version information.
major_version = 0
minor_version = 20
revision = 5
version_number = "{}.{}.{}".format(major_version, minor_version, revision)
# The full version number of this package.
__version__ = "{}.{}.{}{}".format(major_version, minor_version, revision, "-pre.{}".format(prerelease) if development else "")
# Deprecated:
version = __version__

# The blokes who wrote this.
__author__ = "The Silico Dev Team"

# Program date (when we were last updated). This is changed automatically.
_last_updated_string = "12/12/1234"
last_updated = datetime.strptime(_last_updated_string, "%d/%m/%Y")


# Set-up openbabel.
#
# This must be done prior to an obabel import.
#     
# When frozen with pyinstaller we take a version of the openbabel C library with us (along with the relevant python bindings of course).
# This library is split into several .so files corresponding to the various formats obabel supports, and while the main libopenbabel.so file is found automatically, these supplementary library files are not.
# So, when we are frozen, we manually set the location of these library files so openbabel will work.
# If we are not frozen we do not do this as we expect openbabel to be correctly configured.
# The sys attribute 'frozen' is our flag, '_MEIPASS' is the dir location.
# https://pyinstaller.readthedocs.io/en/stable/runtime-information.html#run-time-information
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    frozen = True
else:
    frozen = False

openbabel_version = "3.0.0"
# The sys attribute 'frozen' is our flag, '_MEIPASS' is the dir location.
# https://pyinstaller.readthedocs.io/en/stable/runtime-information.html#run-time-information
if silico.frozen:
    # We need to tell openbabel where its library components are.
    os.environ['BABEL_LIBDIR'] = str(Path(sys._MEIPASS, "openbabel", "lib", openbabel_version))
    
    # And also data.
    os.environ['BABEL_DATADIR'] = str(Path(sys._MEIPASS, "openbabel", "data", openbabel_version))


########################
# Function Definitions #
########################

def default_template_directory():
    return Path(pkg_resources.resource_filename('silico', 'data/templates'))