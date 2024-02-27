"""Computational chemistry management"""

from datetime import datetime
from pathlib import Path
import os
import sys

import silico.log
from silico.datas import get_resource

# Decide on whether we are frozen or not.
# The sys attribute 'frozen' is our flag, '_MEIPASS' is the dir location.
# https://pyinstaller.readthedocs.io/en/stable/runtime-information.html#run-time-information
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    frozen = True
else:
    frozen = False

# Setup openbabel library location.
init_obabel(frozen)

# Silico imports.
import silico.logging


# Version information.
major_version = 6
minor_version = 0
revision = 0
prerelease = 1
# Whether this is a development version.
development = prerelease is not None
demonstration = False
# The full version number of this package.
__version__ = "{}.{}.{}{}{}".format(major_version, minor_version, revision, "-pre.{}".format(prerelease) if development else "", "-demo" if demonstration else "")
# Deprecated:
version = __version__

# The blokes who wrote this.
__author__ = "The Silico Dev Team"

# Program date (when we were last updated). This is changed automatically.
_last_updated_string = "12/12/1234"
last_updated = datetime.strptime(_last_updated_string, "%d/%m/%Y")

if demonstration:
    # Pls no Hack.
    # Date at which the demo will end.
    demo_end_date = date = datetime.strptime("27/02/2024", "%d/%m/%Y")
    
    if datetime.now() > demo_end_date:
        print("The demonstration period has now expired, thank you for trying Digichem!\nFor continued usage, please contact @DigichemProject for a new license.")
        raise Exception("License Expired")

# The sys attribute 'frozen' is our flag, '_MEIPASS' is the dir location.
# https://pyinstaller.readthedocs.io/en/stable/runtime-information.html#run-time-information
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    frozen = True
else:
    frozen = False

import rdkit.RDLogger
# WrapLogs() outputs rdkit logging to python's stderr (which might be redirected to an urwid widget).
# If/when rdkit is further intergrated into silico, this call will likely be moved elsewhere. 
#rdkit.Chem.rdchem.WrapLogs()
# Sadly the behaviour of WrapLogs() is a bit bizzare, although we do get redirection to our custom widgets etc,
# logs are also still dumped to screen...
# for now, disable logging...
rdkit.RDLogger.DisableLog('rdApp.*')

# Setup the logger
silico.log.init_logger()


########################
# Function Definitions #
########################

def default_template_directory():
    return get_resource('data/templates')

# At end to avoid circular imports.
import silico.config
