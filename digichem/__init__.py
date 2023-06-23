"""Computational chemistry management"""

from datetime import datetime
from pathlib import Path
import pkg_resources
import os
import sys

# Deal with obabel early.
from .base import init_obabel

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
major_version = 3
minor_version = 1
revision = 0
prerelease = 2
# Whether this is a development version.
development = False
# Version information.
major_version = 0
minor_version = 20
revision = 6
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
    
    # And also data.
    os.environ['BABEL_DATADIR'] = str(Path(sys._MEIPASS, "openbabel", "data", openbabel_version))
    

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
    return Path(pkg_resources.resource_filename('silico', 'data/templates'))

# At end to avoid circular imports.
import silico.config
