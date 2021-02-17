# Code for result extraction
from .base import Result_extractor
from .base import Result_extractor_group


# These are the class handles for the various section extractors.
SUMMARY_CLASS_HANDLE = ["summary", "all"]    
PDM_CLASS_HANDLE = ["PDM", "permanent dipole moment", "dipole moment"]
TDM_CLASS_HANDLE = ["TDM", "transition dipole moment", "S1 dipole moment"]
SCF_CLASS_HANDLE = ["SCF", "SCF energy", "energy"]
MP_CLASS_HANDLE = ["MP", "MP energy"]
CC_CLASS_HANDLE = ["CC", "CC energy"]
EXCITED_STATE_CLASS_HANDLE = ["ES", "excited", "states", "excited states"]
SOC_CLASS_HANDLE = ["SOC", "spin orbit coupling", "spin-orbit coupling"]
EXCITED_STATE_TRANSITIONS_CLASS_HANDLE = ["TRAN", "transitions", "ES transitions", "excited state transitions"]
GEOM_CLASS_HANDLE = ["GEOM", "geometry", "atoms"]
METADATA_CLASS_HANDLE = ["META", "metadata"]
ORBITALS_CLASS_HANDLE = ["ORB", "ALPHA", "molecular orbitals", "orbitals"]
BETA_CLASS_HANDLE = ["BETA", "beta orbitals"]
VIBRATIONS_CLASS_HANDLE = ["VIB", "FREQ", "vibrations", "frequencies", "vibrational frequencies"]