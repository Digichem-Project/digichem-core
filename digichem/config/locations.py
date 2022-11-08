from pathlib import Path
import pkg_resources
import silico


# A config file included in the silico install directory.
master_config_path = Path(pkg_resources.resource_filename('silico', "data/config/silico.yaml"))

# The location of the system-wide config file (which has precedence over the master).
system_config_location = Path("/etc/silico{}/silico.yaml".format(silico.major_version))

# The location of the user specific config file (which has precedence over the system).
user_config_location = Path(Path.home(), ".config/silico{}/silico.yaml".format(silico.major_version))