# -*- mode: python ; coding: utf-8 -*-

from PyInstaller.utils.hooks import collect_data_files


# Weasprint data needs to go in the root for some reason. Thanks to https://stackoverflow.com/questions/29026827/pyinstaller-does-not-include-dependency-file for the hint.
weasyprint_datas = [(source, dest.split('weasyprint/')[1] if dest != "weasyprint" else ".") for source, dest in collect_data_files('weasyprint')]

datas = weasyprint_datas
datas.extend(collect_data_files('tinycss2'))
datas.extend(collect_data_files('cssselect2'))
datas.extend(collect_data_files('cairocffi'))
datas.extend(collect_data_files('pyphen'))
datas.extend(collect_data_files('pysoc'))
datas.extend(collect_data_files('basis_set_exchange'))

# Add silico data files.
datas.append(('../../silico/data', 'silico/data'))

# Openbabel data too.
datas.append(("/usr/local/share/openbabel/3.1.1/", "openbabel/data/3.1.1"))

# Now add extra binary libraries that we need.
binaries = [
	("/usr/local/lib/openbabel/3.1.1","openbabel/lib/3.1.1"),
    ("/usr/local/lib/libinchi.so.0", "."),
	("/lib64/libz.so.1", ".")
]
