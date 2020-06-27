# -*- mode: python ; coding: utf-8 -*-

from PyInstaller.utils.hooks import collect_data_files

block_cipher = None

# Weasprint data needs to go in the root for some reason. Thanks to https://stackoverflow.com/questions/29026827/pyinstaller-does-not-include-dependency-file for the hint.

weasyprint_datas = [(source, dest.split('weasyprint/')[1] if dest != "weasyprint" else ".") for source, dest in collect_data_files('weasyprint')]

datas = weasyprint_datas
datas.extend(collect_data_files('tinycss2'))
datas.extend(collect_data_files('cssselect2'))
datas.extend(collect_data_files('cairocffi'))
datas.extend(collect_data_files('pyphen'))

# Add silico data files.
datas.append(('../../silico/data', 'silico/data'))

# Now add extra binary libraries that we need.
binaries = [
	("/usr/local/lib/openbabel","openbabel"),
    ("/usr/local/lib/libinchi.so.0", ".")
]

a = Analysis(['../../silico/program/creport.py'],
             binaries=binaries,
             datas=datas,
			 hiddenimports=['cssselect2', 'tinycss2'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='creport',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=True )


import silico
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='creport.{}.CentOS-7.7'.format(silico.version) )
