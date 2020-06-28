# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

# Add silico data files.
datas = [('../../silico/data', 'silico/data')]

# Now add extra binary libraries that we need.
binaries = [
	("/usr/local/lib/openbabel","openbabel"),
    ("/usr/local/lib/libinchi.so.0", ".")
]

a = Analysis(['../../silico/program/cresult.py'],
             binaries=binaries,
             datas=datas,
             hiddenimports=[],
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
          name='cresult',
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
               name='cresult.{}.CentOS-6.10'.format(silico.version) )
