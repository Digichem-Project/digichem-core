# -*- mode: python ; coding: utf-8 -*-

import sys
from pathlib import Path
sys.path.insert(0,str(Path("./").resolve()))

from base import binaries, datas

script = "../../silico/program/_silico_resume.py"
prog_name = "_silico_resume"
package_name = "CentOS-6.10"

a = Analysis([script],
			 binaries=binaries,
			 datas=datas,
			 # 'pkg_resources.py2_warn' see https://github.com/pypa/setuptools/issues/1963
			 hiddenimports=['cssselect2', 'tinycss2', 'pkg_resources.py2_warn'],
			 hookspath=[],
			 runtime_hooks=[],
			 excludes=[],
			 win_no_prefer_redirects=False,
			 win_private_assemblies=False,
			 cipher=None,
			 noarchive=False
)

pyz = PYZ(a.pure, a.zipped_data,
			 cipher=None
)

exe = EXE(pyz,
		a.scripts,
		[],
		exclude_binaries=True,
		name=prog_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True
)

import silico
coll = COLLECT(exe,
				a.binaries,
				a.zipfiles,
				a.datas,
				strip=False,
				upx=True,
				upx_exclude=[],
				name="{}.{}.{}".format(prog_name, silico.version, package_name)
)
