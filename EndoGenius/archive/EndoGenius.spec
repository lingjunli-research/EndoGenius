# -*- mode: python ; coding: utf-8 -*-

import sys ; sys.setrecursionlimit(sys.getrecursionlimit() * 5)

a = Analysis(
    ['GUI.py'],
    pathex=[],
    binaries=[(r'C:\Users\lawashburn\Documents\EndoGeniusDistributions\EndoGenius_v1.1.0\DIA-NN\1.8.1\DIA-NN.exe', '.')],
    datas=[("venv/Lib/site-packages/pyopenms", "./pyopenms/"),('Endogenius_Thumbnail.ico', '.'),('assets', 'assets')],
    hiddenimports=['pkg_resources.py2_warn', 'pkg_resources.extern'],
    hookspath=['C:\\Users\\lawashburn\\Documents\\EndoGeniusDistributions\\EndoGenius_v1.1.0\\hooks'],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['jupyter', 'notebook'],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='EndoGenius',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['EndoGenius_Thumbnail.ico'],
)
