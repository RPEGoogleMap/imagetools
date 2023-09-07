@echo off
python -m pip uninstall -y imagetools
rmdir build /s /q 2>nul
rmdir __pycache__ /s /q 2>nul
rmdir imagetools.egg-info /s /q 2>nul
del imagetools.py 2>nul
del _imagetools* 2>nul
del imagetools_wrap.cpp 2>nul
del laminb.py 2>nul
del _laminb* 2>nul
del laminb_wrap.cpp 2>nul
del mito.py 2>nul
del _mito* 2>nul
del mito_wrap.cpp 2>nul
