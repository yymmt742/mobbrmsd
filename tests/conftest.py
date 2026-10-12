import os
import sys

if sys.platform == "win32":
    import importlib.util
    spec = importlib.util.find_spec("mobbrmsd")
    if spec is not None:
        pkg_dir = os.path.dirname(spec.origin)
        os.add_dll_directory(pkg_dir)
    tdm_gcc = r"C:\TDM-GCC-64\bin"
    if os.path.isdir(tdm_gcc):
        os.add_dll_directory(tdm_gcc)
