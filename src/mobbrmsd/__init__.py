import os as _os
import sys as _sys

if _sys.platform == "win32":
    # Python 3.8+ loads extension modules with LOAD_LIBRARY_SEARCH_DEFAULT_DIRS |
    # LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR (0x1100), which skips PATH.  Python 3.12+
    # also uses these flags in ctypes.CDLL.  Fortran runtime DLLs bundled in the
    # package directory depend on MinGW runtime DLLs (libgcc_s_seh, libwinpthread)
    # that live on PATH but not in the package directory.
    #
    # Fix: use kernel32.LoadLibraryExW with flags=0 (the unrestricted PATH-based
    # search) to pre-load each bundled DLL.  That pulls the full transitive
    # dependency chain into the process's loaded-module cache.  When Python later
    # imports a .pyd with the restricted 0x1100 flags, the already-cached DLLs
    # are satisfied without a new filesystem search.  The add_dll_directory call
    # ensures the package directory itself is in the restricted search set.
    import ctypes as _ctypes
    _k32 = _ctypes.WinDLL("kernel32", use_last_error=True)
    _k32.LoadLibraryExW.restype = _ctypes.c_void_p
    _k32.LoadLibraryExW.argtypes = [_ctypes.c_wchar_p, _ctypes.c_void_p, _ctypes.c_uint32]
    _dll_handles = []  # keep handles alive to prevent FreeLibrary
    _dll_dirs = []     # keep add_dll_directory handles alive
    for _p in __path__:
        try:
            _dlls = [f for f in _os.listdir(_p) if f.lower().endswith(".dll")]
        except OSError:
            continue
        if not _dlls:
            continue
        for _f in _dlls:
            _h = _k32.LoadLibraryExW(_os.path.join(_p, _f), None, 0)
            if _h:
                _dll_handles.append(_ctypes.c_void_p(_h))
        try:
            _dll_dirs.append(_os.add_dll_directory(_p))
        except OSError:
            pass

from .dataclass import molecules, molecular_system
from .mobbrmsd import mobbrmsd
from .demo._coord_generator import coord_generator
from importlib.metadata import version

__version__ = version(__package__)
