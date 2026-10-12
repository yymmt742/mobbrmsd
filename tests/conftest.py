# Windows DLL pre-loading is handled by mobbrmsd/__init__.py using
# kernel32.LoadLibraryExW with flags=0, which resolves Fortran runtime DLLs
# via PATH and caches them for the restricted 0x1100 import flags used by
# Python 3.8+ (and ctypes in 3.12+).  No workaround is needed here.
