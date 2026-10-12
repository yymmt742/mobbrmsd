import glob
import os
import sys

if sys.platform == "win32":
    # Python 3.8+ on Windows does not use PATH for DLL loading.
    # Add directories containing gfortran runtime DLLs explicitly.
    added: set[str] = set()

    for path_dir in os.environ.get("PATH", "").split(os.pathsep):
        if path_dir in added or not os.path.isdir(path_dir):
            continue
        if glob.glob(os.path.join(path_dir, "libgfortran*.dll")):
            try:
                os.add_dll_directory(path_dir)
                added.add(path_dir)
            except (OSError, ValueError):
                pass

    # Also add the package directory for DLLs bundled alongside .pyd files.
    import importlib.util

    spec = importlib.util.find_spec("mobbrmsd")
    if spec is not None:
        pkg_dir = os.path.dirname(spec.origin)
        if pkg_dir not in added:
            try:
                os.add_dll_directory(pkg_dir)
            except (OSError, ValueError):
                pass
