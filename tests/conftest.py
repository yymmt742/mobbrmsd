# Windows DLL pre-loading is handled by mobbrmsd/__init__.py using
# kernel32.LoadLibraryExW with flags=0, which resolves Fortran runtime DLLs
# via PATH and caches them for the restricted 0x1100 import flags used by
# Python 3.8+ (and ctypes in 3.12+).  No workaround is needed here.


def pytest_addoption(parser):
    parser.addoption(
        "--run-perf",
        action="store_true",
        default=False,
        help="Run performance/accuracy comparison tests (slow, requires RDKit and/or GRMSD binary).",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "perf: external-library performance and accuracy benchmarks (skipped unless --run-perf is passed).",
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-perf"):
        return
    skip_perf = __import__("pytest").mark.skip(reason="pass --run-perf to run performance tests")
    for item in items:
        if "perf" in item.keywords:
            item.add_marker(skip_perf)
