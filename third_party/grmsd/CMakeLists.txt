cmake_minimum_required(VERSION 3.14)
project(grmsd_bench NONE)

include(ExternalProject)

# X-GRMSD uses -O3; requires GCC or Clang.
# On Windows, use MinGW (same toolchain as the Fortran core).

ExternalProject_Add(
    xgrmsd
    GIT_REPOSITORY https://github.com/striwata/X-GRMSD
    GIT_TAG        main
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
    CMAKE_GENERATOR "${CMAKE_GENERATOR}"
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
        # X-GRMSD uses cmake_minimum_required(VERSION 3.1), which is below the
        # compatibility floor introduced in cmake 4.x.
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5
    # Build only the GRMSD target (skip Substructure* targets)
    BUILD_COMMAND
        ${CMAKE_COMMAND} --build <BINARY_DIR> --target GRMSD --config Release
    # Copy the finished binary into the top-level build directory so that
    # gen_benchmark.py can find it with --cmake-build-dir.
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy
            <BINARY_DIR>/GRMSD${CMAKE_EXECUTABLE_SUFFIX}
            ${CMAKE_BINARY_DIR}/GRMSD${CMAKE_EXECUTABLE_SUFFIX}
)

message(STATUS
    "X-GRMSD: binary will be installed to "
    "${CMAKE_BINARY_DIR}/GRMSD${CMAKE_EXECUTABLE_SUFFIX}"
)
