find_package(BLAS REQUIRED)

if(BLAS_FOUND)
  add_compile_options(${BLAS_LINKER_FLAGS})
endif()

include(FetchContent)
FetchContent_Declare(
  forbar
  GIT_REPOSITORY https://github.com/yymmt742/forbar
  GIT_TAG origin/main
  GIT_PROGRESS TRUE
  )
FetchContent_MakeAvailable(forbar)

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  FetchContent_Declare(
    unittest
    GIT_REPOSITORY https://github.com/yymmt742/unittest
    GIT_TAG origin/main
    GIT_PROGRESS TRUE
    )
  FetchContent_MakeAvailable(unittest)
endif()

add_subdirectory(mobbrmsd_2d)
add_subdirectory(mobbrmsd_3d)
add_subdirectory(mobbrmsd)

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_subdirectory(test)
endif()
