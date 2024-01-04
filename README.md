# symRMSD
molecular symmetry corrected RMSD by branch and prune

installation
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=[Release, Debug, ...] -DCMAKE_PREFIX_PATH=your_local_cmake_path
make install

It doesn't work with MKL. Please use
```
-DBLA_VENDOR=[OpenBLAS, ATLAS, ...]
```
instead.
