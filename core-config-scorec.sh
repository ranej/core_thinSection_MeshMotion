flags=" -DDEBUG  -O3 -Wall -g -Wextra"

cmake .. \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_C_FLAGS="${flags}" \
  -DCMAKE_CXX_FLAGS="${flags}" \
  -DSIM_MPI="mpich3.3.2" \
  -DENABLE_ZOLTAN=ON \
  -DCMAKE_INSTALL_PREFIX=$PWD/install \
  -DENABLE_SIMMETRIX=True \
  -DENABLE_FIELDSIM=True \
  -DSIM_PARASOLID=True \
  -DPCU_COMPRESS=OFF \
  -DPUMI_FORTRAN_INTERFACE=ON \
  -DCMAKE_Fortran_COMPILER="mpif90" \
  -DCMAKE_BUILD_TYPE=Debug\
..
make -j8
make install
