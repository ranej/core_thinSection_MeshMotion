flags=" -O3 -DNDEBUG -Wall -Wextra"

cmake .. \
   -DCMAKE_C_COMPILER="mpicc" \
   -DCMAKE_CXX_COMPILER="mpicxx" \
   -DCMAKE_C_FLAGS="${flags}"  \
   -DCMAKE_CXX_FLAGS="${flags}" \
   -DSIM_MPI=openmpi4.0.1 \
   -DENABLE_ZOLTAN=ON \
   -DCMAKE_INSTALL_PREFIX=$PWD/install  \
   -DENABLE_SIMMETRIX=True \
   -DENABLE_FIELDSIM=True \
   -DSIM_PARASOLID=True \
   -DPCU_COMPRESS=OFF \
   -DPUMI_FORTRAN_INTERFACE=ON \
   -DCMAKE_Fortran_COMPILER="mpif90" \
..
make -j24
make install










