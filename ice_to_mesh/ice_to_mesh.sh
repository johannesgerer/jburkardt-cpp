#!/bin/bash
#
g++ -c -I /usr/local/include ice_to_mesh.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ice_to_mesh.cpp"
  exit
fi
rm compiler.txt
#
g++ ice_to_mesh.o -L$HOME/lib/$ARCH -lnetcdf -lnetcdf_c++ -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ice_to_mesh.o."
  exit
fi
#
rm ice_to_mesh.o
#
chmod u+x a.out
mv a.out ~/bincpp/$ARCH/ice_to_mesh
#
echo "Executable installed as ~/bincpp/$ARCH/ice_to_mesh"
