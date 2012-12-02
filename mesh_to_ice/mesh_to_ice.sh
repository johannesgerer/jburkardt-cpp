#!/bin/bash
#
g++ -c -g -I /usr/local/include mesh_to_ice.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling mesh_to_ice.cpp"
  exit
fi
rm compiler.txt
#
g++ mesh_to_ice.o -L/usr/local/lib -lnetcdf -lnetcdf_c++ -lm
if [ $? -ne 0 ]; then
  echo "Errors while loading mesh_to_ice.o"
  exit
fi
rm mesh_to_ice.o
#
mv a.out ~/bincpp/$ARCH/mesh_to_ice
#
echo "Executable installed as ~/bincpp/$ARCH/mesh_to_ice"
