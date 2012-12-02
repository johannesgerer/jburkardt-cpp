#!/bin/bash
#
cp netcdf_mpas.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include netcdf_mpas.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling netcdf_mpas.cpp"
  exit
fi
rm compiler.txt
#
mv netcdf_mpas.o ~/libcpp/$ARCH/netcdf_mpas.o
#
echo "Library installed as ~/libcpp/$ARCH/netcdf_mpas.o"
