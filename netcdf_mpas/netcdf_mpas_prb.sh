#!/bin/bash
#
g++ -c -g -I/$HOME/include netcdf_mpas_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling netcdf_mpas_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ netcdf_mpas_prb.o /$HOME/libcpp/$ARCH/netcdf_mpas.o -lnetcdf -lnetcdf_c++ -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading netcdf_mpas_prb.o."
  exit
fi
#
rm netcdf_mpas_prb.o
#
mv a.out netcdf_mpas_prb
./netcdf_mpas_prb > netcdf_mpas_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running netcdf_mpas_prb."
  exit
fi
rm netcdf_mpas_prb
#
echo "Program output written to netcdf_mpas_prb_output.txt"
