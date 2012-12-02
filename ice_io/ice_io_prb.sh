#!/bin/bash
#
g++ -c -g -I /usr/local/include ice_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ice_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ ice_io_prb.o /$HOME/libcpp/$ARCH/ice_io.o -lnetcdf -lnetcdf_c++ -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ice_io_prb.o."
  exit
fi
#
rm ice_io_prb.o
#
mv a.out ice_io_prb
./ice_io_prb > ice_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ice_io_prb."
  exit
fi
rm ice_io_prb
#
echo "Program output written to ice_io_prb_output.txt"
