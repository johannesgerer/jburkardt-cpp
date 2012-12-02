#!/bin/bash
#
g++ -c -g -I$HOME/include sfc_pres_temp_wr.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sfc_pres_temp_wr.cpp"
  exit
fi
rm compiler.txt
#
g++ sfc_pres_temp_wr.o -L$HOME/lib/$ARCH -lnetcdf_c++ -lnetcdf -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sfc_pres_temp_wr.o."
  exit
fi
#
rm sfc_pres_temp_wr.o
#
mv a.out sfc_pres_temp_wr
./sfc_pres_temp_wr > sfc_pres_temp_wr_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sfc_pres_temp_wr."
  exit
fi
rm sfc_pres_temp_wr
#
echo "Program output written to sfc_pres_temp_wr_output.txt"
