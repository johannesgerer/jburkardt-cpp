#!/bin/bash
#
g++ -c -g -I$HOME/include pres_temp_4D_rd.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pres_temp_4D_rd.cpp"
  exit
fi
rm compiler.txt
#
g++ pres_temp_4D_rd.o -L$HOME/lib/$ARCH -lnetcdf_c++ -lnetcdf -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pres_temp_4D_rd.o."
  exit
fi
#
rm pres_temp_4D_rd.o
#
mv a.out pres_temp_4D_rd
./pres_temp_4D_rd > pres_temp_4D_rd_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pres_temp_4D_rd."
  exit
fi
rm pres_temp_4D_rd
#
echo "Program output written to pres_temp_4D_rd_output.txt"
