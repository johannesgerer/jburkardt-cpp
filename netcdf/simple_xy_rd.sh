#!/bin/bash
#
g++ -c -g -I$HOME/include simple_xy_rd.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simple_xy_rd.cpp"
  exit
fi
rm compiler.txt
#
g++ simple_xy_rd.o -L$HOME/lib/$ARCH -lnetcdf_c++ -lnetcdf -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading simple_xy_rd.o."
  exit
fi
#
rm simple_xy_rd.o
#
mv a.out simple_xy_rd
./simple_xy_rd > simple_xy_rd_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running simple_xy_rd."
  exit
fi
rm simple_xy_rd
#
echo "Program output written to simple_xy_rd_output.txt"
