#!/bin/bash
#
g++ -c -g -I$HOME/include simple_xy_wr.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simple_xy_wr.cpp"
  exit
fi
rm compiler.txt
#
g++ simple_xy_wr.o -L$HOME/lib/$ARCH -lnetcdf_c++ -lnetcdf -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading simple_xy_wr.o."
  exit
fi
#
rm simple_xy_wr.o
#
mv a.out simple_xy_wr
./simple_xy_wr > simple_xy_wr_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running simple_xy_wr."
  exit
fi
rm simple_xy_wr
#
echo "Program output written to simple_xy_wr_output.txt"
