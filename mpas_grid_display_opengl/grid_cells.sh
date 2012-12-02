#!/bin/bash
#
g++ -c -I/$HOME/include grid_cells.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grid_cells.cpp"
  exit
fi
rm compiler.txt
#
#  Here is the load statement for Apple's OS X.
#
g++ grid_cells.o -lnetcdf -lnetcdf_c++ -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ grid_cells.o -lnetcdf -lnetcdf_c++ -lm -lm -lGL -lGLU -lglut
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading grid_cells.o"
  exit
fi
#
rm grid_cells.o
mv a.out ~/bincpp/$ARCH/grid_cells
#
echo "Executable installed as ~/bincpp/$ARCH/grid_cells"
