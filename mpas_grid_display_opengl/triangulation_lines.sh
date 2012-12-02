#!/bin/bash
#
g++ -c -I/$HOME/include triangulation_lines.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_lines.cpp"
  exit
fi
rm compiler.txt
#
#  Here is the load statement for Apple's OS X.
#
g++ triangulation_lines.o -lnetcdf -lnetcdf_c++ -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ triangulation_lines.o  -lnetcdf -lnetcdf_c++ -lm -lm -lGL -lGLU -lglut
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_lines.o"
  exit
fi
#
rm triangulation_lines.o
mv a.out ~/bincpp/$ARCH/triangulation_lines
#
echo "Executable installed as ~/bincpp/$ARCH/triangulation_lines"
