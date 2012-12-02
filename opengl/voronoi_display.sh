#!/bin/bash
#
g++ -c voronoi_display.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling voronoi_display.cpp"
  exit
fi
rm compiler.txt
#
#
#  Here is the load statement for Apple's OS X.
#
g++ voronoi_display.o -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ voronoi_display.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Errors linking and loading voronoi_display.o"
  exit
fi
#
rm voronoi_display.o
mv a.out ~/bincpp/$ARCH/voronoi_display
#
echo "Executable installed as ~/bincpp/$ARCH/voronoi_display"
