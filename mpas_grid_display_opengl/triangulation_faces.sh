#!/bin/bash
#
g++ -c -I/$HOME/include triangulation_faces.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_faces.cpp"
  exit
fi
rm compiler.txt
#
#  Here is the load statement for Apple's OS X.
#
g++ triangulation_faces.o -lnetcdf -lnetcdf_c++ -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ triangulation_faces.o -lnetcdf -lnetcdf_c++ -lm -lm -lGL -lGLU -lglut
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_faces.o"
  exit
fi
#
rm triangulation_faces.o
mv a.out ~/bincpp/$ARCH/triangulation_faces
#
echo "Executable installed as ~/bincpp/$ARCH/triangulation_faces"
