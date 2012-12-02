#!/bin/bash
#
g++ -c polygonal_surface_display_opengl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygonal_surface_display_opengl.cpp"
  exit
fi
rm compiler.txt
#
#
#  Here is the load statement for Apple's OS X.
#
g++ polygonal_surface_display_opengl.o -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ polygonal_surface_display_opengl.o -lm -lGL -lGLU -lglut
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading polygonal_surface_display_opengl.o"
  exit
fi
#
rm polygonal_surface_display_opengl.o
mv a.out ~/bincpp/$ARCH/polygonal_surface_display_opengl
#
echo "Executable installed as ~/bincpp/$ARCH/polygonal_surface_display_opengl."
