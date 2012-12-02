#!/bin/bash
#
g++ -c triangulation_display_opengl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_display_opengl.cpp"
  exit
fi
rm compiler.txt
#
#  Here is the load statement for Apple's OS X.
#
g++ triangulation_display_opengl.o -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ triangulation_display_opengl.o -lm -lGL -lGLU -lglut
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_display_opengl.o"
  exit
fi
#
rm triangulation_display_opengl.o
mv a.out ~/bincpp/$ARCH/triangulation_display_opengl
#
echo "Executable installed as ~/bincpp/$ARCH/triangulation_display_opengl"
