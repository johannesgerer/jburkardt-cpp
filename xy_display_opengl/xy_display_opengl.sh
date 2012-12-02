#!/bin/bash
#
g++ -c xy_display_opengl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xy_display_opengl.cppopengl"
  exit
fi
rm compiler.txt
#
#  Here is the load statement for Apple's OS X.
#
g++ xy_display_opengl.o -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ xy_display_opengl.o -lm -lGL -lGLU -lglut
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xy_display_opengl.o"
  exit
fi
#
rm xy_display_opengl.o
mv a.out ~/bincpp/$ARCH/xy_display_opengl
#
echo "Executable installed as ~/bincpp/$ARCH/xy_display_opengl"
