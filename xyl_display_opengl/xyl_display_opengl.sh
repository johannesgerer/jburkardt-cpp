#!/bin/bash
#
g++ -c xyl_display_opengl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xyl_display_opengl.cpp"
  exit
fi
rm compiler.txt
#
#  Here is the load statement for Apple's OS X.
#
g++ xyl_display_opengl.o -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ xyl_display_opengl.o -lm -lGL -lGLU -lglut
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xyl_display_opengl.o"
  exit
fi
#
rm xyl_display_opengl.o
mv a.out ~/bincpp/$ARCH/xyl_display_opengl
#
echo "Executable installed as ~/bincpp/$ARCH/xyl_display_opengl"
