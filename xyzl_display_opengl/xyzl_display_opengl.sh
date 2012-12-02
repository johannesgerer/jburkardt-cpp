#!/bin/bash
#
g++ -c xyzl_display_opengl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xyzl_display_opengl.cpp"
  exit
fi
rm compiler.txt
#
#  Here is the load statement for Apple's OS X.
#
g++ xyzl_display_opengl.o -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ xyzl_display_opengl.o -lm -lGL -lGLU -lglut
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xyzl_display_opengl.o"
  exit
fi
#
rm xyzl_display_opengl.o
mv a.out ~/bincpp/$ARCH/xyzl_display_opengl
#
echo "Executable installed as ~/bincpp/$ARCH/xyzl_display_opengl"
