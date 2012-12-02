#!/bin/bash
#
g++ -c lissajous_opengl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lissajous_opengl.cpp"
  exit
fi
rm compiler.txt
#
#
#  Here is the load statement for Apple's OS X.
#
g++ lissajous_opengl.o -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ lissajous_opengl.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lissajous_opengl.o"
  exit
fi
#
rm lissajous_opengl.o
mv a.out ~/bincpp/$ARCH/lissajous_opengl
#
echo "Executable installed as ~/bincpp/$ARCH/lissajous_opengl"
