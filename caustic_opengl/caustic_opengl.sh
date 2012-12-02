#!/bin/bash
#
g++ -c caustic_opengl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling caustic_opengl.cpp"
  exit
fi
rm compiler.txt
#
#
#  Here is the load statement for Apple's OS X.
#
g++ caustic_opengl.o -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ caustic_opengl.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Errors linking and loading caustic_opengl.o"
  exit
fi
#
rm caustic_opengl.o
mv a.out ~/bincpp/$ARCH/caustic_opengl
#
echo "Executable installed as ~/bincpp/$ARCH/caustic_opengl"
