#!/bin/bash
#
g++ -c shading.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shading.cpp"
  exit
fi
rm compiler.txt
#
g++ shading.o -framework GLUT -framework OpenGL
#g++ shading.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Errors linking and loading shading.o"
  exit
fi
#
rm shading.o
mv a.out ~/bincpp/$ARCH/shading
#
echo "Executable installed as ~/bincpp/$ARCH/shading"
