#!/bin/bash
#
g++ -c life_opengl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling life_opengl.cpp"
  exit
fi
rm compiler.txt
#
g++ life_opengl.o -framework OpenGL -framework GLUT
if [ $? -ne 0 ]; then
  echo "Errors linking life_opengl.o"
  exit
fi
#
rm life_opengl.o
mv a.out ~/bincpp/$ARCH/life_opengl
#
echo "Executable installed as ~/bincpp/$ARCH/life_opengl"
