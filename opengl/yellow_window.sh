#!/bin/bash
#
g++ -c yellow_window.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling yellow_window.cpp"
  exit
fi
rm compiler.txt
#
g++ yellow_window.o -framework GLUT -framework OpenGL
#g++ yellow_window.o -lGL -lGLU -lglut -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading yellow_window.o"
  exit
fi
#
rm yellow_window.o
mv a.out ~/bincpp/$ARCH/yellow_window
#
echo "Executable installed as ~/bincpp/$ARCH/yellow_window"
