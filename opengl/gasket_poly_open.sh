#!/bin/bash
#
g++ -c gasket_poly_open.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gasket_poly_open.cpp"
  exit
fi
rm compiler.txt
#
g++ gasket_poly_open.o -framework GLUT -framework OpenGL
#g++ gasket_poly_open.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gasket_poly_open.o"
  exit
fi
#
rm gasket_poly_open.o
mv a.out ~/bincpp/$ARCH/gasket_poly_open
#
echo "Executable installed as ~/bincpp/$ARCH/gasket_poly_open"
