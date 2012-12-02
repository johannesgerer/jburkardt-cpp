#!/bin/bash
#
g++ -c gasket_poly_filled.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gasket_poly_filled.cpp"
  exit
fi
rm compiler.txt
#
g++ gasket_poly_filled.o -framework GLUT -framework OpenGL
#g++ gasket_poly_filled.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gasket_poly_filled.o"
  exit
fi
#
rm gasket_poly_filled.o
mv a.out ~/bincpp/$ARCH/gasket_poly_filled
#
echo "Executable installed as ~/bincpp/$ARCH/gasket_poly_filled"
