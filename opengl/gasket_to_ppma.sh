#!/bin/bash
#
g++ -c gasket_to_ppma.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gasket_to_ppma.cpp"
  exit
fi
rm compiler.txt
#
g++ gasket_to_ppma.o -framework GLUT -framework OpenGL
#g++ gasket_to_ppma.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gasket_to_ppma.o"
  exit
fi
#
rm gasket_to_ppma.o
mv a.out ~/bincpp/$ARCH/gasket_to_ppma
#
echo "Executable installed as ~/bincpp/$ARCH/gasket_to_ppma"
