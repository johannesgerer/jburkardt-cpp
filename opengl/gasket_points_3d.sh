#!/bin/bash
#
g++ -c gasket_points_3d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gasket_points_3d.cpp"
  exit
fi
rm compiler.txt
#
g++ gasket_points_3d.o -framework GLUT -framework OpenGL
#g++ gasket_points_3d.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gasket_points_3d.o"
  exit
fi
#
rm gasket_points_3d.o
mv a.out ~/bincpp/$ARCH/gasket_points_3d
#
echo "Executable installed as ~/bincpp/$ARCH/gasket_points_3d"
