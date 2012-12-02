#!/bin/bash
#
g++ -c shaded_spheres.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shaded_spheres.cpp"
  exit
fi
rm compiler.txt
#
g++ shaded_spheres.o -framework GLUT -framework OpenGL
#g++ shaded_spheres.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Errors linking and loading shaded_spheres.o"
  exit
fi
#
rm shaded_spheres.o
mv a.out ~/bincpp/$ARCH/shaded_spheres
#
echo "Executable installed as ~/bincpp/$ARCH/shaded_spheres"
