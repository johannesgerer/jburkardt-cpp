#!/bin/bash
#
g++ -c lights_out_opengl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lights_out_opengl.cpp"
  exit
fi
rm compiler.txt
#
g++ lights_out_opengl.o -framework OpenGL -framework GLUT
if [ $? -ne 0 ]; then
  echo "Errors linking lights_out_opengl.o"
  exit
fi
#
rm lights_out_opengl.o
mv a.out ~/bincpp/$ARCH/lights_out_opengl
#
echo "Executable installed as ~/bincpp/$ARCH/lights_out_opengl"
