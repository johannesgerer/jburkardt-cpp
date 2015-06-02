#!/bin/bash
#
g++ -c fem2d_mesh_display_opengl.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_mesh_display_opengl.cpp"
  exit
fi
#
g++ fem2d_mesh_display_opengl.o -lm -lGL -lGLU -lglut
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_mesh_display_opengl.o"
  exit
fi
#
rm fem2d_mesh_display_opengl.o
mv a.out ~/bincpp/$ARCH/fem2d_mesh_display_opengl
#
echo "Executable installed as ~/bincpp/$ARCH/fem2d_mesh_display_opengl"
