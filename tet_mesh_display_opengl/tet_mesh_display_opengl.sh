#!/bin/bash
#
g++ -c tet_mesh_display_opengl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tet_mesh_display_opengl.cpp"
  exit
fi
rm compiler.txt
#
#
#  Here is the load statement for Apple's OS X.
#
g++ tet_mesh_display_opengl.o -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
#g++ tet_mesh_display_opengl.o -lm -lGL -lGLU -lglut
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tet_mesh_display_opengl.o"
  exit
fi
#
rm tet_mesh_display_opengl.o
mv a.out ~/bincpp/$ARCH/tet_mesh_display_opengl
#
echo "Executable installed as ~/bincpp/$ARCH/tet_mesh_display_opengl."
