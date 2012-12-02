#!/bin/bash
#
g++ -c -D _MACOS mpas_draw.C >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mpas_draw.C"
  exit
fi
rm compiler.txt
#
g++ -c netcdf_utils.C >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling netcdf_utils.C"
  exit
fi
rm compiler.txt
#
g++ -c vec_utils.C >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vec_utils.C"
  exit
fi
rm compiler.txt
#
#  Here is the load statement for Apple's OS X.
#
g++ mpas_draw.o netcdf_utils.o vec_utils.o -lm -framework GLUT -framework OpenGL -lnetcdf -lnetcdf_c++
#
#  Here is the load statement for a normal UNIX system!
#
#g++ mpas_draw.o netcdf_utils.o vec_utils.o -lm -lGL -lGLU -lglut -lnetcdf -lnetcdf_c++
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mpas_draw.o + netcdf_utils.o + vec_utils.o"
  exit
fi
#
rm mpas_draw.o
rm netcdf_utils.o
rm vec_utils.o

mv a.out ~/bincpp/$ARCH/mpas_draw
#
echo "Executable installed as ~/bincpp/$ARCH/mpas_draw"
