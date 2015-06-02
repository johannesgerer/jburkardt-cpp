#!/bin/bash
#
g++ -c fem_to_gmsh.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_to_gmsh.cpp"
  exit
fi
#
g++ fem_to_gmsh.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_to_gmsh.o."
  exit
fi
#
rm fem_to_gmsh.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem_to_gmsh
#
echo "Executable installed as ~/bincpp/$ARCH/fem_to_gmsh"
