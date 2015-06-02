#!/bin/bash
#
g++ -c -I/$HOME/include gmsh_io_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling gmsh_io_prb.cpp"
  exit
fi
#
g++ gmsh_io_prb.o /$HOME/libcpp/$ARCH/gmsh_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gmsh_io_prb.o"
  exit
fi
#
rm gmsh_io_prb.o
#
mv a.out gmsh_io_prb
./gmsh_io_prb > gmsh_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running gmsh_io_prb."
  exit
fi
rm gmsh_io_prb
#
echo "Program output written to gmsh_io_prb_output.txt"
