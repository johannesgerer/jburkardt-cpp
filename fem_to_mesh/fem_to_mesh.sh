#!/bin/bash
#
g++ -c fem_to_mesh.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_to_mesh.cpp"
  exit
fi
rm compiler.txt
#
g++ fem_to_mesh.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_to_mesh.o"
  exit
fi
#
rm fem_to_mesh.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem_to_mesh
#
echo "Program installed as ~/bincpp/$ARCH/fem_to_mesh"
