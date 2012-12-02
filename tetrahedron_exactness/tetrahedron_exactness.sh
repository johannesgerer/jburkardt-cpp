#!/bin/bash
#
g++ -c -g tetrahedron_exactness.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_exactness.cpp"
  exit
fi
rm compiler.txt
#
g++ tetrahedron_exactness.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tetrahedron_exactness.o"
  exit
fi
rm tetrahedron_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/tetrahedron_exactness
#
echo "Program installed as ~/bincpp/$ARCH/tetrahedron_exactness"
