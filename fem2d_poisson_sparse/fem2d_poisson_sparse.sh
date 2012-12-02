#!/bin/bash
#
g++ -c -g fem2d_poisson_sparse.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_poisson_sparse.cpp"
  exit
fi
rm compiler.txt
#
mv fem2d_poisson_sparse.o ~/libcpp/$ARCH
#
echo "Partial program installed as ~/libcpp/$ARCH/fem2d_poisson_sparse.o"
