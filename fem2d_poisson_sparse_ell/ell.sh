#!/bin/bash
#
g++ -c -g ell.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ell.cpp"
  exit
fi
rm compiler.txt
#
g++ ~/libcpp/$ARCH/fem2d_poisson_sparse.o ell.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_poisson_sparse.o + ell.o"
  exit
fi
rm ell.o
#
chmod ugo+x a.out
mv a.out fem2d_poisson_sparse_ell
./fem2d_poisson_sparse ell > ell_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ell."
  exit
fi
rm fem2d_poisson_sparse_ell
#
if [ -e ell_elements.eps ]; then
  convert ell_elements.eps ell_elements.png
  rm ell_elements.eps
fi
#
if [ -e ell_nodes.eps ]; then
  convert ell_nodes.eps ell_nodes.png
  rm ell_nodes.eps
fi
#
echo "Program output written to ell_output.txt"
