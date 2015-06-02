#!/bin/bash
#
g++ -c -g ell.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ell.cpp"
  exit
fi
rm compiler.txt
#
g++ ~/libcpp/$ARCH/fem2d_poisson_cg.o ell.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_poisson_cg.o + ell.o"
  exit
fi
rm ell.o
#
chmod ugo+x a.out
mv a.out fem2d_poisson_cg_ell
./fem2d_poisson_cg_ell ell > ell_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ell."
  exit
fi
rm fem2d_poisson_cg_ell
#
#
echo "Program output written to ell_output.txt"
