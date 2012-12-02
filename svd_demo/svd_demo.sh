#!/bin/bash
#
g++ -c -g -I$HOME/include svd_demo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling svd_demo.cpp"
  exit
fi
rm compiler.txt
#
g++ svd_demo.o $HOME/libcpp/$ARCH/linpack_d.o $HOME/libcpp/$ARCH/blas1_d.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading svd_demo.o + linpack_d.o + blas1_d.o."
  exit
fi
#
rm svd_demo.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/svd_demo
#
echo "Executable installed as ~/bincpp/$ARCH/svd_demo"
