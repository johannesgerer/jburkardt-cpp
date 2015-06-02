#!/bin/bash
#
g++ -c -g baffle.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling baffle.cpp"
  exit
fi
rm compiler.txt
#
g++ ~/libcpp/$ARCH/fem2d_poisson_cg.o baffle.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_poisson_cg.o + baffle.o"
  exit
fi
rm baffle.o
#
chmod ugo+x a.out
mv a.out fem2d_poisson_cg_baffle
./fem2d_poisson_cg_baffle baffle > baffle_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running baffle."
  exit
fi
rm fem2d_poisson_cg_baffle
#
if [ -e baffle_elements.eps ]; then
  convert baffle_elements.eps baffle_elements.png
  rm baffle_elements.eps
fi
#
if [ -e baffle_nodes.eps ]; then
  convert baffle_nodes.eps baffle_nodes.png
  rm baffle_nodes.eps
fi
#
echo "Program output written to baffle_output.txt"
