#!/bin/bash
#
g++ -c -g lake.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lake.cpp"
  exit
fi
rm compiler.txt
#
g++ ~/libcpp/$ARCH/fem2d_poisson_sparse.o lake.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_poisson_sparse.o + lake.o"
  exit
fi
rm lake.o
#
chmod ugo+x a.out
mv a.out fem2d_poisson_sparse_lake
./fem2d_poisson_sparse_lake lake > lake_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lake."
  exit
fi
rm fem2d_poisson_sparse_lake
#
if [ -e lake_elements.eps ]; then
  convert lake_elements.eps lake_elements.png
  rm lake_elements.eps
fi
#
if [ -e lake_nodes.eps ]; then
  convert lake_nodes.eps lake_nodes.png
  rm lake_nodes.eps
fi
#
echo "Program output written to lake_output.txt"
