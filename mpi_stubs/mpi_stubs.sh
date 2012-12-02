#!/bin/bash
#
cp mpi_stubs.hpp $HOME/include
#
g++ -c -g mpi_stubs.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mpi_stubs.cpp"
  exit
fi
rm compiler.txt
#
mv mpi_stubs.o ~/libcpp/$ARCH/mpi_stubs.o
#
echo "A new version of mpi_stubs has been created."
