#!/bin/bash
#
#  Compile
#
g++ -c -I/usr/local/include umfpack_simple.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling umfpack_simple.cpp"
  exit
fi
#
#  Link and load
#
#g++ umfpack_simple.o -L/usr/local/lib -L/$HOME/lib/$ARCH -lumfpack -lamd -lcholmod -lcolamd -lm \
#  -lsuitesparseconfig -lblas
gfortran umfpack_simple.o -L/usr/local/lib -L/$HOME/lib/$ARCH -lumfpack -lamd -lcholmod \
  -lcolamd -lm -lsuitesparseconfig -lblas -lstdc++
if [ $? -ne 0 ]; then
  echo "Errors linking and loading umfpack_simple.o"
  exit
fi
rm umfpack_simple.o
mv a.out umfpack_simple
#
#  Run
#
./umfpack_simple > umfpack_simple_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running umfpack_simple"
  exit
fi
rm umfpack_simple
#
#  Terminate.
#
echo "Program output written to umfpack_simple_output.txt"
