#!/bin/bash
#
#  Compile
#
g++ -c -I/usr/local/include umfpack_west.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling umfpack_west.cpp"
  exit
fi
#
#  Link and load
#
#g++ umfpack_west.o -L/usr/local/lib -L/$HOME/lib/$ARCH -lumfpack -lamd -lcholmod -lcolamd -lm \
#  -lsuitesparseconfig -lblas
gfortran umfpack_west.o -L/usr/local/lib -L/$HOME/lib/$ARCH -lumfpack -lamd -lcholmod -lcolamd -lm \
  -lsuitesparseconfig -lblas -lstdc++
if [ $? -ne 0 ]; then
  echo "Errors linking and loading umfpack_west.o"
  exit
fi
rm umfpack_west.o
mv a.out umfpack_west
#
#  Run
#
./umfpack_west > umfpack_west_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running umfpack_west"
  exit
fi
rm umfpack_west
#
#  Terminate.
#
echo "Program output written to umfpack_west_output.txt"
