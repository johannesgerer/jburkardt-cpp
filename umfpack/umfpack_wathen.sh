#!/bin/bash
#
#  Compile
#
g++ -c -I/usr/local/include umfpack_wathen.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling umfpack_wathen.cpp"
  exit
fi
#
#  Link and load
#
#g++ umfpack_wathen.o -L/usr/local/lib -L/$HOME/lib/$ARCH -lumfpack -lamd -lcholmod -lcolamd -lm \
#  -lsuitesparseconfig -lblas
gfortran umfpack_wathen.o -L/usr/local/lib -L/$HOME/lib/$ARCH -lumfpack -lamd -lcholmod -lcolamd -lm \
  -lsuitesparseconfig -lblas -lstdc++
if [ $? -ne 0 ]; then
  echo "Errors linking and loading umfpack_wathen.o"
  exit
fi
rm umfpack_wathen.o
mv a.out umfpack_wathen
#
#  Run
#
./umfpack_wathen > umfpack_wathen_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running umfpack_wathen"
  exit
fi
rm umfpack_wathen
#
#  Terminate.
#
echo "Program output written to umfpack_wathen_output.txt"
