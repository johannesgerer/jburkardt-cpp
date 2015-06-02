#!/bin/bash
#
#  The -I switch allows us to access the include file clapack.h.
#
g++ -c -I/usr/common/clapack clapack_prb3.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling clapack_prb3.cpp"
  exit
fi
#
#  The -L switch allows us to access 4 libraries associated with CLAPACK.
#
g++ clapack_prb3.o -L/usr/common/clapack -lclapack -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading clapack_prb3.o."
  exit
fi
rm clapack_prb3.o
#
mv a.out clapack_prb3
./clapack_prb3 > clapack_prb3_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running clapack_prb3"
  exit
fi
#
rm clapack_prb3
#
echo "clapack_prb3 output written to clapack_prb3_output.txt"

