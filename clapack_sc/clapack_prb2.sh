#!/bin/bash
#
#  The -I switch allows us to access the include file clapack.h.
#
g++ -c -I/usr/common/clapack clapack_prb2.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling clapack_prb2.cpp"
  exit
fi
#
#  The -L switch allows us to access 4 libraries associated with CLAPACK.
#
g++ clapack_prb2.o -L/usr/common/clapack -lclapack -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading clapack_prb2.o."
  exit
fi
rm clapack_prb2.o
#
mv a.out clapack_prb2
./clapack_prb2 > clapack_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running clapack_prb2"
  exit
fi
#
rm clapack_prb2
#
echo "clapack_prb2 output written to clapack_prb2_output.txt"

