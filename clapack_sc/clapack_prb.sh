#!/bin/bash
#
#  The -I switch allows us to access the include file clapack.h.
#
g++ -c -I/usr/common/clapack clapack_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling clapack_prb.cpp"
  exit
fi
#
#  The -L switch allows us to access 4 libraries associated with CLAPACK.
#
g++ clapack_prb.o -L/usr/common/clapack -lclapack -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading clapack_prb.o."
  exit
fi
rm clapack_prb.o
#
mv a.out clapack_prb
./clapack_prb > clapack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running clapack_prb"
  exit
fi
#
rm clapack_prb
#
echo "clapack_prb output written to clapack_prb_output.txt"

