#!/bin/bash
#
g++ -c -g -I/$HOME/include cnf_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cnf_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ cnf_io_prb.o /$HOME/libcpp/$ARCH/cnf_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cnf_io_prb.o."
  exit
fi
#
rm cnf_io_prb.o
#
mv a.out cnf_io_prb
./cnf_io_prb > cnf_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cnf_io_prb."
  exit
fi
rm cnf_io_prb
#
echo "Program output written to cnf_io_prb_output.txt"
