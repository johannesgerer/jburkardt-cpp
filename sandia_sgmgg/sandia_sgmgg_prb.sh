#!/bin/bash
#
g++ -c -g -I/$HOME/include sandia_sgmgg_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_sgmgg_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sandia_sgmgg_prb.o $HOME/libcpp/$ARCH/sandia_sgmgg.o $HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sandia_sgmgg_prb.o."
  exit
fi
#
rm sandia_sgmgg_prb.o
#
mv a.out sandia_sgmgg_prb
./sandia_sgmgg_prb > sandia_sgmgg_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sandia_sgmgg_prb."
  exit
fi
rm sandia_sgmgg_prb
#
echo "Program output written to sandia_sgmgg_prb_output.txt"
