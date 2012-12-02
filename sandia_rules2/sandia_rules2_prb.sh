#!/bin/bash
#
g++ -c -g -I/$HOME/include sandia_rules2_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_rules2_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sandia_rules2_prb.o /$HOME/libcpp/$ARCH/sandia_rules2.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sandia_rules2_prb.o."
  exit
fi
#
rm sandia_rules2_prb.o
#
mv a.out sandia_rules2_prb
./sandia_rules2_prb > sandia_rules2_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sandia_rules2_prb."
  exit
fi
rm sandia_rules2_prb
#
echo "Program output written to sandia_rules2_prb_output.txt"
