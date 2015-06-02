#!/bin/bash
#
g++ -c -I/$HOME/include pyramid_integrals_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_integrals_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ pyramid_integrals_prb.o /$HOME/libcpp/$ARCH/pyramid_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pyramid_integrals_prb.o"
  exit
fi
#
rm pyramid_integrals_prb.o
#
mv a.out pyramid_integrals_prb
./pyramid_integrals_prb > pyramid_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pyramid_integrals_prb."
  exit
fi
rm pyramid_integrals_prb
#
echo "Program output written to pyramid_integrals_prb_output.txt"
