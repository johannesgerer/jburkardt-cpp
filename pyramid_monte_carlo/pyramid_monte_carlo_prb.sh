#!/bin/bash
#
g++ -c -I/$HOME/include pyramid_monte_carlo_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_monte_carlo_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ pyramid_monte_carlo_prb.o /$HOME/libcpp/$ARCH/pyramid_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pyramid_monte_carlo_prb.o"
  exit
fi
#
rm pyramid_monte_carlo_prb.o
#
mv a.out pyramid_monte_carlo_prb
./pyramid_monte_carlo_prb > pyramid_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pyramid_monte_carlo_prb."
  exit
fi
rm pyramid_monte_carlo_prb
#
echo "Program output written to pyramid_monte_carlo_prb_output.txt"
