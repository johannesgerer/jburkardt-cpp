#!/bin/bash
#
g++ -c -g -I/$HOME/include chebyshev_series_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_series_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ chebyshev_series_prb.o /$HOME/libcpp/$ARCH/chebyshev_series.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chebyshev_series_prb.o"
  exit
fi
#
rm chebyshev_series_prb.o
#
mv a.out chebyshev_series_prb
./chebyshev_series_prb > chebyshev_series_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running chebyshev_series_prb."
  exit
fi
rm chebyshev_series_prb
#
echo "Program output written to chebyshev_series_prb_output.txt"
