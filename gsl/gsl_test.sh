#!/bin/bash
#
g++ -c -g gsl_test.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gsl_test.cpp"
  exit
fi
rm compiler.txt
#
g++ gsl_test.o -lm -lgsl -lgslcblas
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gsl_test.o."
  exit
fi
#
rm gsl_test.o
#
mv a.out gsl_test
./gsl_test > gsl_test_output.txt
rm gsl_test
#
echo "Program output written to gsl_test_output.txt"
