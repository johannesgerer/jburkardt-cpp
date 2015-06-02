#!/bin/bash
#
g++ -c -g ellipse.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipse.cpp"
  exit
fi
rm compiler.txt
#
g++ ellipse.o $HOME/libcpp/$ARCH/toms886.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ellipse.o"
  exit
fi
rm ellipse.o
#
mv a.out ellipse
./ellipse > ellipse_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ellipse"
  exit
fi
rm ellipse
#
echo "Test results written to ellipse_output.txt."
