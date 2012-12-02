#!/bin/bash
#
g++ -c -g triangulation_plot.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_plot.cpp"
  exit
fi
rm compiler.txt
#
g++ triangulation_plot.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_plot.o."
  exit
fi
#
rm triangulation_plot.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/triangulation_plot
#
echo "Executable installed as ~/bincpp/$ARCH/triangulation_plot"
