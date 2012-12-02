#!/bin/bash
#
g++ -c -I/usr/local/dislin scatter_plot.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling scatter_plot.cpp"
  exit
fi
rm compiler.txt
#
g++ scatter_plot.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading scatter_plot.o."
  exit
fi
#
rm scatter_plot.o
#
mv a.out scatter_plot
./scatter_plot > scatter_plot_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running scatter_plot."
  exit
fi
rm scatter_plot
#
echo "Program output written to scatter_plot_output.txt"
