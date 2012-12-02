#!/bin/bash
#
g++ -c -I/usr/local/dislin quickplot_pie.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quickplot_pie.cpp"
  exit
fi
rm compiler.txt
#
g++ quickplot_pie.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quickplot_pie.o."
  exit
fi
#
rm quickplot_pie.o
#
mv a.out quickplot_pie
./quickplot_pie > quickplot_pie_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quickplot_pie."
  exit
fi
rm quickplot_pie
#
echo "Program output written to quickplot_pie_output.txt"
