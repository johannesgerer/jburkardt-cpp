#!/bin/bash
#
g++ -c -I/usr/local/dislin quickplot_surface.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quickplot_surface.cpp"
  exit
fi
rm compiler.txt
#
g++ quickplot_surface.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quickplot_surface.o."
  exit
fi
#
rm quickplot_surface.o
#
mv a.out quickplot_surface
./quickplot_surface > quickplot_surface_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quickplot_surface."
  exit
fi
rm quickplot_surface
#
echo "Program output written to quickplot_surface_output.txt"
