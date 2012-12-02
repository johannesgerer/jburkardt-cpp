#!/bin/bash
#
g++ -c dijkstra.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dijkstra.cpp"
  exit
fi
rm compiler.txt
#
g++ dijkstra.o
if [ $? -ne 0 ]; then
  echo "Errors linking dijkstra.o."
  exit
fi
#
rm dijkstra.o
#
mv a.out dijkstra
./dijkstra > dijkstra_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dijkstra."
  exit
fi
rm dijkstra
#
echo "Program output written to dijkstra_output.txt"
