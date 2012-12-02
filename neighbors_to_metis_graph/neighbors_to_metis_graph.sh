#!/bin/bash
#
g++ -c -g -I$HOME/include neighbors_to_metis_graph.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling neighbors_to_metis_graph.cpp"
  exit
fi
rm compiler.txt
#
g++ neighbors_to_metis_graph.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading neighbors_to_metis_graph.o."
  exit
fi
#
rm neighbors_to_metis_graph.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/neighbors_to_metis_graph
#
echo "Executable installed as ~/bincpp/$ARCH/neighbors_to_metis_graph"
