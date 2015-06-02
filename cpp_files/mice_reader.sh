#!/bin/bash
#
g++ -c -g mice_reader.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling mice_reader.cpp"
  exit
fi
rm compiler.txt
#
g++ mice_reader.o
if [ $? -ne 0 ]; then
  echo "Errors while loading mice_reader.o"
  exit
fi
rm mice_reader.o
#
mv a.out mice_reader
./mice_reader > mice_reader_output.txt
rm mice_reader
#
echo "Program output written to mice_reader_output.txt"
