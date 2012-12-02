#!/bin/bash
#
g++ -c -g -I $HOME/include linpack_bench_s.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_bench_s.cpp"
  exit
fi
rm compiler.txt
#
g++ linpack_bench_s.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linpack_bench_s.o."
  exit
fi
#
rm linpack_bench_s.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/linpack_bench_s
#
echo "Executable installed as ~/bincpp/$ARCH/linpack_bench_s."
