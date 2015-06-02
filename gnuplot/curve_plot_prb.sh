#!/bin/bash
#
g++ -o curve_plot_prb curve_plot_prb.cpp curve_plot.cpp
curve_plot_prb
gnuplot < curve_plot_commands.txt

