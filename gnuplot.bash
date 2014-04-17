#!/usr/bin/gnuplot
set terminal png
set xlabel "x"
set ylabel "f(x)"
set title "Example1_plot"  
set term png
set output "Example1_plot.png"
plot "Example1.asc" title "Dan's plot"


