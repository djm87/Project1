#!/usr/bin/gnuplot

reset
set terminal png
set xlabel "x"
set ylabel "u"
set title "Example_plot"  
set term gif animate
set output "Test_plot.gif"
i=1
n=100
d = "test512.asc"
load "animate.gnuplot"


