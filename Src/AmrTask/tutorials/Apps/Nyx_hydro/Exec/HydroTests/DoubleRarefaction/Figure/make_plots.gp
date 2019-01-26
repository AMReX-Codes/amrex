# plot the Test2 results.  Here we assume that we have files 
# called test2x.out, test2y.out, and test2z.out

#set term post eps color
#set output 'figure2.eps'

set term x11

set multiplot;

set size 0.5, 0.5;

set xlabel "x";

set style line 1 lt rgb "black" lw 1

set origin 0.0, 0.5;
set key top center
set ylabel "density";
plot '128/density_x' using 1:2 title   '128' with p pt 1 ps 1.0,\
     '512/density_x' using 1:2 title   '512' with p pt 1 ps 1.0,\
      '../analytic/exact.out' using 1:2 title 'exact' with lines ls 1;
set origin 0.5, 0.5;
set key top left
set ylabel "velocity";
plot '128/velocity_x' using 1:2 title   '128' with p pt 1 ps 1.0,\
     '512/velocity_x' using 1:2 title   '512' with p pt 1 ps 1.0,\
      '../analytic/exact.out' using 1:3 title 'exact' with lines ls 1;

set origin 0.0, 0.0;
set key top center
set ylabel "pressure";
plot '128/pressure_x' using 1:2 title   '128' with p pt 1 ps 1.0,\
     '512/pressure_x' using 1:2 title   '512' with p pt 1 ps 1.0,\
      '../analytic/exact.out' using 1:4 title 'exact' with lines ls 1;

set origin 0.5, 0.0;
set key top center
set ylabel "internal energy";
plot '128/eint_x' using 1:2 title   '128' with p pt 1 ps 1.0,\
     '512/eint_x' using 1:2 title   '512' with p pt 1 ps 1.0,\
      '../analytic/exact.out' using 1:5 title 'exact' with lines ls 1;

unset multiplot;
set term x11;
