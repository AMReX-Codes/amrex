# plot the Sod results.  

#set term post eps color
#set output 'figure1.eps'
set term x11

set multiplot;

set size 0.5, 0.5;

set xlabel "x";

set style line 1 lt rgb "black" lw 1

set origin 0.0, 0.5;
set key top center
set ylabel "density";
plot 'density_x' using 1:2 title   '128' with p pt 1 ps 1.0,\
     'density_z' using 1:2 title   'new' with p pt 1 ps 1.0,\
     '../analytic/sod-exact.out' using 1:2 title 'exact' with lines ls 1;
set origin 0.5, 0.5;
set key top left
set ylabel "velocity";
plot 'velocity_x' using 1:2 title   '128' with p pt 1 ps 1.0,\
     'velocity_z' using 1:2 title   'new' with p pt 1 ps 1.0,\
     '../analytic/sod-exact.out' using 1:3 title 'exact' with lines ls 1;

set origin 0.0, 0.0;
set key top center
set ylabel "pressure";
plot 'pressure_x' using 1:2 title   '128' with p pt 1 ps 1.0,\
     'pressure_z' using 1:2 title   'new' with p pt 1 ps 1.0,\
     '../analytic/sod-exact.out' using 1:4 title 'exact' with lines ls 1;

set origin 0.5, 0.0;
set key top center
set ylabel "internal energy";
plot 'eint_x' using 1:2 title   '128' with p pt 1 ps 1.0,\
     'eint_z' using 1:2 title   'new' with p pt 1 ps 1.0,\
     '../analytic/sod-exact.out' using 1:5 title 'exact' with lines ls 1;

unset multiplot;
set term x11;
