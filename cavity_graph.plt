set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/velocity_vector.png'
set xrange [-0.1:1.1]
set yrange [-0.1:1.1]
adj = 0.1
plot 'outputs/data/velocity.csv' using 1:2:(adj*$3):(adj*$4) with vectors title "Velocity Vector"

set output 'outputs/img/stream_line.png'
set view 0,0
set dgrid3d
set contour base
set nosurface
splot 'outputs/data/phi.csv' using 1:2:3 w lp


pause -1