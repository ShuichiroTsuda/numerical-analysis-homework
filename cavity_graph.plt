set term pngcairo size 1280, 960 font ",24"
set output 'outputs/velocity_vector.png'
set xrange [-0.1:1.1]
set yrange [-0.1:1.1]
adj = 0.1
plot 'cavity_velocity.csv' using 1:2:(adj*$3):(adj*$4) with vectors title "Velocity Vector"
pause -1