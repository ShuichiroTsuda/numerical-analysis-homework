set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/velocity_vector.png'
set xrange [-0.2:1.2]
set yrange [-0.2:1.2]
adj = 0.1
plot "cavity_wall.csv" with lines title "Wall", 'outputs/data/velocity.csv' using 1:2:(adj*$3):(adj*$4) with vectors title "Velocity Vector"

set term pngcairo size 640,480 font ",24"
set output 'outputs/img/stream_line_w_surface.png'
set dgrid3d 50, 50, 5
set contour base
#set cntrparam cubicspline
splot 'outputs/data/phi.csv' using 1:2:3 w l title "Streamline"

set nodgrid3d

set term pngcairo size 640,480 font ",12"
set output 'outputs/img/stream_line_wo_param.png'
set dgrid3d
set view 0,0
set nosurface
set noclabel
set contour base
#set cntrparam cubicspline
splot 'outputs/data/phi.csv' using 1:2:3 w l title "Streamline"

set term pngcairo size 640,480 font ",12"
set output 'outputs/img/stream_line.png'
set view 0,0
set dgrid3d 50, 50, 5
set contour base
set nosurface
set noclabel
#set cntrparam cubicspline
splot 'outputs/data/phi.csv' using 1:2:3 w l title "Streamline"

pause -1