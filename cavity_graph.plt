set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/velocity_vector.png'
set xrange [-0.2:1.2]
set yrange [-0.2:1.2]
adj = 0.15
set size square
plot "cavity_wall.csv" with lines title "Wall", 'outputs/data/velocity.csv' using 1:2:(adj*$3):(adj*$4) with vectors title "Velocity Vector"

set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/velocity_vector_corner1.png'
set xrange[-0.01:0.25]
set yrange[-0.01:0.25]
adj = 50
set size square
set key outside
plot 'cavity_wall.csv' with lines title "Wall", 'outputs/data/velocity_corner1.csv' using 1:2:(adj*$3):(adj*$4) with vectors title "Velocity Vector"

set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/velocity_vector_corner2.png'
set xrange[0.9:1.01]
set yrange[-0.01:0.1]
adj = 50
set size square
plot 'cavity_wall.csv' with lines title "Wall", 'outputs/data/velocity_corner2.csv' using 1:2:(adj*$3):(adj*$4) with vectors title "Velocity Vector"

set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/velocity_abs.png'
set xrange [-0.2:1.2]
set yrange [-0.2:1.2]
adj = 0.15
set size square
plot "cavity_wall.csv" with lines title "Wall", 'outputs/data/velocity_abs.csv' using 1:2:(adj*$3):(adj*$4) with vectors title "Velocity direction "

set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/stream_line.png'
set contour base
set view 0,0
set view equal xy
set key at 2, 1.6
splot 'cavity_wall.csv' with lines title "Wall" nocontour ,'outputs/data/psi.csv' w l title "{/Symbol Y}" nosurface

set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/stream_line_corner1.png'
set contour base
set view 0,0
set view equal xy
set xrange[-0.01:0.30]
set yrange[-0.01:0.30]
set key default
set key at 0.52, 0.33
set cntrparam levels incremental 0.0, -0.00002, -0.1
splot 'cavity_wall.csv' with lines title "Wall" nocontour ,'outputs/data/psi.csv' w l title "{/Symbol Y}" nosurface

set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/stream_line_corner2.png'
set contour base
set view 0,0
set xrange[0.9:1.01]
set yrange[-0.01:0.1]
set view equal xy
set cntrparam levels incremental 0.0, -0.000001, -0.1
set key at 1.09, 0.13
set xtics 0.9, 0.02, 1
splot 'cavity_wall.csv' with lines title "Wall" nocontour ,'outputs/data/psi.csv' w l title "{/Symbol Y}" nosurface

set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/p.png'
set contour base
set view 0,0
set pm3d
set palette rgbformula 22,13,-31
set cntrparam levels 40
set noclabel
set cbrange[-0.2:0.2]
set view equal xy
splot 'cavity_wall.csv' with lines linewidth 2 title "Wall" nocontour ,'outputs/data/p.csv' w l linewidth 2 title "Pressure Contour" nosurface


pause -1