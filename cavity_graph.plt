set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/velocity_vector.png'
set xrange [-0.2:1.2]
set yrange [-0.2:1.2]
adj = 0.15
set size square
plot "cavity_wall.csv" with lines title "Wall", 'outputs/data/velocity.csv' using 1:2:(adj*$3):(adj*$4) with vectors title "Velocity Vector"

set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/stream_line.png'
set contour base
set view 0,0
set view equal xy
set key at 2, 1.6
splot 'cavity_wall.csv' with lines title "Wall" nocontour ,'outputs/data/psi.csv' w l title "{/Symbol Y}" nosurface

set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/p.png'
set contour base
set view 0,0
set pm3d
set palette rgbformula 22,13,-31
set cntrparam levels 10
set noclabel
set view equal xy
splot 'cavity_wall.csv' with lines linewidth 2 title "Wall" nocontour ,'outputs/data/p.csv' w l linewidth 2 title "Pressure Contour" nosurface


set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/p1.png'
set contour base
set view 50,100
set view equal xy
splot 'outputs/data/p.csv' w l title "Streamline"

set view map
set palette rgbformula 22,13,-31
set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/p2.png'
set size square
splot 'outputs/data/p.csv' with pm3d

pause -1