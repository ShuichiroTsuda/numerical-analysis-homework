set term pngcairo size 1280, 960 font ",24"
set output 'outputs/img/velocity_vector.png'
set xrange [-0.2:1.2]
set yrange [-0.2:1.2]
adj = 0.1
plot "cavity_wall.csv" with lines title "Wall", 'outputs/data/velocity.csv' using 1:2:(adj*$3):(adj*$4) with vectors title "Velocity Vector"

set term pngcairo size 640,480 font ",12"
set output 'outputs/img/stream_line.png'
set contour base
set view 0,0
set noclabel
splot 'cavity_wall.csv' with lines title "Wall" nocontour ,'outputs/data/psi.csv' w l title "Streamline" nosurface

set term pngcairo size 640,480 font ",12"
set output 'outputs/img/p.png'
set contour base
set view 0,0
set noclabel
splot 'cavity_wall.csv' with lines title "Wall" nocontour ,'outputs/data/p.csv' w l title "Streamline" nosurface


set term pngcairo size 640,480 font ",12"
set output 'outputs/img/p1.png'
set contour base
set view 50,100
splot 'outputs/data/p.csv' w l title "Streamline"

pause -1