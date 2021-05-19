program cavity
    implicit none

    type :: latticePoint
        double precision u, v
    end type latticePoint

    integer, parameter :: latticeSizeX = 10
    integer, parameter :: latticeSizeY = 10

    double precision, parameter :: reynolds = 10

    double precision :: h = 1.0 / latticeSizeX
    

    type(latticePoint) :: latticePoints(latticeSizeX, latticeSizeY)


    integer :: i, j, k, l, m

    latticePoints = latticePoint(u=0.0, v=0.0)

    open(11, file='lattice.csv')
    do k = 1, latticeSizeX
        write(11, '(*(f0.0:" "))') latticePoints(:, i)%u
    end do
    close(11)

    do m = 1, 100
        do i = 1, latticeSizeX
            do j = 1, latticeSizeY
                if (i == 1) then 
                    latticePoints(i, j)%u = 0
                    latticePoints(i, j)%v = 2.0 * latticePoints(i + 1, j)%u
                else if (i == latticeSizeX) then
                    latticePoints(i, j)%u = 0
                    latticePoints(i, j)%v = 2.0 * latticePoints(i - 1, j)%u
                else if (j == 1) then
                    latticePoints(i, j)%u = 0
                    latticePoints(i, j)%v = 2.0 * latticePoints(i , j + 1)%u  
                else if (j == latticeSizeY) then 
                    latticePoints(i, j)%u = 0
                    latticePoints(i, j)%v = 2.0 * (latticePoints(i , j - 1)%u - h)
                else
                    latticePoints(i, j)%v = &
                        (latticePoints(i-1, j)%v + &
                        latticePoints(i+1, j)%v + &
                        latticePoints(i, j-1)%v + &
                        latticePoints(i, j+1)%v) / 4.0 + &
                        (latticePoints(i-1, j)%v - latticePoints(i+1, j)%v) * &
                        (latticePoints(i, j-1)%u - latticePoints(i, j+1)%u) - &
                        (latticePoints(i-1, j)%u - latticePoints(i+1, j)%u) * &
                        (latticePoints(i, j-1)%v - latticePoints(i, j+1)%v) * &
                        reynolds / 16.0
                    
                    latticePoints(i, j)%u = &
                        (latticePoints(i-1, j)%u + &
                        latticePoints(i+1, j)%u + &
                        latticePoints(i, j+1)%u + &
                        latticePoints(i, j-1)%u - &
                        latticePoints(i, j)%v) / 4.0
                end if
            end do
        end do
    end do
    
    open(10, file='cavity_u.csv')
    do k = 1, latticeSizeY
        write(10, '(*(f0.10:" "))') latticePoints(:, k)%u
    end do
    close(10)

    open(12, file='cavity_v.csv')
    do l = 1, latticeSizeY
        write(12, '(*(f0.10:" "))') latticePoints(:, l)%v
    end do
    close(12)

end program cavity
