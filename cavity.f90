program cavity
    implicit none

    type :: latticePoint
        double precision phi, omega
        double precision u, v
    end type latticePoint

    integer, parameter :: latticeSizeX = 50
    integer, parameter :: latticeSizeY = 50

    double precision, parameter :: reynolds = 10

    double precision :: h = 1.0 / latticeSizeX
    
    type(latticePoint) :: latticePoints(latticeSizeX, latticeSizeY)
    type(latticePoint) :: lastLaticePoint !収束判断のための直前の値

    integer :: i, j, k, l
    integer :: loopCount = 0

    logical :: shouldContinue = .true.

    double precision, parameter :: convergenceThreshold = 0.00001

    latticePoints = latticePoint(phi=0.0, omega=0.0, u=0.0, v=0.0)

    open(11, file='lattice.csv')
    do k = 1, latticeSizeX
        write(11, '(*(f0.0:" "))') latticePoints(:, i)%u
    end do
    close(11)
    
    do while(shouldContinue)
        shouldContinue = .false.
        loopCount = loopCount + 1
        do i = 1, latticeSizeX
            do j = 1, latticeSizeY
                lastLaticePoint = latticePoints(i, j)
                if (i == 1) then 
                    latticePoints(i, j)%phi = 0
                    latticePoints(i, j)%omega = 2.0 * latticePoints(i + 1, j)%phi
                else if (i == latticeSizeX) then
                    latticePoints(i, j)%phi = 0
                    latticePoints(i, j)%omega = 2.0 * latticePoints(i - 1, j)%phi
                else if (j == 1) then
                    latticePoints(i, j)%phi = 0
                    latticePoints(i, j)%omega = 2.0 * latticePoints(i , j + 1)%phi  
                else if (j == latticeSizeY) then 
                    latticePoints(i, j)%phi = 0
                    latticePoints(i, j)%omega = 2.0 * (latticePoints(i , j - 1)%phi - h)
                else
                    latticePoints(i, j)%omega = &
                        (latticePoints(i-1, j)%omega + &
                        latticePoints(i+1, j)%omega + &
                        latticePoints(i, j-1)%omega + &
                        latticePoints(i, j+1)%omega) / 4.0 + &
                        (latticePoints(i-1, j)%omega - latticePoints(i+1, j)%omega) * &
                        (latticePoints(i, j-1)%phi - latticePoints(i, j+1)%phi) - &
                        (latticePoints(i-1, j)%phi - latticePoints(i+1, j)%phi) * &
                        (latticePoints(i, j-1)%omega - latticePoints(i, j+1)%omega) * &
                        reynolds / 16.0
                    
                    latticePoints(i, j)%phi = &
                        (latticePoints(i-1, j)%phi + &
                        latticePoints(i+1, j)%phi + &
                        latticePoints(i, j+1)%phi + &
                        latticePoints(i, j-1)%phi - &
                        latticePoints(i, j)%omega) / 4.0
                end if
                shouldContinue = shouldContinue.or. &
                    isUnconverged(latticePoints(i, j)%omega, lastLaticePoint%omega).or. &
                    isUnconverged(latticePoints(i, j)%phi, lastLaticePoint%phi)
            end do
        end do
    end do

    !v = - h ^ 2 * omegaより変換 
    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            latticePoints(i, j)%omega = - latticePoints(i, j)%omega / h ** 2.0
        end do
    end do


    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            if (i == 1 .or. i == latticeSizeX .or. j == 1) then 
                latticePoints(i, j)%u = 0.0
                latticePoints(i, j)%v = 0.0

            else if (j == latticeSizeY) then 
                latticePoints(i, j)%u = -1.0
                latticePoints(i, j)%v = 0.0
            else
                latticePoints(i, j)%u = &
                    (latticePoints(i, j+1)%phi - latticePoints(i, j-1)%phi) / &
                    h / 2
                latticePoints(i, j)%v = &
                    - (latticePoints(i+1, j)%phi - latticePoints(i-1, j)%phi) / &
                    h / 2
            end if
        end do
    end do

    print *, "loop count: ", loopCount
    
    open(10, file='cavity_phi.csv')
    do k = 1, latticeSizeY
        write(10, '(*(f0.10:" "))') latticePoints(:, k)%phi
    end do
    close(10)

    open(11, file='cavity_omega.csv')
    do l = 1, latticeSizeY
        write(12, '(*(f0.10:" "))') latticePoints(:, l)%omega
    end do
    close(11)

    open(12, file='cavity_u.csv')
    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            write (12,'(3e12.4)') i*h, j*h, latticePoints(i, j)%u
        end do
    end do
    close(12)

    open(13, file='cavity_v.csv')
    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            write (13,'(3e12.4)') i*h, j*h, latticePoints(i, j)%v
        end do
    end do
    close(13)

    open(14, file='cavity_velocity.csv')
    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            write (14,'(4e12.4)') i*h, j*h, latticePoints(i, j)%u, latticePoints(i, j)%v
        end do
    end do
    close(14)

    stop
    contains
    
    logical function isUnconverged(a, b)
        implicit none
        double precision :: a, b 
        isUnconverged = abs(a-b) > convergenceThreshold
        return
    end function isUnconverged

end program cavity
