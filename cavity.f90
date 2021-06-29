program cavity
    implicit none

    type :: latticePoint
        double precision psi, omega
        double precision u, v, D
        double precision p
    end type latticePoint

    integer, parameter :: latticeSizeX = 51
    integer, parameter :: latticeSizeY = 51

    double precision, parameter :: reynolds = 100

    double precision :: h = 1.0 / (latticeSizeX-1)
    double precision :: b = 0.0
    
    type(latticePoint) :: latticePoints(latticeSizeX, latticeSizeY)
    type(latticePoint) :: lastLatticePoints(latticeSizeX, latticeSizeY)

    integer :: i, j, k, l
    integer :: loopCount = 0

    logical :: shouldContinue = .true.

    double precision, parameter :: convergenceThreshold = 0.00001 !収束条件

    latticePoints = latticePoint(psi=0.0, omega=0.0, u=0.0, v=0.0, D=0.0, p=0.0)
    lastLatticePoints = latticePoint(psi=0.0, omega=0.0, u=0.0, v=0.0, D=0.0, p=0.0)


    !lastLatticePoints使わなくていいですというか使うと収束が遅くなります
    do while(shouldContinue)
        shouldContinue = .false.
        loopCount = loopCount + 1
        do i = 1, latticeSizeX
            do j = 1, latticeSizeY
                lastLatticePoints = latticePoints
                if (i == 1) then 
                    latticePoints(i, j)%psi = 0
                    latticePoints(i, j)%omega = 2.0 * lastLatticePoints(i + 1, j)%psi
                else if (i == latticeSizeX) then
                    latticePoints(i, j)%psi = 0
                    latticePoints(i, j)%omega = 2.0 * lastLatticePoints(i - 1, j)%psi
                else if (j == 1) then
                    latticePoints(i, j)%psi = 0
                    latticePoints(i, j)%omega = 2.0 * lastLatticePoints(i , j + 1)%psi  
                else if (j == latticeSizeY) then 
                    latticePoints(i, j)%psi = 0
                    latticePoints(i, j)%omega = 2.0 * (lastLatticePoints(i , j - 1)%psi - h)
                else
                    latticePoints(i, j)%omega = &
                        (lastLatticePoints(i-1, j)%omega + &
                        lastLatticePoints(i+1, j)%omega + &
                        lastLatticePoints(i, j-1)%omega + &
                        lastLatticePoints(i, j+1)%omega) / 4.0 + &
                        ((lastLatticePoints(i+1, j)%psi - lastLatticePoints(i-1, j)%psi) * &
                        (lastLatticePoints(i, j+1)%omega - lastLatticePoints(i, j-1)%omega) - &
                        (lastLatticePoints(i+1, j)%omega - lastLatticePoints(i-1, j)%omega) * &
                        (lastLatticePoints(i, j+1)%psi - lastLatticePoints(i, j-1)%psi)) * &
                        reynolds / 16.0

                    latticePoints(i, j)%psi = &
                        (lastLatticePoints(i-1, j)%psi + &
                        lastLatticePoints(i+1, j)%psi + &
                        lastLatticePoints(i, j+1)%psi + &
                        lastLatticePoints(i, j-1)%psi - &
                        lastLatticePoints(i, j)%omega) / 4.0
                end if
                shouldContinue = shouldContinue.or. &
                    isUnconverged(latticePoints(i, j)%omega, lastLatticePoints(i,j)%omega).or. &
                    isUnconverged(latticePoints(i, j)%psi, lastLatticePoints(i,j)%psi)
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
                    (latticePoints(i, j+1)%psi - latticePoints(i, j-1)%psi) / &
                    h / 2
                latticePoints(i, j)%v = &
                    - (latticePoints(i+1, j)%psi - latticePoints(i-1, j)%psi) / &
                    h / 2
            end if
        end do
    end do

    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            if (i == 1 .or. i == latticeSizeX .or. j == 1) then 
                latticePoints(i, j)%D = 0.0
            else if (j == latticeSizeY) then 
                latticePoints(i, j)%D = 0.0
            else
                latticePoints(i, j)%D = &
                    (latticePoints(i+1, j)%u - latticePoints(i-1, j)%u) / (2 * h) + &
                    (latticePoints(i, j+1)%v - latticePoints(i, j-1)%v) / (2 * h)
            end if
        end do
    end do
    print *, "loop count: ", loopCount

    loopCount = 0
    latticePoints(1, 1)%p = 0
    latticePoints(1, 2)%p = 0
    latticePoints(2, 1)%p = 0

    do j = 2, latticeSizeY-1
        latticePoints(2, j)%p = latticePoints(2,j-1)%p + &
        (latticePoints(3, j)%v +latticePoints(2, j+1)%v + &
        latticePoints(1, j)%v+latticePoints(2, j-1)%v - 4 * latticePoints(2, j)%v) / (reynolds * h) - &
        latticePoints(2,j)%u * (latticePoints(3, j)%v-latticePoints(1,j)%v) /2 - &
        latticePoints(2,j)%v * (latticePoints(2, j+1)%v-latticePoints(2,j-1)%v) /2
    end do

        do j = 2, latticeSizeY-1
            do i = 3, latticeSizeX-1
                    latticePoints(i, j)%p = latticePoints(i-1,j)%p + &
                    (latticePoints(i+1, j)%u +latticePoints(i, j+1)%u + &
                    latticePoints(i-1, j)%u+latticePoints(i, j-1)%u - 4 * latticePoints(i, j)%u) / (reynolds * h) - &
                    latticePoints(i,j)%u * (latticePoints(i+1, j)%u-latticePoints(i-1,j)%u) /2 - &
                    latticePoints(i,j)%v * (latticePoints(i, j+1)%u-latticePoints(i,j-1)%u) /2
            end do
        end do
    
    do j = 1, latticeSizeY
        latticePoints(1, j)%p = latticePoints(2, j)%p
        latticePoints(latticeSizeX, j)%p = latticePoints(latticeSizeX-1, j)%p
    end do
    do i = 1, latticeSizeY
        latticePoints(i, 1)%p = latticePoints(i, 2)%p
        latticePoints(i, latticeSizeY)%p = latticePoints(i, latticeSizeY-1)%p
    end do

    print *, "loop count: ", loopCount
    
    open(10, file='outputs/data/psi.csv')
    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            write (10,'(3e12.4)') (i-1)*h, (j-1)*h, latticePoints(i, j)%psi
        end do
        write(10,'(/)',advance='no')
    end do
    close(10)

    open(11, file='outputs/data/omega.csv')
    do l = 1, latticeSizeY
        write(11, '(*(f0.10:" "))') latticePoints(:, l)%omega
    end do
    close(11)

    open(12, file='outputs/data/u.csv')
    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            write (12,'(3e12.4)') (i-1)*h, (j-1)*h, latticePoints(i, j)%u
        end do
    end do
    close(12)

    open(13, file='outputs/data/v.csv')
    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            write (13,'(3e12.4)') (i-1)*h, (j-1)*h, latticePoints(i, j)%v
        end do
    end do
    close(13)

    open(14, file='outputs/data/velocity.csv')
    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            if (mod(i,2) == 1 .and. mod(j,2) == 1 ) then
                write (14,'(4e12.4)') (i-1)*h, (j-1)*h, latticePoints(i, j)%u, latticePoints(i, j)%v
            end if
        end do
    end do
    close(14)

    open(15, file='outputs/data/p.csv')
    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            write (15,'(3e12.4)') (i-1)*h, (j-1)*h, latticePoints(i, j)%p
        end do
        write(15,'(/)',advance='no')
    end do
    close(15)

    open(16, file='outputs/data/D.csv')
    do i = 1, latticeSizeX
        do j = 1, latticeSizeY
            write (16,'(3e12.4)') (i-1)*h, (j-1)*h, latticePoints(i, j)%D
        end do
        write(16,'(/)',advance='no')
    end do
    close(16)

    stop
    contains
    
    logical function isUnconverged(a, b)
        implicit none
        double precision :: a, b 
        isUnconverged = abs(a-b) > convergenceThreshold
        return
    end function isUnconverged

end program cavity
