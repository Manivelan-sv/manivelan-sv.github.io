! Dynamical instabilities cause extreme events in a theoretical Brusselator model.
! Manivelan, S. V., Sabarathinam, S., Thamilmaran, K., & Manimehan, I.
! Chaos, Solitons & Fractals, 180, 114582. (https://doi.org/10.1016/j.chaos.2024.114582)

program main
    implicit none
    integer :: i, j, nn
    double precision :: t, h, pi

    integer, parameter :: M = 150
    integer, parameter :: MAX_LENGTH = 1000000

    double precision :: omega, A, B, f, y(10)

    nn = 3

    ! Open the file for writing bifurcation data
    open(unit=10, file="bif_b.dat")
    
    pi = 4.0 * atan(1.0)

    ! Loop over different values of A
    A = 0.001
    do while (A < 0.3)
        ! Loop over different values of B
        B = 1.11
        do while (B < 1.13)
            print *, "B =", B

            f = 0.06
            omega = 0.7
            h = (2.0 * pi) / (omega * 1000)

            y(1) = 0.01
            y(2) = 0.02

            double precision :: data(MAX_LENGTH)
            integer :: length

            length = 0

            ! Perform numerical integration using RK4 method
            do j = 1, M
                do i = 1, 1000
                    t = h * dble(i)
                    call RK4(nn, h, t, y, DGL)
                end do

                ! Write data to file and store data for statistical analysis
                if (j > (M - 60)) then
                    write(10, "(4(F20.12))") A, B, y(1), y(2)
                    if (length < MAX_LENGTH) then
                        length = length + 1
                        data(length) = y(2)
                    end if
                end if
            end do

            ! Calculate statistical thresholds
            double precision :: standard_deviation, mean, threshold, threshold1, threshold2
            standard_deviation = calculate_standard_deviation(data, length)
            mean = calculate_mean(data, length)
            threshold = calculate_threshold(mean, standard_deviation)
            threshold1 = calculate_threshold1(mean, standard_deviation)
            threshold2 = calculate_threshold2(mean, standard_deviation)

            character(len=10) :: file_name

            ! Determine file name based on thresholds
            if (y(2) > threshold .and. y(2) < threshold1) then
                write(file_name, "(I0.0,'.dat')") 4
            else if (y(2) > threshold1 .and. y(2) < threshold2) then
                write(file_name, "(I0.0,'.dat')") 8
            else if (y(2) > threshold2) then
                write(file_name, "(I0.0,'.dat')") 12
            end if

            ! Open file for appending and write data
            open(unit=20, file=file_name, action="write", position="append")
            write(20, "(4(F20.12))") A, B, y(1), y(2)
            close(20)

            B = B + 0.01
        end do
        A = A + 0.01
    end do

    ! Close the file
    close(10)

contains

    subroutine RK4(nn, h, t, y, DGL)
        implicit none
        integer, intent(in) :: nn
        double precision, intent(in) :: h, t
        double precision, intent(inout) :: y(nn)
        double precision :: k1(10), k2(10), k3(10), k4(10)
        double precision :: yaux(10)

        ! Perform RK4 integration
        call DGL(t, y, k1)
        do i = 1, nn
            yaux(i) = y(i) + h * k1(i) / 2.0
        end do
        call DGL(t + h / 2.0, yaux, k2)
        do i = 1, nn
            yaux(i) = y(i) + h * k2(i) / 2.0
        end do
        call DGL(t + h / 2.0, yaux, k3)
        do i = 1, nn
            yaux(i) = y(i) + h * k3(i)
        end do
        call DGL(t + h, yaux, k4)
        do i = 1, nn
            y(i) = y(i) + h * (k1(i) / 6.0 + k2(i) / 3.0 + k3(i) / 3.0 + k4(i) / 6.0)
        end do
    end subroutine

    subroutine DGL(t, y, F)
        implicit none
        double precision, intent(in) :: t
        double precision, intent(in) :: y(10)
        double precision, intent(out) :: F(10)
        double precision :: A, B, f, omega

        ! Define the system of differential equations
        A = 1.01
        B = 0.2
        f = 0.06
        omega = 0.7

        F(1) = A - B * y(1) - y(1) + y(1) * y(1) * y(2) + f * sin(y(3))
        F(2) = B * y(1) - y(1) * y(1) * y(2)
        F(3) = omega
    end subroutine

    ! Function to calculate standard deviation
    double precision function calculate_standard_deviation(data, length)
        implicit none
        double precision, intent(in) :: data(length)
        integer, intent(in) :: length
        double precision :: sum, mean, squared_diff_sum, variance

        ! Calculate the mean
        sum = 0.0
        mean = 0.0
        do i = 1, length
            sum = sum + data(i)
        end do
        mean = sum / length

        ! Calculate the variance
        squared_diff_sum = 0.0
        do i = 1, length
            squared_diff_sum = squared_diff_sum + (data(i) - mean)**2
        end do
        variance = squared_diff_sum / length

        ! Return the standard deviation
        calculate_standard_deviation = sqrt(variance)
    end function

    ! Function to calculate mean
    double precision function calculate_mean(data, length)
        implicit none
        double precision, intent(in) :: data(length)
        integer, intent(in) :: length
        double precision :: sum

        ! Calculate the sum
        sum = 0.0
        do i = 1, length
            sum = sum + data(i)
        end do

        ! Return the mean
        calculate_mean = sum / length
    end function

    ! Function to calculate threshold
    double precision function calculate_threshold(mean, standard_deviation)
        implicit none
        double precision, intent(in) :: mean, standard_deviation

        ! Calculate and return the threshold
        calculate_threshold = mean + (4 * standard_deviation)
    end function

    ! Function to calculate threshold1
    double precision function calculate_threshold1(mean, standard_deviation)
        implicit none
        double precision, intent(in) :: mean, standard_deviation

        ! Calculate and return the threshold1
        calculate_threshold1 = mean + (8 * standard_deviation)
    end function

    ! Function to calculate threshold2
    double precision function calculate_threshold2(mean, standard_deviation)
        implicit none
        double precision, intent(in) :: mean, standard_deviation

        ! Calculate and return the threshold2
        calculate_threshold2 = mean + (12 * standard_deviation)
    end function

end program
