program advection
    use utility 
    implicit none
    ! ----- global variables ----- !
    real (fp), dimension(:), allocatable, save :: x, t, u0, uN
    real (fp), dimension(:,:), allocatable, save :: u
    real (fp), save :: x0, xf, t0, tf, a, dx, dt, C
    integer, save :: i, j, N, tN

    ! -------------------------------- !
    !           Driver                 !
    ! -------------------------------- !
    ! General code implementation:
    ! (1) setup initial parameters
    call init_param(N, tN, x0, xf, t0, tf, a)
    ! (2) create grid cells 
    call grid_setup(x, t, dx, dt, C, N, tN, x0, xf, t0, tf, a)
    ! (3) set up boundary conditions
    !    uncomment bc_setup for sine IC
    !    uncomment disc_bc_setup for discontinous IC
    !call bc_setup(u0, uN, u, x)
    call disc_bc_setup(u0, uN, u, x)
    ! (4) perform numerical method
    call advection_FD(u, u0)
    ! (5) Write to output file 
    call output()

    ! -------------------------------------- !
    !             subroutines                !
    ! -------------------------------------- !
    contains 
        ! ------------------------------- !
        ! write output data to file       !
        ! ------------------------------- !
        subroutine output()
            implicit none
            open(20, file = "output.dat", status = "replace")
                do i = 1, tN
                    write(20,*) t(i), (u(i,j), j = 1, N)
                end do
            close(20)
        end subroutine 

        ! ------------------------------- !
        ! initial setup values            !
        ! ------------------------------- !
        subroutine init_param(N, tN, x0, xf, t0, tf, a)
            use utility 
            implicit none
            integer, intent(out) :: N, tN
            real (fp), intent(out) :: x0, xf, t0, tf, a

            ! ---- initial parameters ---- !
            N = 128      ! space resolution 
            tN = 1e4  ! time resolution 
            x0 = -1.0    ! left boundary  
            xf = 1.0    ! right boundary 
            t0 = 0.0    ! initial time 
            tf = 1.0    ! final time 
            a = 1.0     ! advection velocity
        end subroutine 

        ! ------------------------------- !
        ! generate grid points            !
        ! ------------------------------- !
        subroutine grid_setup(x, t, dx, dt, C, N, tN, x0, xf, t0, tf, a)
            use utility 
            implicit none
            integer, intent(in) :: N, tN
            real (fp), intent(in) :: x0, xf, t0, tf, a
            real (fp), dimension(:), allocatable, intent(out) :: x, t
            real (fp), intent(out) :: dx, dt, C

            ! allocate space and time arrays 
            allocate(x(N+4), t(tN))
            x  = 0.0 ! initializing spacial grid 
            t  = 0.0 ! initializing temporal  grid 

            ! setting dx and dt 
            dx = (xf - x0) / N
            dt = (tf - t0) / tN 

            ! spacial grid points 
            x(1) = x0 - 1.5 * dx    ! add ghost cell start
            x(2) = x0 - 0.5 * dx
            x(N+3) = xf + 0.5 * dx
            x(N+4) = xf + 1.5 * dx  ! add ghost cell end
            do i = 3, N + 2             ! interior points 
                x(i) = x0 + (i - 0.5) * dx
            end do 

            !write(*,*) (x(i), i = 1,n)

            ! temporal grid points 
            do i = 1, tN
                t(i) = i * dt
            end do 

            C = a * dt / (2.0 * dx) 
            !C = 0.8
        end subroutine 

        ! ------------------------------- !
        ! Boundary conditions             !
        ! ------------------------------- !
        subroutine bc_setup(u0, uN, u, x)
            use utility
            implicit none
            real (fp), dimension(:), allocatable, intent(in) :: x
            real (fp), dimension(:), allocatable, intent(out) :: u0, uN
            real (fp), dimension(:,:), allocatable, intent(out) :: u 

            ! allocate solution arrays and initialize 
            allocate(u0(size(x)), uN(size(x)), u(tN,size(x)))
            u0 = sin(2.0 * pi * x)  ! initial conditions 
            uN = 0.0                ! initializing timestep 
            u  = 0.0                ! initializing solution
        end subroutine

        subroutine disc_bc_setup(u0, uN, u, x)
            use utility
            implicit none
            real (fp), dimension(:), allocatable, intent(in) :: x
            real (fp), dimension(:), allocatable, intent(out) :: u0, uN
            real (fp), dimension(:,:), allocatable, intent(out) :: u 

            ! allocate solution arrays and initialize 
            allocate(u0(size(x)), uN(size(x)), u(tN,size(x)))

            ! setting discontinious initial condition 
            do i = 1, N
                if (x(i) < 1.0/3.0 .AND. x(i) > -1.0/3.0) then
                    u0(i) = 1.0
                else 
                    u0(i) = 0.0
                end if
            end do 

            !write(*,*) (u0(j), j = 1,N)

            uN = u0              ! initializing timestep 
            u  = 0.0             ! initializing solution
        end subroutine

        ! ------------------------------- !
        ! Solving advection equation      !
        ! u_t - au_x = 0                  !
        ! ------------------------------- !
        subroutine advection_FD(u, u0)  
            use utility
            implicit none
            real (fp), dimension(:), intent(inout) :: u0
            real (fp), dimension(:,:), allocatable, intent(inout) :: u
            ! solution with Lax-Friedrich's method 
            ! update finite central difference in time
            do j = 1, size(t)
                ! from advection finite difference 
                !do i = 2, N - 1
                !    uN(i) = u0(i) - 0.5 * C * (u0(i+1) - u0(i-1))
                !end do
                do i = 3, N + 2
                    uN(i) = 0.5 * (u0(i+1) + u0(i-1)) - 0.5 * C * (u0(i+1) - u0(i-1))
                end do
                u0 = uN             ! updating interior points 
                u0(1) = u0(N+1)     ! periodic boundary condition
                u0(2) = u0(N+2)     ! periodic boundary condition
                u0(N+3) = u0(3)
                u0(N+4) = u0(4)
                u(j,:) = uN       ! save time series solution
            end do
        end subroutine

        subroutine advection_LF(u, u0)  
            use utility
            implicit none
            real (fp), dimension(:), intent(inout) :: u0
            real (fp), dimension(:,:), allocatable, intent(inout) :: u
            ! solution with Lax-Friedrich's method 
            ! update finite central difference in time
            do j = 1, size(t)
                ! from advection finite difference 
                !do i = 2, size(u0) - 1
                !    uN(i) = u0(i) - 0.5 * C * (u0(i+1) - u0(i-1))
                !end do 
                do i = 3, N + 2
                    uN(i) = 0.5 * (u0(i+1) + u0(i-1)) - 0.5 * C * (u0(i+1) - u0(i-1))
                end do
                u0 = uN             ! updating interior points 
                u0(1) = u0(N+1)     ! periodic boundary condition
                u0(2) = u0(N+2)     ! periodic boundary condition
                u0(N+3) = u0(3)
                u0(N+4) = u0(4)
                u(j,:) = uN     ! save time series solution
            end do
        end subroutine

        subroutine advection_LF(u, u0)  
            use utility
            implicit none
            real (fp), dimension(:), intent(inout) :: u0
            real (fp), dimension(:,:), allocatable, intent(inout) :: u
            ! solution with Lax-Friedrich's method 
            do j = 1, size(t)
                do i = 3, N + 2
                    uN(i) = u0(i) - (C / 2.0) * (u0(i+1) - u0(i-1)) + (C ** 2 / 2.0) * (u0(i+1) - 2.0 * u0(i) + u0(i-1))
                end do
                u0 = uN            ! updating interior points 
                u0 = uN             ! updating interior points 
                u0(1) = u0(N+1)     ! periodic boundary condition
                u0(2) = u0(N+2)     ! periodic boundary condition
                u0(N+3) = u0(3)
                u0(N+4) = u0(4)
                u(j,:) = uN        ! save time series solution
            end do
        end subroutine

        subroutine advection_LW(u, u0)  
            use utility
            implicit none
            real (fp), dimension(:), intent(inout) :: u0
            real (fp), dimension(:,:), allocatable, intent(inout) :: u
            ! solution with Lax-Wendroff's method 
            ! update finite central difference in time
            do j = 1, size(t)
                do i = 3, N + 2
                    uN(i) = u0(i) - (C / 2.0) * (u0(i+1) - u0(i-1)) + (C ** 2 / 2.0) * (u0(i+1) - 2.0 * u0(i) + u0(i-1))
                end do
                u0 = uN             ! updating interior points 
                u0(1) = u0(N+1)     ! periodic boundary condition
                u0(2) = u0(N+2)     ! periodic boundary condition
                u0(N+3) = u0(3)
                u0(N+4) = u0(4)
                u(j,:) = uN     ! save time series solution
            end do
        end subroutine
end program 
