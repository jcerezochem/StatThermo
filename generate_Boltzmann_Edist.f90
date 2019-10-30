program generate

    implicit none

    real(8),parameter :: pi = 3.141592d0 , &
                         cm1tokcalmol = 2.8591d-3, &
                         cm1tokJmol   = 1.1963d-2
    ! Allow Boltzmann constant to change (units), so not a parameter
    real(8)           :: kb = 0.695 ! cm-1,

    ! Parameters of the system
    real(8) :: T
    integer :: N

    ! Parameters of the representation
    integer :: Npoints
    real(8) :: Eini, Efin, DE, lnP
    real(8),dimension(:),allocatable :: E,P
    real(8) :: Norm
    logical :: do_gaussian
    logical :: do_ln

    ! Parameters of the distribution
    real(8) :: av, sgm, sgm2
    
    ! Management of units
    character(len=50) :: units

    !Input selections stuff
    logical :: argument_retrieved
    character(len=50) :: arg

    ! Counters
    integer :: i

    ! Defaults
    T       = 298.15d0
    N       = 3
    Eini    = -1.d0
    Efin    = -1.d0
    Npoints = 5000
    do_gaussian = .false.
    do_ln   = .false.
    units   = 'cm-1'

    !Read options from command line
    argument_retrieved=.false.
    do i=1,iargc()
        if (argument_retrieved) then
            argument_retrieved=.false.
            cycle
        endif
        call getarg(i, arg)
        select case (adjustl(arg))
            case ("-T")
                call getarg(i+1, arg)
                read(arg,*) T
                argument_retrieved=.true.
            case ("-N")
                call getarg(i+1, arg)
                read(arg,*) N
                argument_retrieved=.true.
            case ("-Npoints")
                call getarg(i+1, arg)
                read(arg,*) Npoints
                argument_retrieved=.true.
            case ("-Eini")
                call getarg(i+1, arg)
                read(arg,*) Eini
                argument_retrieved=.true.
            case ("-Efin")
                call getarg(i+1, arg)
                read(arg,*) Efin
                argument_retrieved=.true.
            case ("-units")
                call getarg(i+1, units)
                argument_retrieved=.true.
            case ("-gau")
                do_gaussian=.true.
            case ("-h")
                print*, "Program to generate equilibrium distributions of Energy"
                print*, "Options"
                print*, "  -T          Temperature (default: 298.15)"
                print*, "  -N          Number of degrees of freedom (default: 3)"
                print*, "  -Npoints    Number of degrees of points to represent"
!                print*, "  -Eini       Initial energy value (default: 0)"
!                print*, "  -Efin       Final energy value (default: 3*NkBT)"
                print*, "  -gau        Use Gaussian disribution"
                print*, "  -units      Units to output energy (cm-1|kJmol|kcalmol)"
                print*, "              (default: cm-1)"
                print*, ""
                stop
            case default
                write(0,*) "Unknown label ignored:", trim(adjustl(arg))
        end select
    enddo
    
    ! Set units
    write(0,*) "Selected units: "//trim(adjustl(units)) 
    if (units == 'cm-1') then
        write(0,*) "Using: cm-1"
        !kB = kB
    else if (units == 'kJmol') then
        write(0,*) "Using: kJ/mol"
        kB = kB * cm1tokJmol
    else if (units == 'kcalmol') then
        write(0,*) "Using: kcal/mol"
        kB = kB * cm1tokcalmol
    else
        write(0,*) "ERROR: Unknown units: "//trim(adjustl(units)) 
        stop
    endif

    ! Parameters of the distribution
    sgm2= 0.5d0*dfloat(N)*(kB*T)**2
    sgm = dsqrt(sgm2)
    av  = 0.5d0*dfloat(N)*kB*T
    ! Get default Efin/Eini
    if (Efin==-1.d0) then
        Efin = av + 6.d0*sgm
    endif
    if (Eini==-1.d0) then
        Eini = av - 6.d0*sgm
        Eini = max(0.d0,Eini)
    endif
        
    ! Compute DE
    DE = (Efin-Eini)/dfloat(Npoints-1)

    ! Detect errors
    if (N<1) then
        print*, "ERROR: N must be a positive integer"
        stop
    endif
    if (Eini<0.d0) then
        print*, "ERROR: Eini must be >=0.0"
        stop
    endif

    ! Avoid numerical inestabilities
    ! by preventing E=0 with N=1
    if (N<2 .and. Eini == 0.d0 ) then
        Eini = Eini+DE
    endif

    ! Prepare arrays
    allocate(E(Npoints),P(Npoints))

    ! Generate unnormalized distribution and then normalize the data
    ! since analytical normalization factor are numerically challenging
    do i=1,Npoints
        E(i) = Eini + dfloat(i-1) * DE
    enddo
    if (N>200) do_ln=.true.
    if (do_gaussian) then
        do i=1,Npoints
            P(i) = dexp(-(E(i)-av)**2/2.d0/sgm2)
        enddo
    elseif (do_ln) then
        do i=1,Npoints
            lnP  = dfloat(N-2)/2.d0*dlog(E(i)/kB**2/T**2) - E(i)/kB/T
            P(i) = dexp(lnP)
        enddo
    else
        do i=1,Npoints
            P(i) = dsqrt(E(i)/kB/T)**(N-2)/(kB*T) * dexp(-E(i)/(kB*T))
        enddo
    endif
    
    ! Compute normalization (integration with trapezoids)
    Norm=0.d0
    do i=1,Npoints-1
        Norm = Norm + (E(i+1)-E(i))*(P(i+1)+P(i))/2.d0
    enddo

    ! Print
    do i=1,Npoints
        print*, E(i), P(i)/Norm
    enddo

    stop 

end program generate
