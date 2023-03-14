!$UWHPSC/codes/fortran/newton/newton.f90
!http://faculty.washington.edu/rjl/classes/am583s2013/notes/fortran_newton.html
module h2sc_solvers
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_kind_mod , only : r4 => SHR_KIND_R4
  ! module parameters:
  implicit none
  type, public::seepage_argument
     real(kind = r8)::dFlow_downslope
     real(kind = r8)::dFlow_seepage
     real(kind = r8)::dLength_seepage
     real(kind = r8)::dLength_water_table
     real(kind = r8)::dSlope_surface
     real(kind = r8)::dSlope_water_table
     real(kind = r8)::beta_seepage
     real(kind = r8)::m_seepage
     real(kind = r8)::n_seepage
     real(kind = r8)::xBeg
     real(kind = r8)::xEnd
  end type seepage_argument

  !type(seepage_argument):: cSeepage_argu_in
contains
  subroutine newton_method_seepage(x,xNew,cSeepage_argu_in, fNew,residual)
    implicit none
    real (kind = r8), intent(IN) :: x
    real (kind = r8), intent(OUT) :: fNew, residual

    real (kind = r8) :: xNew, f, fprime
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    !! compute function value evaluated at x
    f = gx(x, cSeepage_argu_in)
    !! compute function derivative evaluated at x
    fprime = gx_prime(x, cSeepage_argu_in)
    !! Exit if f' is near or become zero
    if (abs(fprime) < 1.e-12) then
       print *, '[Error: newton_method] Function derivative becomes very close to zero or zero.'
       print *, 'f=',f, 'df/dx =',fprime
       print *, 'Aborting now in order to avoid division by zero.'
       stop
    end if
    !! Algorithm
    xNew = x - f/fprime * 2.0
    fNew = f
    !! Search fails if a newly updated value x is out of the search domain
    if ((xNew < cSeepage_argu_in%xBeg) .or. (xNew > cSeepage_argu_in%xEnd)) then
       print *, '[Error: newton_method] xNew',xNew, 'is out of domain.'
       print *, 'Failed in search. Aborting now.'
       stop
    end if
    !! Calculate a new residual
    residual = abs(xNew - x)
  end subroutine newton_method_seepage

  !===============================================
  !the gx function in equation 29
  !===============================================
  function gx( x , cSeepage_argu_in) result(gx_out)
    implicit none
    !===============================================
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    real(kind = r8)::x
    real(kind = r8)::gx_out
    !===============================================
    real(kind = r8)::dFlow_downslope
    real(kind = r8)::dFlow_seepage
    real(kind = r8)::dLength_seepage
    real(kind = r8)::dLength_water_table
    real(kind = r8)::dSlope_surface
    real(kind = r8)::dSlope_water_table
    real(kind = r8)::beta_seepage
    real(kind = r8)::m_seepage
    real(kind = r8)::n_seepage
    real(kind = r8)::phi2 = 0.0
    real(kind = r8)::dummy0, dummy1, dummy2, dummy3
    real(kind = r8)::dummy4, dummy5, dummy6, dummy7
    !===============================================
    dFlow_seepage= cSeepage_argu_in%dFlow_seepage
    dLength_seepage = cSeepage_argu_in%dLength_seepage
    !===============================================
    dFlow_downslope = cSeepage_argu_in%dFlow_downslope
    dLength_water_table = cSeepage_argu_in%dLength_water_table
    dSlope_surface = cSeepage_argu_in%dSlope_surface
    dSlope_water_table = cSeepage_argu_in%dSlope_water_table
    beta_seepage=cSeepage_argu_in%beta_seepage
    m_seepage=cSeepage_argu_in%m_seepage
    n_seepage=cSeepage_argu_in%n_seepage
    !===============================================
    dummy0 = dLength_water_table * ( x - dLength_seepage)
    dummy1 = dummy0 - (x**2)/2 &
         + (dLength_seepage)**2/2
    dummy2 = tan(dSlope_surface) - tan(dSlope_water_table)
    dummy3 = beta_seepage * ( x - dLength_seepage )/2 &
         *dummy2
    dummy4 = 1 - (1 + dummy3**n_seepage )**(-m_seepage)
    dummy5 = (dFlow_seepage) / (dFlow_seepage + dFlow_downslope)
    dummy6 = ( dLength_water_table - dLength_seepage)**2/2 - phi2
    dummy7 = dummy5 * dummy6
    gx_out = dummy1 * dummy4 - dummy7
    return
  end function gx
  !===============================================
  !the derivative of gx, this is used for Newton's method
  !===============================================
  function gx_prime( x , cSeepage_argu_in) result(gx_prime_out)
    implicit none
    !===============================================
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    real(kind = r8)::x
    real(kind = r8)::gx_prime_out
    !===============================================
    real(kind = r8)::dFlow_downslope
    real(kind = r8)::dFlow_seepage
    real(kind = r8)::dLength_seepage
    real(kind = r8)::dLength_water_table
    real(kind = r8)::dSlope_surface
    real(kind = r8)::dSlope_water_table
    real(kind = r8)::beta_seepage
    real(kind = r8)::m_seepage
    real(kind = r8)::n_seepage
    real(kind = r8)::phi2 = 0.0
    real(kind = r8)::dummy0, dummy1, dummy2, dummy3
    real(kind = r8)::dummy4, dummy5, dummy6, dummy7
    real(kind = r8)::dummy1_prime, dummy3_prime, dummy4_prime
    !===============================================
    dFlow_seepage=cSeepage_argu_in%dFlow_seepage
    dLength_seepage = cSeepage_argu_in%dLength_seepage
    !===============================================
    dFlow_downslope = cSeepage_argu_in%dFlow_downslope
    dLength_water_table = cSeepage_argu_in%dLength_water_table
    dSlope_surface = cSeepage_argu_in%dSlope_surface
    dSlope_water_table = cSeepage_argu_in%dSlope_water_table
    beta_seepage=cSeepage_argu_in%beta_seepage
    m_seepage=cSeepage_argu_in%m_seepage
    n_seepage=cSeepage_argu_in%n_seepage
    !===============================================
    dummy0 = dLength_water_table * ( x - dLength_seepage)
    dummy1 = dummy0 - (x**2)/2 &
         + (dLength_seepage)**2/2
    dummy1_prime = dLength_water_table - x
    dummy2 = tan(dSlope_surface) - tan(dSlope_water_table)
    dummy3 = beta_seepage * ( x - dLength_seepage )/2 &
         * dummy2
    dummy4 = 1 - (1 + dummy3**n_seepage )**(-m_seepage)
    dummy3_prime = beta_seepage / 2 * dummy2
    dummy4_prime = m_seepage * n_seepage * ( dummy3 )**(n_seepage-1) * dummy3_prime
    dummy5 = (dFlow_seepage) / (dFlow_seepage + dFlow_downslope)
    dummy6 = ( dLength_water_table - dLength_seepage)**2/2 - phi2
    dummy7 = dummy5 * dummy6
    gx_prime_out = dummy1_prime * dummy4 + dummy1 * dummy4_prime
    return
  end function gx_prime



  !===============================================
  !a wrapper
  !===============================================
  function f_phi_1( x , cSeepage_argu_in) &
       result(f_1)
    implicit none
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    !===============================================
    real(kind = r4) :: x
    real(kind = r4) :: f_1
    !===============================================
    real(kind = r8)::beta_seepage
    real(kind = r8)::m_seepage
    real(kind = r8)::n_seepage
    !===============================================
    real(kind = r8)::dSlope_surface
    real(kind = r8)::dSlope_water_table
    real(kind = r8)::dummy0, dummy1, dummy2, dummy3
    !===============================================
    beta_seepage=cSeepage_argu_in%beta_seepage
    m_seepage=cSeepage_argu_in%m_seepage
    n_seepage=cSeepage_argu_in%n_seepage
    dSlope_surface=cSeepage_argu_in%dSlope_surface
    dSlope_water_table = cSeepage_argu_in%dSlope_water_table
    !===============================================
    !some quality check
    !===============================================

    if (dSlope_water_table > dSlope_surface) then
       !write(*, *) 'Warning f_phi_1: the water table slope is higher than surface slope!'
    end if
    !===============================================
    dummy0 = tan(dSlope_surface) - tan(dSlope_water_table)
    dummy1 = beta_seepage * dummy0 * x !be careful with the x
    dummy2 = 1 + dummy1 **n_seepage
    dummy3 = dummy2 **(-m_seepage)
    f_1= x * dummy3
    return
  end function f_phi_1


  !===============================================
  !a wrapper
  !===============================================
  function f_phi_2(x, cSeepage_argu_in) &
       result(f_2)
    implicit none
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    !===============================================

    real(kind = r8) :: x
    real(kind = r4) :: f_2
    !===============================================
    real(kind = r8)::beta_seepage
    real(kind = r8)::m_seepage
    real(kind = r8)::n_seepage
    !===============================================
    real(kind = r8)::dLength_seepage
    !===============================================
    real(kind = r8)::dLength_water_table
    real(kind = r8)::dSlope_water_table
    real(kind = r8)::dSlope_surface
    !===============================================
    real(kind = r8)::dummy0, dummy1, dummy2, dummy3, dummy4
    !===============================================
    dLength_seepage = cSeepage_argu_in%dLength_seepage
    !===============================================
    beta_seepage=cSeepage_argu_in%beta_seepage
    m_seepage = cSeepage_argu_in%m_seepage
    n_seepage = cSeepage_argu_in%n_seepage
    dLength_water_table = cSeepage_argu_in%dLength_water_table
    dSlope_surface = cSeepage_argu_in%dSlope_surface
    dSlope_water_table = cSeepage_argu_in%dSlope_water_table

    !===============================================
    !some quality check
    !===============================================

    if (dSlope_water_table > dSlope_surface) then
       !write(*, *) 'Warning f_phi_2: the water table slope is higher than surface slope!'
    end if

    !===============================================
    dummy0 = dLength_water_table - x
    dummy1= tan(dSlope_surface) - tan(dSlope_water_table)
    dummy2 = x - dLength_seepage
    dummy3 = ( beta_seepage * dummy1 * dummy2 )**n_seepage
    dummy4 = (1 + dummy3) **(-m_seepage)
    f_2 = dummy0 * dummy4
    !===============================================
    return
  end function f_phi_2

  subroutine qag_seepage ( f, cSeepage_argu_in, a, b, epsabs, epsrel, key, result, abserr, neval, ier )
    !
    !**************************************************************************
    !
    !! QAG approximates an integral over a finite interval.
    !
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral
    !      I = integral of F over (A,B),
    !    hopefully satisfying
    !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !    QAG is a simple globally adaptive integrator using the strategy of
    !    Aind (Piessens, 1973).  It is possible to choose between 6 pairs of
    !    Gauss-Kronrod quadrature formulae for the rule evaluation component.
    !    The pairs of high degree of precision are suitable for handling
    !    integration difficulties due to a strongly oscillating integrand.
    !
    !  Reference:
    !
    !    R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real F, the name of the function routine, of the form
    !      function f ( x )
    !      real f
    !      real x
    !    which evaluates the integrand function.
    !
    !    Input, real A, B, the limits of integration.
    !
    !    Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Input, integer KEY, chooses the order of the local integration rule:
    !    1,  7 Gauss points, 15 Gauss-Kronrod points,
    !    2, 10 Gauss points, 21 Gauss-Kronrod points,
    !    3, 15 Gauss points, 31 Gauss-Kronrod points,
    !    4, 20 Gauss points, 41 Gauss-Kronrod points,
    !    5, 25 Gauss points, 51 Gauss-Kronrod points,
    !    6, 30 Gauss points, 61 Gauss-Kronrod points.
    !
    !    Output, real RESULT, the estimated value of the integral.
    !
    !    Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !    Output, integer IER, return code.
    !    0, normal and reliable termination of the routine.  It is assumed that the
    !      requested accuracy has been achieved.
    !    1, maximum number of subdivisions allowed has been achieved.  One can
    !      allow more subdivisions by increasing the value of LIMIT in QAG.
    !      However, if this yields no improvement it is advised to analyze the
    !      integrand to determine the integration difficulties.  If the position
    !      of a local difficulty can be determined, such as a singularity or
    !      discontinuity within the interval) one will probably gain from
    !      splitting up the interval at this point and calling the integrator
    !      on the subranges.  If possible, an appropriate special-purpose
    !      integrator should be used which is designed for handling the type
    !      of difficulty involved.
    !    2, the occurrence of roundoff error is detected, which prevents the
    !      requested tolerance from being achieved.
    !    3, extremely bad integrand behavior occurs at some points of the
    !      integration interval.
    !    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
    !
    !  Local parameters:
    !
    !    LIMIT is the maximum number of subintervals allowed in
    !    the subdivision process of QAGE.
    !
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    integer, parameter :: limit = 500
    !
    real a
    real,intent(out):: abserr
    real alist(limit)
    real b
    real blist(limit)
    real elist(limit)
    real epsabs
    real epsrel
    real, external :: f
    integer,intent(out):: ier
    integer iord(limit)
    integer key
    integer last
    !    integer limit
    integer,intent(out):: neval
    real,intent(out):: result
    real rlist(limit)
    !
    call qage_seepage ( f, cSeepage_argu_in, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
         ier, alist, blist, rlist, elist, iord, last )

    return
  end subroutine qag_seepage

  subroutine qage_seepage ( f,cSeepage_argu_in, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
       ier, alist, blist, rlist, elist, iord, last )
    !
    !**************************************************************************
    !
    !! QAGE estimates a definite integral.
    !
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral
    !      I = integral of F over (A,B),
    !    hopefully satisfying
    !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !  Reference:
    !
    !    R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real F, the name of the function routine, of the form
    !      function f ( x )
    !      real f
    !      real x
    !    which evaluates the integrand function.
    !
    !    Input, real A, B, the limits of integration.
    !
    !    Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Input, integer KEY, chooses the order of the local integration rule:
    !    1,  7 Gauss points, 15 Gauss-Kronrod points,
    !    2, 10 Gauss points, 21 Gauss-Kronrod points,
    !    3, 15 Gauss points, 31 Gauss-Kronrod points,
    !    4, 20 Gauss points, 41 Gauss-Kronrod points,
    !    5, 25 Gauss points, 51 Gauss-Kronrod points,
    !    6, 30 Gauss points, 61 Gauss-Kronrod points.
    !
    !    Input, integer LIMIT, the maximum number of subintervals that
    !    can be used.
    !
    !    Output, real RESULT, the estimated value of the integral.
    !
    !    Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !    Output, integer IER, return code.
    !    0, normal and reliable termination of the routine.  It is assumed that the
    !      requested accuracy has been achieved.
    !    1, maximum number of subdivisions allowed has been achieved.  One can
    !      allow more subdivisions by increasing the value of LIMIT in QAG.
    !      However, if this yields no improvement it is advised to analyze the
    !      integrand to determine the integration difficulties.  If the position
    !      of a local difficulty can be determined, such as a singularity or
    !      discontinuity within the interval) one will probably gain from
    !      splitting up the interval at this point and calling the integrator
    !      on the subranges.  If possible, an appropriate special-purpose
    !      integrator should be used which is designed for handling the type
    !      of difficulty involved.
    !    2, the occurrence of roundoff error is detected, which prevents the
    !      requested tolerance from being achieved.
    !    3, extremely bad integrand behavior occurs at some points of the
    !      integration interval.
    !    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
    !
    !    Workspace, real ALIST(LIMIT), BLIST(LIMIT), contains in entries 1
    !    through LAST the left and right ends of the partition subintervals.
    !
    !    Workspace, real RLIST(LIMIT), contains in entries 1 through LAST
    !    the integral approximations on the subintervals.
    !
    !    Workspace, real ELIST(LIMIT), contains in entries 1 through LAST
    !    the absolute error estimates on the subintervals.
    !
    !    Output, integer IORD(LIMIT), the first K elements of which are pointers
    !    to the error estimates over the subintervals, such that
    !    elist(iord(1)), ..., elist(iord(k)) form a decreasing sequence, with
    !    k = last if last <= (limit/2+2), and k = limit+1-last otherwise.
    !
    !    Output, integer LAST, the number of subintervals actually produced
    !    in the subdivision process.
    !
    !  Local parameters:
    !
    !    alist     - list of left end points of all subintervals
    !                       considered up to now
    !    blist     - list of right end points of all subintervals
    !                       considered up to now
    !    elist(i)  - error estimate applying to rlist(i)
    !    maxerr    - pointer to the interval with largest error estimate
    !    errmax    - elist(maxerr)
    !    area      - sum of the integrals over the subintervals
    !    errsum    - sum of the errors over the subintervals
    !    errbnd    - requested accuracy max(epsabs,epsrel*abs(result))
    !    *****1    - variable for the left subinterval
    !    *****2    - variable for the right subinterval
    !    last      - index for subdivision
    !
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    integer limit
    !
    real a
    real,intent(out):: abserr
    real,intent(out):: alist(limit)
    real area
    real area1
    real area12
    real area2
    real a1
    real a2
    real b
    real,intent(out):: blist(limit)
    real b1
    real b2
    real c
    real defabs
    real defab1
    real defab2
    real,intent(out):: elist(limit)
    real epsabs
    real epsrel
    real errbnd
    real errmax
    real error1
    real error2
    real erro12
    real errsum
    real, external :: f
    integer,intent(out):: ier
    integer,intent(out):: iord(limit)
    integer iroff1
    integer iroff2
    !integer k
    integer key
    integer keyf
    integer,intent(out):: last
    integer maxerr
    integer,intent(out):: neval
    integer nrmax
    real resabs
    real,intent(out):: result
    real,intent(out):: rlist(limit)
    !
    !  Test on validity of parameters.
    !
    ier = 0
    neval = 0
    last = 0
    result = 0.0e+00
    abserr = 0.0e+00
    alist(1) = a
    blist(1) = b
    rlist(1) = 0.0e+00
    elist(1) = 0.0e+00
    iord(1) = 0

    if ( epsabs < 0.0e+00 .and. epsrel < 0.0e+00 ) then
       ier = 6
       return
    end if
    !
    !  First approximation to the integral.
    !
    keyf = key
    keyf = max ( keyf, 1 )
    keyf = min ( keyf, 6 )

    c = keyf
    neval = 0

    if ( keyf == 1 ) then
       call qk15_seepage ( f, cSeepage_argu_in, a, b, result, abserr, defabs, resabs )
    else if ( keyf == 2 ) then
       call qk21_seepage ( f, cSeepage_argu_in, a, b, result, abserr, defabs, resabs )
    else if ( keyf == 3 ) then
       call qk31_seepage ( f, cSeepage_argu_in, a, b, result, abserr, defabs, resabs )
    else if ( keyf == 4 ) then
       call qk41_seepage ( f, cSeepage_argu_in, a, b, result, abserr, defabs, resabs )
    else if ( keyf == 5 ) then
       call qk51_seepage ( f, cSeepage_argu_in, a, b, result, abserr, defabs, resabs )
    else if ( keyf == 6 ) then
       call qk61_seepage ( f, cSeepage_argu_in, a, b, result, abserr, defabs, resabs )
    end if

    last = 1
    rlist(1) = result
    elist(1) = abserr
    iord(1) = 1
    !
    !  Test on accuracy.
    !
    errbnd = max ( epsabs, epsrel * abs ( result ) )

    if ( abserr <= 5.0e+01 * epsilon ( defabs ) * defabs .and. &
         abserr > errbnd ) then
       ier = 2
    end if

    if ( limit == 1 ) then
       ier = 1
    end if

    if ( ier /= 0 .or. &
         ( abserr <= errbnd .and. abserr /= resabs ) .or. &
         abserr == 0.0e+00 ) then

       if ( keyf /= 1 ) then
          neval = (10*keyf+1) * (2*neval+1)
       else
          neval = 30 * neval + 15
       end if

       return

    end if
    !
    !  Initialization.
    !
    errmax = abserr
    maxerr = 1
    area = result
    errsum = abserr
    nrmax = 1
    iroff1 = 0
    iroff2 = 0

    do last = 2, limit
       !
       !  Bisect the subinterval with the largest error estimate.
       !
       a1 = alist(maxerr)
       b1 = 0.5E+00 * ( alist(maxerr) + blist(maxerr) )
       a2 = b1
       b2 = blist(maxerr)

       if ( keyf == 1 ) then
          call qk15_seepage ( f, cSeepage_argu_in, a1, b1, area1, error1, resabs, defab1 )
       else if ( keyf == 2 ) then
          call qk21_seepage ( f, cSeepage_argu_in, a1, b1, area1, error1, resabs, defab1 )
       else if ( keyf == 3 ) then
          call qk31_seepage ( f, cSeepage_argu_in, a1, b1, area1, error1, resabs, defab1 )
       else if ( keyf == 4 ) then
          call qk41_seepage ( f, cSeepage_argu_in, a1, b1, area1, error1, resabs, defab1)
       else if ( keyf == 5 ) then
          call qk51_seepage ( f, cSeepage_argu_in, a1, b1, area1, error1, resabs, defab1 )
       else if ( keyf == 6 ) then
          call qk61_seepage ( f, cSeepage_argu_in, a1, b1, area1, error1, resabs, defab1 )
       end if

       if ( keyf == 1 ) then
          call qk15_seepage ( f, cSeepage_argu_in, a2, b2, area2, error2, resabs, defab2 )
       else if ( keyf == 2 ) then
          call qk21_seepage ( f, cSeepage_argu_in, a2, b2, area2, error2, resabs, defab2 )
       else if ( keyf == 3 ) then
          call qk31_seepage ( f, cSeepage_argu_in, a2, b2, area2, error2, resabs, defab2 )
       else if ( keyf == 4 ) then
          call qk41_seepage ( f, cSeepage_argu_in, a2, b2, area2, error2, resabs, defab2 )
       else if ( keyf == 5 ) then
          call qk51_seepage ( f, cSeepage_argu_in, a2, b2, area2, error2, resabs, defab2 )
       else if ( keyf == 6 ) then
          call qk61_seepage ( f, cSeepage_argu_in, a2, b2, area2, error2, resabs, defab2 )
       end if
       !
       !  Improve previous approximations to integral and error and
       !  test for accuracy.
       !
       neval = neval + 1
       area12 = area1 + area2
       erro12 = error1 + error2
       errsum = errsum + erro12 - errmax
       area = area + area12 - rlist(maxerr)

       if ( defab1 /= error1 .and. defab2 /= error2 ) then

          if ( abs ( rlist(maxerr) - area12 ) <= 1.0e-05 * abs ( area12 ) &
               .and. erro12 >= 9.9e-01 * errmax ) then
             iroff1 = iroff1 + 1
          end if

          if ( last > 10 .and. erro12 > errmax ) then
             iroff2 = iroff2 + 1
          end if

       end if

       rlist(maxerr) = area1
       rlist(last) = area2
       errbnd = max ( epsabs, epsrel * abs ( area ) )
       !
       !  Test for roundoff error and eventually set error flag.
       !
       if ( errsum > errbnd ) then

          if ( iroff1 >= 6 .or. iroff2 >= 20 ) then
             ier = 2
          end if
          !
          !  Set error flag in the case that the number of subintervals
          !  equals limit.
          !
          if ( last == limit ) then
             ier = 1
          end if
          !
          !  Set error flag in the case of bad integrand behavior
          !  at a point of the integration range.
          !
          if ( max ( abs ( a1 ), abs ( b2 ) ) <= ( 1.0e+00 + c * 1.0e+03 * &
               epsilon ( a1 ) ) * ( abs ( a2 ) + 1.0e+04 * tiny ( a2 ) ) ) then
             ier = 3
          end if

       end if
       !
       !  Append the newly-created intervals to the list.
       !
       if ( error2 <= error1 ) then
          alist(last) = a2
          blist(maxerr) = b1
          blist(last) = b2
          elist(maxerr) = error1
          elist(last) = error2
       else
          alist(maxerr) = a2
          alist(last) = a1
          blist(last) = b1
          rlist(maxerr) = area2
          rlist(last) = area1
          elist(maxerr) = error2
          elist(last) = error1
       end if
       !
       !  Call QSORT to maintain the descending ordering
       !  in the list of error estimates and select the subinterval
       !  with the largest error estimate (to be bisected next).
       !
       call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

       if ( ier /= 0 .or. errsum <= errbnd ) then
          exit
       end if

    end do
    !
    !  Compute final result.
    !
    result = sum ( rlist(1:last) )

    abserr = errsum

    if ( keyf /= 1 ) then
       neval = ( 10 * keyf + 1 ) * ( 2 * neval + 1 )
    else
       neval = 30 * neval + 15
    end if

    return
  end subroutine qage_seepage

  subroutine qk15_seepage ( f,cSeepage_argu_in, a, b, result, abserr, resabs, resasc )
    !
    !************************************************************************
    !
    !! QK15 carries out a 15 point Gauss-Kronrod quadrature rule.
    !
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Reference:
    !
    !    R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real F, the name of the function routine, of the form
    !      function f ( x )
    !      real f
    !      real x
    !    which evaluates the integrand function.
    !
    !    Input, real A, B, the limits of integration.
    !
    !    Output, real RESULT, the estimated value of the integral.
    !    RESULT is computed by applying the 15-point Kronrod rule (RESK)
    !    obtained by optimal addition of abscissae to the 7-point Gauss rule
    !    (RESG).
    !
    !    Output, real ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real RESASC, approximation to the integral | F-I/(B-A) |
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 15-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 7-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 7-point Gauss rule
    !
    !           wgk    - weights of the 15-point Kronrod rule
    !
    !           wg     - weights of the 7-point Gauss rule
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 7-point Gauss formula
    !           resk   - result of the 15-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    real a
    real absc
    real,intent(out):: abserr
    real b
    real centr
    real dhlgth
    real, external :: f
    real fc
    real fsum
    real fval1
    real fval2
    real fv1(7)
    real fv2(7)
    real hlgth
    integer j
    integer jtw
    integer jtwm1
    real,intent(out):: resabs
    real,intent(out):: resasc
    real resg
    real resk
    real reskh
    real,intent(out):: result
    real wg(4)
    real wgk(8)
    real xgk(8)
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
         9.914553711208126e-01,   9.491079123427585e-01, &
         8.648644233597691e-01,   7.415311855993944e-01, &
         5.860872354676911e-01,   4.058451513773972e-01, &
         2.077849550078985e-01,   0.0e+00              /
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
         2.293532201052922e-02,   6.309209262997855e-02, &
         1.047900103222502e-01,   1.406532597155259e-01, &
         1.690047266392679e-01,   1.903505780647854e-01, &
         2.044329400752989e-01,   2.094821410847278e-01/
    data wg(1),wg(2),wg(3),wg(4)/ &
         1.294849661688697e-01,   2.797053914892767e-01, &
         3.818300505051189e-01,   4.179591836734694e-01/
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 15-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
    fc = f(centr, cSeepage_argu_in)
    resg = fc*wg(4)
    resk = fc*wgk(8)
    resabs = abs(resk)

    do j = 1, 3
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 4
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk * 5.0e-01
    resasc = wgk(8)*abs(fc-reskh)

    do j = 1, 7
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00 ) then
       abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
    end if

    if ( resabs > tiny ( resabs ) / (5.0e+01* epsilon ( resabs ) ) ) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk15_seepage
  subroutine qk21_seepage ( f, cSeepage_argu_in, a, b, result, abserr, resabs, resasc )
    !
    !************************************************************************
    !
    !! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
    !
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Reference:
    !
    !    R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real F, the name of the function routine, of the form
    !      function f ( x )
    !      real f
    !      real x
    !    which evaluates the integrand function.
    !
    !    Input, real A, B, the limits of integration.
    !
    !    Output, real RESULT, the estimated value of the integral.
    !                       result is computed by applying the 21-point
    !                       Kronrod rule (resk) obtained by optimal addition
    !                       of abscissae to the 10-point Gauss rule (resg).
    !
    !    Output, real ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real RESASC, approximation to the integral | F-I/(B-A) |
    !    over [A,B].
    !
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    real a
    real absc
    real,intent(out):: abserr
    real b
    real centr
    real dhlgth
    real, external :: f
    real fc
    real fsum
    real fval1
    real fval2
    real fv1(10)
    real fv2(10)
    real hlgth
    integer j
    integer jtw
    integer jtwm1
    real,intent(out):: resabs
    real,intent(out):: resasc
    real resg
    real resk
    real reskh
    real,intent(out):: result
    real wg(5)
    real wgk(11)
    real xgk(11)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 21-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 10-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 10-point Gauss rule
    !
    !           wgk    - weights of the 21-point Kronrod rule
    !
    !           wg     - weights of the 10-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10),xgk(11)/ &
         9.956571630258081e-01,     9.739065285171717e-01, &
         9.301574913557082e-01,     8.650633666889845e-01, &
         7.808177265864169e-01,     6.794095682990244e-01, &
         5.627571346686047e-01,     4.333953941292472e-01, &
         2.943928627014602e-01,     1.488743389816312e-01, &
         0.000000000000000e+00/
    !
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
         wgk(9),wgk(10),wgk(11)/ &
         1.169463886737187e-02,     3.255816230796473e-02, &
         5.475589657435200e-02,     7.503967481091995e-02, &
         9.312545458369761e-02,     1.093871588022976e-01, &
         1.234919762620659e-01,     1.347092173114733e-01, &
         1.427759385770601e-01,     1.477391049013385e-01, &
         1.494455540029169e-01/
    !
    data wg(1),wg(2),wg(3),wg(4),wg(5)/ &
         6.667134430868814e-02,     1.494513491505806e-01, &
         2.190863625159820e-01,     2.692667193099964e-01, &
         2.955242247147529e-01/
    !
    !
    !           list of major variables
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 10-point Gauss formula
    !           resk   - result of the 21-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 21-point Kronrod approximation to the
    !  integral, and estimate the absolute error.
    !
    resg = 0.0e+00
    fc = f(centr, cSeepage_argu_in)
    resk = wgk(11)*fc
    resabs = abs(resk)

    do j = 1, 5
       jtw = 2*j
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 5
       jtwm1 = 2*j-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk*5.0e-01
    resasc = wgk(11)*abs(fc-reskh)

    do j = 1, 10
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) then
       abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
    end if

    if ( resabs > tiny ( resabs ) /(5.0e+01* epsilon ( resabs ) )) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk21_seepage
  subroutine qk31_seepage ( f, cSeepage_argu_in, a, b, result, abserr, resabs, resasc )
    !
    !************************************************************************
    !
    !! QK31 carries out a 31 point Gauss-Kronrod quadrature rule.
    !
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Reference:
    !
    !    R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real F, the name of the function routine, of the form
    !      function f ( x )
    !      real f
    !      real x
    !    which evaluates the integrand function.
    !
    !    Input, real A, B, the limits of integration.
    !
    !    Output, real RESULT, the estimated value of the integral.
    !                       result is computed by applying the 31-point
    !                       Gauss-Kronrod rule (resk), obtained by optimal
    !                       addition of abscissae to the 15-point Gauss
    !                       rule (resg).
    !
    !    Output, real ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real RESASC, approximation to the integral | F-I/(B-A) |
    !    over [A,B].
    !
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    real a
    real absc
    real abserr
    real b
    real centr
    real dhlgth
    real, external :: f
    real fc
    real fsum
    real fval1
    real fval2
    real fv1(15)
    real fv2(15)
    real hlgth
    integer j
    integer jtw
    integer jtwm1
    real resabs
    real resasc
    real resg
    real resk
    real reskh
    real result
    real wg(8)
    real wgk(16)
    real xgk(16)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 31-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 15-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 15-point Gauss rule
    !
    !           wgk    - weights of the 31-point Kronrod rule
    !
    !           wg     - weights of the 15-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16)/ &
         9.980022986933971e-01,   9.879925180204854e-01, &
         9.677390756791391e-01,   9.372733924007059e-01, &
         8.972645323440819e-01,   8.482065834104272e-01, &
         7.904185014424659e-01,   7.244177313601700e-01, &
         6.509967412974170e-01,   5.709721726085388e-01, &
         4.850818636402397e-01,   3.941513470775634e-01, &
         2.991800071531688e-01,   2.011940939974345e-01, &
         1.011420669187175e-01,   0.0e+00               /
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
         wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16)/ &
         5.377479872923349e-03,   1.500794732931612e-02, &
         2.546084732671532e-02,   3.534636079137585e-02, &
         4.458975132476488e-02,   5.348152469092809e-02, &
         6.200956780067064e-02,   6.985412131872826e-02, &
         7.684968075772038e-02,   8.308050282313302e-02, &
         8.856444305621177e-02,   9.312659817082532e-02, &
         9.664272698362368e-02,   9.917359872179196e-02, &
         1.007698455238756e-01,   1.013300070147915e-01/
    data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
         3.075324199611727e-02,   7.036604748810812e-02, &
         1.071592204671719e-01,   1.395706779261543e-01, &
         1.662692058169939e-01,   1.861610000155622e-01, &
         1.984314853271116e-01,   2.025782419255613e-01/
    !
    !
    !           list of major variables
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 15-point Gauss formula
    !           resk   - result of the 31-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 31-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
    fc = f(centr, cSeepage_argu_in)
    resg = wg(8)*fc
    resk = wgk(16)*fc
    resabs = abs(resk)

    do j = 1, 7
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 8
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk*5.0e-01
    resasc = wgk(16)*abs(fc-reskh)

    do j = 1, 15
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) &
         abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)

    if ( resabs > tiny ( resabs ) /(5.0e+01* epsilon ( resabs ) )) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk31_seepage
  subroutine qk41_seepage ( f,cSeepage_argu_in, a, b, result, abserr, resabs, resasc )
    !
    !************************************************************************
    !
    !! QK41 carries out a 41 point Gauss-Kronrod quadrature rule.
    !
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Reference:
    !
    !    R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real F, the name of the function routine, of the form
    !      function f ( x )
    !      real f
    !      real x
    !    which evaluates the integrand function.
    !
    !    Input, real A, B, the limits of integration.
    !
    !    Output, real RESULT, the estimated value of the integral.
    !                       result is computed by applying the 41-point
    !                       Gauss-Kronrod rule (resk) obtained by optimal
    !                       addition of abscissae to the 20-point Gauss
    !                       rule (resg).
    !
    !    Output, real ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real RESASC, approximation to the integral | F-I/(B-A) |
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 20-point Gauss formula
    !           resk   - result of the 41-point Kronrod formula
    !           reskh  - approximation to mean value of f over (a,b), i.e.
    !                    to i/(b-a)
    !
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    real a
    real absc
    real abserr
    real b
    real centr
    real dhlgth
    real, external :: f
    real fc
    real fsum
    real fval1
    real fval2
    real fv1(20)
    real fv2(20)
    real hlgth
    integer j
    integer jtw
    integer jtwm1
    real resabs
    real resasc
    real resg
    real resk
    real reskh
    real result
    real wg(10)
    real wgk(21)
    real xgk(21)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 41-point Gauss-Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 20-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 20-point Gauss rule
    !
    !           wgk    - weights of the 41-point Gauss-Kronrod rule
    !
    !           wg     - weights of the 20-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16), &
         xgk(17),xgk(18),xgk(19),xgk(20),xgk(21)/ &
         9.988590315882777e-01,   9.931285991850949e-01, &
         9.815078774502503e-01,   9.639719272779138e-01, &
         9.408226338317548e-01,   9.122344282513259e-01, &
         8.782768112522820e-01,   8.391169718222188e-01, &
         7.950414288375512e-01,   7.463319064601508e-01, &
         6.932376563347514e-01,   6.360536807265150e-01, &
         5.751404468197103e-01,   5.108670019508271e-01, &
         4.435931752387251e-01,   3.737060887154196e-01, &
         3.016278681149130e-01,   2.277858511416451e-01, &
         1.526054652409227e-01,   7.652652113349733e-02, &
         0.0e+00               /
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
         wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16), &
         wgk(17),wgk(18),wgk(19),wgk(20),wgk(21)/ &
         3.073583718520532e-03,   8.600269855642942e-03, &
         1.462616925697125e-02,   2.038837346126652e-02, &
         2.588213360495116e-02,   3.128730677703280e-02, &
         3.660016975820080e-02,   4.166887332797369e-02, &
         4.643482186749767e-02,   5.094457392372869e-02, &
         5.519510534828599e-02,   5.911140088063957e-02, &
         6.265323755478117e-02,   6.583459713361842e-02, &
         6.864867292852162e-02,   7.105442355344407e-02, &
         7.303069033278667e-02,   7.458287540049919e-02, &
         7.570449768455667e-02,   7.637786767208074e-02, &
         7.660071191799966e-02/
    data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),wg(10)/ &
         1.761400713915212e-02,   4.060142980038694e-02, &
         6.267204833410906e-02,   8.327674157670475e-02, &
         1.019301198172404e-01,   1.181945319615184e-01, &
         1.316886384491766e-01,   1.420961093183821e-01, &
         1.491729864726037e-01,   1.527533871307259e-01/
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute 41-point Gauss-Kronrod approximation to the
    !  the integral, and estimate the absolute error.
    !
    resg = 0.0e+00
    fc = f(centr, cSeepage_argu_in)
    resk = wgk(21)*fc
    resabs = abs(resk)

    do j = 1, 10
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 10
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk*5.0e-01
    resasc = wgk(21)*abs(fc-reskh)

    do j = 1, 20
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) &
         abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)

    if ( resabs > tiny ( resabs ) /(5.0e+01* epsilon ( resabs ) )) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk41_seepage
  subroutine qk51_seepage ( f,cSeepage_argu_in, a, b, result, abserr, resabs, resasc )
    !
    !************************************************************************
    !
    !! QK51 carries out a 51 point Gauss-Kronrod quadrature rule.
    !
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Reference:
    !
    !    R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real F, the name of the function routine, of the form
    !      function f ( x )
    !      real f
    !      real x
    !    which evaluates the integrand function.
    !
    !    Input, real A, B, the limits of integration.
    !
    !    Output, real RESULT, the estimated value of the integral.
    !                       result is computed by applying the 51-point
    !                       Kronrod rule (resk) obtained by optimal addition
    !                       of abscissae to the 25-point Gauss rule (resg).
    !
    !    Output, real ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real RESASC, approximation to the integral | F-I/(B-A) |
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 25-point Gauss formula
    !           resk   - result of the 51-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    real a
    real absc
    real abserr
    real b
    real centr
    real dhlgth
    real, external :: f
    real fc
    real fsum
    real fval1
    real fval2
    real fv1(25)
    real fv2(25)
    real hlgth
    integer j
    integer jtw
    integer jtwm1
    real resabs
    real resasc
    real resg
    real resk
    real reskh
    real result
    real wg(13)
    real wgk(26)
    real xgk(26)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 51-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 25-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 25-point Gauss rule
    !
    !           wgk    - weights of the 51-point Kronrod rule
    !
    !           wg     - weights of the 25-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14)/ &
         9.992621049926098e-01,   9.955569697904981e-01, &
         9.880357945340772e-01,   9.766639214595175e-01, &
         9.616149864258425e-01,   9.429745712289743e-01, &
         9.207471152817016e-01,   8.949919978782754e-01, &
         8.658470652932756e-01,   8.334426287608340e-01, &
         7.978737979985001e-01,   7.592592630373576e-01, &
         7.177664068130844e-01,   6.735663684734684e-01/
    data xgk(15),xgk(16),xgk(17),xgk(18),xgk(19),xgk(20),xgk(21), &
         xgk(22),xgk(23),xgk(24),xgk(25),xgk(26)/ &
         6.268100990103174e-01,   5.776629302412230e-01, &
         5.263252843347192e-01,   4.730027314457150e-01, &
         4.178853821930377e-01,   3.611723058093878e-01, &
         3.030895389311078e-01,   2.438668837209884e-01, &
         1.837189394210489e-01,   1.228646926107104e-01, &
         6.154448300568508e-02,   0.0e+00               /
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
         wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14)/ &
         1.987383892330316e-03,   5.561932135356714e-03, &
         9.473973386174152e-03,   1.323622919557167e-02, &
         1.684781770912830e-02,   2.043537114588284e-02, &
         2.400994560695322e-02,   2.747531758785174e-02, &
         3.079230016738749e-02,   3.400213027432934e-02, &
         3.711627148341554e-02,   4.008382550403238e-02, &
         4.287284502017005e-02,   4.550291304992179e-02/
    data wgk(15),wgk(16),wgk(17),wgk(18),wgk(19),wgk(20),wgk(21), &
         wgk(22),wgk(23),wgk(24),wgk(25),wgk(26)/ &
         4.798253713883671e-02,   5.027767908071567e-02, &
         5.236288580640748e-02,   5.425112988854549e-02, &
         5.595081122041232e-02,   5.743711636156783e-02, &
         5.868968002239421e-02,   5.972034032417406e-02, &
         6.053945537604586e-02,   6.112850971705305e-02, &
         6.147118987142532e-02,   6.158081806783294e-02/
    data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),wg(10), &
         wg(11),wg(12),wg(13)/ &
         1.139379850102629e-02,   2.635498661503214e-02, &
         4.093915670130631e-02,   5.490469597583519e-02, &
         6.803833381235692e-02,   8.014070033500102e-02, &
         9.102826198296365e-02,   1.005359490670506e-01, &
         1.085196244742637e-01,   1.148582591457116e-01, &
         1.194557635357848e-01,   1.222424429903100e-01, &
         1.231760537267155e-01/
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 51-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
    fc = f(centr, cSeepage_argu_in)
    resg = wg(13)*fc
    resk = wgk(26)*fc
    resabs = abs(resk)

    do j = 1, 12
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 13
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk*5.0e-01
    resasc = wgk(26)*abs(fc-reskh)

    do j = 1, 25
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) then
       abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
    end if

    if ( resabs > tiny ( resabs ) / (5.0e+01* epsilon ( resabs ) ) ) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk51_seepage
  subroutine qk61_seepage ( f,cSeepage_argu_in, a, b, result, abserr, resabs, resasc )
    !
    !************************************************************************
    !
    !! QK61 carries out a 61 point Gauss-Kronrod quadrature rule.
    !
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Reference:
    !
    !    R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real F, the name of the function routine, of the form
    !      function f ( x )
    !      real f
    !      real x
    !    which evaluates the integrand function.
    !
    !    Input, real A, B, the limits of integration.
    !
    !    Output, real RESULT, the estimated value of the integral.
    !                    result is computed by applying the 61-point
    !                    Kronrod rule (resk) obtained by optimal addition of
    !                    abscissae to the 30-point Gauss rule (resg).
    !
    !    Output, real ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real RESASC, approximation to the integral | F-I/(B-A) |
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 30-point Gauss rule
    !           resk   - result of the 61-point Kronrod rule
    !           reskh  - approximation to the mean value of f
    !                    over (a,b), i.e. to i/(b-a)
    !
    type(seepage_argument), intent(IN) :: cSeepage_argu_in
    real a
    real absc
    real abserr
    real b
    real centr
    real dhlgth
    real, external :: f
    real fc
    real fsum
    real fval1
    real fval2
    real fv1(30)
    real fv2(30)
    real hlgth
    integer j
    integer jtw
    integer jtwm1
    real resabs
    real resasc
    real resg
    real resk
    real reskh
    real result
    real wg(15)
    real wgk(31)
    real xgk(31)
    !
    !           the abscissae and weights are given for the
    !           interval (-1,1). because of symmetry only the positive
    !           abscissae and their corresponding weights are given.
    !
    !           xgk   - abscissae of the 61-point Kronrod rule
    !                   xgk(2), xgk(4)  ... abscissae of the 30-point
    !                   Gauss rule
    !                   xgk(1), xgk(3)  ... optimally added abscissae
    !                   to the 30-point Gauss rule
    !
    !           wgk   - weights of the 61-point Kronrod rule
    !
    !           wg    - weigths of the 30-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10)/ &
         9.994844100504906e-01,     9.968934840746495e-01, &
         9.916309968704046e-01,     9.836681232797472e-01, &
         9.731163225011263e-01,     9.600218649683075e-01, &
         9.443744447485600e-01,     9.262000474292743e-01, &
         9.055733076999078e-01,     8.825605357920527e-01/
    data xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16),xgk(17), &
         xgk(18),xgk(19),xgk(20)/ &
         8.572052335460611e-01,     8.295657623827684e-01, &
         7.997278358218391e-01,     7.677774321048262e-01, &
         7.337900624532268e-01,     6.978504947933158e-01, &
         6.600610641266270e-01,     6.205261829892429e-01, &
         5.793452358263617e-01,     5.366241481420199e-01/
    data xgk(21),xgk(22),xgk(23),xgk(24),xgk(25),xgk(26),xgk(27), &
         xgk(28),xgk(29),xgk(30),xgk(31)/ &
         4.924804678617786e-01,     4.470337695380892e-01, &
         4.004012548303944e-01,     3.527047255308781e-01, &
         3.040732022736251e-01,     2.546369261678898e-01, &
         2.045251166823099e-01,     1.538699136085835e-01, &
         1.028069379667370e-01,     5.147184255531770e-02, &
         0.0e+00                   /
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
         wgk(9),wgk(10)/ &
         1.389013698677008e-03,     3.890461127099884e-03, &
         6.630703915931292e-03,     9.273279659517763e-03, &
         1.182301525349634e-02,     1.436972950704580e-02, &
         1.692088918905327e-02,     1.941414119394238e-02, &
         2.182803582160919e-02,     2.419116207808060e-02/
    data wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),wgk(17), &
         wgk(18),wgk(19),wgk(20)/ &
         2.650995488233310e-02,     2.875404876504129e-02, &
         3.090725756238776e-02,     3.298144705748373e-02, &
         3.497933802806002e-02,     3.688236465182123e-02, &
         3.867894562472759e-02,     4.037453895153596e-02, &
         4.196981021516425e-02,     4.345253970135607e-02/
    data wgk(21),wgk(22),wgk(23),wgk(24),wgk(25),wgk(26),wgk(27), &
         wgk(28),wgk(29),wgk(30),wgk(31)/ &
         4.481480013316266e-02,     4.605923827100699e-02, &
         4.718554656929915e-02,     4.818586175708713e-02, &
         4.905543455502978e-02,     4.979568342707421e-02, &
         5.040592140278235e-02,     5.088179589874961e-02, &
         5.122154784925877e-02,     5.142612853745903e-02, &
         5.149472942945157e-02/
    data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
         7.968192496166606e-03,     1.846646831109096e-02, &
         2.878470788332337e-02,     3.879919256962705e-02, &
         4.840267283059405e-02,     5.749315621761907e-02, &
         6.597422988218050e-02,     7.375597473770521e-02/
    data wg(9),wg(10),wg(11),wg(12),wg(13),wg(14),wg(15)/ &
         8.075589522942022e-02,     8.689978720108298e-02, &
         9.212252223778613e-02,     9.636873717464426e-02, &
         9.959342058679527e-02,     1.017623897484055e-01, &
         1.028526528935588e-01/
    !
    centr = 5.0e-01*(b+a)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 61-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
    resg = 0.0e+00
    fc = f(centr)
    resk = wgk(31)*fc
    resabs = abs(resk)

    do j = 1, 15
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 15
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc, cSeepage_argu_in)
       fval2 = f(centr+absc, cSeepage_argu_in)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk * 5.0e-01
    resasc = wgk(31)*abs(fc-reskh)

    do j = 1, 30
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00 .and. abserr /= 0.0e+00) then
       abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
    end if

    if ( resabs > tiny ( resabs ) / (5.0e+01* epsilon ( resabs ) )) then
       abserr = max ( ( epsilon ( resabs ) *5.0e+01)*resabs, abserr )
    end if


    return
  end subroutine qk61_seepage
end module h2sc_solvers
