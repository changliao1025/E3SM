module h2sc_drainage
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_kind_mod , only : r4 => SHR_KIND_R4
  use decompMod , only : bounds_type
  use quadpack
  use h2sc_solvers
  use clm_cpl_indices
  use clm_varctl             , only : iulog
  use shr_sys_mod            , only : shr_sys_flush
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  implicit none
  real(kind = r8), parameter, public :: near_zero =1.0E-10

  !-----------------------------------------------------------------------
  !Added by Chang
  !global h2sc variables


  real(r8) , parameter::pi = 4.D0*DATAN(1.D0)
  real(r8) , parameter::beta_h2sc = 3.6
  real(r8) , parameter::n_h2sc = 2.0
  real(r8) , parameter::m_h2sc = 1.0 - 1/n_h2sc
  type, public::hillslope_class
     integer::nElevation_interval
     integer::iFlag_hillslope

     real(kind = r8)::dChannel_depth
     real(kind = r8)::dArea_hillslope
     real(kind = r8)::dElevation_max
     real(kind = r8)::dElevation_min
     real(kind = r8)::dWidth_hillslope
     real(kind = r8)::dLength_hillslope
     real(kind = r8)::dHeight_below_river
     real(kind = r8)::dSlope_bedrock  !radian
     real(kind = r8)::dSlope_surface
     real(kind = r8)::dThickness_critical_zone
     real(kind = r8), dimension(:), allocatable:: aElevation_profile
     !===============================================
     !may be moved to other data type
     !===============================================
     real(kind = r8)::dSaturate_hydraulic_conductivity
     real(kind = r8)::beta_hillslope
     real(kind = r8)::m_hillslope
     real(kind = r8)::n_hillslope
   contains
     procedure::initialize_hillslope_class
     procedure::calculate_surface_slope
     procedure::calculate_bedrock_slope
  end type hillslope_class
  type, public::groundwater_class
     real(kind = r8)::dElevation_water_table !the groundwater level
   contains
     procedure::calculate_groundwater_level
  end type groundwater_class
  type, public::seepage_class
     integer:: iSeepage_type !the type of seepage face
     integer:: iSeepage_count
     integer:: iFlag_intersect
     real(kind = r8)::dFlow_seepage
     real(kind = r8)::dFlow_downslope
     real(kind = r8)::dLength_seepage
     real(kind = r8)::dLength_seepage_dominant
     real(kind = r8)::dLength_water_table
     real(kind = r8)::phi
     real(kind = r8)::dSlope_water_table
     integer, dimension(:), allocatable:: aSeepage_mask
     type(seepage_argument):: cSeepage_argu
   contains
     procedure::calculate_seepage_type
     procedure::calculate_water_table_properties_without_seepage
     procedure::calculate_water_table_properties_with_seepage
     procedure::calculate_seepage_extent
     procedure::calculate_seepage_length
     procedure::calculate_downslope_flow
     procedure::calculate_seepage_flow
     procedure::calculate_seepage_dominant_distance
     procedure::setup_seepage_solver_argument
     procedure::calculate_phi_1
     procedure::calculate_phi_2
  end type seepage_class
  type, public::soil_class
     real(kind = r8)::dDrainage_downslope
     real(kind = r8)::dDrainage_seepage
     real(kind = r8)::psi
     real(kind = r8)::dSaturate_hydraulic_conductivity
     !===============================================
     !parameters from van Genuchten equation
     !===============================================
     real(kind = r8)::beta_soil
     real(kind = r8)::m_soil
     real(kind = r8)::n_soil
   contains
     procedure::run_soil
     procedure::calculate_drainage_downslope
     procedure::calculate_drainage_seepage
     procedure::calculate_psi
  end type soil_class
  type, public :: soil_column_class
     integer :: iSoil_type
     real(kind = r8):: dDelta_h !the height difference
     real(kind = r8):: dDistance !distance from stream channel
     real(kind = r8):: dElevation
     real(kind = r8):: dFactor_seepage
     real(kind = r8):: dDrainage_et
     real(kind = r8):: dDrainage !the total drainage from the soil?
     real(kind = r8):: dThickness_critical_zone
     !===============================================
     type(soil_class):: cSoil
     type(groundwater_class):: cGroundwater
     !===============================================
   contains
     procedure::run_soil_column
     procedure::calculate_delta_h
     procedure::calculate_soil_column_drainage
     procedure::calculate_seepage_factor
  end type soil_column_class
  !===============================================
  !interfaces
  !===============================================
  interface hillslope_class
     module procedure initialize_hillslope_class
  end interface hillslope_class
  interface soil_column_class
     module procedure initialize_soil_column_class
  end interface soil_column_class
  interface soil_class
     module procedure initialize_soil_class
  end interface soil_class
  interface groundwater_class
     module procedure initialize_groundwater_class
  end interface groundwater_class
  type(hillslope_class), DIMENSION(:), allocatable:: cHillslope
  type(seepage_class), DIMENSION(:), allocatable:: cSeepage
  !be careful here
  !type(seepage_argument), allocatable:: cSeepage_argu (:)
contains

  subroutine initialize_hillslope_class(this)
    implicit none
    class(hillslope_class), intent(inout) :: this
    this%dElevation_max = MAXVAL(this%aElevation_profile)
    this%dElevation_min = MINVAL(this%aElevation_profile)
    if(this%dElevation_max > 8848.0) then
       this%dElevation_max = 8848.0
    end if
    if(this%dElevation_min < -413.0) then
       this%dElevation_min = -413.0
    end if
    return
  end subroutine initialize_hillslope_class
  !===============================================
  !calculate surface slope
  !===============================================
  function calculate_surface_slope(this) &
       result(dSlope_surface_out)
    !===============================================
    implicit none
    class(hillslope_class), intent(inout) :: this
    real(kind = r8) :: dDistance
    real(kind = r8) :: dElevation_difference
    real(kind = r8) :: dElevation_max
    real(kind = r8) :: dElevation_min
    real(kind = r8) :: dSlope_surface_out
    real(kind = r8) :: dDummy
    dElevation_max = this%dElevation_max
    dElevation_min = this%dElevation_min
    dElevation_difference = dElevation_max - dElevation_min
    dDistance = this%dLength_hillslope
    !the equation used to calculate slope
    !tan(slope) = elevation / distance
    dDummy = dElevation_difference / dDistance
    !tan/atan function:
    !radian is used, so dSlope_surface_out is radian
    !to convert it to angle, dSlope_surface_angle = dSlope_surface_out / pi * 180.0
    !for example, the maximum slope is 0.2. so the angle is 0.2 /pi * 180 = 11.45 degree
    !and tan(11.45 / 180 * pi) = 0.2
    dSlope_surface_out = atan( dDummy )
    !===============================================
    this%dSlope_surface = dSlope_surface_out
    !===============================================
    return
  end function calculate_surface_slope
  function calculate_bedrock_slope(this) &
       result(dSlope_bedrock_out)
    !===============================================
    implicit none
    class(hillslope_class), intent(inout) :: this
    real(kind = r8) :: dElevation_difference
    real(kind = r8) :: dDistance
    real(kind = r8) :: dSlope_bedrock_out
    real(kind = r8) :: dElevation_max
    real(kind = r8) :: dElevation_min
    real(kind = r8) :: dThickness_critical_zone_min
    real(kind = r8) :: dThickness_critical_zone_max
    real(kind = r8) :: dDummy
    dElevation_max = this%dElevation_max
    dElevation_min = this%dElevation_min
    !we assume that the soil thickness decreases at the elevation peak
    dThickness_critical_zone_min = this%dThickness_critical_zone * 0.95
    dThickness_critical_zone_max = this%dThickness_critical_zone
    !we assume the bedrock is almost paralle with surface
    !as a result, the bedrock slope is the same with surface slope
    dElevation_difference = ( dElevation_max - dThickness_critical_zone_min) &
         - ( dElevation_min - dThickness_critical_zone_max )
    dDistance = this%dLength_hillslope
    dDummy = dElevation_difference / dDistance
    dSlope_bedrock_out = atan( dDummy )
    !===============================================
    this%dSlope_bedrock = dSlope_bedrock_out
    !===============================================
    return
  end function calculate_bedrock_slope

  !===============================================
  !calculate the downslope drainage flow
  !===============================================
  function calculate_downslope_flow( this , &
       dArea_hillslope_in, &
       dHeight_below_river_in, &
       dLength_hillslope_in, &
       dSaturate_hydraulic_conductivity_in, &
       dSlope_surface_in, &
       dWidth_hillslope_in ) &
       result(dFlow_downslope_out)
    implicit none
    !===============================================
    class(seepage_class), intent(inout):: this
    real(kind = r8)::dFlow_downslope_out
    real(kind = r8)::dArea_hillslope_in
    real(kind = r8)::dHeight_below_river_in
    real(kind = r8)::dLength_hillslope_in
    real(kind = r8)::dSaturate_hydraulic_conductivity_in !mm/s
    real(kind = r8)::dSlope_surface_in !the slope of hillslope surface
    real(kind = r8)::dWidth_hillslope_in
    !===============================================
    real(kind = r8)::dLength_seepage !the total length of seepage
    real(kind = r8)::dLength_water_table
    real(kind = r8)::dSlope_water_table !the slope of water table
    real(kind = r8)::dDummy0, dDummy1, dDummy2, dDummy3, dDummy4, dDummy5, dDummy6
    real(kind = r8)::dDummy_slope
    !===============================================
    !retrieve value from seepage object
    !===============================================
    dLength_seepage = this%dLength_seepage
    dLength_water_table = this%dLength_water_table
    dSlope_water_table = this%dSlope_water_table
    !===============================================
    if( dLength_seepage <= near_zero ) then
       dDummy_slope = dSlope_water_table
    else
       dDummy_slope = dSlope_surface_in
    endif
    !length L from meter to mm
    dDummy0 =  dLength_water_table * 1000.0
    !head pressure delta H, average
    dDummy1 = dDummy0 * tan(dDummy_slope)
    !cross area A
    dDummy2 = dHeight_below_river_in * 1000.0 * dWidth_hillslope_in  * 1000.0
    !(delta H / L)
    dDummy3 = dDummy1 / dDummy0
    !q = k A (delta H / L)
    dDummy4 = dSaturate_hydraulic_conductivity_in * dDummy2 * dDummy3

    !drainage area, normalization
    dDummy5 = dArea_hillslope_in * 1000.0 * 1000.0
    !dFlow_downslope_out = dDummy4 / dDummy5, which is exactly the same with the following:

    dDummy6 =  (dLength_seepage + dLength_water_table)/dLength_hillslope_in
    if(dDummy6 > 1)then
       dDummy6 = 1
    end if
    dFlow_downslope_out = dSaturate_hydraulic_conductivity_in * dHeight_below_river_in &
         * tan(dDummy_slope) / dLength_hillslope_in

    !===============================================
    this%dFlow_downslope = dFlow_downslope_out
    !===============================================
    return
  end function calculate_downslope_flow
  function calculate_groundwater_level(this, &
       iSeepage_type_in, &
       dElevation_in, &
       dElevation_max_in,&
       dElevation_min_in, &
       dElevation_water_table_in, &
       dThickness_critical_zone_in) &
       result(dElevation_water_table_out)
    implicit none
    !===============================================
    class(groundwater_class), intent(inout) :: this
    integer ::iSeepage_type_in
    real(kind = r8)::dElevation_in
    real(kind = r8)::dElevation_water_table_in
    real(kind = r8)::dThickness_critical_zone_in
    real(kind = r8)::dElevation_water_table_out
    real(kind = r8):: a, b, x, y
    real(kind = r8):: dElevation_max_in
    real(kind = r8):: dElevation_min_in
    !===============================================
    !iSeepage_type = this%cSeepage%iSeepage_type
    ! dThickness_soil = 50.0 !unit: meter
    if (iSeepage_type_in == 0 ) then
       !lower than outlet
       !useing the equation sets
       !calculate the cooefficients in y = ax + b
       a = (dElevation_max_in - dThickness_critical_zone_in - dElevation_water_table_in) &
            / ( dElevation_max_in - dElevation_min_in )
       b = dElevation_water_table_in &
            - a * ( dElevation_min_in - dElevation_water_table_in )
    else if( iSeepage_type_in == 1 ) then
       !in between
       b = dElevation_water_table_in
       a = 1 - ( dThickness_critical_zone_in ) &
            / ( dElevation_max_in - dElevation_water_table_in)
    else
       ! higher than peakb
       a = 0
       b = dElevation_water_table_in
    end if
    x = dElevation_in - dElevation_water_table_in
    y = a * x + b
    dElevation_water_table_out = y
    this%dElevation_water_table = dElevation_water_table_out
    return
  end function calculate_groundwater_level
  !===============================================
  !calculate the seepage type
  !===============================================
  function calculate_seepage_type (this , &
       dElevation_max_in, &
       dElevation_min_in, &
       dElevation_water_table_in) &
       result(iSeepage_type_out)
    implicit none
    !===============================================
    class(seepage_class), intent(inout):: this
    integer::iSeepage_type_out
    real(kind = r8)::dElevation_max_in
    real(kind = r8)::dElevation_min_in
    real(kind = r8)::dElevation_water_table_in
    !===============================================
    iSeepage_type_out =1
    if (dElevation_water_table_in <= dElevation_min_in) then
       !lower than outlet
       iSeepage_type_out = 1
    else if( dElevation_water_table_in >= dElevation_max_in) then
       ! higher than peak
       iSeepage_type_out = 3
    else
       !in between
       iSeepage_type_out = 2
    end if
    this%iSeepage_type = iSeepage_type_out
    return
  end function calculate_seepage_type
  !===============================================
  !calculate the water table slope when there is no
  !seepage face
  !in this case, obvious the outlet is the lowest point
  !===============================================
  function calculate_water_table_properties_without_seepage(this, &
       dElevation_max_in, &
       dElevation_min_in, &
       dElevation_water_table_in, &
       dLength_hillslope_in, &
       dSlope_bedrock_in, &
       dSlope_surface_in, &
       dThickness_critical_zone_in ) &
       result(dSlope_water_table_out)
    !===============================================
    implicit none
    class(seepage_class), intent(inout) :: this
    real(kind = r8) :: dElevation_min_in
    real(kind = r8) :: dElevation_max_in
    real(kind = r8) :: dElevation_water_table_in
    real(kind = r8) :: dLength_hillslope_in
    real(kind = r8) :: dSlope_bedrock_in
    real(kind = r8) :: dSlope_surface_in
    real(kind = r8) :: dThickness_critical_zone_in
    !===============================================
    real(kind = r8) :: dSlope_water_table_reference
    real(kind = r8) :: dSlope_water_table_out
    real(kind = r8) :: dDummy1, dDummy2, dDummy3, dDummy4
    real(kind = r8) :: a1, b1, a2, b2, a3, b3, a5, b5
    real(kind = r8) :: x23, y23
    real(kind = r8) ::  y3_1, y3_2, x5_2, y5_2
    real(kind = r8) :: dRange0, dRange1
    !we need to guarantee whether the bottowm of critical zone at summit is higher than watertable
    a5 = tan(dSlope_surface_in) * 0.9
    dSlope_water_table_reference = a5
    b5 = dElevation_min_in
    x5_2 = dLength_hillslope_in
    y5_2 = a5 * x5_2 + b5
    dRange0 = dThickness_critical_zone_in
    dRange1 = dThickness_critical_zone_in + a5 * x5_2

    dDummy1 = dElevation_min_in - dElevation_water_table_in
    dDummy2 = dDummy1 / dThickness_critical_zone_in
    dDummy3 = dDummy2**1.0
    dDummy4 = dRange1 * dDummy3
    y3_2 = y5_2 - dDummy4
    y3_1 = dElevation_water_table_in
    a3 = (y3_2 - y3_1) / dLength_hillslope_in
    dSlope_water_table_out = atan( a3 )
    !with the slope known, we can also check the intersect length because it might intersect with the bed rock
    b3 = dElevation_water_table_in
    !x = (b3-b2) / (a2-a3)
    b2 = dElevation_min_in - dThickness_critical_zone_in
    a2 = tan(dSlope_bedrock_in)
    x23 = (b3-b2) / (a2-a3)
    y23 = a2 * x23 + b2
    if( x23 < 0 ) then
       write(iulog ,*) 'x23: ',a2, b2, a3, b3, y3_1, y3_2
    end if
    if( x23  <= dLength_hillslope_in) then
       !intersect place is less than the whole hillslope
       this%iFlag_intersect = 1
       this%dLength_water_table = x23
    else
       !not intersect
       this%iFlag_intersect = 0
       this%dLength_water_table = dLength_hillslope_in
    end if
    if( dSlope_water_table_out < 0.0 ) then
       write(iulog ,*) dElevation_max_in, dElevation_min_in, dDummy1, dDummy2, dDummy3, dDummy4, y5_2, y3_1, y3_2, dRange0, dRange1
       call shr_sys_flush(iulog)
    end if
    !===============================================
    this%dSlope_water_table = dSlope_water_table_out
    !===============================================
    return
  end function calculate_water_table_properties_without_seepage

  !this is a function used to location the index satisfying some matrix operations
  function idl_where(aArray_in, m, dValue, iFlag) result (aIndex_out)
    implicit NONE
    integer :: m , i, j, iFlag
    real(kind = r4), dimension (m), intent(in) :: aArray_in
    real(kind = r4):: dValue
    integer, dimension(m) :: aIndex_out
    integer:: iCount
    integer:: iFlag_allocation
    aIndex_out(:)=0
    iCount = 0
    select case (iFlag)
    case (1) !equal
       do j = 1, m, 1
          if ( aArray_in(j) == dValue ) then
             iCount = iCount + 1
             aIndex_out(iCount) = j
          end if
       end do
    case (2) !not equal
       do j = 1, m, 1
          if ( aArray_in(j) /= dValue ) then
             iCount = iCount + 1
             aIndex_out(iCount) = j
          end if
       end do
    case (3) !greater
       do j = 1, m, 1
          if ( aArray_in(j) > dValue ) then
             iCount = iCount + 1
             aIndex_out(iCount) = j
          end if
       end do
    case (4) !less
       do j = 1, m, 1
          if ( aArray_in(j) < dValue ) then
             iCount = iCount + 1
             aIndex_out(iCount) = j
          end if
       end do
    end select
    !check count
    if (iCount > 0) then
       !allocate the result
       return
    else
       !there is no match
       aIndex_out = -1
    end if
    return
  end function idl_where
  !===============================================
  !calculate the slope when there is seepage
  !===============================================
  function calculate_water_table_properties_with_seepage(this, &
       nElevation_interval_in, &
       dElevation_max_in, &
       dElevation_min_in, &
       dElevation_water_table_in , &
       dLength_hillslope_in, &
       dSlope_bedrock_in, &
       dSlope_surface_in, &
       dThickness_critical_zone_in) &
       result(dSlope_water_table_out)
    !===============================================
    implicit none
    class(seepage_class), intent(inout) :: this
    real(kind = r8) :: dSlope_water_table_out  !unit radian
    integer, intent(in) :: nElevation_interval_in
    real(kind = r8), intent(in) :: dElevation_min_in
    real(kind = r8), intent(in) :: dElevation_max_in
    real(kind = r8), intent(in) :: dElevation_water_table_in
    real(kind = r8) , intent(in) :: dLength_hillslope_in
    real(kind = r8) , intent(in) :: dSlope_bedrock_in
    real(kind = r8) , intent(in) :: dSlope_surface_in
    real(kind = r8) , intent(in) :: dThickness_critical_zone_in
    !===============================================
    real(kind = r8) :: dSlope_water_table_reference
    real(kind = r8) :: dThickness_critical_zone_max
    real(kind = r8) :: dDummy1, dDummy2, dDummy3, dDummy4, dDummy5, dDummy6
    real(kind = r8) :: a1, b1, a2, b2, a3, b3, a4, b4, a5, b5
    real(kind = r8) :: dRange2, dRange3
    real(kind = r8) :: x14, y14, x24, y24
    real(kind = r8) :: x4_2, y4_2, x5_2, y5_2
    integer, dimension(nElevation_interval_in):: aSeepage_mask
    integer, dimension(:), allocatable :: aIndex
    aSeepage_mask = this%aSeepage_mask

    !lowest: watertable
    a1 = tan(dSlope_surface_in)

    b1 = dElevation_min_in
    a2 = tan(dSlope_bedrock_in)
    b2 = dElevation_min_in - dThickness_critical_zone_in

    a5 = tan(dSlope_surface_in) * 0.9
    dSlope_water_table_reference = a5
    b5 = dElevation_min_in
    x5_2 = dLength_hillslope_in
    y5_2 = a5 * x5_2 + b5

    dRange2 =   dElevation_max_in - dElevation_min_in
    dRange3 = dElevation_max_in - y5_2

    dDummy1 = dElevation_water_table_in - dElevation_min_in
    dDummy2 = dElevation_max_in - dElevation_water_table_in
    if ( dDummy2  < near_zero ) then
       !the water table is almost at peak elevation
       !in this case, we have to adjust the water table

       dSlope_water_table_out = 0.0
       this%iFlag_intersect = 0
       this%dLength_water_table = dLength_hillslope_in
       this%dSlope_water_table = dSlope_water_table_out
    else
       !water table is at hillslope
       dDummy3 = dDummy1 / dRange2
       dDummy4 =  dDummy3**0.5
       dDummy5 = dRange3 * dDummy4
       y4_2 = y5_2 + dDummy5
       y14 = dElevation_water_table_in
       x14 = (y14 - b1)  /  a1
       x4_2 = dLength_hillslope_in
       a4 = ( y4_2 - y14)/( x4_2 - x14 )
       dSlope_water_table_out = atan( a4  )
       b4 = y14 - a4 * x14
       x24 = (b4-b2) / (a2-a4)
       y24 = a2 * x24 + b2

       if( x24 < 0 ) then
          write(iulog ,*) 'x24: ',a2, b2, a4, b4
       end if
       if( dSlope_water_table_out < 0.0 ) then
          write(*,*) dElevation_max_in, dElevation_min_in, dDummy3, dDummy4, dDummy5,  a4,  b1, x14, y14
       end if
       if( x24 < dLength_hillslope_in ) then
          !intersect
          this%iFlag_intersect = 0
          this%dLength_water_table = x24 - x14
       else
          !not intersect
          this%iFlag_intersect = 1
          this%dLength_water_table = dLength_hillslope_in - x14
       end if
       this%dSlope_water_table = dSlope_water_table_out
       !===============================================
    end if
    return
  end function calculate_water_table_properties_with_seepage
  !===============================================
  !calculate the seepage face extent
  !===============================================
  function calculate_seepage_extent (this , &
       nElevation_interval_in, &
       dElevation_water_table_in, &
       aElevation_profile_in ) &
       result(aSeepage_mask_out)
    implicit none
    !===============================================
    class(seepage_class), intent(inout):: this
    integer, intent(in) ::nElevation_interval_in
    real(kind = r8), intent(in) ::dElevation_water_table_in
    real(kind = r8), dimension(nElevation_interval_in), intent(in) ::aElevation_profile_in
    integer, dimension(nElevation_interval_in) :: aSeepage_mask_out
    integer::iSeepage_type
    !===============================================
    aSeepage_mask_out(:) = 0
    iSeepage_type = this%iSeepage_type
    !===============================================
    select case (iSeepage_type)
    case (1)
       !no seepage at all
    case (2)
       ! in this case, we need to calculate the spatial extent based on dem and stream networks
       !use where to location seepage pixels
       !be careful with where
       where (aElevation_profile_in <= dElevation_water_table_in )
          !seepage
          aSeepage_mask_out = 1
       elsewhere
          !upland
       end where
       where (aElevation_profile_in > dElevation_water_table_in)
          !seepage
          aSeepage_mask_out = 0
       elsewhere
          !upland
       end where
    case (3)
       !completely flooded
       aSeepage_mask_out(:) = 1
    end select

    this%aSeepage_mask = aSeepage_mask_out

    this%iSeepage_count =  sum( this%aSeepage_mask )
    return
  end function calculate_seepage_extent


  function calculate_seepage_length(this,&
       nElevation_interval_in,&
       dLength_hillslope_in) &
       result(dLength_seepage_out)
    implicit none
    !===============================================
    class(seepage_class), intent(inout):: this
    integer, intent(in)  :: nElevation_interval_in
    real(kind = r8), intent(in) ::dLength_hillslope_in
    real(kind = r8) ::dLength_seepage_out
    real(kind = r8)::dummy
    integer :: iIndex
    iIndex = sum( this%aSeepage_mask )

    dummy = real(iIndex)/ nElevation_interval_in
    dLength_seepage_out = dummy * dLength_hillslope_in
    this%dLength_seepage = dLength_seepage_out
    return
  end function calculate_seepage_length
  !===============================================
  !calculate the seepage flow
  !===============================================
  function calculate_seepage_flow(this, &
       dLength_hillslope_in,&
       dSaturate_hydraulic_conductivity_in, &
       dSlope_surface_in) &
       result(dFlow_seepage_out)
    implicit none
    !===============================================
    class(seepage_class), intent(inout):: this
    real(kind = r8), intent(in) ::dLength_hillslope_in
    real(kind = r8), intent(in) ::dSaturate_hydraulic_conductivity_in
    real(kind = r8), intent(in) ::dSlope_surface_in
    real(kind = r8) ::dFlow_seepage_out
    real(kind = r8) ::dLength_seepage
    real(kind = r8) ::dLength_water_table
    real(kind = r8) ::dDummy6
    !===============================================
    dLength_seepage = this%dLength_seepage
    dLength_water_table = this%dLength_water_table
    dDummy6 =  (dLength_seepage + dLength_water_table)/dLength_hillslope_in
    if(dDummy6 > 1)then
       dDummy6 = 1
    end if
    dFlow_seepage_out = dSaturate_hydraulic_conductivity_in &
         * dLength_seepage * (tan(dSlope_surface_in) **2) / dLength_hillslope_in
    this%dFlow_seepage = dFlow_seepage_out

    !write(iulog ,*) 'seepage flow: ', dFlow_seepage_out, dSaturate_hydraulic_conductivity_in, dLength_seepage, dSlope_surface_in
    !call shr_sys_flush(iulog)

    return
  end function calculate_seepage_flow
  !===============================================
  !calculate the dominant seepage distance, equation 29
  !===============================================
  function calculate_seepage_dominant_distance(this) &
       result(dLength_seepage_dominant_out)
    implicit none
    !===============================================
    class(seepage_class), intent(inout):: this
    real(kind = r8) :: dLength_seepage_dominant_out
    real(kind = r8) :: dLength_water_table
    type(seepage_argument):: cSeepage_argu
    real(kind = r8):: dummy0
    real(kind = r8):: x, f, xNew
    real(kind = r8):: residual, threshold
    integer :: nIter , maxIter
    dLength_water_table = this%dLength_water_table
    cSeepage_argu = this%cSeepage_argu
    !===============================================
    nIter = 1
    maxIter = 10000
    !set the condition every time step
    residual = gx(dLength_water_table, cSeepage_argu)
    if (residual <= 0.0) then
       dLength_seepage_dominant_out =  dLength_water_table
    else
       x = dLength_water_table
       !we need to find out the unique solution to the function
       !call newton_method(gx, gx_prime, x0, x, iters, debug)
       do while ((residual > near_zero) .and. (nIter < maxIter))
          !! Search using conventional Newton's method
          call newton_method_seepage(x,xNew, cSeepage_argu, f,residual)
          !! Store a newly updated root into an array
          !x_history(nIter) = xNew
          !! Save for the next search iteration
          x = xNew
          !! Store a newly updated residual into an array
          !res_history(nIter) = residual
          !! Update iteration number
          nIter = nIter + 1
          write(*, *) 'The x is: ', x , 'at iteration: ' , nIter , 'and residual is: ', residual
       end do
       !===============================================
       dLength_seepage_dominant_out = x
    end if

    this%dLength_seepage_dominant = dLength_seepage_dominant_out
    !===============================================
    return
  end function calculate_seepage_dominant_distance
  !===============================================
  !set up arguments for solver
  !===============================================
  subroutine setup_seepage_solver_argument(this, &
       dFlow_downslope_in, &
       dLength_hillslope_in, &
       dSlope_surface_in,&
       beta_seepage_in, m_seepage_in, n_seepage_in )
    !use h2sc_solvers
    implicit none
    class(seepage_class), intent(inout):: this
    real(kind = r8), intent(in)::dFlow_downslope_in
    real(kind = r8), intent(in)::dLength_hillslope_in
    real(kind = r8), intent(in)::dSlope_surface_in
    real(kind = r8), intent(in)::beta_seepage_in
    real(kind = r8), intent(in)::m_seepage_in
    real(kind = r8), intent(in)::n_seepage_in
    this%cSeepage_argu%dFlow_seepage=this%dFlow_seepage
    this%cSeepage_argu%dLength_seepage = this%dLength_seepage
    this%cSeepage_argu%xBeg = this%dLength_seepage
    this%cSeepage_argu%dLength_water_table = this%dLength_water_table
    this%cSeepage_argu%dSlope_water_table = this%dSlope_water_table
    !===============================================
    this%cSeepage_argu%dFlow_downslope = dFlow_downslope_in
    this%cSeepage_argu%xEnd = dLength_hillslope_in
    this%cSeepage_argu%dSlope_surface = dSlope_surface_in
    this%cSeepage_argu%beta_seepage=beta_seepage_in
    this%cSeepage_argu%m_seepage=m_seepage_in
    this%cSeepage_argu%n_seepage=n_seepage_in
    return
  end subroutine setup_seepage_solver_argument
  !===============================================
  !calculate phi using integration
  !===============================================
  function calculate_phi_1(this, &
       dLength_water_table_in    ) &
       result(phi_1_out)
    implicit none
    !===============================================
    class(seepage_class), intent(inout):: this
    type(seepage_argument) :: cSeepage_argu
    real(kind = r4) :: phi_1_out
    !===============================================
    real(kind = r4) :: abserr
    real(kind = r4) , parameter :: epsabs = 0.0E+00
    real(kind = r4) , parameter :: epsrel = 0.001E+00
    integer :: ier
    integer , parameter :: key = 6
    integer :: neval
    real (kind = r4) :: dummy0
    !===============================================
    real(kind = r8)::dLength_water_table_in
    !===============================================
    cSeepage_argu = this%cSeepage_argu
    dummy0 = real(dLength_water_table_in, 4)
    !===============================================
    call qag_seepage ( f_phi_1, cSeepage_argu, 0.0, dummy0, epsabs, epsrel, key, phi_1_out, abserr, neval, ier )
    !===============================================
    this%phi=phi_1_out
    return
  end function calculate_phi_1
  !===============================================
  !calculate phi using integration
  !===============================================
  function calculate_phi_2(this,&
       dLength_seepage_in,&
       dLength_water_table_in ) &
       result(phi_2_out)
    implicit none
    !===============================================
    class(seepage_class), intent(inout):: this
    real(kind = r4)::phi_2_out
    type(seepage_argument) :: cSeepage_argu
    !===============================================
    real(kind = r4) :: abserr
    real(kind = r4) , parameter :: epsabs = 0.0E+00
    real(kind = r4) , parameter :: epsrel = 0.001E+00
    integer :: ier
    integer , parameter :: key = 6
    integer :: neval
    real(kind = r4) ::dummy0, dummy1
    !===============================================
    real(kind = r8)::dLength_seepage_in
    real(kind = r8)::dLength_water_table_in
    !===============================================
    cSeepage_argu = this%cSeepage_argu
    dummy0 = REAL(dLength_seepage_in , 4)
    dummy1 = REAL(dLength_water_table_in, 4)
    if (dummy0 > dummy1)then
       write(*,*) 'Warning: the seepage length is larger than water table length'
       stop
    end if
    !===============================================
    call qag_seepage ( f_phi_2, cSeepage_argu, dummy0, dummy1, epsabs, epsrel, key, phi_2_out, abserr, neval, ier )
    !===============================================
    this%phi = phi_2_out
    !===============================================
    return
  end function calculate_phi_2
  function run_soil(this, &
       dDistance_in,&
       dElevation_water_table_in ,&
       dFlow_downslope_in, &
       dFlow_seepage_in, &
       dHeight_below_river_in, &
       dLength_hillslope_in, &
       dLength_seepage_in,&
       dLength_water_table_in, &
       dSlope_intersect_in, &
       dSlope_surface_in, &
       dSlope_water_table_in,&
       phi_in) &
       result(error_code)
    implicit none
    class(soil_class), intent(inout):: this
    real(kind = r8)::dDistance_in
    real(kind = r8)::dElevation_water_table_in
    real(kind = r8)::dFlow_downslope_in
    real(kind = r8)::dFlow_seepage_in
    real(kind = r8)::dHeight_below_river_in
    real(kind = r8)::dLength_hillslope_in
    real(kind = r8)::dLength_seepage_in
    real(kind = r8)::dLength_water_table_in
    real(kind = r8)::dSlope_intersect_in
    real(kind = r8)::dSlope_surface_in
    real(kind = r8)::dSlope_water_table_in
    real(kind = r8)::phi_in
    real(kind = r8)::dummy0 , dummy1, dummy2
    integer::error_code
    dummy0 = this%calculate_psi( dDistance_in, &
         dElevation_water_table_in,&
         dHeight_below_river_in, &
         dSlope_surface_in)
    dummy1 = this%calculate_drainage_downslope(dDistance_in, &
         dFlow_downslope_in, &
         dLength_hillslope_in , &
         dLength_water_table_in,&
         dSlope_intersect_in, &
         dSlope_surface_in, &
         dSlope_water_table_in, &
         phi_in)
    dummy2 = this%calculate_drainage_seepage(dDistance_in ,&
         dFlow_downslope_in,&
         dLength_hillslope_in,&
         dLength_water_table_in,&
         dSlope_intersect_in,&
         dSlope_surface_in,&
         dSlope_water_table_in, &
         phi_in,&
         dLength_seepage_in , &
         dFlow_seepage_in)
    error_code = 1
    return
  end function run_soil
  !===============================================
  !calculate psi
  !===============================================
  function calculate_psi(this, &
       dDistance_in, &
       dElevation_water_table_in,&
       dHeight_below_river_in, &
       dSlope_surface_in) &
       result(psi_out)
    implicit none
    !===============================================
    class(soil_class), intent(inout):: this
    real(kind = r8)::dDistance_in
    real(kind = r8)::dElevation_water_table_in
    real(kind = r8)::psi_out
    !===============================================
    real(kind = r8):: beta_soil
    real(kind = r8):: m_soil
    real(kind = r8):: n_soil
    !===============================================
    real(kind = r8):: dHeight_below_river_in
    real(kind = r8):: dSlope_surface_in
    !===============================================
    real(kind = r8):: dummy0, dummy1, dummy2
    !===============================================
    m_soil=this%m_soil
    n_soil=this%n_soil
    beta_soil=this%beta_soil
    !===============================================
    !===============================================
    dummy0 = dHeight_below_river_in &
         + dDistance_in * tan(dSlope_surface_in) &
         - dElevation_water_table_in
    dummy1 = dummy0 ** n_soil
    dummy2 = 1 + beta_soil**n_soil * dummy1
    psi_out = 1 - dummy2 **(-1*m_soil)
    !===============================================
    this%psi=psi_out
    return
  end function calculate_psi
  !===============================================
  !calculate d1 in function
  !===============================================
  function calculate_drainage_downslope(this, &
       dDistance_in,&
       dFlow_downslope_in,&
       dLength_hillslope_in , &
       dLength_water_table_in,&
       dSlope_intersect_in,&
       dSlope_surface_in, &
       dSlope_water_table_in, &
       phi_in) &
       result(dDrainage_downslope_out)
    implicit none
    !===============================================
    class(soil_class), intent(inout):: this
    real(kind = r8)::dDistance_in
    real(kind = r8)::dDrainage_downslope_out
    !===============================================
    real(kind = r8)::psi
    real(kind = r8)::n_soil
    !===============================================
    real(kind = r8)::dFlow_downslope_in !from hillslope
    real(kind = r8)::dLength_hillslope_in
    real(kind = r8)::dLength_water_table_in
    real(kind = r8)::dSlope_intersect_in
    real(kind = r8)::dSlope_surface_in
    real(kind = r8)::dSlope_water_table_in
    !===============================================
    real(kind = r8)::phi_in !from seepage
    !===============================================
    real(kind = r8)::dummy0, dummy1
    !===============================================
    psi = this%psi
    n_soil = this%n_soil
    !===============================================
    if ( abs(dSlope_surface_in - dSlope_water_table_in) < near_zero) then
       dDrainage_downslope_out= -1 * dFlow_downslope_in * (n_soil+2) &
            * ( dLength_hillslope_in ** (n_soil+1) )/( dLength_water_table_in **(n_soil+2) )
    else if ( dSlope_water_table_in < dSlope_intersect_in &
         .and. dDistance_in < dLength_water_table_in) then
       dummy0 = dLength_water_table_in **2 / 2.0 - phi_in
       dummy1 = psi / dummy0
       dDrainage_downslope_out = -1 * dFlow_downslope_in * dDistance_in * dummy1
    else
       dDrainage_downslope_out = 0.0
    endif
    if (dDrainage_downslope_out > 0.0) then
       !write( fID_log, *) 'Something is wrong here!'
    end if
    this%dDrainage_downslope = dDrainage_downslope_out
    return
  end function calculate_drainage_downslope
  !===============================================
  !calculate seepage drainage
  !===============================================
  function calculate_drainage_seepage(this, &
       dDistance_in ,&
       dFlow_downslope_in,&
       dLength_hillslope_in,&
       dLength_water_table_in,&
       dSlope_intersect_in,&
       dSlope_surface_in,&
       dSlope_water_table_in, &
       phi_in,&
       dLength_seepage_in , &
       dFlow_seepage_in) &
       result(dDrainage_seepage_out)
    implicit none
    !===============================================
    class(soil_class), intent(inout):: this
    real(kind = r8)::dDistance_in
    real(kind = r8)::dDrainage_seepage_out
    !===============================================
    real(kind = r8)::n_soil
    real(kind = r8)::psi
    !===============================================
    real(kind = r8)::dFlow_downslope_in
    real(kind = r8)::dLength_hillslope_in
    real(kind = r8)::dLength_water_table_in
    real(kind = r8)::dSlope_intersect_in
    real(kind = r8)::dSlope_surface_in
    real(kind = r8)::dSlope_water_table_in
    !===============================================
    real(kind = r8)::phi_in
    real(kind = r8)::dLength_seepage_in
    real(kind = r8)::dFlow_seepage_in
    !===============================================
    real(kind = r8)::dummy0, dummy1, dummy2, dummy3
    !===============================================
    n_soil=this%n_soil
    psi=this%psi
    !===============================================
    !===============================================
    if(abs(dSlope_surface_in - dSlope_water_table_in) < near_zero) then
       dummy0 = -1 * (dFlow_downslope_in + dFlow_seepage_in) * (n_soil+1) * (n_soil+2)
       dummy1 = (dLength_hillslope_in - dDistance_in) * dDistance_in**n_soil
       dummy2 = (dLength_hillslope_in - dLength_seepage_in)**(n_soil+2)
       dDrainage_seepage_out = dummy0 * dummy1 / dummy2
    else if (dSlope_water_table_in < dSlope_intersect_in &
         .and. dDistance_in < dLength_water_table_in) then
       dummy0 = -1 * (dFlow_downslope_in + dFlow_seepage_in)
       dummy1 = (dLength_hillslope_in - dDistance_in) * psi
       dummy2 = (dLength_hillslope_in- dLength_seepage_in)**2 / 2.0
       dummy3 = dummy2 - phi_in
       dDrainage_seepage_out = dummy0 * dummy1 /dummy3
    endif
    if (dDrainage_seepage_out > 0.0) then
       !write( fID_log, *) 'Something is wrong here!'
    end if
    this%dDrainage_seepage = dDrainage_seepage_out
    return
  end function calculate_drainage_seepage
  !===============================================
  !this function is based on equation 1
  !===============================================
  function calculate_delta_h(this, &
       dHeight_below_river_in, &
       dSlope_bedrock_in,&
       dSlope_surface_in) &
       result(dDelta_h_out)
    implicit none
    !===============================================
    class(soil_column_class), intent(inout) :: this
    real(kind = r8) :: dDelta_h_out
    !===============================================
    real(kind = r8) :: dDistance
    real(kind = r8) :: dHeight_below_river_in
    real(kind = r8) :: dSlope_bedrock_in
    real(kind = r8) :: dSlope_surface_in
    real(kind = r8) :: dummy0, dummy1
    !===============================================
    dDistance = this%dDistance
    dummy0 = tan(dSlope_surface_in)
    dummy1 = tan(dSlope_bedrock_in)
    dDelta_h_out = dHeight_below_river_in + dDistance * (dummy0 - dummy1)
    this%dDelta_h = dDelta_h_out
    return
  end function calculate_delta_h
  !===============================================
  !calculate the total drainage
  !===============================================
  function calculate_soil_column_drainage(this,&
       iSeepage_type_in,&
       dDrainage_downslope_in,&
       dDrainage_et_in,&
       dDrainage_seepage_in,&
       dLength_seepage_in,&
       dLength_seepage_dominant_in,&
       dLength_water_table_in) &
       result(dDrainage_out)
    implicit none
    !===============================================
    class(soil_column_class), intent(inout):: this
    integer :: iSeepage_type_in
    real(kind = r8)::dDistance
    real(kind = r8)::dFactor_seepage
    real(kind = r8)::dDistance_in !the distance to the stream
    real(kind = r8)::dDrainage_out !the total drainage
    !===============================================
    !retrieve value from soil object
    !===============================================
    real(kind = r8)::dDrainage_downslope_in !the drainage from downslope
    real(kind = r8)::dDrainage_et_in !drainage from ET
    real(kind = r8)::dDrainage_seepage_in !the drainage from seepage
    real(kind = r8)::dLength_seepage_in !the total length of seepage
    real(kind = r8)::dLength_seepage_dominant_in
    real(kind = r8):: dLength_water_table_in
    real(kind = r8)::dummy0
    !===============================================
    dDistance = this%dDistance
    !mixture function
    select case(iSeepage_type_in)
    case (1) !no seepage
       dDrainage_out = dDrainage_downslope_in + dDrainage_et_in
    case (2) !partially seepage
       if (dDistance <= dLength_seepage_in ) then
          dDrainage_out = dDrainage_seepage_in + dDrainage_et_in
       else
          !calulate transition factor
          dummy0 = this%calculate_seepage_factor( dLength_seepage_dominant_in,&
               dLength_water_table_in)
          dFactor_seepage = this%dFactor_seepage
          dDrainage_out = dDrainage_downslope_in * (1.0 - dFactor_seepage) &
               + dDrainage_seepage_in* (dFactor_seepage) &
               + dDrainage_et_in
       endif
    case (3)
       !i have no idea yet
    end select
    if (dDrainage_out > 0.0) then
       !write( fID_log, *) 'Something is wrong here!'
    end if
    ! we decide to change the sign here
  end function calculate_soil_column_drainage

  function calculate_seepage_factor(this,&
       x,&
       dLength_water_table_in) &
       result(dFactor_seepage_out)
    implicit none
    !===============================================
    class(soil_column_class), intent(inout):: this
    real(kind = r8) ::x,c
    real(kind = r8) ::dFactor_seepage_out
    !===============================================
    real(kind = r8) :: dLength_water_table_in
    real(kind = r8) :: x0
    real(kind = r8) :: dummy0, dummy1, dummy2, dummy3, dummy4
    x0 = 3 * dLength_water_table_in /4.0
    dummy0 = -tan( 0.8 * pi /2 )
    dummy1 = (1.0/x0 + 1) / ( x0 - dLength_water_table_in )
    c = dummy0 / dummy1
    dummy2 = 1/x + 1 / ( x -dLength_water_table_in)
    dummy3 = -c * dummy2
    dummy4 = 1 + 2.0 / pi * atan( dummy3 )
    dFactor_seepage_out = 0.5 * dummy4
    !===============================================
    this%dFactor_seepage = dFactor_seepage_out
    !===============================================
  end function calculate_seepage_factor
  function run_soil_column(this, &
       iSeepage_type_in,&
       dDrainage_downslope_in,&
       dDrainage_et_in,&
       dDrainage_seepage_in,&
       dElevation_max_in,&
       dElevation_min_in,&
       dElevation_water_table_in,&
       dFlow_downslope_in,&
       dFlow_seepage_in, &
       dHeight_below_river_in, &
       dLength_hillslope_in, &
       dLength_seepage_dominant_in,&
       dLength_seepage_in,&
       dLength_water_table_in, &
       dSlope_intersect_in, &
       dSlope_surface_in, &
       dSlope_water_table_in,&
       dThickness_critical_zone_in, &
       phi_in) &
       result(error_code)
    implicit none
    class(soil_column_class), intent(inout):: this
    integer:: error_code
    integer:: iSeepage_type_in
    real(kind = r8)::dDistance
    real(kind = r8)::dElevation
    real(kind = r8):: dDrainage_downslope_in
    real(kind = r8):: dDrainage_et_in
    real(kind = r8):: dDrainage_seepage_in
    real(kind = r8):: dElevation_max_in
    real(kind = r8):: dElevation_min_in
    real(kind = r8):: dElevation_water_table_in
    real(kind = r8):: dFlow_downslope_in
    real(kind = r8):: dFlow_seepage_in
    real(kind = r8):: dHeight_below_river_in
    real(kind = r8):: dLength_hillslope_in
    real(kind = r8):: dLength_seepage_dominant_in
    real(kind = r8):: dLength_seepage_in
    real(kind = r8):: dLength_water_table_in
    real(kind = r8):: dSlope_intersect_in
    real(kind = r8):: dSlope_surface_in
    real(kind = r8):: dSlope_water_table_in
    real(kind = r8):: dThickness_critical_zone_in
    real(kind = r8):: phi_in
    real(kind = r8)::dummy0 , dummy1, dummy2, dummy3
    !===============================================
    error_code = 1
    dDistance= this%dDistance
    dElevation = this%dElevation
    dummy2 = this%cGroundwater%calculate_groundwater_level( iSeepage_type_in, &
         dElevation, &
         dElevation_max_in, &
         dElevation_min_in, &
         dElevation_water_table_in, &
         dThickness_critical_zone_in)
    error_code = this%cSoil%run_soil(    dDistance,&
         dElevation_water_table_in ,&
         dFlow_downslope_in, &
         dFlow_seepage_in, &
         dHeight_below_river_in, &
         dLength_hillslope_in, &
         dLength_seepage_in,&
         dLength_water_table_in, &
         dSlope_intersect_in, &
         dSlope_surface_in, &
         dSlope_water_table_in,&
         phi_in)
    dummy1 = this%calculate_soil_column_drainage( iSeepage_type_in,&
         dDrainage_downslope_in,&
         dDrainage_et_in,&
         dDrainage_seepage_in,&
         dLength_seepage_in,&
         dLength_seepage_dominant_in,&
         dLength_water_table_in)
    return
  end function run_soil_column
  !initialize the h2sc module

  subroutine elm_Init_h2sc(bounds)
    use clm_instMod, only: rof2lnd_vars
    use clm_varpar , only : nlevgrnd, nlayer, nlayert, nlevsoi
    use clm_varcon , only : zsoi, dzsoi, zisoi, spval
    use decompmod , only: bounds_type
    use ColumnType , only : col_pp
    type(bounds_type) , intent(in) :: bounds ! bounds
    integer :: g, i ,l ,c ! indices
    integer :: ier
    integer :: nElevation_interval
    real(r8) :: dThickness_critical_zone
    nElevation_interval=11
    dThickness_critical_zone = 0.0
    do i=1, nlevgrnd
       dThickness_critical_zone = dThickness_critical_zone + dzsoi(i)
    end do
    allocate(cHillslope(bounds%begg:bounds%endg),stat=ier)
    allocate(cSeepage(bounds%begc:bounds%endc),stat=ier)
    do g = bounds%begg,bounds%endg
       cHillslope(g)%nElevation_interval = nElevation_interval
       cHillslope(g)%dThickness_critical_zone = dThickness_critical_zone
       cHillslope(g)%beta_hillslope = beta_h2sc
       cHillslope(g)%m_hillslope = m_h2sc
       cHillslope(g)%n_hillslope = n_h2sc
       allocate( cHillslope(g)%aElevation_profile( nElevation_interval ) )
    end do
  end subroutine elm_Init_h2sc

  subroutine InitHillslopeAllocate(this)
    implicit none
    class(hillslope_class) :: this
  end subroutine InitHillslopeAllocate
  subroutine InitHillslopeCold(this)
    implicit none
    class(hillslope_class) :: this
    return
  end subroutine InitHillslopeCold

  function initialize_soil_column_class( dDistance_in, &
       dElevation_in) &
       result(cSoil_column_out)
    implicit none
    !===============================================
    real(kind = r8) :: dDistance_in
    real(kind = r8) :: dElevation_in
    type(soil_column_class) :: cSoil_column_out
    !===============================================
    cSoil_column_out%dDistance = dDistance_in
    cSoil_column_out%dElevation = dElevation_in
    return
  end function initialize_soil_column_class
  function initialize_soil_class( dSaturate_hydraulic_conductivity_in , &
       beta_soil_in, &
       m_soil_in, &
       n_soil_in ) &
       result(cSoil_out)
    implicit none
    !===============================================
    real(kind = r8)::dSaturate_hydraulic_conductivity_in
    real(kind = r8)::beta_soil_in
    real(kind = r8)::m_soil_in
    real(kind = r8)::n_soil_in
    !===============================================
    type(soil_class)::cSoil_out
    !===============================================
    cSoil_out%dSaturate_hydraulic_conductivity = dSaturate_hydraulic_conductivity_in
    cSoil_out%beta_soil = beta_soil_in
    cSoil_out%m_soil = m_soil_in
    cSoil_out%n_soil = n_soil_in
    return
  end function initialize_soil_class
  function initialize_groundwater_class( dElevation_water_table_in) &
       result(cGroundwater_out)
    implicit none
    real(kind = r8)::dElevation_water_table_in
    type(groundwater_class)::cGroundwater_out
    cGroundwater_out%dElevation_water_table = dElevation_water_table_in
    return
  end function initialize_groundwater_class

end module H2SC_drainage
