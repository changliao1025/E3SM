module SoilTextureMod
implicit none
 !--------------------------------------------------------------------------
 !                            INPOLY
 !   Function to tell if a point is inside a polygon or not.
 !--------------------------------------------------------------------------
 !   Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.
 !
 !   Please feel free to use this source code for any purpose, commercial
 !   or otherwise, as long as you don't restrict anyone else's use of
 !   this source code.  Please give credit where credit is due.
 !
 !   Point-in-polygon algorithm, created especially for World-Wide Web
 !   servers to process image maps with mouse-clickable regions.
 !
 !   Home for this file:  http://www.gcomm.com/develop/inpoly.c
 !
 !                                       6/19/95 - Bob Stein & Craig Yap
 !                                       stein@gcomm.com
 !                                       craig@cse.fau.edu
 !--------------------------------------------------------------------------
 !   Modified by:
 !   Aris Gerakis, apr. 1998: 1.  translated to Fortran
 !                            2.  made it work with real coordinates
 !                            3.  now resolves the case where point falls
 !                                on polygon border.
 !   Aris Gerakis, nov. 1998: Fixed error caused by hardware arithmetic
 !   Aris Gerakis, july 1999: Now all borderline cases are valid
 !--------------------------------------------------------------------------
 !   Glossary:
 !   function inpoly: true=inside, false=outside (is target point inside
 !                    a 2D polygon?)
 !   poly(*,2):  polygon points, [0]=x, [1]=y
 !   npoints: number of points in polygon
 !   xt: x (horizontal) of target point
 !   yt: y (vertical) of target point
 !--------------------------------------------------------------------------
contains

logical function inpoly (poly, npoints, xt, yt)

implicit none

! Declare arguments:

integer :: npoints
real    :: poly(7, 2), xt, yt

! Declare local variables:

real    :: xnew, ynew, xold, yold, x1, y1, x2, y2
integer :: i
logical :: inside, on_border

inside = .false.
on_border = .false.

if (npoints < 3)  then
  inpoly = .false.
  return
end if

xold = poly(npoints,1)
yold = poly(npoints,2)

do i = 1 , npoints
  xnew = poly(i,1)
  ynew = poly(i,2)

  if (xnew > xold)  then
    x1 = xold
    x2 = xnew
    y1 = yold
    y2 = ynew
  else
    x1 = xnew
    x2 = xold
    y1 = ynew
    y2 = yold
  end if

  ! The outer IF is the 'straddle' test and the 'vertical border' test.
  ! The inner IF is the 'non-vertical border' test and the 'north' test.  

  ! The first statement checks whether a north pointing vector crosses  
  ! (stradles) the straight segment.  There are two possibilities, depe-
  ! nding on whether xnew < xold or xnew > xold.  The '<' is because edge 
  ! must be "open" at left, which is necessary to keep correct count when 
  ! vector 'licks' a vertix of a polygon.  

  if ((xnew < xt .and. xt <= xold) .or. (.not. xnew < xt .and. &
     .not. xt <= xold)) then
    ! The test point lies on a non-vertical border:
    if ((yt-y1)*(x2-x1) == (y2-y1)*(xt-x1)) then
      on_border = .true. 
    ! Check if segment is north of test point.  If yes, reverse the 
    ! value of INSIDE.  The +0.001 was necessary to avoid errors due   
    ! arithmetic (e.g., when clay = 98.87 and sand = 1.13):   
    elseif ((yt-y1)*(x2-x1) < (y2-y1)*(xt-x1) + 0.001) then
      inside = .not.inside ! cross a segment
    endif
  ! This is the rare case when test point falls on vertical border or  
  ! left edge of non-vertical border. The left x-coordinate must be  
  ! common.  The slope requirement must be met, but also point must be
  ! between the lower and upper y-coordinate of border segment.  There 
  ! are two possibilities,  depending on whether ynew < yold or ynew > 
  ! yold:
  elseif ((xnew == xt .or. xold == xt) .and. (yt-y1)*(x2-x1) == &
    (y2-y1)*(xt-x1) .and. ((ynew <= yt .and. yt <= yold) .or. &
    (.not. ynew < yt .and. .not. yt < yold))) then
    on_border = .true. 
  endif

  xold = xnew
  yold = ynew

enddo

! If test point is not on a border, the function result is the last state 
! of INSIDE variable.  Otherwise, INSIDE doesn't matter.  The point is
! inside the polygon if it falls on any of its borders:

if (.not. on_border) then
   inpoly = inside
else
   inpoly = .true.
endif

return

end function inpoly


!* +-----------------------------------------------------------------------
!* | WHAT TEXTURE?
!* | Function to classify a soil in the triangle based on sand and clay %
!* +-----------------------------------------------------------------------
!* | Created by: aris gerakis, apr. 98
!* | Modified by: aris gerakis, june 99.  Now check all polygons instead of
!* | stopping when a right solution is found.  This to cover all borderline 
!* | cases.
!* +-----------------------------------------------------------------------

character*42 function what_soil_texture (sand, clay)

implicit none

! Declare arguments:

real, intent(in) :: clay, sand

! Declare local variables:

!logical   :: inpoly 
real      :: silty_loam(7,2), sandy(7,2), silty_clay_loam(7,2), loam(7,2), &
             clay_loam(7,2), sandy_loam(7,2), silty_clay(7,2), &
             sandy_clay_loam(7,2), loamy_sand(7,2), clayey(7,2), silt(7,2), &
             sandy_clay(7,2)
character :: texture*42

!Initalize polygon coordinates:
! First write all x's, pad with 0s if less than 7 vertices, then
! write all y's, pad with 0s if less than 7 vertices:

data silty_loam/0, 0, 23, 50, 20, 8, 0, 12, 27, 27, 0, 0, 12, 0/
data sandy/85, 90, 100, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0/
data silty_clay_loam/0, 0, 20, 20, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
data loam/43, 23, 45, 52, 52, 0, 0, 7, 27, 27, 20, 7, 0, 0/
data clay_loam/20, 20, 45, 45, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
data sandy_loam/50, 43, 52, 52, 80, 85, 70, 0, 7, 7, 20, 20, 15, 0/
data silty_clay/0, 0, 20, 0, 0, 0, 0, 40, 60, 40, 0, 0, 0, 0/
data sandy_clay_loam/52, 45, 45, 65, 80, 0, 0, 20, 27, 35, 35, 20, 0, 0/
data loamy_sand/70, 85, 90, 85, 0, 0, 0, 0, 15, 10, 0, 0, 0, 0/
data clayey/20, 0, 0, 45, 45, 0, 0, 40, 60, 100, 55, 40, 0, 0/
data silt/0, 0, 8, 20, 0, 0, 0, 0, 12, 12, 0, 0, 0, 0/
data sandy_clay/45, 45, 65, 0, 0, 0, 0, 35, 55, 35, 0, 0, 0, 0/

! Find polygon(s) where the point is.

texture = ' '

if (sand .ge. 0.0 .and. clay .ge. 0.0) then
   if (inpoly(silty_loam, 6, sand, clay)) then
      texture = trim(texture)//'/silt loam'
   endif
   if (inpoly(sandy, 3, sand, clay)) then
      texture = trim(texture)//'/sand'
   endif
   if (inpoly(silty_clay_loam, 4, sand, clay)) then
      texture = trim(texture)//'/silty clay loam'
   endif
   if (inpoly(loam, 5, sand, clay)) then
      texture = trim(texture)//'/loam'
   endif
   if (inpoly(clay_loam, 4, sand, clay)) then
      texture = trim(texture)//'/clay loam'
   endif
   if (inpoly(sandy_loam, 7, sand, clay)) then
      texture = trim(texture)//'/sandy loam'
   endif
   if (inpoly(silty_clay, 3, sand, clay)) then
      texture = trim(texture)//'/silty clay'
   endif
   if (inpoly(sandy_clay_loam, 5, sand, clay)) then
      texture = trim(texture)//'/sandy clay loam'
   endif
   if (inpoly(loamy_sand, 4, sand, clay)) then
      texture = trim(texture)//'/loamy sand'
   endif
   if (inpoly(clayey, 5, sand, clay)) then
      texture = trim(texture)//'/clay'
   endif
   if (inpoly(silt, 4, sand, clay)) then
      texture = trim(texture)//'/silt'
   endif
   if (inpoly(sandy_clay, 3, sand, clay)) then
      texture = trim(texture)//'/sandy clay'
   endif
endif

if (texture == ' ') then
   texture = '/unknown texture'
   !write (unit = 6, fmt = 1000) sand, clay
   !1000 format (/, 1x, 'Texture not found for ', f5.1, ' sand and ', f5.1, ' clay')

   !write(*, *) 'Texture not found for: ',   sand, clay
else 
   !write(*, *) 'Texture is : ',texture
end if

what_soil_texture = texture

return
end function what_soil_texture

character*42 function what_soil_texture_single (sand, clay)

implicit none

! Declare arguments:

real, intent(in) :: clay, sand

! Declare local variables:

!logical   :: inpoly 
real      :: silty_loam(7,2), sandy(7,2), silty_clay_loam(7,2), loam(7,2), &
             clay_loam(7,2), sandy_loam(7,2), silty_clay(7,2), &
             sandy_clay_loam(7,2), loamy_sand(7,2), clayey(7,2), silt(7,2), &
             sandy_clay(7,2)
character :: texture*42

!Initalize polygon coordinates:
! First write all x's, pad with 0s if less than 7 vertices, then
! write all y's, pad with 0s if less than 7 vertices:

data silty_loam/0, 0, 23, 50, 20, 8, 0, 12, 27, 27, 0, 0, 12, 0/
data sandy/85, 90, 100, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0/
data silty_clay_loam/0, 0, 20, 20, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
data loam/43, 23, 45, 52, 52, 0, 0, 7, 27, 27, 20, 7, 0, 0/
data clay_loam/20, 20, 45, 45, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
data sandy_loam/50, 43, 52, 52, 80, 85, 70, 0, 7, 7, 20, 20, 15, 0/
data silty_clay/0, 0, 20, 0, 0, 0, 0, 40, 60, 40, 0, 0, 0, 0/
data sandy_clay_loam/52, 45, 45, 65, 80, 0, 0, 20, 27, 35, 35, 20, 0, 0/
data loamy_sand/70, 85, 90, 85, 0, 0, 0, 0, 15, 10, 0, 0, 0, 0/
data clayey/20, 0, 0, 45, 45, 0, 0, 40, 60, 100, 55, 40, 0, 0/
data silt/0, 0, 8, 20, 0, 0, 0, 0, 12, 12, 0, 0, 0, 0/
data sandy_clay/45, 45, 65, 0, 0, 0, 0, 35, 55, 35, 0, 0, 0, 0/

! Find polygon(s) where the point is.

texture = ' '

if (sand .ge. 0.0 .and. clay .ge. 0.0) then
   if (inpoly(silty_loam, 6, sand, clay)) then
      texture = 'silt loam'
   else 
      if (inpoly(sandy, 3, sand, clay)) then
         texture = 'sand'
      else
         if (inpoly(silty_clay_loam, 4, sand, clay)) then
            texture = 'silty clay loam'
         else
            if (inpoly(loam, 5, sand, clay)) then
               texture = 'loam'
            else
               if (inpoly(clay_loam, 4, sand, clay)) then
                  texture = 'clay loam'
               else
                  if (inpoly(sandy_loam, 7, sand, clay)) then
                     texture = 'sandy loam'
                  else
                     if (inpoly(silty_clay, 3, sand, clay)) then
                        texture = 'silty clay'
                        else
                           if (inpoly(sandy_clay_loam, 5, sand, clay)) then
                              texture = 'sandy clay loam'
                           else
                              if (inpoly(loamy_sand, 4, sand, clay)) then
                                 texture = 'loamy sand'
                              else
                                 if (inpoly(clayey, 5, sand, clay)) then
                                    texture = 'clay'
                                 else
                                    if (inpoly(silt, 4, sand, clay)) then
                                       texture = 'silt'
                                    else
                                       if (inpoly(sandy_clay, 3, sand, clay)) then
                                          texture = 'sandy clay'
                                       else 
                                          texture = 'silt'
                                       endif
                                    endif
                                 endif
                              endif
                           endif
                        endif
                  endif
               endif
            endif
         endif
      endif
   endif
   
endif

!if (texture == ' ') then
!   texture = '/unknown texture'
!   !write (unit = 6, fmt = 1000) sand, clay
!   !1000 format (/, 1x, 'Texture not found for ', f5.1, ' sand and ', f5.1, ' clay')
!
!   write(*, *) 'Texture not found for: ',   sand, clay
!else 
!   write(*, *) 'Texture is : ',texture
!end if

what_soil_texture_single = texture

return
end function what_soil_texture_single


end module SoilTextureMod