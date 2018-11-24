!-----------------------------------------------------------------------------------------
! AtmosRefraction.f95 :compute atmospheric refraction
!
!          Copyright (C) 2018 by DING Junsheng, All rights reserved.
!
! version :$Revision: 1.1 $ $Date: 2018/10/30 21:48:06 $
! history : 2018/10/30  1.0  new
!-----------------------------------------------------------------------------------------
module constant
  implicit none
    real(8),parameter::pi     = 3.141592653589793238462643383280
    real(8),parameter::epsilo = 1e-6
end module constant

module AtmosRefraction
  use constant
  implicit none

contains
! degrees to radians-----------------------------------------------------------------------
real(8) function deg2rad(angle)
  implicit none
  real*8::angle
    deg2rad = pi*angle/180
end function deg2rad

! radians to degrees------------------------------------------------------------------------
real(8) function rad2deg(angle)
  implicit none
  real*8::angle
    rad2deg=180*angle/pi
end function rad2deg

! azimuth and altitude angle----------------------------------------------------------------
! compute azimuth and altitude angle in the observed direction
! args    : double press       i   temperature at the station(degrees centigrade)
!           double temp        i   pressure at the station(hpa)
!           double azimuth     i   azimuth at station in horizontal coordinate system
!           double altangle    i   altitude angle at station in horizontal coordinate system
! return  : double atmref      o   atmospheric refraction
! ------------------------------------------------------------------------------------------
real(8) function getobsangle(press,temp,azimuth,altangle)
  implicit none
  integer::n=0
  real*8 ::r0,r1,tag,press,temp,azimuth,altangle,a=59.92,b=0.0665
    r0=altangle
    r1=r0-deg2rad((a*tan(r0)-b*tan(r0)**3)/3600)

    do while(abs(r1-r0)>deg2rad(epsilo))
        n=n+1
		write(*,*) abs(r0-r1),n
		  tag=(r0+r1)/2
		if((tag+deg2rad((a*tan(r0)-b*tan(r0)**3)/3600))>altangle)then
          r0=tag
        else
          r1=tag
        endif
    end do

    getobsangle=deg2rad((r0*(press/1013.25)*(273.15/(273.15+temp)))/3600)

end function getobsangle

end module AtmosRefraction

program main
!implicit none
use AtmosRefraction
    real*8::get_obsangle,azimuth,altangle,press=980.0,temp=20.0
	real*8::a_zimuth=30.0,a_ltangle=50.0
	real*8::aa,bb,cc,dd
	  azimuth=deg2rad(a_zimuth)
	  altangle=deg2rad(90-a_ltangle)
	write(*,*) "     gap of iterations     number of iterations"
      get_obsangle=getobsangle(press,temp,azimuth,altangle)

    write(*,'("AtmosRefraction：",F12.9," μrad = ",F6.3," ”")') 1e6*get_obsangle,3600*rad2deg(get_obsangle)
    write(*,'("obsDirection：",F10.7," ° ","obsHeight：",F10.7," °")') &
    INT(3.6E6*a_zimuth+0.5)/3.6E6,INT(3.6E6*(a_ltangle+rad2deg(get_obsangle))+0.5)/3600000.0D0

!--------------------------------------------------------------------------------------------------------------
!   write(*,'(F20.9)') 180000629.0D0 / 3.6E6
!   write(*,'(F20.9)') 180.0006290D0 / 3.6
!   write(*,'(F20.9)') 18000062.90D0 / 360000

!   write(*,'(F20.9)') 629.0       / 3.6E6
!   write(*,'(F20.9)') 180000000.0 / 3.6E6
!   write(*,'(F20.9)') 629.0/3.6E6 + 180000000.0D0/3.6E6

!   aa=629.0/3.6E6
!   bb=180000000.0/3.6E6
!   write(*,'(F20.9)') aa
!   write(*,'(F20.9)') bb
!   write(*,'(F20.9)') aa+bb

!   cc=180000629.0D0
!   dd=3.6E6
!   write(*,*) kind(cc),huge(cc),tiny(cc)
!   write(*,'(F20.9)') cc
!   write(*,'(F20.9)') dd
!   write(*,'(F20.9)') cc/dd
!--------------------------------------------------------------------------------------------------------------

end
