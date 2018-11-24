/*--------------------------------------------------------------------------------------
* AtmosRefraction.c : compute atmospheric refraction
*
*          Copyright (C) 2018 by DING Junsheng, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2018/10/28 21:48:06 $
* history : 2018/10/28  1.0  new
*--------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>

/* constants --------------------------------------------------------------------------*/
#define epsilon     1E-6            /* Iteration exit threshold(degree) */
#define PI          3.141592653589793238462643383280

/* Degrees to Radians -----------------------------------------------------------------*/
double deg2rad(double angle)
{
    return PI*angle/180;
}
/* Radians to Degrees -----------------------------------------------------------------*/
double rad2deg(double angle)
{
    return angle*180/PI;
}
/* azimuth and altitude Angle-----------------------------------------------------------
* compute azimuth and altitude Angle in the observed direction
* args    : double PRESS       I   temperature at the station(degrees centigrade)
*           double TEMP        I   pressure at the station(hpa)
*           double azimuth     I   azimuth at station in horizontal coordinate system
*           double AltAngle    I   altitude Angle at station in horizontal coordinate system
* return  : double AtmRef      O   atmospheric refraction
*--------------------------------------------------------------------------------------*/
double GetObsAngle(double PRESS,double TEMP,double azimuth,double AltAngle)
{
    double A=59.92,B=0.0665;
    double R0=AltAngle;
    double R1=R0-deg2rad((A*tan(R0)-B*pow(tan(R0),3))/3600);
	double tag;
	int    n=0;
    while(fabs(R1-R0)>deg2rad(epsilon))
    {
		printf("%18.12f%10d\n",fabs(R0-R1),++n);		
		tag=(R0+R1)/2;
		(tag+deg2rad((A*tan(tag)-B*pow(tan(tag),3))/3600))>AltAngle ? R0=tag : R1=tag;
    }
    return deg2rad((R0*(PRESS/1013.25)*(273.15/(273.15+TEMP)))/3600);
}

int main()
{
	double azimuth_d=30,AltAngle_d=50;
	double azimuth=deg2rad(azimuth_d);
	double AltAngle=deg2rad(90-AltAngle_d);
	double PRESS=980;
	double TEMP=20;
	printf("     两次迭代变化     迭代次数\n");
	double AtRef=GetObsAngle(PRESS,TEMP,azimuth,AltAngle);
    printf("\n蒙气差:%10.7f μrad = %7.3f \"\n",1E6*AtRef,3600*rad2deg(AtRef));
	printf("\n观测方向方位:%15.7f ° 观测方向高度角:%15.7f °\n\n",
		int(3.6E6*azimuth_d+0.5)/3.6E6,int(3.6E6*(AltAngle_d+rad2deg(AtRef))+0.5)/3.6E6);
	return 0;
}
