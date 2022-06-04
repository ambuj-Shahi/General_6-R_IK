/** 
 

 @file theta2.c
 @brief This file is used to calculate the value of theta2 with the values of theta1, theta3, theta4, theta5, and the D-H parameters  as input.
       .
 
 */


# include <stdio.h>
# include <math.h>
# include <complex.h>
# include "cx.h"



struct Theta theta2(double* coeff, double _Complex *thlist)
{
    
    struct Theta sctheta2;

    double a1=coeff[0];
    double a2=coeff[1];
    double a3=coeff[2];
    double a4=coeff[3];
    double a5=coeff[4];
    double a6=coeff[5];

    double alpha1=coeff[6];
    double alpha2=coeff[7];
    double alpha3=coeff[8];
    double alpha4=coeff[9];
    double alpha5=coeff[10];
    double alpha6=coeff[11];

    double d1=coeff[12];
    double d2=coeff[13];
    double d3=coeff[14];
    double d4=coeff[15];
    double d5=coeff[16];
    double d6=coeff[17];

    double lambda1=coeff[18];
    double lambda2=coeff[19];
    double lambda3=coeff[20];
    double lambda4=coeff[21];
    double lambda5=coeff[22];
    double lambda6=coeff[23];
    double lx=coeff[24];
    double ly=coeff[25];
    double lz=coeff[26];
    double mu1=coeff[27];
    double mu2=coeff[28];
    double mu3=coeff[29];
    double mu4=coeff[30];
    double mu5=coeff[31];
    double mu6=coeff[32];
    double mx=coeff[33];
    double my=coeff[34];
    double mz=coeff[35];
    double nx=coeff[36];
    double ny=coeff[37];
    double nz=coeff[38];
    double px=coeff[39];
    double py=coeff[40];
    double pz=coeff[41];
	
	  double _Complex theta1r = thlist[0];
    double _Complex theta3r = thlist[1];
    double _Complex theta4r = thlist[2];
    double _Complex theta5r = thlist[3];




  double _Complex costheta2 = -(cpow(cpow(-a1 - (a6*lx + d6*mu6*mx + d6*lambda6*nx - px)*ccos(theta1r) - (a6*ly + d6*mu6*my + d6*lambda6*ny - py)*csin(theta1r),2) - 
       (-(mu1*(d1 + a6*lz + d6*mu6*mz + d6*lambda6*nz - pz)) - 
          lambda1*((a6*ly + d6*mu6*my + d6*lambda6*ny - py)*ccos(theta1r) - (a6*lx + d6*mu6*mx + d6*lambda6*nx - px)*csin(theta1r)))*
        (mu1*(d1 + a6*lz + d6*mu6*mz + d6*lambda6*nz - pz) + lambda1*
           ((a6*ly + d6*mu6*my + d6*lambda6*ny - py)*ccos(theta1r) - (a6*lx + d6*mu6*mx + d6*lambda6*nx - px)*csin(theta1r))),-1)*
     ((-a1 - (a6*lx + d6*mu6*mx + d6*lambda6*nx - px)*ccos(theta1r) - (a6*ly + d6*mu6*my + d6*lambda6*ny - py)*csin(theta1r))*
        (-a2 - ccos(theta3r)*(a3 + ccos(theta4r)*(a4 + a5*ccos(theta5r)) + d5*mu4*csin(theta4r) - a5*lambda4*csin(theta4r)*csin(theta5r)) - 
          csin(theta3r)*(mu3*(d4 + d5*lambda4 + a5*mu4*csin(theta5r)) + 
             lambda3*(-((a4 + a5*ccos(theta5r))*csin(theta4r)) + ccos(theta4r)*(d5*mu4 - a5*lambda4*csin(theta5r))))) - 
       (-(mu1*(d1 + a6*lz + d6*mu6*mz + d6*lambda6*nz - pz)) - 
          lambda1*((a6*ly + d6*mu6*my + d6*lambda6*ny - py)*ccos(theta1r) - (a6*lx + d6*mu6*mx + d6*lambda6*nx - px)*csin(theta1r)))*
        (-(mu2*(d3 + a4*mu3*csin(theta4r) + a5*mu3*ccos(theta5r)*csin(theta4r) + mu3*ccos(theta4r)*(-(d5*mu4) + a5*lambda4*csin(theta5r)) + 
               lambda3*(d4 + d5*lambda4 + a5*mu4*csin(theta5r)))) - 
          lambda2*(-(csin(theta3r)*(a3 + ccos(theta4r)*(a4 + a5*ccos(theta5r)) + d5*mu4*csin(theta4r) - a5*lambda4*csin(theta4r)*csin(theta5r))) + 
             ccos(theta3r)*(mu3*(d4 + d5*lambda4 + a5*mu4*csin(theta5r)) + 
                lambda3*(-((a4 + a5*ccos(theta5r))*csin(theta4r)) + ccos(theta4r)*(d5*mu4 - a5*lambda4*csin(theta5r))))))));

double _Complex sintheta2 = -(cpow(2*a1*a6*lx*ccos(theta1r) + 2*a6*d1*lambda1*ly*mu1*ccos(theta1r) + 2*a1*d6*mu6*mx*ccos(theta1r) + 2*d1*d6*lambda1*mu1*mu6*my*ccos(theta1r) + 
       2*a6*d6*lambda1*lz*mu1*mu6*my*ccos(theta1r) + 2*a6*d6*lambda1*ly*mu1*mu6*mz*ccos(theta1r) + 2*a1*d6*lambda6*nx*ccos(theta1r) + 
       2*d1*d6*lambda1*lambda6*mu1*ny*ccos(theta1r) + 2*a6*d6*lambda1*lambda6*lz*mu1*ny*ccos(theta1r) + 2*a6*d6*lambda1*lambda6*ly*mu1*nz*ccos(theta1r) - 
       2*a1*px*ccos(theta1r) - 2*d1*lambda1*mu1*py*ccos(theta1r) - 2*a6*lambda1*lz*mu1*py*ccos(theta1r) - 2*d6*lambda1*mu1*mu6*mz*py*ccos(theta1r) - 
       2*d6*lambda1*lambda6*mu1*nz*py*ccos(theta1r) - 2*a6*lambda1*ly*mu1*pz*ccos(theta1r) - 2*d6*lambda1*mu1*mu6*my*pz*ccos(theta1r) - 
       2*d6*lambda1*lambda6*mu1*ny*pz*ccos(theta1r) + 2*lambda1*mu1*py*pz*ccos(theta1r) + cpow(a1,2) + 2*lambda1*ly*lz*mu1*ccos(theta1r)*cpow(a6,2) + 
       2*lambda1*lambda6*mu1*mu6*mz*ny*ccos(theta1r)*cpow(d6,2) + 2*lambda1*lambda6*mu1*mu6*my*nz*ccos(theta1r)*cpow(d6,2) + 
       2*lambda1*mu1*ny*nz*ccos(theta1r)*cpow(d6,2)*cpow(lambda6,2) + 2*a6*d1*lz*cpow(mu1,2) + 2*d1*d6*mu6*mz*cpow(mu1,2) + 2*a6*d6*lz*mu6*mz*cpow(mu1,2) + 
       2*d1*d6*lambda6*nz*cpow(mu1,2) + 2*a6*d6*lambda6*lz*nz*cpow(mu1,2) - 2*d1*pz*cpow(mu1,2) - 2*a6*lz*pz*cpow(mu1,2) - 2*d6*mu6*mz*pz*cpow(mu1,2) - 
       2*d6*lambda6*nz*pz*cpow(mu1,2) + cpow(d1,2)*cpow(mu1,2) + 2*lambda6*mu6*mz*nz*cpow(d6,2)*cpow(mu1,2) + cpow(a6,2)*cpow(lz,2)*cpow(mu1,2) + 
       2*lambda1*mu1*my*mz*ccos(theta1r)*cpow(d6,2)*cpow(mu6,2) + cpow(d6,2)*cpow(mu1,2)*cpow(mu6,2)*cpow(mz,2) + cpow(d6,2)*cpow(lambda6,2)*cpow(mu1,2)*cpow(nz,2) + 
       cpow(mu1,2)*cpow(pz,2) + 2*a6*d6*lx*mu6*mx*cpow(ccos(theta1r),2) + 2*a6*d6*lambda6*lx*nx*cpow(ccos(theta1r),2) - 2*a6*lx*px*cpow(ccos(theta1r),2) - 
       2*d6*mu6*mx*px*cpow(ccos(theta1r),2) - 2*d6*lambda6*nx*px*cpow(ccos(theta1r),2) + 2*lambda6*mu6*mx*nx*cpow(d6,2)*cpow(ccos(theta1r),2) + 
       2*a6*d6*ly*mu6*my*cpow(lambda1,2)*cpow(ccos(theta1r),2) + 2*a6*d6*lambda6*ly*ny*cpow(lambda1,2)*cpow(ccos(theta1r),2) - 
       2*a6*ly*py*cpow(lambda1,2)*cpow(ccos(theta1r),2) - 2*d6*mu6*my*py*cpow(lambda1,2)*cpow(ccos(theta1r),2) - 
       2*d6*lambda6*ny*py*cpow(lambda1,2)*cpow(ccos(theta1r),2) + 2*lambda6*mu6*my*ny*cpow(d6,2)*cpow(lambda1,2)*cpow(ccos(theta1r),2) + 
       cpow(a6,2)*cpow(lx,2)*cpow(ccos(theta1r),2) + cpow(a6,2)*cpow(lambda1,2)*cpow(ly,2)*cpow(ccos(theta1r),2) + 
       cpow(d6,2)*cpow(mu6,2)*cpow(mx,2)*cpow(ccos(theta1r),2) + cpow(d6,2)*cpow(lambda1,2)*cpow(mu6,2)*cpow(my,2)*cpow(ccos(theta1r),2) + 
       cpow(d6,2)*cpow(lambda6,2)*cpow(nx,2)*cpow(ccos(theta1r),2) + cpow(d6,2)*cpow(lambda1,2)*cpow(lambda6,2)*cpow(ny,2)*cpow(ccos(theta1r),2) + 
       cpow(px,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(py,2)*cpow(ccos(theta1r),2) + 2*a6*d6*ly*mu6*my*cpow(csin(theta1r),2) + 
       2*a6*d6*lambda6*ly*ny*cpow(csin(theta1r),2) - 2*a6*ly*py*cpow(csin(theta1r),2) - 2*d6*mu6*my*py*cpow(csin(theta1r),2) - 
       2*d6*lambda6*ny*py*cpow(csin(theta1r),2) + 2*lambda6*mu6*my*ny*cpow(d6,2)*cpow(csin(theta1r),2) + 2*a6*d6*lx*mu6*mx*cpow(lambda1,2)*cpow(csin(theta1r),2) + 
       2*a6*d6*lambda6*lx*nx*cpow(lambda1,2)*cpow(csin(theta1r),2) - 2*a6*lx*px*cpow(lambda1,2)*cpow(csin(theta1r),2) - 
       2*d6*mu6*mx*px*cpow(lambda1,2)*cpow(csin(theta1r),2) - 2*d6*lambda6*nx*px*cpow(lambda1,2)*cpow(csin(theta1r),2) + 
       2*lambda6*mu6*mx*nx*cpow(d6,2)*cpow(lambda1,2)*cpow(csin(theta1r),2) + cpow(a6,2)*cpow(lambda1,2)*cpow(lx,2)*cpow(csin(theta1r),2) + 
       cpow(a6,2)*cpow(ly,2)*cpow(csin(theta1r),2) + cpow(d6,2)*cpow(lambda1,2)*cpow(mu6,2)*cpow(mx,2)*cpow(csin(theta1r),2) + 
       cpow(d6,2)*cpow(mu6,2)*cpow(my,2)*cpow(csin(theta1r),2) + cpow(d6,2)*cpow(lambda1,2)*cpow(lambda6,2)*cpow(nx,2)*cpow(csin(theta1r),2) + 
       cpow(d6,2)*cpow(lambda6,2)*cpow(ny,2)*cpow(csin(theta1r),2) + cpow(lambda1,2)*cpow(px,2)*cpow(csin(theta1r),2) + cpow(py,2)*cpow(csin(theta1r),2) + 
       2*a1*a6*ly*csin(theta1r) - 2*a6*d1*lambda1*lx*mu1*csin(theta1r) - 2*d1*d6*lambda1*mu1*mu6*mx*csin(theta1r) - 
       2*a6*d6*lambda1*lz*mu1*mu6*mx*csin(theta1r) + 2*a1*d6*mu6*my*csin(theta1r) - 2*a6*d6*lambda1*lx*mu1*mu6*mz*csin(theta1r) - 
       2*d1*d6*lambda1*lambda6*mu1*nx*csin(theta1r) - 2*a6*d6*lambda1*lambda6*lz*mu1*nx*csin(theta1r) + 2*a1*d6*lambda6*ny*csin(theta1r) - 
       2*a6*d6*lambda1*lambda6*lx*mu1*nz*csin(theta1r) + 2*d1*lambda1*mu1*px*csin(theta1r) + 2*a6*lambda1*lz*mu1*px*csin(theta1r) + 
       2*d6*lambda1*mu1*mu6*mz*px*csin(theta1r) + 2*d6*lambda1*lambda6*mu1*nz*px*csin(theta1r) - 2*a1*py*csin(theta1r) + 
       2*a6*lambda1*lx*mu1*pz*csin(theta1r) + 2*d6*lambda1*mu1*mu6*mx*pz*csin(theta1r) + 2*d6*lambda1*lambda6*mu1*nx*pz*csin(theta1r) - 
       2*lambda1*mu1*px*pz*csin(theta1r) + 2*a6*d6*ly*mu6*mx*ccos(theta1r)*csin(theta1r) + 2*a6*d6*lx*mu6*my*ccos(theta1r)*csin(theta1r) + 
       2*a6*d6*lambda6*ly*nx*ccos(theta1r)*csin(theta1r) + 2*a6*d6*lambda6*lx*ny*ccos(theta1r)*csin(theta1r) - 2*a6*ly*px*ccos(theta1r)*csin(theta1r) - 
       2*d6*mu6*my*px*ccos(theta1r)*csin(theta1r) - 2*d6*lambda6*ny*px*ccos(theta1r)*csin(theta1r) - 2*a6*lx*py*ccos(theta1r)*csin(theta1r) - 
       2*d6*mu6*mx*py*ccos(theta1r)*csin(theta1r) - 2*d6*lambda6*nx*py*ccos(theta1r)*csin(theta1r) + 2*px*py*ccos(theta1r)*csin(theta1r) - 
       2*lambda1*lx*lz*mu1*cpow(a6,2)*csin(theta1r) + 2*lx*ly*ccos(theta1r)*cpow(a6,2)*csin(theta1r) - 2*lambda1*lambda6*mu1*mu6*mz*nx*cpow(d6,2)*csin(theta1r) - 
       2*lambda1*lambda6*mu1*mu6*mx*nz*cpow(d6,2)*csin(theta1r) + 2*lambda6*mu6*my*nx*ccos(theta1r)*cpow(d6,2)*csin(theta1r) + 
       2*lambda6*mu6*mx*ny*ccos(theta1r)*cpow(d6,2)*csin(theta1r) - 2*a6*d6*ly*mu6*mx*ccos(theta1r)*cpow(lambda1,2)*csin(theta1r) - 
       2*a6*d6*lx*mu6*my*ccos(theta1r)*cpow(lambda1,2)*csin(theta1r) - 2*a6*d6*lambda6*ly*nx*ccos(theta1r)*cpow(lambda1,2)*csin(theta1r) - 
       2*a6*d6*lambda6*lx*ny*ccos(theta1r)*cpow(lambda1,2)*csin(theta1r) + 2*a6*ly*px*ccos(theta1r)*cpow(lambda1,2)*csin(theta1r) + 
       2*d6*mu6*my*px*ccos(theta1r)*cpow(lambda1,2)*csin(theta1r) + 2*d6*lambda6*ny*px*ccos(theta1r)*cpow(lambda1,2)*csin(theta1r) + 
       2*a6*lx*py*ccos(theta1r)*cpow(lambda1,2)*csin(theta1r) + 2*d6*mu6*mx*py*ccos(theta1r)*cpow(lambda1,2)*csin(theta1r) + 
       2*d6*lambda6*nx*py*ccos(theta1r)*cpow(lambda1,2)*csin(theta1r) - 2*px*py*ccos(theta1r)*cpow(lambda1,2)*csin(theta1r) - 
       2*lx*ly*ccos(theta1r)*cpow(a6,2)*cpow(lambda1,2)*csin(theta1r) - 2*lambda6*mu6*my*nx*ccos(theta1r)*cpow(d6,2)*cpow(lambda1,2)*csin(theta1r) - 
       2*lambda6*mu6*mx*ny*ccos(theta1r)*cpow(d6,2)*cpow(lambda1,2)*csin(theta1r) - 2*lambda1*mu1*nx*nz*cpow(d6,2)*cpow(lambda6,2)*csin(theta1r) + 
       2*nx*ny*ccos(theta1r)*cpow(d6,2)*cpow(lambda6,2)*csin(theta1r) - 2*nx*ny*ccos(theta1r)*cpow(d6,2)*cpow(lambda1,2)*cpow(lambda6,2)*csin(theta1r) - 
       2*lambda1*mu1*mx*mz*cpow(d6,2)*cpow(mu6,2)*csin(theta1r) + 2*mx*my*ccos(theta1r)*cpow(d6,2)*cpow(mu6,2)*csin(theta1r) - 
       2*mx*my*ccos(theta1r)*cpow(d6,2)*cpow(lambda1,2)*cpow(mu6,2)*csin(theta1r),-1)*
     (a2*d1*mu1 + a2*a6*lz*mu1 + a1*d3*mu2 + a1*d4*lambda3*mu2 + a1*d5*lambda3*lambda4*mu2 + a2*d6*mu1*mu6*mz + a2*d6*lambda6*mu1*nz - a2*mu1*pz + 
       a2*a6*lambda1*ly*ccos(theta1r) + a6*d3*lx*mu2*ccos(theta1r) + a6*d4*lambda3*lx*mu2*ccos(theta1r) + a6*d5*lambda3*lambda4*lx*mu2*ccos(theta1r) + 
       d3*d6*mu2*mu6*mx*ccos(theta1r) + d4*d6*lambda3*mu2*mu6*mx*ccos(theta1r) + d5*d6*lambda3*lambda4*mu2*mu6*mx*ccos(theta1r) + 
       a2*d6*lambda1*mu6*my*ccos(theta1r) + d3*d6*lambda6*mu2*nx*ccos(theta1r) + d4*d6*lambda3*lambda6*mu2*nx*ccos(theta1r) + 
       d5*d6*lambda3*lambda4*lambda6*mu2*nx*ccos(theta1r) + a2*d6*lambda1*lambda6*ny*ccos(theta1r) - d3*mu2*px*ccos(theta1r) - 
       d4*lambda3*mu2*px*ccos(theta1r) - d5*lambda3*lambda4*mu2*px*ccos(theta1r) - a2*lambda1*py*ccos(theta1r) + a3*d1*mu1*ccos(theta3r) + 
       a3*a6*lz*mu1*ccos(theta3r) + a1*d4*lambda2*mu3*ccos(theta3r) + a1*d5*lambda2*lambda4*mu3*ccos(theta3r) + a3*d6*mu1*mu6*mz*ccos(theta3r) + 
       a3*d6*lambda6*mu1*nz*ccos(theta3r) - a3*mu1*pz*ccos(theta3r) + a3*a6*lambda1*ly*ccos(theta1r)*ccos(theta3r) + 
       a6*d4*lambda2*lx*mu3*ccos(theta1r)*ccos(theta3r) + a6*d5*lambda2*lambda4*lx*mu3*ccos(theta1r)*ccos(theta3r) + 
       d4*d6*lambda2*mu3*mu6*mx*ccos(theta1r)*ccos(theta3r) + d5*d6*lambda2*lambda4*mu3*mu6*mx*ccos(theta1r)*ccos(theta3r) + 
       a3*d6*lambda1*mu6*my*ccos(theta1r)*ccos(theta3r) + d4*d6*lambda2*lambda6*mu3*nx*ccos(theta1r)*ccos(theta3r) + 
       d5*d6*lambda2*lambda4*lambda6*mu3*nx*ccos(theta1r)*ccos(theta3r) + a3*d6*lambda1*lambda6*ny*ccos(theta1r)*ccos(theta3r) - 
       d4*lambda2*mu3*px*ccos(theta1r)*ccos(theta3r) - d5*lambda2*lambda4*mu3*px*ccos(theta1r)*ccos(theta3r) - a3*lambda1*py*ccos(theta1r)*ccos(theta3r) - 
       a1*d5*mu2*mu3*mu4*ccos(theta4r) - a6*d5*lx*mu2*mu3*mu4*ccos(theta1r)*ccos(theta4r) - d5*d6*mu2*mu3*mu4*mu6*mx*ccos(theta1r)*ccos(theta4r) - 
       d5*d6*lambda6*mu2*mu3*mu4*nx*ccos(theta1r)*ccos(theta4r) + d5*mu2*mu3*mu4*px*ccos(theta1r)*ccos(theta4r) + a4*d1*mu1*ccos(theta3r)*ccos(theta4r) + 
       a4*a6*lz*mu1*ccos(theta3r)*ccos(theta4r) + a1*d5*lambda2*lambda3*mu4*ccos(theta3r)*ccos(theta4r) + a4*d6*mu1*mu6*mz*ccos(theta3r)*ccos(theta4r) + 
       a4*d6*lambda6*mu1*nz*ccos(theta3r)*ccos(theta4r) - a4*mu1*pz*ccos(theta3r)*ccos(theta4r) + a4*a6*lambda1*ly*ccos(theta1r)*ccos(theta3r)*ccos(theta4r) + 
       a6*d5*lambda2*lambda3*lx*mu4*ccos(theta1r)*ccos(theta3r)*ccos(theta4r) + d5*d6*lambda2*lambda3*mu4*mu6*mx*ccos(theta1r)*ccos(theta3r)*ccos(theta4r) + 
       a4*d6*lambda1*mu6*my*ccos(theta1r)*ccos(theta3r)*ccos(theta4r) + d5*d6*lambda2*lambda3*lambda6*mu4*nx*ccos(theta1r)*ccos(theta3r)*ccos(theta4r) + 
       a4*d6*lambda1*lambda6*ny*ccos(theta1r)*ccos(theta3r)*ccos(theta4r) - d5*lambda2*lambda3*mu4*px*ccos(theta1r)*ccos(theta3r)*ccos(theta4r) - 
       a4*lambda1*py*ccos(theta1r)*ccos(theta3r)*ccos(theta4r) + a5*d1*mu1*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) + 
       a5*a6*lz*mu1*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) + a5*d6*mu1*mu6*mz*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) + 
       a5*d6*lambda6*mu1*nz*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) - a5*mu1*pz*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) + 
       a5*a6*lambda1*ly*ccos(theta1r)*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) + a5*d6*lambda1*mu6*my*ccos(theta1r)*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) + 
       a5*d6*lambda1*lambda6*ny*ccos(theta1r)*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) - a5*lambda1*py*ccos(theta1r)*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) - 
       a2*a6*lambda1*lx*csin(theta1r) + a6*d3*ly*mu2*csin(theta1r) + a6*d4*lambda3*ly*mu2*csin(theta1r) + a6*d5*lambda3*lambda4*ly*mu2*csin(theta1r) - 
       a2*d6*lambda1*mu6*mx*csin(theta1r) + d3*d6*mu2*mu6*my*csin(theta1r) + d4*d6*lambda3*mu2*mu6*my*csin(theta1r) + 
       d5*d6*lambda3*lambda4*mu2*mu6*my*csin(theta1r) - a2*d6*lambda1*lambda6*nx*csin(theta1r) + d3*d6*lambda6*mu2*ny*csin(theta1r) + 
       d4*d6*lambda3*lambda6*mu2*ny*csin(theta1r) + d5*d6*lambda3*lambda4*lambda6*mu2*ny*csin(theta1r) + a2*lambda1*px*csin(theta1r) - 
       d3*mu2*py*csin(theta1r) - d4*lambda3*mu2*py*csin(theta1r) - d5*lambda3*lambda4*mu2*py*csin(theta1r) - a3*a6*lambda1*lx*ccos(theta3r)*csin(theta1r) + 
       a6*d4*lambda2*ly*mu3*ccos(theta3r)*csin(theta1r) + a6*d5*lambda2*lambda4*ly*mu3*ccos(theta3r)*csin(theta1r) - 
       a3*d6*lambda1*mu6*mx*ccos(theta3r)*csin(theta1r) + d4*d6*lambda2*mu3*mu6*my*ccos(theta3r)*csin(theta1r) + 
       d5*d6*lambda2*lambda4*mu3*mu6*my*ccos(theta3r)*csin(theta1r) - a3*d6*lambda1*lambda6*nx*ccos(theta3r)*csin(theta1r) + 
       d4*d6*lambda2*lambda6*mu3*ny*ccos(theta3r)*csin(theta1r) + d5*d6*lambda2*lambda4*lambda6*mu3*ny*ccos(theta3r)*csin(theta1r) + 
       a3*lambda1*px*ccos(theta3r)*csin(theta1r) - d4*lambda2*mu3*py*ccos(theta3r)*csin(theta1r) - d5*lambda2*lambda4*mu3*py*ccos(theta3r)*csin(theta1r) - 
       a6*d5*ly*mu2*mu3*mu4*ccos(theta4r)*csin(theta1r) - d5*d6*mu2*mu3*mu4*mu6*my*ccos(theta4r)*csin(theta1r) - 
       d5*d6*lambda6*mu2*mu3*mu4*ny*ccos(theta4r)*csin(theta1r) + d5*mu2*mu3*mu4*py*ccos(theta4r)*csin(theta1r) - 
       a4*a6*lambda1*lx*ccos(theta3r)*ccos(theta4r)*csin(theta1r) + a6*d5*lambda2*lambda3*ly*mu4*ccos(theta3r)*ccos(theta4r)*csin(theta1r) - 
       a4*d6*lambda1*mu6*mx*ccos(theta3r)*ccos(theta4r)*csin(theta1r) + d5*d6*lambda2*lambda3*mu4*mu6*my*ccos(theta3r)*ccos(theta4r)*csin(theta1r) - 
       a4*d6*lambda1*lambda6*nx*ccos(theta3r)*ccos(theta4r)*csin(theta1r) + d5*d6*lambda2*lambda3*lambda6*mu4*ny*ccos(theta3r)*ccos(theta4r)*csin(theta1r) + 
       a4*lambda1*px*ccos(theta3r)*ccos(theta4r)*csin(theta1r) - d5*lambda2*lambda3*mu4*py*ccos(theta3r)*ccos(theta4r)*csin(theta1r) - 
       a5*a6*lambda1*lx*ccos(theta3r)*ccos(theta4r)*ccos(theta5r)*csin(theta1r) - a5*d6*lambda1*mu6*mx*ccos(theta3r)*ccos(theta4r)*ccos(theta5r)*csin(theta1r) - 
       a5*d6*lambda1*lambda6*nx*ccos(theta3r)*ccos(theta4r)*ccos(theta5r)*csin(theta1r) + a5*lambda1*px*ccos(theta3r)*ccos(theta4r)*ccos(theta5r)*csin(theta1r) - 
       a1*a3*lambda2*csin(theta3r) + d1*d4*mu1*mu3*csin(theta3r) + d1*d5*lambda4*mu1*mu3*csin(theta3r) + a6*d4*lz*mu1*mu3*csin(theta3r) + 
       a6*d5*lambda4*lz*mu1*mu3*csin(theta3r) + d4*d6*mu1*mu3*mu6*mz*csin(theta3r) + d5*d6*lambda4*mu1*mu3*mu6*mz*csin(theta3r) + 
       d4*d6*lambda6*mu1*mu3*nz*csin(theta3r) + d5*d6*lambda4*lambda6*mu1*mu3*nz*csin(theta3r) - d4*mu1*mu3*pz*csin(theta3r) - 
       d5*lambda4*mu1*mu3*pz*csin(theta3r) - a3*a6*lambda2*lx*ccos(theta1r)*csin(theta3r) + a6*d4*lambda1*ly*mu3*ccos(theta1r)*csin(theta3r) + 
       a6*d5*lambda1*lambda4*ly*mu3*ccos(theta1r)*csin(theta3r) - a3*d6*lambda2*mu6*mx*ccos(theta1r)*csin(theta3r) + 
       d4*d6*lambda1*mu3*mu6*my*ccos(theta1r)*csin(theta3r) + d5*d6*lambda1*lambda4*mu3*mu6*my*ccos(theta1r)*csin(theta3r) - 
       a3*d6*lambda2*lambda6*nx*ccos(theta1r)*csin(theta3r) + d4*d6*lambda1*lambda6*mu3*ny*ccos(theta1r)*csin(theta3r) + 
       d5*d6*lambda1*lambda4*lambda6*mu3*ny*ccos(theta1r)*csin(theta3r) + a3*lambda2*px*ccos(theta1r)*csin(theta3r) - 
       d4*lambda1*mu3*py*ccos(theta1r)*csin(theta3r) - d5*lambda1*lambda4*mu3*py*ccos(theta1r)*csin(theta3r) - a1*a4*lambda2*ccos(theta4r)*csin(theta3r) + 
       d1*d5*lambda3*mu1*mu4*ccos(theta4r)*csin(theta3r) + a6*d5*lambda3*lz*mu1*mu4*ccos(theta4r)*csin(theta3r) + 
       d5*d6*lambda3*mu1*mu4*mu6*mz*ccos(theta4r)*csin(theta3r) + d5*d6*lambda3*lambda6*mu1*mu4*nz*ccos(theta4r)*csin(theta3r) - 
       d5*lambda3*mu1*mu4*pz*ccos(theta4r)*csin(theta3r) - a4*a6*lambda2*lx*ccos(theta1r)*ccos(theta4r)*csin(theta3r) + 
       a6*d5*lambda1*lambda3*ly*mu4*ccos(theta1r)*ccos(theta4r)*csin(theta3r) - a4*d6*lambda2*mu6*mx*ccos(theta1r)*ccos(theta4r)*csin(theta3r) + 
       d5*d6*lambda1*lambda3*mu4*mu6*my*ccos(theta1r)*ccos(theta4r)*csin(theta3r) - a4*d6*lambda2*lambda6*nx*ccos(theta1r)*ccos(theta4r)*csin(theta3r) + 
       d5*d6*lambda1*lambda3*lambda6*mu4*ny*ccos(theta1r)*ccos(theta4r)*csin(theta3r) + a4*lambda2*px*ccos(theta1r)*ccos(theta4r)*csin(theta3r) - 
       d5*lambda1*lambda3*mu4*py*ccos(theta1r)*ccos(theta4r)*csin(theta3r) - a1*a5*lambda2*ccos(theta4r)*ccos(theta5r)*csin(theta3r) - 
       a5*a6*lambda2*lx*ccos(theta1r)*ccos(theta4r)*ccos(theta5r)*csin(theta3r) - a5*d6*lambda2*mu6*mx*ccos(theta1r)*ccos(theta4r)*ccos(theta5r)*csin(theta3r) - 
       a5*d6*lambda2*lambda6*nx*ccos(theta1r)*ccos(theta4r)*ccos(theta5r)*csin(theta3r) + a5*lambda2*px*ccos(theta1r)*ccos(theta4r)*ccos(theta5r)*csin(theta3r) - 
       a3*a6*lambda2*ly*csin(theta1r)*csin(theta3r) - a6*d4*lambda1*lx*mu3*csin(theta1r)*csin(theta3r) - 
       a6*d5*lambda1*lambda4*lx*mu3*csin(theta1r)*csin(theta3r) - d4*d6*lambda1*mu3*mu6*mx*csin(theta1r)*csin(theta3r) - 
       d5*d6*lambda1*lambda4*mu3*mu6*mx*csin(theta1r)*csin(theta3r) - a3*d6*lambda2*mu6*my*csin(theta1r)*csin(theta3r) - 
       d4*d6*lambda1*lambda6*mu3*nx*csin(theta1r)*csin(theta3r) - d5*d6*lambda1*lambda4*lambda6*mu3*nx*csin(theta1r)*csin(theta3r) - 
       a3*d6*lambda2*lambda6*ny*csin(theta1r)*csin(theta3r) + d4*lambda1*mu3*px*csin(theta1r)*csin(theta3r) + 
       d5*lambda1*lambda4*mu3*px*csin(theta1r)*csin(theta3r) + a3*lambda2*py*csin(theta1r)*csin(theta3r) - 
       a4*a6*lambda2*ly*ccos(theta4r)*csin(theta1r)*csin(theta3r) - a6*d5*lambda1*lambda3*lx*mu4*ccos(theta4r)*csin(theta1r)*csin(theta3r) - 
       d5*d6*lambda1*lambda3*mu4*mu6*mx*ccos(theta4r)*csin(theta1r)*csin(theta3r) - a4*d6*lambda2*mu6*my*ccos(theta4r)*csin(theta1r)*csin(theta3r) - 
       d5*d6*lambda1*lambda3*lambda6*mu4*nx*ccos(theta4r)*csin(theta1r)*csin(theta3r) - a4*d6*lambda2*lambda6*ny*ccos(theta4r)*csin(theta1r)*csin(theta3r) + 
       d5*lambda1*lambda3*mu4*px*ccos(theta4r)*csin(theta1r)*csin(theta3r) + a4*lambda2*py*ccos(theta4r)*csin(theta1r)*csin(theta3r) - 
       a5*a6*lambda2*ly*ccos(theta4r)*ccos(theta5r)*csin(theta1r)*csin(theta3r) - a5*d6*lambda2*mu6*my*ccos(theta4r)*ccos(theta5r)*csin(theta1r)*csin(theta3r) - 
       a5*d6*lambda2*lambda6*ny*ccos(theta4r)*ccos(theta5r)*csin(theta1r)*csin(theta3r) + a5*lambda2*py*ccos(theta4r)*ccos(theta5r)*csin(theta1r)*csin(theta3r) + 
       a1*a4*mu2*mu3*csin(theta4r) + a4*a6*lx*mu2*mu3*ccos(theta1r)*csin(theta4r) + a4*d6*mu2*mu3*mu6*mx*ccos(theta1r)*csin(theta4r) + 
       a4*d6*lambda6*mu2*mu3*nx*ccos(theta1r)*csin(theta4r) - a4*mu2*mu3*px*ccos(theta1r)*csin(theta4r) - a1*a4*lambda2*lambda3*ccos(theta3r)*csin(theta4r) + 
       d1*d5*mu1*mu4*ccos(theta3r)*csin(theta4r) + a6*d5*lz*mu1*mu4*ccos(theta3r)*csin(theta4r) + d5*d6*mu1*mu4*mu6*mz*ccos(theta3r)*csin(theta4r) + 
       d5*d6*lambda6*mu1*mu4*nz*ccos(theta3r)*csin(theta4r) - d5*mu1*mu4*pz*ccos(theta3r)*csin(theta4r) - 
       a4*a6*lambda2*lambda3*lx*ccos(theta1r)*ccos(theta3r)*csin(theta4r) + a6*d5*lambda1*ly*mu4*ccos(theta1r)*ccos(theta3r)*csin(theta4r) - 
       a4*d6*lambda2*lambda3*mu6*mx*ccos(theta1r)*ccos(theta3r)*csin(theta4r) + d5*d6*lambda1*mu4*mu6*my*ccos(theta1r)*ccos(theta3r)*csin(theta4r) - 
       a4*d6*lambda2*lambda3*lambda6*nx*ccos(theta1r)*ccos(theta3r)*csin(theta4r) + d5*d6*lambda1*lambda6*mu4*ny*ccos(theta1r)*ccos(theta3r)*csin(theta4r) + 
       a4*lambda2*lambda3*px*ccos(theta1r)*ccos(theta3r)*csin(theta4r) - d5*lambda1*mu4*py*ccos(theta1r)*ccos(theta3r)*csin(theta4r) + 
       a1*a5*mu2*mu3*ccos(theta5r)*csin(theta4r) + a5*a6*lx*mu2*mu3*ccos(theta1r)*ccos(theta5r)*csin(theta4r) + 
       a5*d6*mu2*mu3*mu6*mx*ccos(theta1r)*ccos(theta5r)*csin(theta4r) + a5*d6*lambda6*mu2*mu3*nx*ccos(theta1r)*ccos(theta5r)*csin(theta4r) - 
       a5*mu2*mu3*px*ccos(theta1r)*ccos(theta5r)*csin(theta4r) - a1*a5*lambda2*lambda3*ccos(theta3r)*ccos(theta5r)*csin(theta4r) - 
       a5*a6*lambda2*lambda3*lx*ccos(theta1r)*ccos(theta3r)*ccos(theta5r)*csin(theta4r) - 
       a5*d6*lambda2*lambda3*mu6*mx*ccos(theta1r)*ccos(theta3r)*ccos(theta5r)*csin(theta4r) - 
       a5*d6*lambda2*lambda3*lambda6*nx*ccos(theta1r)*ccos(theta3r)*ccos(theta5r)*csin(theta4r) + 
       a5*lambda2*lambda3*px*ccos(theta1r)*ccos(theta3r)*ccos(theta5r)*csin(theta4r) + a4*a6*ly*mu2*mu3*csin(theta1r)*csin(theta4r) + 
       a4*d6*mu2*mu3*mu6*my*csin(theta1r)*csin(theta4r) + a4*d6*lambda6*mu2*mu3*ny*csin(theta1r)*csin(theta4r) - a4*mu2*mu3*py*csin(theta1r)*csin(theta4r) - 
       a4*a6*lambda2*lambda3*ly*ccos(theta3r)*csin(theta1r)*csin(theta4r) - a6*d5*lambda1*lx*mu4*ccos(theta3r)*csin(theta1r)*csin(theta4r) - 
       d5*d6*lambda1*mu4*mu6*mx*ccos(theta3r)*csin(theta1r)*csin(theta4r) - a4*d6*lambda2*lambda3*mu6*my*ccos(theta3r)*csin(theta1r)*csin(theta4r) - 
       d5*d6*lambda1*lambda6*mu4*nx*ccos(theta3r)*csin(theta1r)*csin(theta4r) - a4*d6*lambda2*lambda3*lambda6*ny*ccos(theta3r)*csin(theta1r)*csin(theta4r) + 
       d5*lambda1*mu4*px*ccos(theta3r)*csin(theta1r)*csin(theta4r) + a4*lambda2*lambda3*py*ccos(theta3r)*csin(theta1r)*csin(theta4r) + 
       a5*a6*ly*mu2*mu3*ccos(theta5r)*csin(theta1r)*csin(theta4r) + a5*d6*mu2*mu3*mu6*my*ccos(theta5r)*csin(theta1r)*csin(theta4r) + 
       a5*d6*lambda6*mu2*mu3*ny*ccos(theta5r)*csin(theta1r)*csin(theta4r) - a5*mu2*mu3*py*ccos(theta5r)*csin(theta1r)*csin(theta4r) - 
       a5*a6*lambda2*lambda3*ly*ccos(theta3r)*ccos(theta5r)*csin(theta1r)*csin(theta4r) - 
       a5*d6*lambda2*lambda3*mu6*my*ccos(theta3r)*ccos(theta5r)*csin(theta1r)*csin(theta4r) - 
       a5*d6*lambda2*lambda3*lambda6*ny*ccos(theta3r)*ccos(theta5r)*csin(theta1r)*csin(theta4r) + 
       a5*lambda2*lambda3*py*ccos(theta3r)*ccos(theta5r)*csin(theta1r)*csin(theta4r) - a4*d1*lambda3*mu1*csin(theta3r)*csin(theta4r) - 
       a4*a6*lambda3*lz*mu1*csin(theta3r)*csin(theta4r) - a1*d5*lambda2*mu4*csin(theta3r)*csin(theta4r) - a4*d6*lambda3*mu1*mu6*mz*csin(theta3r)*csin(theta4r) - 
       a4*d6*lambda3*lambda6*mu1*nz*csin(theta3r)*csin(theta4r) + a4*lambda3*mu1*pz*csin(theta3r)*csin(theta4r) - 
       a4*a6*lambda1*lambda3*ly*ccos(theta1r)*csin(theta3r)*csin(theta4r) - a6*d5*lambda2*lx*mu4*ccos(theta1r)*csin(theta3r)*csin(theta4r) - 
       d5*d6*lambda2*mu4*mu6*mx*ccos(theta1r)*csin(theta3r)*csin(theta4r) - a4*d6*lambda1*lambda3*mu6*my*ccos(theta1r)*csin(theta3r)*csin(theta4r) - 
       d5*d6*lambda2*lambda6*mu4*nx*ccos(theta1r)*csin(theta3r)*csin(theta4r) - a4*d6*lambda1*lambda3*lambda6*ny*ccos(theta1r)*csin(theta3r)*csin(theta4r) + 
       d5*lambda2*mu4*px*ccos(theta1r)*csin(theta3r)*csin(theta4r) + a4*lambda1*lambda3*py*ccos(theta1r)*csin(theta3r)*csin(theta4r) - 
       a5*d1*lambda3*mu1*ccos(theta5r)*csin(theta3r)*csin(theta4r) - a5*a6*lambda3*lz*mu1*ccos(theta5r)*csin(theta3r)*csin(theta4r) - 
       a5*d6*lambda3*mu1*mu6*mz*ccos(theta5r)*csin(theta3r)*csin(theta4r) - a5*d6*lambda3*lambda6*mu1*nz*ccos(theta5r)*csin(theta3r)*csin(theta4r) + 
       a5*lambda3*mu1*pz*ccos(theta5r)*csin(theta3r)*csin(theta4r) - a5*a6*lambda1*lambda3*ly*ccos(theta1r)*ccos(theta5r)*csin(theta3r)*csin(theta4r) - 
       a5*d6*lambda1*lambda3*mu6*my*ccos(theta1r)*ccos(theta5r)*csin(theta3r)*csin(theta4r) - 
       a5*d6*lambda1*lambda3*lambda6*ny*ccos(theta1r)*ccos(theta5r)*csin(theta3r)*csin(theta4r) + 
       a5*lambda1*lambda3*py*ccos(theta1r)*ccos(theta5r)*csin(theta3r)*csin(theta4r) + a4*a6*lambda1*lambda3*lx*csin(theta1r)*csin(theta3r)*csin(theta4r) - 
       a6*d5*lambda2*ly*mu4*csin(theta1r)*csin(theta3r)*csin(theta4r) + a4*d6*lambda1*lambda3*mu6*mx*csin(theta1r)*csin(theta3r)*csin(theta4r) - 
       d5*d6*lambda2*mu4*mu6*my*csin(theta1r)*csin(theta3r)*csin(theta4r) + a4*d6*lambda1*lambda3*lambda6*nx*csin(theta1r)*csin(theta3r)*csin(theta4r) - 
       d5*d6*lambda2*lambda6*mu4*ny*csin(theta1r)*csin(theta3r)*csin(theta4r) - a4*lambda1*lambda3*px*csin(theta1r)*csin(theta3r)*csin(theta4r) + 
       d5*lambda2*mu4*py*csin(theta1r)*csin(theta3r)*csin(theta4r) + a5*a6*lambda1*lambda3*lx*ccos(theta5r)*csin(theta1r)*csin(theta3r)*csin(theta4r) + 
       a5*d6*lambda1*lambda3*mu6*mx*ccos(theta5r)*csin(theta1r)*csin(theta3r)*csin(theta4r) + 
       a5*d6*lambda1*lambda3*lambda6*nx*ccos(theta5r)*csin(theta1r)*csin(theta3r)*csin(theta4r) - 
       a5*lambda1*lambda3*px*ccos(theta5r)*csin(theta1r)*csin(theta3r)*csin(theta4r) + a1*a5*lambda3*mu2*mu4*csin(theta5r) + 
       a5*a6*lambda3*lx*mu2*mu4*ccos(theta1r)*csin(theta5r) + a5*d6*lambda3*mu2*mu4*mu6*mx*ccos(theta1r)*csin(theta5r) + 
       a5*d6*lambda3*lambda6*mu2*mu4*nx*ccos(theta1r)*csin(theta5r) - a5*lambda3*mu2*mu4*px*ccos(theta1r)*csin(theta5r) + 
       a1*a5*lambda2*mu3*mu4*ccos(theta3r)*csin(theta5r) + a5*a6*lambda2*lx*mu3*mu4*ccos(theta1r)*ccos(theta3r)*csin(theta5r) + 
       a5*d6*lambda2*mu3*mu4*mu6*mx*ccos(theta1r)*ccos(theta3r)*csin(theta5r) + a5*d6*lambda2*lambda6*mu3*mu4*nx*ccos(theta1r)*ccos(theta3r)*csin(theta5r) - 
       a5*lambda2*mu3*mu4*px*ccos(theta1r)*ccos(theta3r)*csin(theta5r) + a1*a5*lambda4*mu2*mu3*ccos(theta4r)*csin(theta5r) + 
       a5*a6*lambda4*lx*mu2*mu3*ccos(theta1r)*ccos(theta4r)*csin(theta5r) + a5*d6*lambda4*mu2*mu3*mu6*mx*ccos(theta1r)*ccos(theta4r)*csin(theta5r) + 
       a5*d6*lambda4*lambda6*mu2*mu3*nx*ccos(theta1r)*ccos(theta4r)*csin(theta5r) - a5*lambda4*mu2*mu3*px*ccos(theta1r)*ccos(theta4r)*csin(theta5r) - 
       a1*a5*lambda2*lambda3*lambda4*ccos(theta3r)*ccos(theta4r)*csin(theta5r) - 
       a5*a6*lambda2*lambda3*lambda4*lx*ccos(theta1r)*ccos(theta3r)*ccos(theta4r)*csin(theta5r) - 
       a5*d6*lambda2*lambda3*lambda4*mu6*mx*ccos(theta1r)*ccos(theta3r)*ccos(theta4r)*csin(theta5r) - 
       a5*d6*lambda2*lambda3*lambda4*lambda6*nx*ccos(theta1r)*ccos(theta3r)*ccos(theta4r)*csin(theta5r) + 
       a5*lambda2*lambda3*lambda4*px*ccos(theta1r)*ccos(theta3r)*ccos(theta4r)*csin(theta5r) + a5*a6*lambda3*ly*mu2*mu4*csin(theta1r)*csin(theta5r) + 
       a5*d6*lambda3*mu2*mu4*mu6*my*csin(theta1r)*csin(theta5r) + a5*d6*lambda3*lambda6*mu2*mu4*ny*csin(theta1r)*csin(theta5r) - 
       a5*lambda3*mu2*mu4*py*csin(theta1r)*csin(theta5r) + a5*a6*lambda2*ly*mu3*mu4*ccos(theta3r)*csin(theta1r)*csin(theta5r) + 
       a5*d6*lambda2*mu3*mu4*mu6*my*ccos(theta3r)*csin(theta1r)*csin(theta5r) + a5*d6*lambda2*lambda6*mu3*mu4*ny*ccos(theta3r)*csin(theta1r)*csin(theta5r) - 
       a5*lambda2*mu3*mu4*py*ccos(theta3r)*csin(theta1r)*csin(theta5r) + a5*a6*lambda4*ly*mu2*mu3*ccos(theta4r)*csin(theta1r)*csin(theta5r) + 
       a5*d6*lambda4*mu2*mu3*mu6*my*ccos(theta4r)*csin(theta1r)*csin(theta5r) + a5*d6*lambda4*lambda6*mu2*mu3*ny*ccos(theta4r)*csin(theta1r)*csin(theta5r) - 
       a5*lambda4*mu2*mu3*py*ccos(theta4r)*csin(theta1r)*csin(theta5r) - 
       a5*a6*lambda2*lambda3*lambda4*ly*ccos(theta3r)*ccos(theta4r)*csin(theta1r)*csin(theta5r) - 
       a5*d6*lambda2*lambda3*lambda4*mu6*my*ccos(theta3r)*ccos(theta4r)*csin(theta1r)*csin(theta5r) - 
       a5*d6*lambda2*lambda3*lambda4*lambda6*ny*ccos(theta3r)*ccos(theta4r)*csin(theta1r)*csin(theta5r) + 
       a5*lambda2*lambda3*lambda4*py*ccos(theta3r)*ccos(theta4r)*csin(theta1r)*csin(theta5r) + a5*d1*mu1*mu3*mu4*csin(theta3r)*csin(theta5r) + 
       a5*a6*lz*mu1*mu3*mu4*csin(theta3r)*csin(theta5r) + a5*d6*mu1*mu3*mu4*mu6*mz*csin(theta3r)*csin(theta5r) + 
       a5*d6*lambda6*mu1*mu3*mu4*nz*csin(theta3r)*csin(theta5r) - a5*mu1*mu3*mu4*pz*csin(theta3r)*csin(theta5r) + 
       a5*a6*lambda1*ly*mu3*mu4*ccos(theta1r)*csin(theta3r)*csin(theta5r) + a5*d6*lambda1*mu3*mu4*mu6*my*ccos(theta1r)*csin(theta3r)*csin(theta5r) + 
       a5*d6*lambda1*lambda6*mu3*mu4*ny*ccos(theta1r)*csin(theta3r)*csin(theta5r) - a5*lambda1*mu3*mu4*py*ccos(theta1r)*csin(theta3r)*csin(theta5r) - 
       a5*d1*lambda3*lambda4*mu1*ccos(theta4r)*csin(theta3r)*csin(theta5r) - a5*a6*lambda3*lambda4*lz*mu1*ccos(theta4r)*csin(theta3r)*csin(theta5r) - 
       a5*d6*lambda3*lambda4*mu1*mu6*mz*ccos(theta4r)*csin(theta3r)*csin(theta5r) - 
       a5*d6*lambda3*lambda4*lambda6*mu1*nz*ccos(theta4r)*csin(theta3r)*csin(theta5r) + a5*lambda3*lambda4*mu1*pz*ccos(theta4r)*csin(theta3r)*csin(theta5r) - 
       a5*a6*lambda1*lambda3*lambda4*ly*ccos(theta1r)*ccos(theta4r)*csin(theta3r)*csin(theta5r) - 
       a5*d6*lambda1*lambda3*lambda4*mu6*my*ccos(theta1r)*ccos(theta4r)*csin(theta3r)*csin(theta5r) - 
       a5*d6*lambda1*lambda3*lambda4*lambda6*ny*ccos(theta1r)*ccos(theta4r)*csin(theta3r)*csin(theta5r) + 
       a5*lambda1*lambda3*lambda4*py*ccos(theta1r)*ccos(theta4r)*csin(theta3r)*csin(theta5r) - 
       a5*a6*lambda1*lx*mu3*mu4*csin(theta1r)*csin(theta3r)*csin(theta5r) - a5*d6*lambda1*mu3*mu4*mu6*mx*csin(theta1r)*csin(theta3r)*csin(theta5r) - 
       a5*d6*lambda1*lambda6*mu3*mu4*nx*csin(theta1r)*csin(theta3r)*csin(theta5r) + a5*lambda1*mu3*mu4*px*csin(theta1r)*csin(theta3r)*csin(theta5r) + 
       a5*a6*lambda1*lambda3*lambda4*lx*ccos(theta4r)*csin(theta1r)*csin(theta3r)*csin(theta5r) + 
       a5*d6*lambda1*lambda3*lambda4*mu6*mx*ccos(theta4r)*csin(theta1r)*csin(theta3r)*csin(theta5r) + 
       a5*d6*lambda1*lambda3*lambda4*lambda6*nx*ccos(theta4r)*csin(theta1r)*csin(theta3r)*csin(theta5r) - 
       a5*lambda1*lambda3*lambda4*px*ccos(theta4r)*csin(theta1r)*csin(theta3r)*csin(theta5r) - a5*d1*lambda4*mu1*ccos(theta3r)*csin(theta4r)*csin(theta5r) - 
       a5*a6*lambda4*lz*mu1*ccos(theta3r)*csin(theta4r)*csin(theta5r) - a5*d6*lambda4*mu1*mu6*mz*ccos(theta3r)*csin(theta4r)*csin(theta5r) - 
       a5*d6*lambda4*lambda6*mu1*nz*ccos(theta3r)*csin(theta4r)*csin(theta5r) + a5*lambda4*mu1*pz*ccos(theta3r)*csin(theta4r)*csin(theta5r) - 
       a5*a6*lambda1*lambda4*ly*ccos(theta1r)*ccos(theta3r)*csin(theta4r)*csin(theta5r) - 
       a5*d6*lambda1*lambda4*mu6*my*ccos(theta1r)*ccos(theta3r)*csin(theta4r)*csin(theta5r) - 
       a5*d6*lambda1*lambda4*lambda6*ny*ccos(theta1r)*ccos(theta3r)*csin(theta4r)*csin(theta5r) + 
       a5*lambda1*lambda4*py*ccos(theta1r)*ccos(theta3r)*csin(theta4r)*csin(theta5r) + 
       a5*a6*lambda1*lambda4*lx*ccos(theta3r)*csin(theta1r)*csin(theta4r)*csin(theta5r) + 
       a5*d6*lambda1*lambda4*mu6*mx*ccos(theta3r)*csin(theta1r)*csin(theta4r)*csin(theta5r) + 
       a5*d6*lambda1*lambda4*lambda6*nx*ccos(theta3r)*csin(theta1r)*csin(theta4r)*csin(theta5r) - 
       a5*lambda1*lambda4*px*ccos(theta3r)*csin(theta1r)*csin(theta4r)*csin(theta5r) + a1*a5*lambda2*lambda4*csin(theta3r)*csin(theta4r)*csin(theta5r) + 
       a5*a6*lambda2*lambda4*lx*ccos(theta1r)*csin(theta3r)*csin(theta4r)*csin(theta5r) + 
       a5*d6*lambda2*lambda4*mu6*mx*ccos(theta1r)*csin(theta3r)*csin(theta4r)*csin(theta5r) + 
       a5*d6*lambda2*lambda4*lambda6*nx*ccos(theta1r)*csin(theta3r)*csin(theta4r)*csin(theta5r) - 
       a5*lambda2*lambda4*px*ccos(theta1r)*csin(theta3r)*csin(theta4r)*csin(theta5r) + 
       a5*a6*lambda2*lambda4*ly*csin(theta1r)*csin(theta3r)*csin(theta4r)*csin(theta5r) + 
       a5*d6*lambda2*lambda4*mu6*my*csin(theta1r)*csin(theta3r)*csin(theta4r)*csin(theta5r) + 
       a5*d6*lambda2*lambda4*lambda6*ny*csin(theta1r)*csin(theta3r)*csin(theta4r)*csin(theta5r) - 
       a5*lambda2*lambda4*py*csin(theta1r)*csin(theta3r)*csin(theta4r)*csin(theta5r)));
	  
    sctheta2.sintheta = sintheta2;
    sctheta2.costheta = costheta2;

    return(sctheta2);

}
	   
