/** 
 

 @file theta1.c
 @brief This file is used to calculate the value of theta1 with the values of theta3, theta4, theta5, and the D-H parameters  as input.
       .
 
 */


# include <stdio.h>
# include <math.h>
# include <complex.h>
# include "cx.h"



struct Theta theta1(double* coeff, double _Complex *thlist)
{
    
    struct Theta sctheta1;

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

  	double _Complex theta3r = thlist[0];
    double _Complex theta4r = thlist[1];
    double _Complex theta5r = thlist[2];




    double _Complex sintheta1 = -(cpow(-(a6*ly*mu6*mx*cpow(mu1,2)) + a6*lx*mu6*my*cpow(mu1,2) - a6*lambda6*ly*nx*cpow(mu1,2) + a6*lambda6*lx*ny*cpow(mu1,2) - mu6*my*px*cpow(mu1,2) - lambda6*ny*px*cpow(mu1,2) + mu6*mx*py*cpow(mu1,2) + lambda6*nx*py*cpow(mu1,2),-1)*((a6*ly*mu1 + d6*mu1*mu6*my + d6*lambda6*mu1*ny - mu1*py)*(lambda2*lambda3*lambda4*lambda5 - lambda1*(mu6*mz + lambda6*nz) - lambda4*lambda5*mu2*mu3*ccos(theta3r) + (-(lambda2*lambda5*mu3*mu4) - lambda3*lambda5*mu2*mu4*ccos(theta3r))*ccos(theta4r) + (-(lambda2*lambda3*mu4*mu5) + mu2*mu3*mu4*mu5*ccos(theta3r))*ccos(theta5r) + (-(lambda2*lambda4*mu3*mu5) - lambda3*lambda4*mu2*mu5*ccos(theta3r))*ccos(theta4r)*ccos(theta5r) + lambda5*mu2*mu4*csin(theta3r)*csin(theta4r) + lambda4*mu2*mu5*ccos(theta5r)*csin(theta3r)*csin(theta4r) + mu2*mu5*ccos(theta4r)*csin(theta3r)*csin(theta5r) + (lambda2*mu3*mu5 + lambda3*mu2*mu5*ccos(theta3r))*csin(theta4r)*csin(theta5r)) + (mu1*mu6*my + lambda6*mu1*ny)*(d2 + lambda2*(d3 + lambda3*(d4 + d5*lambda4)) + lambda1*(d1 + a6*lz + d6*mu6*mz + d6*lambda6*nz - pz) - 
                        (d4 + d5*lambda4)*mu2*mu3*ccos(theta3r) + a3*mu2*csin(theta3r) + a5*mu2*ccos(theta4r)*ccos(theta5r)*csin(theta3r) + ccos(theta4r)*(-(d5*lambda2*mu3*mu4) - d5*lambda3*mu2*mu4*ccos(theta3r) + a4*mu2*csin(theta3r)) + (a5*lambda2*mu3 + a5*lambda3*mu2*ccos(theta3r))*ccos(theta5r)*csin(theta4r) + (a4*lambda2*mu3 + a4*lambda3*mu2*ccos(theta3r) + d5*mu2*mu4*csin(theta3r))*csin(theta4r) + (a5*lambda2*lambda3*mu4 - a5*mu2*mu3*mu4*ccos(theta3r))*csin(theta5r) + (a5*lambda2*lambda4*mu3 + a5*lambda3*lambda4*mu2*ccos(theta3r))*ccos(theta4r)*csin(theta5r) - a5*lambda4*mu2*csin(theta3r)*csin(theta4r)*csin(theta5r))));

    double _Complex costheta1 = -(cpow(mu1,-1)*cpow(-(a6*ly*mu6*mx) + a6*lx*mu6*my - a6*lambda6*ly*nx + a6*lambda6*lx*ny - mu6*my*px - lambda6*ny*px + mu6*mx*py + lambda6*nx*py,-1)*
     (a6*lambda2*lambda3*lambda4*lambda5*lx + d2*mu6*mx + d1*lambda1*mu6*mx + d3*lambda2*mu6*mx + d4*lambda2*lambda3*mu6*mx + 
       d5*lambda2*lambda3*lambda4*mu6*mx + d6*lambda2*lambda3*lambda4*lambda5*mu6*mx + a6*lambda1*lz*mu6*mx - a6*lambda1*lx*mu6*mz + d2*lambda6*nx + 
       d1*lambda1*lambda6*nx + d3*lambda2*lambda6*nx + d4*lambda2*lambda3*lambda6*nx + d5*lambda2*lambda3*lambda4*lambda6*nx + 
       d6*lambda2*lambda3*lambda4*lambda5*lambda6*nx + a6*lambda1*lambda6*lz*nx - a6*lambda1*lambda6*lx*nz - lambda2*lambda3*lambda4*lambda5*px + 
       lambda1*mu6*mz*px + lambda1*lambda6*nz*px - lambda1*mu6*mx*pz - lambda1*lambda6*nx*pz - a6*lambda4*lambda5*lx*mu2*mu3*ccos(theta3r) - 
       d4*mu2*mu3*mu6*mx*ccos(theta3r) - d5*lambda4*mu2*mu3*mu6*mx*ccos(theta3r) - d6*lambda4*lambda5*mu2*mu3*mu6*mx*ccos(theta3r) - 
       d4*lambda6*mu2*mu3*nx*ccos(theta3r) - d5*lambda4*lambda6*mu2*mu3*nx*ccos(theta3r) - d6*lambda4*lambda5*lambda6*mu2*mu3*nx*ccos(theta3r) + 
       lambda4*lambda5*mu2*mu3*px*ccos(theta3r) - a6*lambda2*lambda5*lx*mu3*mu4*ccos(theta4r) - d5*lambda2*mu3*mu4*mu6*mx*ccos(theta4r) - 
       d6*lambda2*lambda5*mu3*mu4*mu6*mx*ccos(theta4r) - d5*lambda2*lambda6*mu3*mu4*nx*ccos(theta4r) - d6*lambda2*lambda5*lambda6*mu3*mu4*nx*ccos(theta4r) + 
       lambda2*lambda5*mu3*mu4*px*ccos(theta4r) - a6*lambda3*lambda5*lx*mu2*mu4*ccos(theta3r)*ccos(theta4r) - 
       d5*lambda3*mu2*mu4*mu6*mx*ccos(theta3r)*ccos(theta4r) - d6*lambda3*lambda5*mu2*mu4*mu6*mx*ccos(theta3r)*ccos(theta4r) - 
       d5*lambda3*lambda6*mu2*mu4*nx*ccos(theta3r)*ccos(theta4r) - d6*lambda3*lambda5*lambda6*mu2*mu4*nx*ccos(theta3r)*ccos(theta4r) + 
       lambda3*lambda5*mu2*mu4*px*ccos(theta3r)*ccos(theta4r) - a6*lambda2*lambda3*lx*mu4*mu5*ccos(theta5r) - d6*lambda2*lambda3*mu4*mu5*mu6*mx*ccos(theta5r) - 
       d6*lambda2*lambda3*lambda6*mu4*mu5*nx*ccos(theta5r) + lambda2*lambda3*mu4*mu5*px*ccos(theta5r) + a6*lx*mu2*mu3*mu4*mu5*ccos(theta3r)*ccos(theta5r) + 
       d6*mu2*mu3*mu4*mu5*mu6*mx*ccos(theta3r)*ccos(theta5r) + d6*lambda6*mu2*mu3*mu4*mu5*nx*ccos(theta3r)*ccos(theta5r) - 
       mu2*mu3*mu4*mu5*px*ccos(theta3r)*ccos(theta5r) - a6*lambda2*lambda4*lx*mu3*mu5*ccos(theta4r)*ccos(theta5r) - 
       d6*lambda2*lambda4*mu3*mu5*mu6*mx*ccos(theta4r)*ccos(theta5r) - d6*lambda2*lambda4*lambda6*mu3*mu5*nx*ccos(theta4r)*ccos(theta5r) + 
       lambda2*lambda4*mu3*mu5*px*ccos(theta4r)*ccos(theta5r) - a6*lambda3*lambda4*lx*mu2*mu5*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) - 
       d6*lambda3*lambda4*mu2*mu5*mu6*mx*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) - 
       d6*lambda3*lambda4*lambda6*mu2*mu5*nx*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) + lambda3*lambda4*mu2*mu5*px*ccos(theta3r)*ccos(theta4r)*ccos(theta5r) + 
       a3*mu2*mu6*mx*csin(theta3r) + a3*lambda6*mu2*nx*csin(theta3r) + a4*mu2*mu6*mx*ccos(theta4r)*csin(theta3r) + a4*lambda6*mu2*nx*ccos(theta4r)*csin(theta3r) + 
       a5*mu2*mu6*mx*ccos(theta4r)*ccos(theta5r)*csin(theta3r) + a5*lambda6*mu2*nx*ccos(theta4r)*ccos(theta5r)*csin(theta3r) + 
       a4*lambda2*mu3*mu6*mx*csin(theta4r) + a4*lambda2*lambda6*mu3*nx*csin(theta4r) + a4*lambda3*mu2*mu6*mx*ccos(theta3r)*csin(theta4r) + 
       a4*lambda3*lambda6*mu2*nx*ccos(theta3r)*csin(theta4r) + a5*lambda2*mu3*mu6*mx*ccos(theta5r)*csin(theta4r) + 
       a5*lambda2*lambda6*mu3*nx*ccos(theta5r)*csin(theta4r) + a5*lambda3*mu2*mu6*mx*ccos(theta3r)*ccos(theta5r)*csin(theta4r) + 
       a5*lambda3*lambda6*mu2*nx*ccos(theta3r)*ccos(theta5r)*csin(theta4r) + a6*lambda5*lx*mu2*mu4*csin(theta3r)*csin(theta4r) + 
       d5*mu2*mu4*mu6*mx*csin(theta3r)*csin(theta4r) + d6*lambda5*mu2*mu4*mu6*mx*csin(theta3r)*csin(theta4r) + d5*lambda6*mu2*mu4*nx*csin(theta3r)*csin(theta4r) + 
       d6*lambda5*lambda6*mu2*mu4*nx*csin(theta3r)*csin(theta4r) - lambda5*mu2*mu4*px*csin(theta3r)*csin(theta4r) + 
       a6*lambda4*lx*mu2*mu5*ccos(theta5r)*csin(theta3r)*csin(theta4r) + d6*lambda4*mu2*mu5*mu6*mx*ccos(theta5r)*csin(theta3r)*csin(theta4r) + 
       d6*lambda4*lambda6*mu2*mu5*nx*ccos(theta5r)*csin(theta3r)*csin(theta4r) - lambda4*mu2*mu5*px*ccos(theta5r)*csin(theta3r)*csin(theta4r) + 
       a5*lambda2*lambda3*mu4*mu6*mx*csin(theta5r) + a5*lambda2*lambda3*lambda6*mu4*nx*csin(theta5r) - a5*mu2*mu3*mu4*mu6*mx*ccos(theta3r)*csin(theta5r) - 
       a5*lambda6*mu2*mu3*mu4*nx*ccos(theta3r)*csin(theta5r) + a5*lambda2*lambda4*mu3*mu6*mx*ccos(theta4r)*csin(theta5r) + 
       a5*lambda2*lambda4*lambda6*mu3*nx*ccos(theta4r)*csin(theta5r) + a5*lambda3*lambda4*mu2*mu6*mx*ccos(theta3r)*ccos(theta4r)*csin(theta5r) + 
       a5*lambda3*lambda4*lambda6*mu2*nx*ccos(theta3r)*ccos(theta4r)*csin(theta5r) + a6*lx*mu2*mu5*ccos(theta4r)*csin(theta3r)*csin(theta5r) + 
       d6*mu2*mu5*mu6*mx*ccos(theta4r)*csin(theta3r)*csin(theta5r) + d6*lambda6*mu2*mu5*nx*ccos(theta4r)*csin(theta3r)*csin(theta5r) - 
       mu2*mu5*px*ccos(theta4r)*csin(theta3r)*csin(theta5r) + a6*lambda2*lx*mu3*mu5*csin(theta4r)*csin(theta5r) + 
       d6*lambda2*mu3*mu5*mu6*mx*csin(theta4r)*csin(theta5r) + d6*lambda2*lambda6*mu3*mu5*nx*csin(theta4r)*csin(theta5r) - 
       lambda2*mu3*mu5*px*csin(theta4r)*csin(theta5r) + a6*lambda3*lx*mu2*mu5*ccos(theta3r)*csin(theta4r)*csin(theta5r) + 
       d6*lambda3*mu2*mu5*mu6*mx*ccos(theta3r)*csin(theta4r)*csin(theta5r) + d6*lambda3*lambda6*mu2*mu5*nx*ccos(theta3r)*csin(theta4r)*csin(theta5r) - 
       lambda3*mu2*mu5*px*ccos(theta3r)*csin(theta4r)*csin(theta5r) - a5*lambda4*mu2*mu6*mx*csin(theta3r)*csin(theta4r)*csin(theta5r) - 
       a5*lambda4*lambda6*mu2*nx*csin(theta3r)*csin(theta4r)*csin(theta5r)));
   
   sctheta1.sintheta = sintheta1;
   sctheta1.costheta = costheta1;

    return(sctheta1);

}
	   
