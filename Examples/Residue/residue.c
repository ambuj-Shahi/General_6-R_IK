/** 
 

 @file residue.c
 @brief This file is used to calculate residue1 to residue6 obtained by back-substituting the IK solutions into Eq.1-Eq.6 of the paper by Manocha and Canny
        and then taking the difference of the L.H.S. and R.H.S. of the aforementioned equations.
       .
 
 */


# include <stdio.h>
# include <math.h>
# include <complex.h>
# include "cx.h"

struct Residue residue(double *coeff, double _Complex *thlist)
{

    struct Residue resd;

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
    double _Complex theta2r = thlist[1];
    double _Complex theta3r = thlist[2];
    double _Complex theta4r = thlist[3];
    double _Complex theta5r = thlist[4];




    double _Complex res1 = a2 + ccos(theta2r)*(a1 + (a6*lx + d6*mu6*mx + d6*lambda6*nx - px)*ccos(theta1r) + (a6*ly + d6*mu6*my + d6*lambda6*ny - py)*csin(theta1r)) + 
                (mu1*(d1 + a6*lz + d6*mu6*mz + d6*lambda6*nz - pz) + lambda1*
                ((a6*ly + d6*mu6*my + d6*lambda6*ny - py)*ccos(theta1r) - (a6*lx + d6*mu6*mx + d6*lambda6*nx - px)*csin(theta1r)))*csin(theta2r) + 
                ccos(theta3r)*(a3 + ccos(theta4r)*(a4 + a5*ccos(theta5r)) + d5*mu4*csin(theta4r) - a5*lambda4*csin(theta4r)*csin(theta5r)) + 
                csin(theta3r)*(mu3*(d4 + d5*lambda4 + a5*mu4*csin(theta5r)) + 
                lambda3*(-((a4 + a5*ccos(theta5r))*csin(theta4r)) + ccos(theta4r)*(d5*mu4 - a5*lambda4*csin(theta5r))));

   double _Complex res2 = -(mu1*(d1 + a6*lz + d6*mu6*mz + d6*lambda6*nz - pz)*ccos(theta2r)) - 
                lambda1*ccos(theta2r)*((a6*ly + d6*mu6*my + d6*lambda6*ny - py)*ccos(theta1r) - (a6*lx + d6*mu6*mx + d6*lambda6*nx - px)*csin(theta1r)) + 
                (a1 + (a6*lx + d6*mu6*mx + d6*lambda6*nx - px)*ccos(theta1r) + (a6*ly + d6*mu6*my + d6*lambda6*ny - py)*csin(theta1r))*csin(theta2r) + 
                mu2*(d3 + a4*mu3*csin(theta4r) + a5*mu3*ccos(theta5r)*csin(theta4r) + mu3*ccos(theta4r)*(-(d5*mu4) + a5*lambda4*csin(theta5r)) + 
                lambda3*(d4 + d5*lambda4 + a5*mu4*csin(theta5r))) + lambda2*
                (-(csin(theta3r)*(a3 + ccos(theta4r)*(a4 + a5*ccos(theta5r)) + d5*mu4*csin(theta4r) - a5*lambda4*csin(theta4r)*csin(theta5r))) + 
                ccos(theta3r)*(mu3*(d4 + d5*lambda4 + a5*mu4*csin(theta5r)) + 
                lambda3*(-((a4 + a5*ccos(theta5r))*csin(theta4r)) + ccos(theta4r)*(d5*mu4 - a5*lambda4*csin(theta5r)))));

    double _Complex res3 = d2 + lambda1*(d1 + a6*lz + d6*mu6*mz + d6*lambda6*nz - pz) - mu1*
                ((a6*ly + d6*mu6*my + d6*lambda6*ny - py)*ccos(theta1r) - (a6*lx + d6*mu6*mx + d6*lambda6*nx - px)*csin(theta1r)) + 
                lambda2*(d3 + a4*mu3*csin(theta4r) + a5*mu3*ccos(theta5r)*csin(theta4r) + mu3*ccos(theta4r)*(-(d5*mu4) + a5*lambda4*csin(theta5r)) + 
                lambda3*(d4 + d5*lambda4 + a5*mu4*csin(theta5r))) + mu2*(csin(theta3r)*
                (a3 + ccos(theta4r)*(a4 + a5*ccos(theta5r)) + d5*mu4*csin(theta4r) - a5*lambda4*csin(theta4r)*csin(theta5r)) - 
                ccos(theta3r)*(mu3*(d4 + d5*lambda4 + a5*mu4*csin(theta5r)) + 
                lambda3*(-((a4 + a5*ccos(theta5r))*csin(theta4r)) + ccos(theta4r)*(d5*mu4 - a5*lambda4*csin(theta5r)))));

    double _Complex res4 = -(ccos(theta2r)*(mu6*(mx*ccos(theta1r) + my*csin(theta1r)) + lambda6*(nx*ccos(theta1r) + ny*csin(theta1r)))) - 
                (mu1*(mu6*mz + lambda6*nz) + lambda1*(mu6*(my*ccos(theta1r) - mx*csin(theta1r)) + lambda6*(ny*ccos(theta1r) - nx*csin(theta1r))))*csin(theta2r) + 
                ccos(theta3r)*(lambda5*mu4*csin(theta4r) + mu5*(lambda4*ccos(theta5r)*csin(theta4r) + ccos(theta4r)*csin(theta5r))) + 
                csin(theta3r)*(-(mu3*mu4*mu5*ccos(theta5r)) + lambda4*(lambda5*mu3 + lambda3*mu5*ccos(theta4r)*ccos(theta5r)) + 
                lambda3*(lambda5*mu4*ccos(theta4r) - mu5*csin(theta4r)*csin(theta5r)));

   double _Complex res5 = ccos(theta2r)*(mu1*(mu6*mz + lambda6*nz) + lambda1*(mu6*(my*ccos(theta1r) - mx*csin(theta1r)) + lambda6*(ny*ccos(theta1r) - nx*csin(theta1r)))) - 
                (mu6*(mx*ccos(theta1r) + my*csin(theta1r)) + lambda6*(nx*ccos(theta1r) + ny*csin(theta1r)))*csin(theta2r) + 
                mu2*(lambda5*(lambda3*lambda4 - mu3*mu4*ccos(theta4r)) - mu5*(lambda3*mu4 + lambda4*mu3*ccos(theta4r))*ccos(theta5r) + 
                mu3*mu5*csin(theta4r)*csin(theta5r)) - lambda2*(csin(theta3r)*
                (lambda5*mu4*csin(theta4r) + mu5*(lambda4*ccos(theta5r)*csin(theta4r) + ccos(theta4r)*csin(theta5r))) - 
                ccos(theta3r)*(-(mu3*mu4*mu5*ccos(theta5r)) + lambda4*(lambda5*mu3 + lambda3*mu5*ccos(theta4r)*ccos(theta5r)) + 
                lambda3*(lambda5*mu4*ccos(theta4r) - mu5*csin(theta4r)*csin(theta5r))));

    double _Complex res6 = -(lambda1*(mu6*mz + lambda6*nz)) - mu1*(mu6*(-(my*ccos(theta1r)) + mx*csin(theta1r)) + lambda6*(-(ny*ccos(theta1r)) + nx*csin(theta1r))) + 
                lambda2*(lambda5*(lambda3*lambda4 - mu3*mu4*ccos(theta4r)) - mu5*(lambda3*mu4 + lambda4*mu3*ccos(theta4r))*ccos(theta5r) + 
                mu3*mu5*csin(theta4r)*csin(theta5r)) + mu2*(csin(theta3r)*(lambda5*mu4*csin(theta4r) + 
                mu5*(lambda4*ccos(theta5r)*csin(theta4r) + ccos(theta4r)*csin(theta5r))) - 
                ccos(theta3r)*(-(mu3*mu4*mu5*ccos(theta5r)) + lambda4*(lambda5*mu3 + lambda3*mu5*ccos(theta4r)*ccos(theta5r)) + 
                lambda3*(lambda5*mu4*ccos(theta4r) - mu5*csin(theta4r)*csin(theta5r))));

    resd.residue1 = res1;
    resd.residue2 = res2;
    resd.residue3 = res3;
    resd.residue4 = res4;
    resd.residue5 = res5;
    resd.residue6 = res6;
		 

    return(resd);

}
	   
