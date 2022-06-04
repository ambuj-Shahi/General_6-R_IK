/** 
 

 @file theta6.c
 @brief This file is used to calculate the value of theta6 with the values of theta1, theta2, theta3, theta4, theta5, and the D-H parameters  as input.
       .
 
 */


# include <stdio.h>
# include <math.h>
# include <complex.h>
# include "cx.h"

struct Theta theta6(double* coeff, double _Complex *thlist)
{

    struct Theta sctheta6;

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





double _Complex sintheta6 = ly*(-(mu1*ccos(theta1r)*cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + 
           cpow(mu1,2)*cpow(csin(theta1r),2),-1)*(cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + 
              cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(mu2*cpow(ccos(theta2r),2) + mu2*cpow(csin(theta2r),2))*
            (-(mu3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                   cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda5*ccos(theta5r)*
                    cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                      cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                      cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
                   cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                      cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                      cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                    (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))) + 
              cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
               (-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                      cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                    cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                      cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
                 mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r)\
                  - lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) + 
              lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                 cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
                 mu4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
                 lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                  csin(theta5r))) + cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
              cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(lambda2*cpow(ccos(theta2r),2) + lambda2*cpow(csin(theta2r),2))*
            (cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*(lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
               (lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
                 cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                  (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))) + 
              cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
               (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                    cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                  cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                    cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
                 lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                  csin(theta5r))))) + cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + 
         cpow(mu1,2)*cpow(csin(theta1r),2),-1)*(cpow(lambda1,2)*csin(theta1r) + cpow(mu1,2)*csin(theta1r))*
       (-(lambda2*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
              cpow(mu2,2)*cpow(csin(theta2r),2),-1)*csin(theta2r)*(-(mu3*ccos(theta3r)*
                 cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                   cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda5*ccos(theta5r)*
                    cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                      cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                      cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
                   cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                      cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                      cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                    (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))) + 
              cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
               (-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                      cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                    cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                      cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
                 mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r)\
                  - lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) + 
              lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                 cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
                 mu4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
                 lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                  csin(theta5r)))) + mu2*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*csin(theta2r)*(cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + 
               cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
             (lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))) + (ccos(theta2r)*cpow(lambda2,2) + ccos(theta2r)*cpow(mu2,2))*
          cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (mu3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))*csin(theta3r) + 
            (ccos(theta3r)*cpow(lambda3,2) + ccos(theta3r)*cpow(mu3,2))*
             cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                    cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                  cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                    cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
               mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r) - 
               lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) - 
            lambda3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*csin(theta3r)*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r)))) + lambda1*ccos(theta1r)*cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + 
         cpow(lambda1,2)*cpow(csin(theta1r),2) + cpow(mu1,2)*cpow(csin(theta1r),2),-1)*
       (lambda2*ccos(theta2r)*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(-(mu3*ccos(theta3r)*
               cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*(lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
                 cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                  (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             (-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
               mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r) - 
               lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               mu4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))) - mu2*ccos(theta2r)*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + 
            cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),-1)*
             (lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
             (lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))) + cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(cpow(lambda2,2)*csin(theta2r) + cpow(mu2,2)*csin(theta2r))*
          (mu3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))*csin(theta3r) + 
            (ccos(theta3r)*cpow(lambda3,2) + ccos(theta3r)*cpow(mu3,2))*
             cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                    cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                  cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                    cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
               mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r) - 
               lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) - 
            lambda3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*csin(theta3r)*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))))) + lz*(cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + 
         cpow(mu1,2)*cpow(csin(theta1r),2),-1)*(lambda1*cpow(ccos(theta1r),2) + lambda1*cpow(csin(theta1r),2))*
       (cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (mu2*cpow(ccos(theta2r),2) + mu2*cpow(csin(theta2r),2))*(-(mu3*ccos(theta3r)*
               cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*(lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
                 cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                  (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             (-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
               mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r) - 
               lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               mu4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))) + cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(lambda2*cpow(ccos(theta2r),2) + lambda2*cpow(csin(theta2r),2))*
          (cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),-1)*
             (lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
             (lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r)))) + cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + 
         cpow(mu1,2)*cpow(csin(theta1r),2),-1)*(mu1*cpow(ccos(theta1r),2) + mu1*cpow(csin(theta1r),2))*
       (lambda2*ccos(theta2r)*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(-(mu3*ccos(theta3r)*
               cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*(lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
                 cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                  (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             (-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
               mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r) - 
               lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               mu4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))) - mu2*ccos(theta2r)*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + 
            cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),-1)*
             (lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
             (lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))) + cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(cpow(lambda2,2)*csin(theta2r) + cpow(mu2,2)*csin(theta2r))*
          (mu3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))*csin(theta3r) + 
            (ccos(theta3r)*cpow(lambda3,2) + ccos(theta3r)*cpow(mu3,2))*
             cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                    cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                  cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                    cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
               mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r) - 
               lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) - 
            lambda3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*csin(theta3r)*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))))) + lx*(mu1*cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + 
         cpow(mu1,2)*cpow(csin(theta1r),2),-1)*csin(theta1r)*(cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + 
            cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(mu2*cpow(ccos(theta2r),2) + mu2*cpow(csin(theta2r),2))*
          (-(mu3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                 cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda5*ccos(theta5r)*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
                 cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                  (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             (-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
               mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r) - 
               lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               mu4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))) + cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(lambda2*cpow(ccos(theta2r),2) + lambda2*cpow(csin(theta2r),2))*
          (cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),-1)*
             (lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
             (lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r)))) + (ccos(theta1r)*cpow(lambda1,2) + ccos(theta1r)*cpow(mu1,2))*
       cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + cpow(mu1,2)*cpow(csin(theta1r),2),-1)*
       (-(lambda2*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
              cpow(mu2,2)*cpow(csin(theta2r),2),-1)*csin(theta2r)*(-(mu3*ccos(theta3r)*
                 cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                   cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda5*ccos(theta5r)*
                    cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                      cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                      cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
                   cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                      cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                      cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                    (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))) + 
              cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
               (-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                      cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                    cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                      cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
                 mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r)\
                  - lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) + 
              lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                 cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
                 mu4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
                 lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                  csin(theta5r)))) + mu2*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*csin(theta2r)*(cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + 
               cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
             (lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))) + (ccos(theta2r)*cpow(lambda2,2) + ccos(theta2r)*cpow(mu2,2))*
          cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (mu3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))*csin(theta3r) + 
            (ccos(theta3r)*cpow(lambda3,2) + ccos(theta3r)*cpow(mu3,2))*
             cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                    cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                  cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                    cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
               mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r) - 
               lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) - 
            lambda3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*csin(theta3r)*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r)))) - lambda1*cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + 
         cpow(mu1,2)*cpow(csin(theta1r),2),-1)*csin(theta1r)*(lambda2*ccos(theta2r)*
          cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (-(mu3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                 cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda5*ccos(theta5r)*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
                 cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                  (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             (-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
               mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r) - 
               lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               mu4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))) - mu2*ccos(theta2r)*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + 
            cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),-1)*
             (lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
             (lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r))) + cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(cpow(lambda2,2)*csin(theta2r) + cpow(mu2,2)*csin(theta2r))*
          (mu3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2)) + 
               cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
                 -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda4*cpow(ccos(theta4r),2) + lambda4*cpow(csin(theta4r),2))*
                (mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)))*csin(theta3r) + 
            (ccos(theta3r)*cpow(lambda3,2) + ccos(theta3r)*cpow(mu3,2))*
             cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(-(lambda4*lambda5*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                    cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                  cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                    cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)) + 
               mu4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2))*csin(theta4r) - 
               lambda5*(ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta5r)) - 
            lambda3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*csin(theta3r)*
             (lambda4*lambda5*ccos(theta4r)*ccos(theta5r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                  cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                  cpow(mu5,2)*cpow(csin(theta5r),2),-1) - mu4*ccos(theta4r)*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu5*cpow(ccos(theta5r),2) + mu5*cpow(csin(theta5r),2)) - 
               lambda5*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r))*
                csin(theta5r)))));

double _Complex costheta6 = lz*(cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + cpow(mu1,2)*cpow(csin(theta1r),2),-1)*
       (mu1*cpow(ccos(theta1r),2) + mu1*cpow(csin(theta1r),2))*(-(mu2*ccos(theta2r)*
            cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),
             -1)*(cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                 cpow(mu3,2)*cpow(csin(theta3r),2),-1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                 cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
               cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),
                -1)*(lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
               (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) + 
              cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
               ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
                 lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))))) + 
         cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (cpow(lambda2,2)*csin(theta2r) + cpow(mu2,2)*csin(theta2r))*
          (mu3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
               cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
               cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*csin(theta3r)*
             (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) - 
            lambda3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*csin(theta3r)*
             ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            (ccos(theta3r)*cpow(lambda3,2) + ccos(theta3r)*cpow(mu3,2))*
             cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))) + 
         lambda2*ccos(theta2r)*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(-(mu3*ccos(theta3r)*
               cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                 cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                 cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
               (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             ((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))))) + 
      cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + cpow(mu1,2)*cpow(csin(theta1r),2),-1)*
       (lambda1*cpow(ccos(theta1r),2) + lambda1*cpow(csin(theta1r),2))*
       (cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (lambda2*cpow(ccos(theta2r),2) + lambda2*cpow(csin(theta2r),2))*
          (cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),-1)*
             cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
              -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
               cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
             (mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))) + 
         cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (mu2*cpow(ccos(theta2r),2) + mu2*cpow(csin(theta2r),2))*(-(mu3*ccos(theta3r)*
               cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                 cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                 cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
               (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             ((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))))) + 
   lx*(-(lambda1*cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + 
           cpow(mu1,2)*cpow(csin(theta1r),2),-1)*csin(theta1r)*(-(mu2*ccos(theta2r)*
              cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),
               -1)*(cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                   cpow(mu3,2)*cpow(csin(theta3r),2),-1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                   cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                 cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                   cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
                 (mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) + 
                cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                   cpow(mu3,2)*cpow(csin(theta3r),2),-1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
                 ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                    cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                      cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                      cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
                   lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                      cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                      cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))))) + 
           cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
            (cpow(lambda2,2)*csin(theta2r) + cpow(mu2,2)*csin(theta2r))*
            (mu3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                 cpow(mu3,2)*cpow(csin(theta3r),2),-1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                 cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
               cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),
                -1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*csin(theta3r)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) - 
              lambda3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                 cpow(mu3,2)*cpow(csin(theta3r),2),-1)*csin(theta3r)*
               ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
                 lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
              (ccos(theta3r)*cpow(lambda3,2) + ccos(theta3r)*cpow(mu3,2))*
               cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
                 lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                  (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))) + 
           lambda2*ccos(theta2r)*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
              cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(-(mu3*ccos(theta3r)*
                 cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                   cpow(mu3,2)*cpow(csin(theta3r),2),-1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                   cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
                 cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
                   cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
                 (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
              lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                 cpow(mu3,2)*cpow(csin(theta3r),2),-1)*((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
                 lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
              cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
               ((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
                 lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                  (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))))) + 
      mu1*cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + cpow(mu1,2)*cpow(csin(theta1r),2),-1)*
       csin(theta1r)*(cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(lambda2*cpow(ccos(theta2r),2) + lambda2*cpow(csin(theta2r),2))*
          (cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),-1)*
             cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
              -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
               cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
             (mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))) + 
         cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (mu2*cpow(ccos(theta2r),2) + mu2*cpow(csin(theta2r),2))*(-(mu3*ccos(theta3r)*
               cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                 cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                 cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
               (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             ((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))))) + 
      (ccos(theta1r)*cpow(lambda1,2) + ccos(theta1r)*cpow(mu1,2))*cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + 
         cpow(lambda1,2)*cpow(csin(theta1r),2) + cpow(mu1,2)*cpow(csin(theta1r),2),-1)*
       (mu2*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),
           -1)*csin(theta2r)*(cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
               cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
             cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),
              -1)*(lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
             (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))) + 
         (ccos(theta2r)*cpow(lambda2,2) + ccos(theta2r)*cpow(mu2,2))*
          cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (mu3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
               cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
               cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*csin(theta3r)*
             (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) - 
            lambda3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*csin(theta3r)*
             ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            (ccos(theta3r)*cpow(lambda3,2) + ccos(theta3r)*cpow(mu3,2))*
             cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))) - 
         lambda2*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*csin(theta2r)*(-(mu3*ccos(theta3r)*
               cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                 cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                 cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
               (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             ((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))))) + 
   ly*(lambda1*ccos(theta1r)*cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + 
         cpow(mu1,2)*cpow(csin(theta1r),2),-1)*(-(mu2*ccos(theta2r)*
            cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),
             -1)*(cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
                 cpow(mu3,2)*cpow(csin(theta3r),2),-1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
                 cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
               cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),
                -1)*(lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
               (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) + 
              cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
               ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                  cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
                 lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                    cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                    cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))))) + 
         cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (cpow(lambda2,2)*csin(theta2r) + cpow(mu2,2)*csin(theta2r))*
          (mu3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
               cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
               cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*csin(theta3r)*
             (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) - 
            lambda3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*csin(theta3r)*
             ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            (ccos(theta3r)*cpow(lambda3,2) + ccos(theta3r)*cpow(mu3,2))*
             cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))) + 
         lambda2*ccos(theta2r)*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(-(mu3*ccos(theta3r)*
               cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                 cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                 cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
               (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             ((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))))) - 
      mu1*ccos(theta1r)*cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + 
         cpow(mu1,2)*cpow(csin(theta1r),2),-1)*(cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*(lambda2*cpow(ccos(theta2r),2) + lambda2*cpow(csin(theta2r),2))*
          (cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),-1)*
             cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),
              -1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + 
               cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*
             (mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))) + 
         cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (mu2*cpow(ccos(theta2r),2) + mu2*cpow(csin(theta2r),2))*(-(mu3*ccos(theta3r)*
               cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                 cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                 cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
               (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             ((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))))) + 
      cpow(cpow(lambda1,2)*cpow(ccos(theta1r),2) + cpow(mu1,2)*cpow(ccos(theta1r),2) + cpow(lambda1,2)*cpow(csin(theta1r),2) + cpow(mu1,2)*cpow(csin(theta1r),2),-1)*
       (cpow(lambda1,2)*csin(theta1r) + cpow(mu1,2)*csin(theta1r))*
       (mu2*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),
           -1)*csin(theta2r)*(cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + 
               cpow(lambda4,2)*cpow(csin(theta4r),2) + cpow(mu4,2)*cpow(csin(theta4r),2),-1)*
             cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),
              -1)*(lambda3*cpow(ccos(theta3r),2) + lambda3*cpow(csin(theta3r),2))*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
             (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(mu3*cpow(ccos(theta3r),2) + mu3*cpow(csin(theta3r),2))*
             ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))) + 
         (ccos(theta2r)*cpow(lambda2,2) + ccos(theta2r)*cpow(mu2,2))*
          cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + cpow(mu2,2)*cpow(csin(theta2r),2),-1)*
          (mu3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
               cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
               cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*csin(theta3r)*
             (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)) - 
            lambda3*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*csin(theta3r)*
             ((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            (ccos(theta3r)*cpow(lambda3,2) + ccos(theta3r)*cpow(mu3,2))*
             cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r)))) - 
         lambda2*cpow(cpow(lambda2,2)*cpow(ccos(theta2r),2) + cpow(mu2,2)*cpow(ccos(theta2r),2) + cpow(lambda2,2)*cpow(csin(theta2r),2) + 
            cpow(mu2,2)*cpow(csin(theta2r),2),-1)*csin(theta2r)*(-(mu3*ccos(theta3r)*
               cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
                -1)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                 cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                 cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(mu4*cpow(ccos(theta4r),2) + mu4*cpow(csin(theta4r),2))*
               (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            lambda3*ccos(theta3r)*cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + 
               cpow(mu3,2)*cpow(csin(theta3r),2),-1)*((ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda4,2)*csin(theta4r) + cpow(mu4,2)*csin(theta4r)) + 
               lambda4*ccos(theta4r)*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*(cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))) + 
            cpow(cpow(lambda3,2)*cpow(ccos(theta3r),2) + cpow(mu3,2)*cpow(ccos(theta3r),2) + cpow(lambda3,2)*cpow(csin(theta3r),2) + cpow(mu3,2)*cpow(csin(theta3r),2),
              -1)*(cpow(lambda3,2)*csin(theta3r) + cpow(mu3,2)*csin(theta3r))*
             ((ccos(theta4r)*cpow(lambda4,2) + ccos(theta4r)*cpow(mu4,2))*(ccos(theta5r)*cpow(lambda5,2) + ccos(theta5r)*cpow(mu5,2))*
                cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1) - 
               lambda4*cpow(cpow(lambda4,2)*cpow(ccos(theta4r),2) + cpow(mu4,2)*cpow(ccos(theta4r),2) + cpow(lambda4,2)*cpow(csin(theta4r),2) + 
                  cpow(mu4,2)*cpow(csin(theta4r),2),-1)*cpow(cpow(lambda5,2)*cpow(ccos(theta5r),2) + cpow(mu5,2)*cpow(ccos(theta5r),2) + 
                  cpow(lambda5,2)*cpow(csin(theta5r),2) + cpow(mu5,2)*cpow(csin(theta5r),2),-1)*csin(theta4r)*
                (cpow(lambda5,2)*csin(theta5r) + cpow(mu5,2)*csin(theta5r))))));
	  
    sctheta6.sintheta = sintheta6;
    sctheta6.costheta = costheta6;

    return(sctheta6);

}
	   
