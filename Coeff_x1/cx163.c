
# include <stdio.h>
# include <math.h>

double cx163(double* coeff)
{

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



double dummyexp = -(mu1*pow(a1,-1)*(lambda2*(-2*a1*lambda1*(d5*mu5 - d3*sin(alpha3 - alpha4 - alpha5) + d4*sin(alpha4 + alpha5)) + (a3 - a4 - a5)*mu1*(-2*d5*lambda5 - 2*a6*lx*mu6*mx - 2*a6*ly*mu6*my - 2*d1*mu6*mz - 2*a6*lz*mu6*mz - 2*a6*lambda6*lx*nx - 2*a6*lambda6*ly*ny - 2*d1*lambda6*nz - 2*a6*lambda6*lz*nz + 2*mu6*mx*px + 2*lambda6*nx*px + 2*mu6*my*py + 2*lambda6*ny*py + 2*mu6*mz*pz + 2*lambda6*nz*pz - 2*d3*cos(alpha3 - alpha4 - alpha5) - 2*d4*cos(alpha4 + alpha5) - d6*pow(mx,2) + d6*cos(2*alpha6)*pow(mx,2) - d6*pow(my,2) + d6*cos(2*alpha6)*pow(my,2) - d6*pow(mz,2) + d6*cos(2*alpha6)*pow(mz,2) - d6*pow(nx,2) - d6*cos(2*alpha6)*pow(nx,2) - d6*pow(ny,2) - d6*cos(2*alpha6)*pow(ny,2) - d6*pow(nz,2) - d6*cos(2*alpha6)*pow(nz,2) - 2*d6*mx*nx*sin(2*alpha6) - 2*d6*my*ny*sin(2*alpha6) - 2*d6*mz*nz*sin(2*alpha6))) + (2*a1*d4*mz*cos(alpha3 - alpha6) + 2*a1*d5*mz*cos(alpha3 - alpha4 - alpha6) + 2*a1*d6*mz*cos(alpha3 - alpha4 - alpha5 - alpha6) - 2*a1*d4*mz*cos(alpha3 + alpha6) - 2*a1*d5*mz*cos(alpha3 - alpha4 + alpha6) - 2*a1*d6*mz*cos(alpha3 - alpha4 - alpha5 + alpha6) + a2*d5*sin(alpha1 - alpha2 - alpha5) - a2*d5*sin(alpha1 + alpha2 - alpha5) + a2*d4*sin(alpha1 - alpha2 - alpha4 - alpha5) - a2*d4*sin(alpha1 + alpha2 - alpha4 - alpha5) + 4*a1*d1*sin(alpha3 - alpha4 - alpha5) + 4*a1*a6*lz*sin(alpha3 - alpha4 - alpha5) - 4*a1*pz*sin(alpha3 - alpha4 - alpha5) + 2*a1*d2*sin(alpha1 + alpha3 - alpha4 - alpha5) - 2*a3*d2*sin(alpha1 + alpha3 - alpha4 - alpha5) + 2*a4*d2*sin(alpha1 + alpha3 - alpha4 - alpha5) + 2*a5*d2*sin(alpha1 + alpha3 - alpha4 - alpha5) + a2*d3*sin(alpha1 - alpha2 + alpha3 - alpha4 - alpha5) - a2*d3*sin(alpha1 + alpha2 + alpha3 - alpha4 - alpha5) - a2*d5*sin(alpha1 - alpha2 + alpha5) + a2*d5*sin(alpha1 + alpha2 + alpha5) - a2*d4*sin(alpha1 - alpha2 + alpha4 + alpha5) + a2*d4*sin(alpha1 + alpha2 + alpha4 + alpha5) - 2*a1*d2*sin(alpha1 - alpha3 + alpha4 + alpha5) - 2*a3*d2*sin(alpha1 - alpha3 + alpha4 + alpha5) + 2*a4*d2*sin(alpha1 - alpha3 + alpha4 + alpha5) + 2*a5*d2*sin(alpha1 - alpha3 + alpha4 + alpha5) - a2*d3*sin(alpha1 - alpha2 - alpha3 + alpha4 + alpha5) + a2*d3*sin(alpha1 + alpha2 - alpha3 + alpha4 + alpha5) + 2*a1*d4*nz*sin(alpha3 - alpha6) + 2*a1*d5*nz*sin(alpha3 - alpha4 - alpha6) + 2*a1*d6*nz*sin(alpha3 - alpha4 - alpha5 - alpha6) + 2*a1*d4*nz*sin(alpha3 + alpha6) + 2*a1*d5*nz*sin(alpha3 - alpha4 + alpha6) + 2*a1*d6*nz*sin(alpha3 - alpha4 - alpha5 + alpha6))/2.));

return(dummyexp);

}