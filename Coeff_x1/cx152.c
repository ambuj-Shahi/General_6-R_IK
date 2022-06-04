
# include <stdio.h>
# include <math.h>

double cx152(double* coeff)
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



double dummyexp = -2*mu1*pow(a1,-1)*(2*a5*lambda5*mu1*(d3 + d2*lambda2 + d5*mu3*mu4) - 2*a1*d3*lambda1*mu5 - 2*a1*d1*lambda2*mu5 - 2*a1*d2*lambda1*lambda2*mu5 - 2*a1*a6*lambda2*lz*mu5 - 2*a2*d2*mu1*mu2*mu5 - 2*a3*d4*mu1*mu3*mu5 + 2*a4*d4*mu1*mu3*mu5 - 2*a3*d5*lambda4*mu1*mu3*mu5 + 2*a4*d5*lambda4*mu1*mu3*mu5 - 2*a1*d5*lambda1*mu3*mu4*mu5 + 2*a5*a6*lx*mu1*mu3*mu4*mu6*mx + 2*a5*a6*ly*mu1*mu3*mu4*mu6*my + 2*a5*d1*mu1*mu3*mu4*mu6*mz + 2*a5*a6*lz*mu1*mu3*mu4*mu6*mz - 2*a1*d6*lambda2*mu5*mu6*mz + 2*a5*a6*lambda6*lx*mu1*mu3*mu4*nx + 4*a5*d6*lambda6*mu1*mu3*mu4*mu6*mx*nx + 2*a5*a6*lambda6*ly*mu1*mu3*mu4*ny + 4*a5*d6*lambda6*mu1*mu3*mu4*mu6*my*ny + 2*a5*d1*lambda6*mu1*mu3*mu4*nz + 2*a5*a6*lambda6*lz*mu1*mu3*mu4*nz - 2*a1*d6*lambda2*lambda6*mu5*nz + 4*a5*d6*lambda6*mu1*mu3*mu4*mu6*mz*nz - 2*a5*mu1*mu3*mu4*mu6*mx*px - 2*a5*lambda6*mu1*mu3*mu4*nx*px - 2*a5*mu1*mu3*mu4*mu6*my*py - 2*a5*lambda6*mu1*mu3*mu4*ny*py + 2*a1*lambda2*mu5*pz - 2*a5*mu1*mu3*mu4*mu6*mz*pz - 2*a5*lambda6*mu1*mu3*mu4*nz*pz + a5*d6*mu1*mu3*mu4*pow(mx,2) - a5*d6*mu1*mu3*mu4*pow(lambda6,2)*pow(mx,2) + a5*d6*mu1*mu3*mu4*pow(mu6,2)*pow(mx,2) + a5*d6*mu1*mu3*mu4*pow(my,2) - a5*d6*mu1*mu3*mu4*pow(lambda6,2)*pow(my,2) + a5*d6*mu1*mu3*mu4*pow(mu6,2)*pow(my,2) + a5*d6*mu1*mu3*mu4*pow(mz,2) - a5*d6*mu1*mu3*mu4*pow(lambda6,2)*pow(mz,2) + a5*d6*mu1*mu3*mu4*pow(mu6,2)*pow(mz,2) + a5*d6*mu1*mu3*mu4*pow(nx,2) + a5*d6*mu1*mu3*mu4*pow(lambda6,2)*pow(nx,2) - a5*d6*mu1*mu3*mu4*pow(mu6,2)*pow(nx,2) + a5*d6*mu1*mu3*mu4*pow(ny,2) + a5*d6*mu1*mu3*mu4*pow(lambda6,2)*pow(ny,2) - a5*d6*mu1*mu3*mu4*pow(mu6,2)*pow(ny,2) + a5*d6*mu1*mu3*mu4*pow(nz,2) + a5*d6*mu1*mu3*mu4*pow(lambda6,2)*pow(nz,2) - a5*d6*mu1*mu3*mu4*pow(mu6,2)*pow(nz,2) + lambda3*(2*(a5*d4*lambda5*mu1 - (a1*d4*lambda1 + (-a3 + a4)*d5*mu1*mu4)*mu5) + lambda4*(-2*a1*d5*lambda1*mu5 - a5*mu1*(-2*d5*lambda5 - 2*a6*lx*mu6*mx - 2*a6*ly*mu6*my - 2*d1*mu6*mz - 2*a6*lz*mu6*mz + 2*mu6*mx*px + 2*mu6*my*py + 2*mu6*mz*pz - 2*lambda6*(a6*lx*nx + a6*ly*ny + d1*nz + a6*lz*nz - nx*px - ny*py - nz*pz) - d6*pow(mx,2) + d6*cos(2*alpha6)*pow(mx,2) - d6*pow(my,2) + d6*cos(2*alpha6)*pow(my,2) - d6*pow(mz,2) + d6*cos(2*alpha6)*pow(mz,2) - d6*pow(nx,2) - d6*cos(2*alpha6)*pow(nx,2) - d6*pow(ny,2) - d6*cos(2*alpha6)*pow(ny,2) - d6*pow(nz,2) - d6*cos(2*alpha6)*pow(nz,2) - 2*d6*mx*nx*sin(2*alpha6) - 2*d6*my*ny*sin(2*alpha6) - 2*d6*mz*nz*sin(2*alpha6)))));

return(dummyexp);

}