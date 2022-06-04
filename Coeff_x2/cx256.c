
# include <stdio.h>
# include <math.h>

double cx256(double* coeff)
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



double dummyexp = (mu1*pow(a1,-1)*(-4*a2*a3*lambda5*mu1*mu4 - 4*a4*a5*lambda5*mu1*mu4 - 4*d3*d4*lambda3*lambda5*mu1*mu4 + 4*a6*d1*lambda5*lz*mu1*mu4 - 4*d2*d4*lambda5*mu1*mu2*mu3*mu4 - 4*d4*d5*mu1*mu5 - 4*d3*d5*lambda3*mu1*mu5 - 4*a2*a3*lambda4*mu1*mu5 - 4*a4*a5*lambda4*mu1*mu5 - 4*d3*d4*lambda3*lambda4*mu1*mu5 + 4*a6*d1*lambda4*lz*mu1*mu5 - 4*d2*d5*mu1*mu2*mu3*mu5 - 4*d2*d4*lambda4*mu1*mu2*mu3*mu5 + 4*a6*d5*lx*mu1*mu4*mu6*mx + 4*a6*d6*lambda5*lx*mu1*mu4*mu6*mx + 4*a6*d6*lambda4*lx*mu1*mu5*mu6*mx + 4*a6*d5*ly*mu1*mu4*mu6*my + 4*a6*d6*lambda5*ly*mu1*mu4*mu6*my + 4*a6*d6*lambda4*ly*mu1*mu5*mu6*my + 4*a1*a4*mu2*mu3*mu6*mz + 4*a1*a5*mu2*mu3*mu6*mz + 4*d1*d5*mu1*mu4*mu6*mz + 4*d1*d6*lambda5*mu1*mu4*mu6*mz + 4*a6*d5*lz*mu1*mu4*mu6*mz + 4*a6*d6*lambda5*lz*mu1*mu4*mu6*mz + 4*d1*d6*lambda4*mu1*mu5*mu6*mz + 4*a6*d6*lambda4*lz*mu1*mu5*mu6*mz + 4*a6*d5*lambda6*lx*mu1*mu4*nx + 4*a6*d6*lambda5*lambda6*lx*mu1*mu4*nx + 4*a6*d6*lambda4*lambda6*lx*mu1*mu5*nx + 8*d5*d6*lambda6*mu1*mu4*mu6*mx*nx + 4*a6*d5*lambda6*ly*mu1*mu4*ny + 4*a6*d6*lambda5*lambda6*ly*mu1*mu4*ny + 4*a6*d6*lambda4*lambda6*ly*mu1*mu5*ny + 8*d5*d6*lambda6*mu1*mu4*mu6*my*ny + 4*a1*a4*lambda6*mu2*mu3*nz + 4*a1*a5*lambda6*mu2*mu3*nz + 4*d1*d5*lambda6*mu1*mu4*nz + 4*d1*d6*lambda5*lambda6*mu1*mu4*nz + 4*a6*d5*lambda6*lz*mu1*mu4*nz + 4*a6*d6*lambda5*lambda6*lz*mu1*mu4*nz + 4*d1*d6*lambda4*lambda6*mu1*mu5*nz + 4*a6*d6*lambda4*lambda6*lz*mu1*mu5*nz + 8*d5*d6*lambda6*mu1*mu4*mu6*mz*nz - 4*a6*lambda5*lx*mu1*mu4*px - 4*a6*lambda4*lx*mu1*mu5*px - 4*d5*mu1*mu4*mu6*mx*px - 4*d6*lambda5*mu1*mu4*mu6*mx*px - 4*d6*lambda4*mu1*mu5*mu6*mx*px - 4*d5*lambda6*mu1*mu4*nx*px - 4*d6*lambda5*lambda6*mu1*mu4*nx*px - 4*d6*lambda4*lambda6*mu1*mu5*nx*px - 4*a6*lambda5*ly*mu1*mu4*py - 4*a6*lambda4*ly*mu1*mu5*py - 4*d5*mu1*mu4*mu6*my*py - 4*d6*lambda5*mu1*mu4*mu6*my*py - 4*d6*lambda4*mu1*mu5*mu6*my*py - 4*d5*lambda6*mu1*mu4*ny*py - 4*d6*lambda5*lambda6*mu1*mu4*ny*py - 4*d6*lambda4*lambda6*mu1*mu5*ny*py - 4*d1*lambda5*mu1*mu4*pz - 4*a6*lambda5*lz*mu1*mu4*pz - 4*d1*lambda4*mu1*mu5*pz - 4*a6*lambda4*lz*mu1*mu5*pz - 4*d5*mu1*mu4*mu6*mz*pz - 4*d6*lambda5*mu1*mu4*mu6*mz*pz - 4*d6*lambda4*mu1*mu5*mu6*mz*pz - 4*d5*lambda6*mu1*mu4*nz*pz - 4*d6*lambda5*lambda6*mu1*mu4*nz*pz - 4*d6*lambda4*lambda6*mu1*mu5*nz*pz - 4*a1*(a4 + a5)*lambda1*cos(alpha4 + alpha5) - 2*lambda5*mu1*mu4*pow(a1,2) - 2*lambda4*mu1*mu5*pow(a1,2) + 2*lambda5*mu1*mu4*pow(a2,2) + 2*lambda4*mu1*mu5*pow(a2,2) + 2*lambda5*mu1*mu4*pow(a3,2) + 2*lambda4*mu1*mu5*pow(a3,2) - 2*lambda5*mu1*mu4*pow(a4,2) - 2*lambda4*mu1*mu5*pow(a4,2) - 2*lambda5*mu1*mu4*pow(a5,2) - 2*lambda4*mu1*mu5*pow(a5,2) + 2*lambda5*mu1*mu4*pow(d1,2) + 2*lambda4*mu1*mu5*pow(d1,2) - 2*lambda5*mu1*mu4*pow(d2,2) - 2*lambda4*mu1*mu5*pow(d2,2) - 2*lambda5*mu1*mu4*pow(d3,2) - 2*lambda4*mu1*mu5*pow(d3,2) - 2*lambda5*mu1*mu4*pow(d4,2) - 2*lambda4*mu1*mu5*pow(d4,2) + 2*lambda5*mu1*mu4*pow(d5,2) - 2*lambda4*mu1*mu5*pow(d5,2) + 4*lambda5*lambda6*mu1*mu4*mu6*mx*nx*pow(d6,2) + 4*lambda4*lambda6*mu1*mu5*mu6*mx*nx*pow(d6,2) + 4*lambda5*lambda6*mu1*mu4*mu6*my*ny*pow(d6,2) + 4*lambda4*lambda6*mu1*mu5*mu6*my*ny*pow(d6,2) + 4*lambda5*lambda6*mu1*mu4*mu6*mz*nz*pow(d6,2) + 4*lambda4*lambda6*mu1*mu5*mu6*mz*nz*pow(d6,2) + 2*lambda5*mu1*mu4*pow(a6,2)*pow(lx,2) + 2*lambda4*mu1*mu5*pow(a6,2)*pow(lx,2) + 2*lambda5*mu1*mu4*pow(a6,2)*pow(ly,2) + 2*lambda4*mu1*mu5*pow(a6,2)*pow(ly,2) + 2*lambda5*mu1*mu4*pow(a6,2)*pow(lz,2) + 2*lambda4*mu1*mu5*pow(a6,2)*pow(lz,2) + 2*d5*d6*mu1*mu4*pow(mx,2) + lambda5*mu1*mu4*pow(d6,2)*pow(mx,2) + lambda4*mu1*mu5*pow(d6,2)*pow(mx,2) - 2*d5*d6*mu1*mu4*pow(lambda6,2)*pow(mx,2) - lambda5*mu1*mu4*pow(d6,2)*pow(lambda6,2)*pow(mx,2) - lambda4*mu1*mu5*pow(d6,2)*pow(lambda6,2)*pow(mx,2) + 2*d5*d6*mu1*mu4*pow(mu6,2)*pow(mx,2) + lambda5*mu1*mu4*pow(d6,2)*pow(mu6,2)*pow(mx,2) + lambda4*mu1*mu5*pow(d6,2)*pow(mu6,2)*pow(mx,2) + 2*d5*d6*mu1*mu4*pow(my,2) + lambda5*mu1*mu4*pow(d6,2)*pow(my,2) + lambda4*mu1*mu5*pow(d6,2)*pow(my,2) - 2*d5*d6*mu1*mu4*pow(lambda6,2)*pow(my,2) - lambda5*mu1*mu4*pow(d6,2)*pow(lambda6,2)*pow(my,2) - lambda4*mu1*mu5*pow(d6,2)*pow(lambda6,2)*pow(my,2) + 2*d5*d6*mu1*mu4*pow(mu6,2)*pow(my,2) + lambda5*mu1*mu4*pow(d6,2)*pow(mu6,2)*pow(my,2) + lambda4*mu1*mu5*pow(d6,2)*pow(mu6,2)*pow(my,2) + 2*d5*d6*mu1*mu4*pow(mz,2) + lambda5*mu1*mu4*pow(d6,2)*pow(mz,2) + lambda4*mu1*mu5*pow(d6,2)*pow(mz,2) - 2*d5*d6*mu1*mu4*pow(lambda6,2)*pow(mz,2) - lambda5*mu1*mu4*pow(d6,2)*pow(lambda6,2)*pow(mz,2) - lambda4*mu1*mu5*pow(d6,2)*pow(lambda6,2)*pow(mz,2) + 2*d5*d6*mu1*mu4*pow(mu6,2)*pow(mz,2) + lambda5*mu1*mu4*pow(d6,2)*pow(mu6,2)*pow(mz,2) + lambda4*mu1*mu5*pow(d6,2)*pow(mu6,2)*pow(mz,2) + 2*d5*d6*mu1*mu4*pow(nx,2) + lambda5*mu1*mu4*pow(d6,2)*pow(nx,2) + lambda4*mu1*mu5*pow(d6,2)*pow(nx,2) + 2*d5*d6*mu1*mu4*pow(lambda6,2)*pow(nx,2) + lambda5*mu1*mu4*pow(d6,2)*pow(lambda6,2)*pow(nx,2) + lambda4*mu1*mu5*pow(d6,2)*pow(lambda6,2)*pow(nx,2) - 2*d5*d6*mu1*mu4*pow(mu6,2)*pow(nx,2) - lambda5*mu1*mu4*pow(d6,2)*pow(mu6,2)*pow(nx,2) - lambda4*mu1*mu5*pow(d6,2)*pow(mu6,2)*pow(nx,2) + 2*d5*d6*mu1*mu4*pow(ny,2) + lambda5*mu1*mu4*pow(d6,2)*pow(ny,2) + lambda4*mu1*mu5*pow(d6,2)*pow(ny,2) + 2*d5*d6*mu1*mu4*pow(lambda6,2)*pow(ny,2) + lambda5*mu1*mu4*pow(d6,2)*pow(lambda6,2)*pow(ny,2) + lambda4*mu1*mu5*pow(d6,2)*pow(lambda6,2)*pow(ny,2) - 2*d5*d6*mu1*mu4*pow(mu6,2)*pow(ny,2) - lambda5*mu1*mu4*pow(d6,2)*pow(mu6,2)*pow(ny,2) - lambda4*mu1*mu5*pow(d6,2)*pow(mu6,2)*pow(ny,2) + 2*d5*d6*mu1*mu4*pow(nz,2) + lambda5*mu1*mu4*pow(d6,2)*pow(nz,2) + lambda4*mu1*mu5*pow(d6,2)*pow(nz,2) + 2*d5*d6*mu1*mu4*pow(lambda6,2)*pow(nz,2) + lambda5*mu1*mu4*pow(d6,2)*pow(lambda6,2)*pow(nz,2) + lambda4*mu1*mu5*pow(d6,2)*pow(lambda6,2)*pow(nz,2) - 2*d5*d6*mu1*mu4*pow(mu6,2)*pow(nz,2) - lambda5*mu1*mu4*pow(d6,2)*pow(mu6,2)*pow(nz,2) - lambda4*mu1*mu5*pow(d6,2)*pow(mu6,2)*pow(nz,2) + 2*lambda5*mu1*mu4*pow(px,2) + 2*lambda4*mu1*mu5*pow(px,2) + 2*lambda5*mu1*mu4*pow(py,2) + 2*lambda4*mu1*mu5*pow(py,2) + 2*lambda5*mu1*mu4*pow(pz,2) + 2*lambda4*mu1*mu5*pow(pz,2) + 4*lambda2*(lambda3*(-(d2*d4*lambda5*mu1*mu4) - d2*d5*mu1*mu5 - d2*d4*lambda4*mu1*mu5 + a1*a4*mu6*mz + a1*a5*mu6*mz + a1*(a4 + a5)*lambda6*nz) - d2*d3*mu1*sin(alpha4 + alpha5))))/2.;

return(dummyexp);

}