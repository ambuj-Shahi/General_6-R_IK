
# include <stdio.h>
# include <math.h>

double cx237(double* coeff)
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



double dummyexp = pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1)*(a2*a6*lambda3*lambda4*lambda5*ly*mu2*mu6*mx - a3*a6*lambda3*lambda4*lambda5*ly*mu2*mu6*mx - a4*a6*lambda3*lambda4*lambda5*ly*mu2*mu6*mx + a5*a6*lambda3*lambda4*lambda5*ly*mu2*mu6*mx + a6*d1*lambda4*lambda5*lx*mu2*mu3*mu6*mx + a6*d1*lambda3*lambda5*lx*mu2*mu4*mu6*mx - a2*a6*lambda5*ly*mu2*mu3*mu4*mu6*mx + a3*a6*lambda5*ly*mu2*mu3*mu4*mu6*mx + a4*a6*lambda5*ly*mu2*mu3*mu4*mu6*mx - a5*a6*lambda5*ly*mu2*mu3*mu4*mu6*mx - a6*d1*lambda3*lambda4*lx*mu2*mu5*mu6*mx + a2*a6*lambda4*ly*mu2*mu3*mu5*mu6*mx - a3*a6*lambda4*ly*mu2*mu3*mu5*mu6*mx - a4*a6*lambda4*ly*mu2*mu3*mu5*mu6*mx + a5*a6*lambda4*ly*mu2*mu3*mu5*mu6*mx + a2*a6*lambda3*ly*mu2*mu4*mu5*mu6*mx - a3*a6*lambda3*ly*mu2*mu4*mu5*mu6*mx - a4*a6*lambda3*ly*mu2*mu4*mu5*mu6*mx + a5*a6*lambda3*ly*mu2*mu4*mu5*mu6*mx + a6*d1*lx*mu2*mu3*mu4*mu5*mu6*mx - a2*a6*lambda3*lambda4*lambda5*lx*mu2*mu6*my + a3*a6*lambda3*lambda4*lambda5*lx*mu2*mu6*my + a4*a6*lambda3*lambda4*lambda5*lx*mu2*mu6*my - a5*a6*lambda3*lambda4*lambda5*lx*mu2*mu6*my + a6*d1*lambda4*lambda5*ly*mu2*mu3*mu6*my + a6*d1*lambda3*lambda5*ly*mu2*mu4*mu6*my + a2*a6*lambda5*lx*mu2*mu3*mu4*mu6*my - a3*a6*lambda5*lx*mu2*mu3*mu4*mu6*my - a4*a6*lambda5*lx*mu2*mu3*mu4*mu6*my + a5*a6*lambda5*lx*mu2*mu3*mu4*mu6*my - a6*d1*lambda3*lambda4*ly*mu2*mu5*mu6*my - a2*a6*lambda4*lx*mu2*mu3*mu5*mu6*my + a3*a6*lambda4*lx*mu2*mu3*mu5*mu6*my + a4*a6*lambda4*lx*mu2*mu3*mu5*mu6*my - a5*a6*lambda4*lx*mu2*mu3*mu5*mu6*my - a2*a6*lambda3*lx*mu2*mu4*mu5*mu6*my + a3*a6*lambda3*lx*mu2*mu4*mu5*mu6*my + a4*a6*lambda3*lx*mu2*mu4*mu5*mu6*my - a5*a6*lambda3*lx*mu2*mu4*mu5*mu6*my + a6*d1*ly*mu2*mu3*mu4*mu5*mu6*my - d1*lambda4*lambda5*mu2*mu3*mu6*mx*px - a6*lambda4*lambda5*lz*mu2*mu3*mu6*mx*px - d1*lambda3*lambda5*mu2*mu4*mu6*mx*px - a6*lambda3*lambda5*lz*mu2*mu4*mu6*mx*px + d1*lambda3*lambda4*mu2*mu5*mu6*mx*px + a6*lambda3*lambda4*lz*mu2*mu5*mu6*mx*px - d1*mu2*mu3*mu4*mu5*mu6*mx*px - a6*lz*mu2*mu3*mu4*mu5*mu6*mx*px + a2*lambda3*lambda4*lambda5*mu2*mu6*my*px - a3*lambda3*lambda4*lambda5*mu2*mu6*my*px - a4*lambda3*lambda4*lambda5*mu2*mu6*my*px + a5*lambda3*lambda4*lambda5*mu2*mu6*my*px - a2*lambda5*mu2*mu3*mu4*mu6*my*px + a3*lambda5*mu2*mu3*mu4*mu6*my*px + a4*lambda5*mu2*mu3*mu4*mu6*my*px - a5*lambda5*mu2*mu3*mu4*mu6*my*px + a2*lambda4*mu2*mu3*mu5*mu6*my*px - a3*lambda4*mu2*mu3*mu5*mu6*my*px - a4*lambda4*mu2*mu3*mu5*mu6*my*px + a5*lambda4*mu2*mu3*mu5*mu6*my*px + a2*lambda3*mu2*mu4*mu5*mu6*my*px - a3*lambda3*mu2*mu4*mu5*mu6*my*px - a4*lambda3*mu2*mu4*mu5*mu6*my*px + a5*lambda3*mu2*mu4*mu5*mu6*my*px + 2*a6*lambda4*lambda5*lx*mu2*mu3*mu6*mz*px + 2*a6*lambda3*lambda5*lx*mu2*mu4*mu6*mz*px - 2*a6*lambda3*lambda4*lx*mu2*mu5*mu6*mz*px + 2*a6*lx*mu2*mu3*mu4*mu5*mu6*mz*px - a2*lambda3*lambda4*lambda5*mu2*mu6*mx*py + a3*lambda3*lambda4*lambda5*mu2*mu6*mx*py + a4*lambda3*lambda4*lambda5*mu2*mu6*mx*py - a5*lambda3*lambda4*lambda5*mu2*mu6*mx*py + a2*lambda5*mu2*mu3*mu4*mu6*mx*py - a3*lambda5*mu2*mu3*mu4*mu6*mx*py - a4*lambda5*mu2*mu3*mu4*mu6*mx*py + a5*lambda5*mu2*mu3*mu4*mu6*mx*py - a2*lambda4*mu2*mu3*mu5*mu6*mx*py + a3*lambda4*mu2*mu3*mu5*mu6*mx*py + a4*lambda4*mu2*mu3*mu5*mu6*mx*py - a5*lambda4*mu2*mu3*mu5*mu6*mx*py - a2*lambda3*mu2*mu4*mu5*mu6*mx*py + a3*lambda3*mu2*mu4*mu5*mu6*mx*py + a4*lambda3*mu2*mu4*mu5*mu6*mx*py - a5*lambda3*mu2*mu4*mu5*mu6*mx*py - d1*lambda4*lambda5*mu2*mu3*mu6*my*py - a6*lambda4*lambda5*lz*mu2*mu3*mu6*my*py - d1*lambda3*lambda5*mu2*mu4*mu6*my*py - a6*lambda3*lambda5*lz*mu2*mu4*mu6*my*py + d1*lambda3*lambda4*mu2*mu5*mu6*my*py + a6*lambda3*lambda4*lz*mu2*mu5*mu6*my*py - d1*mu2*mu3*mu4*mu5*mu6*my*py - a6*lz*mu2*mu3*mu4*mu5*mu6*my*py + 2*a6*lambda4*lambda5*ly*mu2*mu3*mu6*mz*py + 2*a6*lambda3*lambda5*ly*mu2*mu4*mu6*mz*py - 2*a6*lambda3*lambda4*ly*mu2*mu5*mu6*mz*py + 2*a6*ly*mu2*mu3*mu4*mu5*mu6*mz*py - a6*lambda4*lambda5*lx*mu2*mu3*mu6*mx*pz - a6*lambda3*lambda5*lx*mu2*mu4*mu6*mx*pz + a6*lambda3*lambda4*lx*mu2*mu5*mu6*mx*pz - a6*lx*mu2*mu3*mu4*mu5*mu6*mx*pz - a6*lambda4*lambda5*ly*mu2*mu3*mu6*my*pz - a6*lambda3*lambda5*ly*mu2*mu4*mu6*my*pz + a6*lambda3*lambda4*ly*mu2*mu5*mu6*my*pz - a6*ly*mu2*mu3*mu4*mu5*mu6*my*pz + lambda4*lambda5*mu2*mu3*mu6*mx*px*pz + lambda3*lambda5*mu2*mu4*mu6*mx*px*pz - lambda3*lambda4*mu2*mu5*mu6*mx*px*pz + mu2*mu3*mu4*mu5*mu6*mx*px*pz + lambda4*lambda5*mu2*mu3*mu6*my*py*pz + lambda3*lambda5*mu2*mu4*mu6*my*py*pz - lambda3*lambda4*mu2*mu5*mu6*my*py*pz + mu2*mu3*mu4*mu5*mu6*my*py*pz + lambda4*lambda5*lx*lz*mu2*mu3*mu6*mx*pow(a6,2) + lambda3*lambda5*lx*lz*mu2*mu4*mu6*mx*pow(a6,2) - lambda3*lambda4*lx*lz*mu2*mu5*mu6*mx*pow(a6,2) + lx*lz*mu2*mu3*mu4*mu5*mu6*mx*pow(a6,2) + lambda4*lambda5*ly*lz*mu2*mu3*mu6*my*pow(a6,2) + lambda3*lambda5*ly*lz*mu2*mu4*mu6*my*pow(a6,2) - lambda3*lambda4*ly*lz*mu2*mu5*mu6*my*pow(a6,2) + ly*lz*mu2*mu3*mu4*mu5*mu6*my*pow(a6,2) - lambda4*lambda5*mu2*mu3*mu6*mz*pow(a6,2)*pow(lx,2) - lambda3*lambda5*mu2*mu4*mu6*mz*pow(a6,2)*pow(lx,2) + lambda3*lambda4*mu2*mu5*mu6*mz*pow(a6,2)*pow(lx,2) - mu2*mu3*mu4*mu5*mu6*mz*pow(a6,2)*pow(lx,2) - lambda4*lambda5*mu2*mu3*mu6*mz*pow(a6,2)*pow(ly,2) - lambda3*lambda5*mu2*mu4*mu6*mz*pow(a6,2)*pow(ly,2) + lambda3*lambda4*mu2*mu5*mu6*mz*pow(a6,2)*pow(ly,2) - mu2*mu3*mu4*mu5*mu6*mz*pow(a6,2)*pow(ly,2) + a1*a6*lambda1*lambda4*lambda5*ly*mu2*mu3*mu6*mx*pow(mu1,-1) + a1*a6*lambda1*lambda3*lambda5*ly*mu2*mu4*mu6*mx*pow(mu1,-1) - a1*a6*lambda1*lambda3*lambda4*ly*mu2*mu5*mu6*mx*pow(mu1,-1) + a1*a6*lambda1*ly*mu2*mu3*mu4*mu5*mu6*mx*pow(mu1,-1) - a1*a6*lambda1*lambda4*lambda5*lx*mu2*mu3*mu6*my*pow(mu1,-1) - a1*a6*lambda1*lambda3*lambda5*lx*mu2*mu4*mu6*my*pow(mu1,-1) + a1*a6*lambda1*lambda3*lambda4*lx*mu2*mu5*mu6*my*pow(mu1,-1) - a1*a6*lambda1*lx*mu2*mu3*mu4*mu5*mu6*my*pow(mu1,-1) + a1*lambda1*lambda4*lambda5*mu2*mu3*mu6*my*px*pow(mu1,-1) + a1*lambda1*lambda3*lambda5*mu2*mu4*mu6*my*px*pow(mu1,-1) - a1*lambda1*lambda3*lambda4*mu2*mu5*mu6*my*px*pow(mu1,-1) + a1*lambda1*mu2*mu3*mu4*mu5*mu6*my*px*pow(mu1,-1) - a1*lambda1*lambda4*lambda5*mu2*mu3*mu6*mx*py*pow(mu1,-1) - a1*lambda1*lambda3*lambda5*mu2*mu4*mu6*mx*py*pow(mu1,-1) + a1*lambda1*lambda3*lambda4*mu2*mu5*mu6*mx*py*pow(mu1,-1) - a1*lambda1*mu2*mu3*mu4*mu5*mu6*mx*py*pow(mu1,-1) - a6*d2*lx*mx*mz*pow(mu6,2) - 2*a6*d1*lambda1*lx*mx*mz*pow(mu6,2) - a1*a6*ly*mu1*mx*mz*pow(mu6,2) - a6*d4*lx*mu2*mu3*mx*mz*pow(mu6,2) - a6*d5*lambda4*lx*mu2*mu3*mx*mz*pow(mu6,2) - a6*d6*lambda4*lambda5*lx*mu2*mu3*mx*mz*pow(mu6,2) - a6*d5*lambda3*lx*mu2*mu4*mx*mz*pow(mu6,2) - a6*d6*lambda3*lambda5*lx*mu2*mu4*mx*mz*pow(mu6,2) + a6*d6*lambda3*lambda4*lx*mu2*mu5*mx*mz*pow(mu6,2) - a6*d6*lx*mu2*mu3*mu4*mu5*mx*mz*pow(mu6,2) - a6*d2*ly*my*mz*pow(mu6,2) - 2*a6*d1*lambda1*ly*my*mz*pow(mu6,2) + a1*a6*lx*mu1*my*mz*pow(mu6,2) - a6*d4*ly*mu2*mu3*my*mz*pow(mu6,2) - a6*d5*lambda4*ly*mu2*mu3*my*mz*pow(mu6,2) - a6*d6*lambda4*lambda5*ly*mu2*mu3*my*mz*pow(mu6,2) - a6*d5*lambda3*ly*mu2*mu4*my*mz*pow(mu6,2) - a6*d6*lambda3*lambda5*ly*mu2*mu4*my*mz*pow(mu6,2) + a6*d6*lambda3*lambda4*ly*mu2*mu5*my*mz*pow(mu6,2) - a6*d6*ly*mu2*mu3*mu4*mu5*my*mz*pow(mu6,2) + 2*a6*lambda1*ly*mx*my*px*pow(mu6,2) + d2*mx*mz*px*pow(mu6,2) + 2*d1*lambda1*mx*mz*px*pow(mu6,2) + 2*a6*lambda1*lz*mx*mz*px*pow(mu6,2) + d4*mu2*mu3*mx*mz*px*pow(mu6,2) + d5*lambda4*mu2*mu3*mx*mz*px*pow(mu6,2) + d6*lambda4*lambda5*mu2*mu3*mx*mz*px*pow(mu6,2) + d5*lambda3*mu2*mu4*mx*mz*px*pow(mu6,2) + d6*lambda3*lambda5*mu2*mu4*mx*mz*px*pow(mu6,2) - d6*lambda3*lambda4*mu2*mu5*mx*mz*px*pow(mu6,2) + d6*mu2*mu3*mu4*mu5*mx*mz*px*pow(mu6,2) - a1*mu1*my*mz*px*pow(mu6,2) + 2*a6*lambda1*lx*mx*my*py*pow(mu6,2) + a1*mu1*mx*mz*py*pow(mu6,2) + d2*my*mz*py*pow(mu6,2) + 2*d1*lambda1*my*mz*py*pow(mu6,2) + 2*a6*lambda1*lz*my*mz*py*pow(mu6,2) + d4*mu2*mu3*my*mz*py*pow(mu6,2) + d5*lambda4*mu2*mu3*my*mz*py*pow(mu6,2) + d6*lambda4*lambda5*mu2*mu3*my*mz*py*pow(mu6,2) + d5*lambda3*mu2*mu4*my*mz*py*pow(mu6,2) + d6*lambda3*lambda5*mu2*mu4*my*mz*py*pow(mu6,2) - d6*lambda3*lambda4*mu2*mu5*my*mz*py*pow(mu6,2) + d6*mu2*mu3*mu4*mu5*my*mz*py*pow(mu6,2) - 2*lambda1*mx*my*px*py*pow(mu6,2) + 2*a6*lambda1*lx*mx*mz*pz*pow(mu6,2) + 2*a6*lambda1*ly*my*mz*pz*pow(mu6,2) - 2*lambda1*mx*mz*px*pz*pow(mu6,2) - 2*lambda1*my*mz*py*pz*pow(mu6,2) - 2*lambda1*lx*ly*mx*my*pow(a6,2)*pow(mu6,2) - 2*lambda1*lx*lz*mx*mz*pow(a6,2)*pow(mu6,2) - 2*lambda1*ly*lz*my*mz*pow(a6,2)*pow(mu6,2) - a1*a6*ly*mx*mz*pow(lambda1,2)*pow(mu1,-1)*pow(mu6,2) + a1*a6*lx*my*mz*pow(lambda1,2)*pow(mu1,-1)*pow(mu6,2) - a1*my*mz*px*pow(lambda1,2)*pow(mu1,-1)*pow(mu6,2) + a1*mx*mz*py*pow(lambda1,2)*pow(mu1,-1)*pow(mu6,2) + d1*d2*pow(mu6,2)*pow(mx,2) + a6*d2*lz*pow(mu6,2)*pow(mx,2) + 2*a6*d1*lambda1*lz*pow(mu6,2)*pow(mx,2) + d1*d4*mu2*mu3*pow(mu6,2)*pow(mx,2) + d1*d5*lambda4*mu2*mu3*pow(mu6,2)*pow(mx,2) + d1*d6*lambda4*lambda5*mu2*mu3*pow(mu6,2)*pow(mx,2) + a6*d4*lz*mu2*mu3*pow(mu6,2)*pow(mx,2) + a6*d5*lambda4*lz*mu2*mu3*pow(mu6,2)*pow(mx,2) + a6*d6*lambda4*lambda5*lz*mu2*mu3*pow(mu6,2)*pow(mx,2) + d1*d5*lambda3*mu2*mu4*pow(mu6,2)*pow(mx,2) + d1*d6*lambda3*lambda5*mu2*mu4*pow(mu6,2)*pow(mx,2) + a6*d5*lambda3*lz*mu2*mu4*pow(mu6,2)*pow(mx,2) + a6*d6*lambda3*lambda5*lz*mu2*mu4*pow(mu6,2)*pow(mx,2) - d1*d6*lambda3*lambda4*mu2*mu5*pow(mu6,2)*pow(mx,2) - a6*d6*lambda3*lambda4*lz*mu2*mu5*pow(mu6,2)*pow(mx,2) + d1*d6*mu2*mu3*mu4*mu5*pow(mu6,2)*pow(mx,2) + a6*d6*lz*mu2*mu3*mu4*mu5*pow(mu6,2)*pow(mx,2) - 2*a6*lambda1*ly*py*pow(mu6,2)*pow(mx,2) - d2*pz*pow(mu6,2)*pow(mx,2) - 2*d1*lambda1*pz*pow(mu6,2)*pow(mx,2) - 2*a6*lambda1*lz*pz*pow(mu6,2)*pow(mx,2) - d4*mu2*mu3*pz*pow(mu6,2)*pow(mx,2) - d5*lambda4*mu2*mu3*pz*pow(mu6,2)*pow(mx,2) - d6*lambda4*lambda5*mu2*mu3*pz*pow(mu6,2)*pow(mx,2) - d5*lambda3*mu2*mu4*pz*pow(mu6,2)*pow(mx,2) - d6*lambda3*lambda5*mu2*mu4*pz*pow(mu6,2)*pow(mx,2) + d6*lambda3*lambda4*mu2*mu5*pz*pow(mu6,2)*pow(mx,2) - d6*mu2*mu3*mu4*mu5*pz*pow(mu6,2)*pow(mx,2) + lambda1*pow(d1,2)*pow(mu6,2)*pow(mx,2) + lambda1*pow(a6,2)*pow(ly,2)*pow(mu6,2)*pow(mx,2) + lambda1*pow(a6,2)*pow(lz,2)*pow(mu6,2)*pow(mx,2) + d1*d2*pow(mu6,2)*pow(my,2) + a6*d2*lz*pow(mu6,2)*pow(my,2) + 2*a6*d1*lambda1*lz*pow(mu6,2)*pow(my,2) + d1*d4*mu2*mu3*pow(mu6,2)*pow(my,2) + d1*d5*lambda4*mu2*mu3*pow(mu6,2)*pow(my,2) + d1*d6*lambda4*lambda5*mu2*mu3*pow(mu6,2)*pow(my,2) + a6*d4*lz*mu2*mu3*pow(mu6,2)*pow(my,2) + a6*d5*lambda4*lz*mu2*mu3*pow(mu6,2)*pow(my,2) + a6*d6*lambda4*lambda5*lz*mu2*mu3*pow(mu6,2)*pow(my,2) + d1*d5*lambda3*mu2*mu4*pow(mu6,2)*pow(my,2) + d1*d6*lambda3*lambda5*mu2*mu4*pow(mu6,2)*pow(my,2) + a6*d5*lambda3*lz*mu2*mu4*pow(mu6,2)*pow(my,2) + a6*d6*lambda3*lambda5*lz*mu2*mu4*pow(mu6,2)*pow(my,2) - d1*d6*lambda3*lambda4*mu2*mu5*pow(mu6,2)*pow(my,2) - a6*d6*lambda3*lambda4*lz*mu2*mu5*pow(mu6,2)*pow(my,2) + d1*d6*mu2*mu3*mu4*mu5*pow(mu6,2)*pow(my,2) + a6*d6*lz*mu2*mu3*mu4*mu5*pow(mu6,2)*pow(my,2) - 2*a6*lambda1*lx*px*pow(mu6,2)*pow(my,2) - d2*pz*pow(mu6,2)*pow(my,2) - 2*d1*lambda1*pz*pow(mu6,2)*pow(my,2) - 2*a6*lambda1*lz*pz*pow(mu6,2)*pow(my,2) - d4*mu2*mu3*pz*pow(mu6,2)*pow(my,2) - d5*lambda4*mu2*mu3*pz*pow(mu6,2)*pow(my,2) - d6*lambda4*lambda5*mu2*mu3*pz*pow(mu6,2)*pow(my,2) - d5*lambda3*mu2*mu4*pz*pow(mu6,2)*pow(my,2) - d6*lambda3*lambda5*mu2*mu4*pz*pow(mu6,2)*pow(my,2) + d6*lambda3*lambda4*mu2*mu5*pz*pow(mu6,2)*pow(my,2) - d6*mu2*mu3*mu4*mu5*pz*pow(mu6,2)*pow(my,2) + lambda1*pow(d1,2)*pow(mu6,2)*pow(my,2) + lambda1*pow(a6,2)*pow(lx,2)*pow(mu6,2)*pow(my,2) + lambda1*pow(a6,2)*pow(lz,2)*pow(mu6,2)*pow(my,2) - 2*a6*lambda1*lx*px*pow(mu6,2)*pow(mz,2) - 2*a6*lambda1*ly*py*pow(mu6,2)*pow(mz,2) + lambda1*pow(a6,2)*pow(lx,2)*pow(mu6,2)*pow(mz,2) + lambda1*pow(a6,2)*pow(ly,2)*pow(mu6,2)*pow(mz,2) - lambda4*lambda5*mu2*mu3*mu6*mz*pow(px,2) - lambda3*lambda5*mu2*mu4*mu6*mz*pow(px,2) + lambda3*lambda4*mu2*mu5*mu6*mz*pow(px,2) - mu2*mu3*mu4*mu5*mu6*mz*pow(px,2) + lambda1*pow(mu6,2)*pow(my,2)*pow(px,2) + lambda1*pow(mu6,2)*pow(mz,2)*pow(px,2) - lambda4*lambda5*mu2*mu3*mu6*mz*pow(py,2) - lambda3*lambda5*mu2*mu4*mu6*mz*pow(py,2) + lambda3*lambda4*mu2*mu5*mu6*mz*pow(py,2) - mu2*mu3*mu4*mu5*mu6*mz*pow(py,2) + lambda1*pow(mu6,2)*pow(mx,2)*pow(py,2) + lambda1*pow(mu6,2)*pow(mz,2)*pow(py,2) + lambda6*(-(a2*a6*lambda5*ly*mu2*mu3*mu4*nx) + a3*a6*lambda5*ly*mu2*mu3*mu4*nx + a4*a6*lambda5*ly*mu2*mu3*mu4*nx - a5*a6*lambda5*ly*mu2*mu3*mu4*nx + a6*d1*lx*mu2*mu3*mu4*mu5*nx + 4*a6*d1*lambda1*lz*mu6*mx*nx - a6*d2*lx*mu6*mz*nx - 2*a6*d1*lambda1*lx*mu6*mz*nx - a1*a6*ly*mu1*mu6*mz*nx - a6*d4*lx*mu2*mu3*mu6*mz*nx - a6*d6*lx*mu2*mu3*mu4*mu5*mu6*mz*nx + a2*a6*lambda5*lx*mu2*mu3*mu4*ny - a3*a6*lambda5*lx*mu2*mu3*mu4*ny - a4*a6*lambda5*lx*mu2*mu3*mu4*ny + a5*a6*lambda5*lx*mu2*mu3*mu4*ny + a6*d1*ly*mu2*mu3*mu4*mu5*ny + 4*a6*d1*lambda1*lz*mu6*my*ny - a6*d2*ly*mu6*mz*ny - 2*a6*d1*lambda1*ly*mu6*mz*ny + a1*a6*lx*mu1*mu6*mz*ny - a6*d4*ly*mu2*mu3*mu6*mz*ny - a6*d6*ly*mu2*mu3*mu4*mu5*mu6*mz*ny - a6*d2*lx*mu6*mx*nz - 2*a6*d1*lambda1*lx*mu6*mx*nz - a1*a6*ly*mu1*mu6*mx*nz - a6*d4*lx*mu2*mu3*mu6*mx*nz - a6*d6*lx*mu2*mu3*mu4*mu5*mu6*mx*nz - a6*d2*ly*mu6*my*nz - 2*a6*d1*lambda1*ly*mu6*my*nz + a1*a6*lx*mu1*mu6*my*nz - a6*d4*ly*mu2*mu3*mu6*my*nz - a6*d6*ly*mu2*mu3*mu4*mu5*mu6*my*nz - d1*mu2*mu3*mu4*mu5*nx*px - a6*lz*mu2*mu3*mu4*mu5*nx*px + d2*mu6*mz*nx*px + d4*mu2*mu3*mu6*mz*nx*px + d6*mu2*mu3*mu4*mu5*mu6*mz*nx*px - a2*lambda5*mu2*mu3*mu4*ny*px + a3*lambda5*mu2*mu3*mu4*ny*px + a4*lambda5*mu2*mu3*mu4*ny*px - a5*lambda5*mu2*mu3*mu4*ny*px - 4*a6*lambda1*lx*mu6*my*ny*px - a1*mu1*mu6*mz*ny*px + 2*a6*lx*mu2*mu3*mu4*mu5*nz*px + d2*mu6*mx*nz*px + d4*mu2*mu3*mu6*mx*nz*px + d6*mu2*mu3*mu4*mu5*mu6*mx*nz*px - a1*mu1*mu6*my*nz*px - 4*a6*lambda1*lx*mu6*mz*nz*px + a2*lambda5*mu2*mu3*mu4*nx*py - a3*lambda5*mu2*mu3*mu4*nx*py - a4*lambda5*mu2*mu3*mu4*nx*py + a5*lambda5*mu2*mu3*mu4*nx*py - 4*a6*lambda1*ly*mu6*mx*nx*py + a1*mu1*mu6*mz*nx*py - d1*mu2*mu3*mu4*mu5*ny*py - a6*lz*mu2*mu3*mu4*mu5*ny*py + d2*mu6*mz*ny*py + d4*mu2*mu3*mu6*mz*ny*py + d6*mu2*mu3*mu4*mu5*mu6*mz*ny*py + 2*a6*ly*mu2*mu3*mu4*mu5*nz*py + a1*mu1*mu6*mx*nz*py + d2*mu6*my*nz*py + d4*mu2*mu3*mu6*my*nz*py + d6*mu2*mu3*mu4*mu5*mu6*my*nz*py - 4*a6*lambda1*ly*mu6*mz*nz*py - 2*lambda1*mu6*my*nx*px*py - 2*lambda1*mu6*mx*ny*px*py - a6*lx*mu2*mu3*mu4*mu5*nx*pz - 2*d2*mu6*mx*nx*pz - 4*d1*lambda1*mu6*mx*nx*pz - 4*a6*lambda1*lz*mu6*mx*nx*pz - 2*d4*mu2*mu3*mu6*mx*nx*pz - 2*d6*mu2*mu3*mu4*mu5*mu6*mx*nx*pz - a6*ly*mu2*mu3*mu4*mu5*ny*pz - 2*d2*mu6*my*ny*pz - 4*d1*lambda1*mu6*my*ny*pz - 4*a6*lambda1*lz*mu6*my*ny*pz - 2*d4*mu2*mu3*mu6*my*ny*pz - 2*d6*mu2*mu3*mu4*mu5*mu6*my*ny*pz + mu2*mu3*mu4*mu5*nx*px*pz - 2*lambda1*mu6*mz*nx*px*pz - 2*lambda1*mu6*mx*nz*px*pz + mu2*mu3*mu4*mu5*ny*py*pz - 2*lambda1*mu6*mz*ny*py*pz - 2*lambda1*mu6*my*nz*py*pz + lx*lz*mu2*mu3*mu4*mu5*nx*pow(a6,2) - 2*lambda1*lx*ly*mu6*my*nx*pow(a6,2) - 2*lambda1*lx*lz*mu6*mz*nx*pow(a6,2) + ly*lz*mu2*mu3*mu4*mu5*ny*pow(a6,2) - 2*lambda1*lx*ly*mu6*mx*ny*pow(a6,2) - 2*lambda1*ly*lz*mu6*mz*ny*pow(a6,2) - 2*lambda1*lx*lz*mu6*mx*nz*pow(a6,2) - 2*lambda1*ly*lz*mu6*my*nz*pow(a6,2) - mu2*mu3*mu4*mu5*nz*pow(a6,2)*pow(lx,2) - mu2*mu3*mu4*mu5*nz*pow(a6,2)*pow(ly,2) + a1*a6*lambda1*ly*mu2*mu3*mu4*mu5*nx*pow(mu1,-1) - a1*a6*lambda1*lx*mu2*mu3*mu4*mu5*ny*pow(mu1,-1) + a1*lambda1*mu2*mu3*mu4*mu5*ny*px*pow(mu1,-1) - a1*lambda1*mu2*mu3*mu4*mu5*nx*py*pow(mu1,-1) - a1*a6*ly*mu6*mz*nx*pow(lambda1,2)*pow(mu1,-1) + a1*a6*lx*mu6*mz*ny*pow(lambda1,2)*pow(mu1,-1) - a1*a6*ly*mu6*mx*nz*pow(lambda1,2)*pow(mu1,-1) + a1*a6*lx*mu6*my*nz*pow(lambda1,2)*pow(mu1,-1) - a1*mu6*mz*ny*px*pow(lambda1,2)*pow(mu1,-1) - a1*mu6*my*nz*px*pow(lambda1,2)*pow(mu1,-1) + a1*mu6*mz*nx*py*pow(lambda1,2)*pow(mu1,-1) + a1*mu6*mx*nz*py*pow(lambda1,2)*pow(mu1,-1) - mu2*mu3*mu4*mu5*nz*pow(px,2) - mu2*mu3*mu4*mu5*nz*pow(py,2) - lambda4*mu2*mu3*(-((a2 - a3 - a4 + a5)*mu5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)) + d5*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz) + lambda5*(-(a6*d1*lx*nx) - a6*d1*ly*ny + d1*nx*px + a6*lz*nx*px - 2*a6*lx*nz*px + d1*ny*py + a6*lz*ny*py - 2*a6*ly*nz*py + a6*lx*nx*pz + a6*ly*ny*pz - nx*px*pz - ny*py*pz + d6*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz) - lx*lz*nx*pow(a6,2) - ly*lz*ny*pow(a6,2) + nz*pow(a6,2)*pow(lx,2) + nz*pow(a6,2)*pow(ly,2) - a1*lambda1*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*pow(mu1,-1) + nz*pow(px,2) + nz*pow(py,2))) + lambda3*mu2*(lambda4*((a2 - a3 - a4 + a5)*lambda5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) + mu5*(-(a6*d1*lx*nx) - a6*d1*ly*ny + d1*nx*px + a6*lz*nx*px - 2*a6*lx*nz*px + d1*ny*py + a6*lz*ny*py - 2*a6*ly*nz*py - d6*mu6*(2*d1*(mx*nx + my*ny) - a6*(lx*mz*nx + ly*mz*ny - 2*lz*(mx*nx + my*ny) + lx*mx*nz + ly*my*nz) + mz*nx*px + mx*nz*px + mz*ny*py + my*nz*py) + a6*lx*nx*pz + a6*ly*ny*pz - nx*px*pz - ny*py*pz - lx*lz*nx*pow(a6,2) - ly*lz*ny*pow(a6,2) + nz*pow(a6,2)*pow(lx,2) + nz*pow(a6,2)*pow(ly,2) - a1*lambda1*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*pow(mu1,-1) + nz*pow(px,2) + nz*pow(py,2))) - mu4*(-((a2 - a3 - a4 + a5)*mu5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)) + d5*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz) + lambda5*(-(a6*d1*lx*nx) - a6*d1*ly*ny + d1*nx*px + a6*lz*nx*px - 2*a6*lx*nz*px + d1*ny*py + a6*lz*ny*py - 2*a6*ly*nz*py + a6*lx*nx*pz + a6*ly*ny*pz - nx*px*pz - ny*py*pz + d6*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz) - lx*lz*nx*pow(a6,2) - ly*lz*ny*pow(a6,2) + nz*pow(a6,2)*pow(lx,2) + nz*pow(a6,2)*pow(ly,2) - a1*lambda1*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*pow(mu1,-1) + nz*pow(px,2) + nz*pow(py,2))))) + lambda1*pow(mu6,2)*pow(mx,2)*pow(pz,2) + lambda1*pow(mu6,2)*pow(my,2)*pow(pz,2) + (pow(lambda6,2)*(-2*a6*d2*lx*nx*nz - 2*a1*a6*ly*mu1*nx*nz - 2*a6*d2*ly*ny*nz + 2*a1*a6*lx*mu1*ny*nz + 2*d2*nx*nz*px - 2*a1*mu1*ny*nz*px + 2*a1*mu1*nx*nz*py + 2*d2*ny*nz*py + a6*d4*lx*nx*nz*cos(alpha2 + alpha3) + a6*d4*ly*ny*nz*cos(alpha2 + alpha3) - d4*nx*nz*px*cos(alpha2 + alpha3) - d4*ny*nz*py*cos(alpha2 + alpha3) - a6*d5*lx*nx*nz*cos(alpha2 - alpha3 - alpha4) - a6*d5*ly*ny*nz*cos(alpha2 - alpha3 - alpha4) + d5*nx*nz*px*cos(alpha2 - alpha3 - alpha4) + d5*ny*nz*py*cos(alpha2 - alpha3 - alpha4) + a6*d5*lx*nx*nz*cos(alpha2 + alpha3 + alpha4) + a6*d5*ly*ny*nz*cos(alpha2 + alpha3 + alpha4) - d5*nx*nz*px*cos(alpha2 + alpha3 + alpha4) - d5*ny*nz*py*cos(alpha2 + alpha3 + alpha4) + a6*d6*lx*nx*nz*cos(alpha2 + alpha3 + alpha4 - alpha5) + a6*d6*ly*ny*nz*cos(alpha2 + alpha3 + alpha4 - alpha5) - d6*nx*nz*px*cos(alpha2 + alpha3 + alpha4 - alpha5) - d6*ny*nz*py*cos(alpha2 + alpha3 + alpha4 - alpha5) - a6*d6*lx*nx*nz*cos(alpha2 - alpha3 - alpha4 + alpha5) - a6*d6*ly*ny*nz*cos(alpha2 - alpha3 - alpha4 + alpha5) + d6*nx*nz*px*cos(alpha2 - alpha3 - alpha4 + alpha5) + d6*ny*nz*py*cos(alpha2 - alpha3 - alpha4 + alpha5) + 2*d1*d2*pow(nx,2) + 2*a6*d2*lz*pow(nx,2) - 2*d2*pz*pow(nx,2) - d1*d4*cos(alpha2 + alpha3)*pow(nx,2) - a6*d4*lz*cos(alpha2 + alpha3)*pow(nx,2) + d4*pz*cos(alpha2 + alpha3)*pow(nx,2) + d1*d5*cos(alpha2 - alpha3 - alpha4)*pow(nx,2) + a6*d5*lz*cos(alpha2 - alpha3 - alpha4)*pow(nx,2) - d5*pz*cos(alpha2 - alpha3 - alpha4)*pow(nx,2) - d1*d5*cos(alpha2 + alpha3 + alpha4)*pow(nx,2) - a6*d5*lz*cos(alpha2 + alpha3 + alpha4)*pow(nx,2) + d5*pz*cos(alpha2 + alpha3 + alpha4)*pow(nx,2) - d1*d6*cos(alpha2 + alpha3 + alpha4 - alpha5)*pow(nx,2) - a6*d6*lz*cos(alpha2 + alpha3 + alpha4 - alpha5)*pow(nx,2) + d6*pz*cos(alpha2 + alpha3 + alpha4 - alpha5)*pow(nx,2) + d1*d6*cos(alpha2 - alpha3 - alpha4 + alpha5)*pow(nx,2) + a6*d6*lz*cos(alpha2 - alpha3 - alpha4 + alpha5)*pow(nx,2) - d6*pz*cos(alpha2 - alpha3 - alpha4 + alpha5)*pow(nx,2) + 2*d1*d2*pow(ny,2) + 2*a6*d2*lz*pow(ny,2) - 2*d2*pz*pow(ny,2) - d1*d4*cos(alpha2 + alpha3)*pow(ny,2) - a6*d4*lz*cos(alpha2 + alpha3)*pow(ny,2) + d4*pz*cos(alpha2 + alpha3)*pow(ny,2) + d1*d5*cos(alpha2 - alpha3 - alpha4)*pow(ny,2) + a6*d5*lz*cos(alpha2 - alpha3 - alpha4)*pow(ny,2) - d5*pz*cos(alpha2 - alpha3 - alpha4)*pow(ny,2) - d1*d5*cos(alpha2 + alpha3 + alpha4)*pow(ny,2) - a6*d5*lz*cos(alpha2 + alpha3 + alpha4)*pow(ny,2) + d5*pz*cos(alpha2 + alpha3 + alpha4)*pow(ny,2) - d1*d6*cos(alpha2 + alpha3 + alpha4 - alpha5)*pow(ny,2) - a6*d6*lz*cos(alpha2 + alpha3 + alpha4 - alpha5)*pow(ny,2) + d6*pz*cos(alpha2 + alpha3 + alpha4 - alpha5)*pow(ny,2) + d1*d6*cos(alpha2 - alpha3 - alpha4 + alpha5)*pow(ny,2) + a6*d6*lz*cos(alpha2 - alpha3 - alpha4 + alpha5)*pow(ny,2) - d6*pz*cos(alpha2 - alpha3 - alpha4 + alpha5)*pow(ny,2) + d4*cos(alpha2 - alpha3)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))) + 2*lambda1*(2*a6*ly*nx*ny*px + 2*a6*lz*nx*nz*px + 2*a6*lx*nx*ny*py + 2*a6*lz*ny*nz*py - 2*nx*ny*px*py + 2*a6*lx*nx*nz*pz + 2*a6*ly*ny*nz*pz - 2*nx*nz*px*pz - 2*ny*nz*py*pz - 2*lx*ly*nx*ny*pow(a6,2) - 2*lx*lz*nx*nz*pow(a6,2) - 2*ly*lz*ny*nz*pow(a6,2) - a1*lambda1*nz*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*pow(mu1,-1) - 2*a6*ly*py*pow(nx,2) - 2*a6*lz*pz*pow(nx,2) + pow(a6,2)*pow(ly,2)*pow(nx,2) + pow(a6,2)*pow(lz,2)*pow(nx,2) - 2*a6*lx*px*pow(ny,2) - 2*a6*lz*pz*pow(ny,2) + pow(a6,2)*pow(lx,2)*pow(ny,2) + pow(a6,2)*pow(lz,2)*pow(ny,2) + pow(d1,2)*(pow(nx,2) + pow(ny,2)) + 2*d1*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))) - 2*a6*lx*px*pow(nz,2) - 2*a6*ly*py*pow(nz,2) + pow(a6,2)*pow(lx,2)*pow(nz,2) + pow(a6,2)*pow(ly,2)*pow(nz,2) + pow(ny,2)*pow(px,2) + pow(nz,2)*pow(px,2) + pow(nx,2)*pow(py,2) + pow(nz,2)*pow(py,2) + pow(nx,2)*pow(pz,2) + pow(ny,2)*pow(pz,2))))/2. + d1*d2*mx*nx*sin(2*alpha6) + a6*d2*lz*mx*nx*sin(2*alpha6) + d1*d4*mu2*mu3*mx*nx*sin(2*alpha6) + d1*d5*lambda4*mu2*mu3*mx*nx*sin(2*alpha6) + d1*d6*lambda4*lambda5*mu2*mu3*mx*nx*sin(2*alpha6) + a6*d4*lz*mu2*mu3*mx*nx*sin(2*alpha6) + a6*d5*lambda4*lz*mu2*mu3*mx*nx*sin(2*alpha6) + a6*d6*lambda4*lambda5*lz*mu2*mu3*mx*nx*sin(2*alpha6) + d1*d5*lambda3*mu2*mu4*mx*nx*sin(2*alpha6) + d1*d6*lambda3*lambda5*mu2*mu4*mx*nx*sin(2*alpha6) + a6*d5*lambda3*lz*mu2*mu4*mx*nx*sin(2*alpha6) + a6*d6*lambda3*lambda5*lz*mu2*mu4*mx*nx*sin(2*alpha6) + d1*d6*mu2*mu3*mu4*mu5*mx*nx*sin(2*alpha6) + a6*d6*lz*mu2*mu3*mu4*mu5*mx*nx*sin(2*alpha6) + d1*d2*my*ny*sin(2*alpha6) + a6*d2*lz*my*ny*sin(2*alpha6) + d1*d4*mu2*mu3*my*ny*sin(2*alpha6) + d1*d5*lambda4*mu2*mu3*my*ny*sin(2*alpha6) + d1*d6*lambda4*lambda5*mu2*mu3*my*ny*sin(2*alpha6) + a6*d4*lz*mu2*mu3*my*ny*sin(2*alpha6) + a6*d5*lambda4*lz*mu2*mu3*my*ny*sin(2*alpha6) + a6*d6*lambda4*lambda5*lz*mu2*mu3*my*ny*sin(2*alpha6) + d1*d5*lambda3*mu2*mu4*my*ny*sin(2*alpha6) + d1*d6*lambda3*lambda5*mu2*mu4*my*ny*sin(2*alpha6) + a6*d5*lambda3*lz*mu2*mu4*my*ny*sin(2*alpha6) + a6*d6*lambda3*lambda5*lz*mu2*mu4*my*ny*sin(2*alpha6) + d1*d6*mu2*mu3*mu4*mu5*my*ny*sin(2*alpha6) + a6*d6*lz*mu2*mu3*mu4*mu5*my*ny*sin(2*alpha6) + a6*lambda1*ly*my*nx*px*sin(2*alpha6) + d1*lambda1*mz*nx*px*sin(2*alpha6) + a6*lambda1*lz*mz*nx*px*sin(2*alpha6) + a6*lambda1*ly*mx*ny*px*sin(2*alpha6) + d1*lambda1*mx*nz*px*sin(2*alpha6) + a6*lambda1*lz*mx*nz*px*sin(2*alpha6) + a6*lambda1*lx*my*nx*py*sin(2*alpha6) + a6*lambda1*lx*mx*ny*py*sin(2*alpha6) + d1*lambda1*mz*ny*py*sin(2*alpha6) + a6*lambda1*lz*mz*ny*py*sin(2*alpha6) + d1*lambda1*my*nz*py*sin(2*alpha6) + a6*lambda1*lz*my*nz*py*sin(2*alpha6) + d6*lambda3*lambda4*mu2*mu5*mx*nx*pz*sin(2*alpha6) + a6*lambda1*lx*mz*nx*pz*sin(2*alpha6) + d6*lambda3*lambda4*mu2*mu5*my*ny*pz*sin(2*alpha6) + a6*lambda1*ly*mz*ny*pz*sin(2*alpha6) + a6*lambda1*lx*mx*nz*pz*sin(2*alpha6) + a6*lambda1*ly*my*nz*pz*sin(2*alpha6) + lambda1*mx*nx*pow(d1,2)*sin(2*alpha6) + lambda1*my*ny*pow(d1,2)*sin(2*alpha6) + lambda1*my*ny*pow(a6,2)*pow(lx,2)*sin(2*alpha6) + lambda1*mz*nz*pow(a6,2)*pow(lx,2)*sin(2*alpha6) + lambda1*mx*nx*pow(a6,2)*pow(ly,2)*sin(2*alpha6) + lambda1*mz*nz*pow(a6,2)*pow(ly,2)*sin(2*alpha6) + lambda1*mx*nx*pow(a6,2)*pow(lz,2)*sin(2*alpha6) + lambda1*my*ny*pow(a6,2)*pow(lz,2)*sin(2*alpha6) + lambda1*my*ny*pow(px,2)*sin(2*alpha6) + lambda1*mz*nz*pow(px,2)*sin(2*alpha6) + lambda1*mx*nx*pow(py,2)*sin(2*alpha6) + lambda1*mz*nz*pow(py,2)*sin(2*alpha6) + lambda1*mx*nx*pow(pz,2)*sin(2*alpha6) + lambda1*my*ny*pow(pz,2)*sin(2*alpha6) + lambda2*(-(a2*a6*lambda4*lambda5*ly*mu3*mu6*mx) + a3*a6*lambda4*lambda5*ly*mu3*mu6*mx + a4*a6*lambda4*lambda5*ly*mu3*mu6*mx - a5*a6*lambda4*lambda5*ly*mu3*mu6*mx - a6*d1*lambda5*lx*mu3*mu4*mu6*mx + a6*d1*lambda4*lx*mu3*mu5*mu6*mx - a2*a6*ly*mu3*mu4*mu5*mu6*mx + a3*a6*ly*mu3*mu4*mu5*mu6*mx + a4*a6*ly*mu3*mu4*mu5*mu6*mx - a5*a6*ly*mu3*mu4*mu5*mu6*mx + a2*a6*lambda4*lambda5*lx*mu3*mu6*my - a3*a6*lambda4*lambda5*lx*mu3*mu6*my - a4*a6*lambda4*lambda5*lx*mu3*mu6*my + a5*a6*lambda4*lambda5*lx*mu3*mu6*my - a6*d1*lambda5*ly*mu3*mu4*mu6*my + a6*d1*lambda4*ly*mu3*mu5*mu6*my + a2*a6*lx*mu3*mu4*mu5*mu6*my - a3*a6*lx*mu3*mu4*mu5*mu6*my - a4*a6*lx*mu3*mu4*mu5*mu6*my + a5*a6*lx*mu3*mu4*mu5*mu6*my + d1*lambda5*mu3*mu4*mu6*mx*px + a6*lambda5*lz*mu3*mu4*mu6*mx*px - d1*lambda4*mu3*mu5*mu6*mx*px - a6*lambda4*lz*mu3*mu5*mu6*mx*px - a2*lambda4*lambda5*mu3*mu6*my*px + a3*lambda4*lambda5*mu3*mu6*my*px + a4*lambda4*lambda5*mu3*mu6*my*px - a5*lambda4*lambda5*mu3*mu6*my*px - a2*mu3*mu4*mu5*mu6*my*px + a3*mu3*mu4*mu5*mu6*my*px + a4*mu3*mu4*mu5*mu6*my*px - a5*mu3*mu4*mu5*mu6*my*px - 2*a6*lambda5*lx*mu3*mu4*mu6*mz*px + 2*a6*lambda4*lx*mu3*mu5*mu6*mz*px + a2*lambda4*lambda5*mu3*mu6*mx*py - a3*lambda4*lambda5*mu3*mu6*mx*py - a4*lambda4*lambda5*mu3*mu6*mx*py + a5*lambda4*lambda5*mu3*mu6*mx*py + a2*mu3*mu4*mu5*mu6*mx*py - a3*mu3*mu4*mu5*mu6*mx*py - a4*mu3*mu4*mu5*mu6*mx*py + a5*mu3*mu4*mu5*mu6*mx*py + d1*lambda5*mu3*mu4*mu6*my*py + a6*lambda5*lz*mu3*mu4*mu6*my*py - d1*lambda4*mu3*mu5*mu6*my*py - a6*lambda4*lz*mu3*mu5*mu6*my*py - 2*a6*lambda5*ly*mu3*mu4*mu6*mz*py + 2*a6*lambda4*ly*mu3*mu5*mu6*mz*py + a6*lambda5*lx*mu3*mu4*mu6*mx*pz - a6*lambda4*lx*mu3*mu5*mu6*mx*pz + a6*lambda5*ly*mu3*mu4*mu6*my*pz - a6*lambda4*ly*mu3*mu5*mu6*my*pz - lambda5*mu3*mu4*mu6*mx*px*pz + lambda4*mu3*mu5*mu6*mx*px*pz - lambda5*mu3*mu4*mu6*my*py*pz + lambda4*mu3*mu5*mu6*my*py*pz - lambda5*lx*lz*mu3*mu4*mu6*mx*pow(a6,2) + lambda4*lx*lz*mu3*mu5*mu6*mx*pow(a6,2) - lambda5*ly*lz*mu3*mu4*mu6*my*pow(a6,2) + lambda4*ly*lz*mu3*mu5*mu6*my*pow(a6,2) + lambda5*mu3*mu4*mu6*mz*pow(a6,2)*pow(lx,2) - lambda4*mu3*mu5*mu6*mz*pow(a6,2)*pow(lx,2) + lambda5*mu3*mu4*mu6*mz*pow(a6,2)*pow(ly,2) - lambda4*mu3*mu5*mu6*mz*pow(a6,2)*pow(ly,2) - a1*a6*lambda1*lambda5*ly*mu3*mu4*mu6*mx*pow(mu1,-1) + a1*a6*lambda1*lambda4*ly*mu3*mu5*mu6*mx*pow(mu1,-1) + a1*a6*lambda1*lambda5*lx*mu3*mu4*mu6*my*pow(mu1,-1) - a1*a6*lambda1*lambda4*lx*mu3*mu5*mu6*my*pow(mu1,-1) - a1*lambda1*lambda5*mu3*mu4*mu6*my*px*pow(mu1,-1) + a1*lambda1*lambda4*mu3*mu5*mu6*my*px*pow(mu1,-1) + a1*lambda1*lambda5*mu3*mu4*mu6*mx*py*pow(mu1,-1) - a1*lambda1*lambda4*mu3*mu5*mu6*mx*py*pow(mu1,-1) - a6*d3*lx*mx*mz*pow(mu6,2) + a6*d5*lx*mu3*mu4*mx*mz*pow(mu6,2) + a6*d6*lambda5*lx*mu3*mu4*mx*mz*pow(mu6,2) - a6*d6*lambda4*lx*mu3*mu5*mx*mz*pow(mu6,2) - a6*d3*ly*my*mz*pow(mu6,2) + a6*d5*ly*mu3*mu4*my*mz*pow(mu6,2) + a6*d6*lambda5*ly*mu3*mu4*my*mz*pow(mu6,2) - a6*d6*lambda4*ly*mu3*mu5*my*mz*pow(mu6,2) + d3*mx*mz*px*pow(mu6,2) - d5*mu3*mu4*mx*mz*px*pow(mu6,2) - d6*lambda5*mu3*mu4*mx*mz*px*pow(mu6,2) + d6*lambda4*mu3*mu5*mx*mz*px*pow(mu6,2) + d3*my*mz*py*pow(mu6,2) - d5*mu3*mu4*my*mz*py*pow(mu6,2) - d6*lambda5*mu3*mu4*my*mz*py*pow(mu6,2) + d6*lambda4*mu3*mu5*my*mz*py*pow(mu6,2) + d1*d3*pow(mu6,2)*pow(mx,2) + a6*d3*lz*pow(mu6,2)*pow(mx,2) - d1*d5*mu3*mu4*pow(mu6,2)*pow(mx,2) - d1*d6*lambda5*mu3*mu4*pow(mu6,2)*pow(mx,2) - a6*d5*lz*mu3*mu4*pow(mu6,2)*pow(mx,2) - a6*d6*lambda5*lz*mu3*mu4*pow(mu6,2)*pow(mx,2) + d1*d6*lambda4*mu3*mu5*pow(mu6,2)*pow(mx,2) + a6*d6*lambda4*lz*mu3*mu5*pow(mu6,2)*pow(mx,2) - d3*pz*pow(mu6,2)*pow(mx,2) + d5*mu3*mu4*pz*pow(mu6,2)*pow(mx,2) + d6*lambda5*mu3*mu4*pz*pow(mu6,2)*pow(mx,2) - d6*lambda4*mu3*mu5*pz*pow(mu6,2)*pow(mx,2) + d1*d3*pow(mu6,2)*pow(my,2) + a6*d3*lz*pow(mu6,2)*pow(my,2) - d1*d5*mu3*mu4*pow(mu6,2)*pow(my,2) - d1*d6*lambda5*mu3*mu4*pow(mu6,2)*pow(my,2) - a6*d5*lz*mu3*mu4*pow(mu6,2)*pow(my,2) - a6*d6*lambda5*lz*mu3*mu4*pow(mu6,2)*pow(my,2) + d1*d6*lambda4*mu3*mu5*pow(mu6,2)*pow(my,2) + a6*d6*lambda4*lz*mu3*mu5*pow(mu6,2)*pow(my,2) - d3*pz*pow(mu6,2)*pow(my,2) + d5*mu3*mu4*pz*pow(mu6,2)*pow(my,2) + d6*lambda5*mu3*mu4*pz*pow(mu6,2)*pow(my,2) - d6*lambda4*mu3*mu5*pz*pow(mu6,2)*pow(my,2) + lambda5*mu3*mu4*mu6*mz*pow(px,2) - lambda4*mu3*mu5*mu6*mz*pow(px,2) + lambda5*mu3*mu4*mu6*mz*pow(py,2) - lambda4*mu3*mu5*mu6*mz*pow(py,2) - lambda6*(a2*a6*ly*mu3*mu4*mu5*nx - a3*a6*ly*mu3*mu4*mu5*nx - a4*a6*ly*mu3*mu4*mu5*nx + a5*a6*ly*mu3*mu4*mu5*nx + 2*d1*d5*mu3*mu4*mu6*mx*nx + 2*a6*d5*lz*mu3*mu4*mu6*mx*nx + a6*d3*lx*mu6*mz*nx - a6*d5*lx*mu3*mu4*mu6*mz*nx - a2*a6*lx*mu3*mu4*mu5*ny + a3*a6*lx*mu3*mu4*mu5*ny + a4*a6*lx*mu3*mu4*mu5*ny - a5*a6*lx*mu3*mu4*mu5*ny + 2*d1*d5*mu3*mu4*mu6*my*ny + 2*a6*d5*lz*mu3*mu4*mu6*my*ny + a6*d3*ly*mu6*mz*ny - a6*d5*ly*mu3*mu4*mu6*mz*ny + a6*d3*lx*mu6*mx*nz - a6*d5*lx*mu3*mu4*mu6*mx*nz + a6*d3*ly*mu6*my*nz - a6*d5*ly*mu3*mu4*mu6*my*nz - d3*mu6*mz*nx*px + d5*mu3*mu4*mu6*mz*nx*px + a2*mu3*mu4*mu5*ny*px - a3*mu3*mu4*mu5*ny*px - a4*mu3*mu4*mu5*ny*px + a5*mu3*mu4*mu5*ny*px - d3*mu6*mx*nz*px + d5*mu3*mu4*mu6*mx*nz*px - a2*mu3*mu4*mu5*nx*py + a3*mu3*mu4*mu5*nx*py + a4*mu3*mu4*mu5*nx*py - a5*mu3*mu4*mu5*nx*py - d3*mu6*mz*ny*py + d5*mu3*mu4*mu6*mz*ny*py - d3*mu6*my*nz*py + d5*mu3*mu4*mu6*my*nz*py + 2*d3*mu6*mx*nx*pz + 2*d3*mu6*my*ny*pz - lambda5*mu3*mu4*(-(a6*d1*lx*nx) - a6*d1*ly*ny + d1*nx*px + a6*lz*nx*px - 2*a6*lx*nz*px + d1*ny*py + a6*lz*ny*py - 2*a6*ly*nz*py - d6*mu6*(2*d1*(mx*nx + my*ny) - a6*(lx*mz*nx + ly*mz*ny - 2*lz*(mx*nx + my*ny) + lx*mx*nz + ly*my*nz) + mz*nx*px + mx*nz*px + mz*ny*py + my*nz*py) + a6*lx*nx*pz + a6*ly*ny*pz - nx*px*pz - ny*py*pz - lx*lz*nx*pow(a6,2) - ly*lz*ny*pow(a6,2) + nz*pow(a6,2)*pow(lx,2) + nz*pow(a6,2)*pow(ly,2) - a1*lambda1*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*pow(mu1,-1) + nz*pow(px,2) + nz*pow(py,2)) + lambda4*mu3*((a2 - a3 - a4 + a5)*lambda5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) + mu5*(-(a6*d1*lx*nx) - a6*d1*ly*ny + d1*nx*px + a6*lz*nx*px - 2*a6*lx*nz*px + d1*ny*py + a6*lz*ny*py - 2*a6*ly*nz*py + a6*lx*nx*pz + a6*ly*ny*pz - nx*px*pz - ny*py*pz + d6*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz) - lx*lz*nx*pow(a6,2) - ly*lz*ny*pow(a6,2) + nz*pow(a6,2)*pow(lx,2) + nz*pow(a6,2)*pow(ly,2) - a1*lambda1*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*pow(mu1,-1) + nz*pow(px,2) + nz*pow(py,2)))) + pow(lambda6,2)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2))))*(d3 - mu3*(d5*mu4 + d6*sin(alpha4 - alpha5))) + d1*d3*mx*nx*sin(2*alpha6) + a6*d3*lz*mx*nx*sin(2*alpha6) + d1*d6*lambda4*mu3*mu5*mx*nx*sin(2*alpha6) + a6*d6*lambda4*lz*mu3*mu5*mx*nx*sin(2*alpha6) + d1*d3*my*ny*sin(2*alpha6) + a6*d3*lz*my*ny*sin(2*alpha6) + d1*d6*lambda4*mu3*mu5*my*ny*sin(2*alpha6) + a6*d6*lambda4*lz*mu3*mu5*my*ny*sin(2*alpha6) + d5*mu3*mu4*mx*nx*pz*sin(2*alpha6) + d6*lambda5*mu3*mu4*mx*nx*pz*sin(2*alpha6) + d5*mu3*mu4*my*ny*pz*sin(2*alpha6) + d6*lambda5*mu3*mu4*my*ny*pz*sin(2*alpha6) + lambda3*(-(a2*a6*lambda5*ly*mu4*mu6*mx) + a3*a6*lambda5*ly*mu4*mu6*mx + a4*a6*lambda5*ly*mu4*mu6*mx - a5*a6*lambda5*ly*mu4*mu6*mx + a6*d1*lx*mu4*mu5*mu6*mx + a2*a6*lambda5*lx*mu4*mu6*my - a3*a6*lambda5*lx*mu4*mu6*my - a4*a6*lambda5*lx*mu4*mu6*my + a5*a6*lambda5*lx*mu4*mu6*my + a6*d1*ly*mu4*mu5*mu6*my - a2*a6*lambda5*lambda6*ly*mu4*nx + a3*a6*lambda5*lambda6*ly*mu4*nx + a4*a6*lambda5*lambda6*ly*mu4*nx - a5*a6*lambda5*lambda6*ly*mu4*nx + a6*d1*lambda6*lx*mu4*mu5*nx - a6*d4*lambda6*lx*mu6*mz*nx - a6*d6*lambda6*lx*mu4*mu5*mu6*mz*nx + a2*a6*lambda5*lambda6*lx*mu4*ny - a3*a6*lambda5*lambda6*lx*mu4*ny - a4*a6*lambda5*lambda6*lx*mu4*ny + a5*a6*lambda5*lambda6*lx*mu4*ny + a6*d1*lambda6*ly*mu4*mu5*ny - a6*d4*lambda6*ly*mu6*mz*ny - a6*d6*lambda6*ly*mu4*mu5*mu6*mz*ny - a6*d4*lambda6*lx*mu6*mx*nz - a6*d6*lambda6*lx*mu4*mu5*mu6*mx*nz - a6*d4*lambda6*ly*mu6*my*nz - a6*d6*lambda6*ly*mu4*mu5*mu6*my*nz - d1*mu4*mu5*mu6*mx*px - a6*lz*mu4*mu5*mu6*mx*px - a2*lambda5*mu4*mu6*my*px + a3*lambda5*mu4*mu6*my*px + a4*lambda5*mu4*mu6*my*px - a5*lambda5*mu4*mu6*my*px + 2*a6*lx*mu4*mu5*mu6*mz*px - d1*lambda6*mu4*mu5*nx*px - a6*lambda6*lz*mu4*mu5*nx*px + d4*lambda6*mu6*mz*nx*px + d6*lambda6*mu4*mu5*mu6*mz*nx*px - a2*lambda5*lambda6*mu4*ny*px + a3*lambda5*lambda6*mu4*ny*px + a4*lambda5*lambda6*mu4*ny*px - a5*lambda5*lambda6*mu4*ny*px + 2*a6*lambda6*lx*mu4*mu5*nz*px + d4*lambda6*mu6*mx*nz*px + d6*lambda6*mu4*mu5*mu6*mx*nz*px + a2*lambda5*mu4*mu6*mx*py - a3*lambda5*mu4*mu6*mx*py - a4*lambda5*mu4*mu6*mx*py + a5*lambda5*mu4*mu6*mx*py - d1*mu4*mu5*mu6*my*py - a6*lz*mu4*mu5*mu6*my*py + 2*a6*ly*mu4*mu5*mu6*mz*py + a2*lambda5*lambda6*mu4*nx*py - a3*lambda5*lambda6*mu4*nx*py - a4*lambda5*lambda6*mu4*nx*py + a5*lambda5*lambda6*mu4*nx*py - d1*lambda6*mu4*mu5*ny*py - a6*lambda6*lz*mu4*mu5*ny*py + d4*lambda6*mu6*mz*ny*py + d6*lambda6*mu4*mu5*mu6*mz*ny*py + 2*a6*lambda6*ly*mu4*mu5*nz*py + d4*lambda6*mu6*my*nz*py + d6*lambda6*mu4*mu5*mu6*my*nz*py - a6*lx*mu4*mu5*mu6*mx*pz - a6*ly*mu4*mu5*mu6*my*pz - a6*lambda6*lx*mu4*mu5*nx*pz - 2*d4*lambda6*mu6*mx*nx*pz - 2*d6*lambda6*mu4*mu5*mu6*mx*nx*pz - a6*lambda6*ly*mu4*mu5*ny*pz - 2*d4*lambda6*mu6*my*ny*pz - 2*d6*lambda6*mu4*mu5*mu6*my*ny*pz + mu4*mu5*mu6*mx*px*pz + lambda6*mu4*mu5*nx*px*pz + mu4*mu5*mu6*my*py*pz + lambda6*mu4*mu5*ny*py*pz + lx*lz*mu4*mu5*mu6*mx*pow(a6,2) + ly*lz*mu4*mu5*mu6*my*pow(a6,2) + lambda6*lx*lz*mu4*mu5*nx*pow(a6,2) + lambda6*ly*lz*mu4*mu5*ny*pow(a6,2) - a6*d4*lx*nx*nz*pow(lambda6,2) - a6*d6*lx*mu4*mu5*nx*nz*pow(lambda6,2) - a6*d4*ly*ny*nz*pow(lambda6,2) - a6*d6*ly*mu4*mu5*ny*nz*pow(lambda6,2) + d4*nx*nz*px*pow(lambda6,2) + d6*mu4*mu5*nx*nz*px*pow(lambda6,2) + d4*ny*nz*py*pow(lambda6,2) + d6*mu4*mu5*ny*nz*py*pow(lambda6,2) - mu4*mu5*mu6*mz*pow(a6,2)*pow(lx,2) - lambda6*mu4*mu5*nz*pow(a6,2)*pow(lx,2) - mu4*mu5*mu6*mz*pow(a6,2)*pow(ly,2) - lambda6*mu4*mu5*nz*pow(a6,2)*pow(ly,2) + a1*a6*lambda1*ly*mu4*mu5*mu6*mx*pow(mu1,-1) - a1*a6*lambda1*lx*mu4*mu5*mu6*my*pow(mu1,-1) + a1*a6*lambda1*lambda6*ly*mu4*mu5*nx*pow(mu1,-1) - a1*a6*lambda1*lambda6*lx*mu4*mu5*ny*pow(mu1,-1) + a1*lambda1*mu4*mu5*mu6*my*px*pow(mu1,-1) + a1*lambda1*lambda6*mu4*mu5*ny*px*pow(mu1,-1) - a1*lambda1*mu4*mu5*mu6*mx*py*pow(mu1,-1) - a1*lambda1*lambda6*mu4*mu5*nx*py*pow(mu1,-1) - a6*d4*lx*mx*mz*pow(mu6,2) - a6*d6*lx*mu4*mu5*mx*mz*pow(mu6,2) - a6*d4*ly*my*mz*pow(mu6,2) - a6*d6*ly*mu4*mu5*my*mz*pow(mu6,2) + d4*mx*mz*px*pow(mu6,2) + d6*mu4*mu5*mx*mz*px*pow(mu6,2) + d4*my*mz*py*pow(mu6,2) + d6*mu4*mu5*my*mz*py*pow(mu6,2) + d1*d4*pow(mu6,2)*pow(mx,2) + a6*d4*lz*pow(mu6,2)*pow(mx,2) + d1*d6*mu4*mu5*pow(mu6,2)*pow(mx,2) + a6*d6*lz*mu4*mu5*pow(mu6,2)*pow(mx,2) - d4*pz*pow(mu6,2)*pow(mx,2) - d6*mu4*mu5*pz*pow(mu6,2)*pow(mx,2) + d1*d4*pow(mu6,2)*pow(my,2) + a6*d4*lz*pow(mu6,2)*pow(my,2) + d1*d6*mu4*mu5*pow(mu6,2)*pow(my,2) + a6*d6*lz*mu4*mu5*pow(mu6,2)*pow(my,2) - d4*pz*pow(mu6,2)*pow(my,2) - d6*mu4*mu5*pz*pow(mu6,2)*pow(my,2) + d1*d4*pow(lambda6,2)*pow(nx,2) + a6*d4*lz*pow(lambda6,2)*pow(nx,2) + d1*d6*mu4*mu5*pow(lambda6,2)*pow(nx,2) + a6*d6*lz*mu4*mu5*pow(lambda6,2)*pow(nx,2) - d4*pz*pow(lambda6,2)*pow(nx,2) - d6*mu4*mu5*pz*pow(lambda6,2)*pow(nx,2) + d1*d4*pow(lambda6,2)*pow(ny,2) + a6*d4*lz*pow(lambda6,2)*pow(ny,2) + d1*d6*mu4*mu5*pow(lambda6,2)*pow(ny,2) + a6*d6*lz*mu4*mu5*pow(lambda6,2)*pow(ny,2) - d4*pz*pow(lambda6,2)*pow(ny,2) - d6*mu4*mu5*pz*pow(lambda6,2)*pow(ny,2) - mu4*mu5*mu6*mz*pow(px,2) - lambda6*mu4*mu5*nz*pow(px,2) - mu4*mu5*mu6*mz*pow(py,2) - lambda6*mu4*mu5*nz*pow(py,2) + d1*d4*mx*nx*sin(2*alpha6) + a6*d4*lz*mx*nx*sin(2*alpha6) + d1*d6*mu4*mu5*mx*nx*sin(2*alpha6) + a6*d6*lz*mu4*mu5*mx*nx*sin(2*alpha6) + d1*d4*my*ny*sin(2*alpha6) + a6*d4*lz*my*ny*sin(2*alpha6) + d1*d6*mu4*mu5*my*ny*sin(2*alpha6) + a6*d6*lz*mu4*mu5*my*ny*sin(2*alpha6) + lambda4*(-(lambda6*(-((a2 - a3 - a4 + a5)*mu5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)) + d5*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz))) + mu6*(2*d5*lambda6*(d1 + a6*lz)*(mx*nx + my*ny) + (a2 - a3 - a4 + a5)*mu5*(a6*ly*mx - a6*lx*my + my*px - mx*py) + d5*mu6*(mx*mz*px + my*mz*py - pz*pow(mx,2) - pz*pow(my,2) + d1*(pow(mx,2) + pow(my,2)) + a6*(-((lx*mx + ly*my)*mz) + lz*(pow(mx,2) + pow(my,2))))) + d5*pow(lambda6,2)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))) + lambda5*(d6*pow(mu6,2)*(mx*mz*px + my*mz*py - pz*pow(mx,2) - pz*pow(my,2) + d1*(pow(mx,2) + pow(my,2)) + a6*(-((lx*mx + ly*my)*mz) + lz*(pow(mx,2) + pow(my,2)))) + d6*pow(lambda6,2)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))) - mu6*(d1*mx*px + d1*my*py - mx*px*pz - my*py*pz + a6*(-(d1*(lx*mx + ly*my)) + lz*mx*px - 2*lx*mz*px + lz*my*py - 2*ly*mz*py + lx*mx*pz + ly*my*pz) + pow(a6,2)*(-(lx*lz*mx) + ly*(-(lz*my) + ly*mz) + mz*pow(lx,2)) - a1*lambda1*(a6*ly*mx - a6*lx*my + my*px - mx*py)*pow(mu1,-1) + mz*pow(px,2) + mz*pow(py,2)) - lambda6*(-(a6*d1*lx*nx) - a6*d1*ly*ny + d1*nx*px + a6*lz*nx*px - 2*a6*lx*nz*px + d1*ny*py + a6*lz*ny*py - 2*a6*ly*nz*py + a6*lx*nx*pz + a6*ly*ny*pz - nx*px*pz - ny*py*pz + d6*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz) - lx*lz*nx*pow(a6,2) - ly*lz*ny*pow(a6,2) + nz*pow(a6,2)*pow(lx,2) + nz*pow(a6,2)*pow(ly,2) - a1*lambda1*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*pow(mu1,-1) + nz*pow(px,2) + nz*pow(py,2)) + d6*(d1 + a6*lz)*(mx*nx + my*ny)*sin(2*alpha6))))));

return(dummyexp);

}