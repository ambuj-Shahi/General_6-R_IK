
# include <stdio.h>
# include <math.h>

double cx011(double* coeff)
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



double dummyexp = (pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1)*(-8*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*pow(d6,2)*pow(lambda6,3)*(pow(nx,2) + pow(ny,2) + pow(nz,2)) + 8*d6*pow(lambda6,2)*(2*a6*d6*lx*mu6*mx*nx*ny - 2*a6*d6*ly*mu6*my*nx*ny - 2*a6*d1*ly*nx*nz - 2*a6*d6*ly*mu6*mz*nx*nz + 2*a6*d1*lx*ny*nz + 2*a6*d6*lx*mu6*mz*ny*nz - 4*a6*lx*nx*ny*px - 2*d6*mu6*mx*nx*ny*px - 2*d1*ny*nz*px - 2*a6*lz*ny*nz*px - 2*d6*mu6*mz*ny*nz*px + 4*a6*ly*nx*ny*py + 2*d6*mu6*my*nx*ny*py + 2*d1*nx*nz*py + 2*a6*lz*nx*nz*py + 2*d6*mu6*mz*nx*nz*py + 2*a6*ly*nx*nz*pz - 2*a6*lx*ny*nz*pz + 2*ny*nz*px*pz - 2*nx*nz*py*pz - 2*ly*lz*nx*nz*pow(a6,2) + 2*lx*lz*ny*nz*pow(a6,2) + 2*nx*ny*pow(a6,2)*pow(lx,2) - 2*nx*ny*pow(a6,2)*pow(ly,2) - 3*a6*d6*ly*mu6*mx*pow(nx,2) + a6*d6*lx*mu6*my*pow(nx,2) + 2*a6*ly*px*pow(nx,2) - d6*mu6*my*px*pow(nx,2) + 2*a6*lx*py*pow(nx,2) + 3*d6*mu6*mx*py*pow(nx,2) - 2*px*py*pow(nx,2) - 2*lx*ly*pow(a6,2)*pow(nx,2) - a6*d6*ly*mu6*mx*pow(ny,2) + 3*a6*d6*lx*mu6*my*pow(ny,2) - 2*a6*ly*px*pow(ny,2) - 3*d6*mu6*my*px*pow(ny,2) - 2*a6*lx*py*pow(ny,2) + d6*mu6*mx*py*pow(ny,2) + 2*px*py*pow(ny,2) + 2*lx*ly*pow(a6,2)*pow(ny,2) - 2*a1*(d2 + d3*lambda2 + d4*cos(alpha2 + alpha3) + d5*cos(alpha2 + alpha3 - alpha4) + d6*cos(alpha2 + alpha3 - alpha4 + alpha5))*pow(mu1,-1)*(pow(nx,2) + pow(ny,2)) - a6*d6*ly*mu6*mx*pow(nz,2) + a6*d6*lx*mu6*my*pow(nz,2) - d6*mu6*my*px*pow(nz,2) + d6*mu6*mx*py*pow(nz,2) + 2*nx*ny*pow(px,2) - 2*nx*ny*pow(py,2) - lambda1*pow(mu1,-1)*(d1*mz*ny*px + d1*my*nz*px + 2*mx*nx*px*py - 2*my*ny*px*py + mz*nx*py*pz + mx*nz*py*pz + ly*(2*lx*mx*nx + ly*my*nx + lz*mz*nx + ly*mx*ny - 2*lx*my*ny + lz*mx*nz)*pow(a6,2) + 2*a1*(nx*nz*px + ny*(nz*py - ny*pz) - pz*pow(nx,2) + d1*(pow(nx,2) + pow(ny,2))) + a6*(d1*ly*(mz*nx + mx*nz) - 2*ly*mx*nx*px + 2*lx*my*nx*px + 2*lx*mx*ny*px + 2*ly*my*ny*px + lz*mz*ny*px + lz*my*nz*px - 2*lx*mx*nx*py - 2*ly*my*nx*py - 2*ly*mx*ny*py + 2*lx*my*ny*py + lx*mz*ny*pz + lx*my*nz*pz + 2*a1*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))) + my*nx*pow(py,2) + mx*ny*pow(py,2))) + 8*(2*a2*a3*a6*ly*mu6*mx - 2*a2*a4*a6*ly*mu6*mx - 2*a3*a4*a6*ly*mu6*mx + 2*a2*a5*a6*ly*mu6*mx + 2*a3*a5*a6*ly*mu6*mx - 2*a4*a5*a6*ly*mu6*mx + 2*a6*d3*d4*lambda3*ly*mu6*mx + 2*a6*d4*d5*lambda4*ly*mu6*mx + 2*a6*d3*d5*lambda3*lambda4*ly*mu6*mx - 2*a6*d2*d4*ly*mu2*mu3*mu6*mx - 2*a6*d2*d5*lambda4*ly*mu2*mu3*mu6*mx + 2*a6*d2*d5*lambda3*ly*mu2*mu4*mu6*mx + 2*a6*d3*d5*ly*mu3*mu4*mu6*mx - 2*a2*a3*a6*lx*mu6*my + 2*a2*a4*a6*lx*mu6*my + 2*a3*a4*a6*lx*mu6*my - 2*a2*a5*a6*lx*mu6*my - 2*a3*a5*a6*lx*mu6*my + 2*a4*a5*a6*lx*mu6*my - 2*a6*d3*d4*lambda3*lx*mu6*my - 2*a6*d4*d5*lambda4*lx*mu6*my - 2*a6*d3*d5*lambda3*lambda4*lx*mu6*my + 2*a6*d2*d4*lx*mu2*mu3*mu6*my + 2*a6*d2*d5*lambda4*lx*mu2*mu3*mu6*my - 2*a6*d2*d5*lambda3*lx*mu2*mu4*mu6*my - 2*a6*d3*d5*lx*mu3*mu4*mu6*my + 2*a2*a3*mu6*my*px - 2*a2*a4*mu6*my*px - 2*a3*a4*mu6*my*px + 2*a2*a5*mu6*my*px + 2*a3*a5*mu6*my*px - 2*a4*a5*mu6*my*px + 2*d3*d4*lambda3*mu6*my*px + 2*d4*d5*lambda4*mu6*my*px + 2*d3*d5*lambda3*lambda4*mu6*my*px - 2*d2*d4*mu2*mu3*mu6*my*px - 2*d2*d5*lambda4*mu2*mu3*mu6*my*px + 2*d2*d5*lambda3*mu2*mu4*mu6*my*px + 2*d3*d5*mu3*mu4*mu6*my*px - 2*a2*a3*mu6*mx*py + 2*a2*a4*mu6*mx*py + 2*a3*a4*mu6*mx*py - 2*a2*a5*mu6*mx*py - 2*a3*a5*mu6*mx*py + 2*a4*a5*mu6*mx*py - 2*d3*d4*lambda3*mu6*mx*py - 2*d4*d5*lambda4*mu6*mx*py - 2*d3*d5*lambda3*lambda4*mu6*mx*py + 2*d2*d4*mu2*mu3*mu6*mx*py + 2*d2*d5*lambda4*mu2*mu3*mu6*mx*py - 2*d2*d5*lambda3*mu2*mu4*mu6*mx*py - 2*d3*d5*mu3*mu4*mu6*mx*py - 2*a6*lx*mu6*mx*px*py - a6*ly*mu6*mx*pow(a1,2) + a6*lx*mu6*my*pow(a1,2) - mu6*my*px*pow(a1,2) + mu6*mx*py*pow(a1,2) + a6*ly*mu6*mx*pow(a2,2) - a6*lx*mu6*my*pow(a2,2) + mu6*my*px*pow(a2,2) - mu6*mx*py*pow(a2,2) + a6*ly*mu6*mx*pow(a3,2) - a6*lx*mu6*my*pow(a3,2) + mu6*my*px*pow(a3,2) - mu6*mx*py*pow(a3,2) + a6*ly*mu6*mx*pow(a4,2) - a6*lx*mu6*my*pow(a4,2) + mu6*my*px*pow(a4,2) - mu6*mx*py*pow(a4,2) + a6*ly*mu6*mx*pow(a5,2) - a6*lx*mu6*my*pow(a5,2) + mu6*my*px*pow(a5,2) - mu6*mx*py*pow(a5,2) + 2*lx*ly*mu6*mx*px*pow(a6,2) + a6*ly*mu6*mx*pow(d2,2) - a6*lx*mu6*my*pow(d2,2) + mu6*my*px*pow(d2,2) - mu6*mx*py*pow(d2,2) + a6*ly*mu6*mx*pow(d3,2) - a6*lx*mu6*my*pow(d3,2) + mu6*my*px*pow(d3,2) - mu6*mx*py*pow(d3,2) + a6*ly*mu6*mx*pow(d4,2) - a6*lx*mu6*my*pow(d4,2) + mu6*my*px*pow(d4,2) - mu6*mx*py*pow(d4,2) + a6*ly*mu6*mx*pow(d5,2) - a6*lx*mu6*my*pow(d5,2) + mu6*my*px*pow(d5,2) - mu6*mx*py*pow(d5,2) - 2*a6*d1*lz*mu6*my*px*pow(lambda1,2) + 2*a6*d1*lz*mu6*mx*py*pow(lambda1,2) + 2*a6*ly*mu6*my*px*py*pow(lambda1,2) + 2*a6*d1*ly*mu6*mx*pz*pow(lambda1,2) - 2*a6*d1*lx*mu6*my*pz*pow(lambda1,2) + 2*d1*mu6*my*px*pz*pow(lambda1,2) + 2*a6*lz*mu6*my*px*pz*pow(lambda1,2) - 2*d1*mu6*mx*py*pz*pow(lambda1,2) - 2*a6*lz*mu6*mx*py*pz*pow(lambda1,2) - 2*d1*ly*lz*mu6*mx*pow(a6,2)*pow(lambda1,2) + 2*d1*lx*lz*mu6*my*pow(a6,2)*pow(lambda1,2) - 2*lx*ly*mu6*my*py*pow(a6,2)*pow(lambda1,2) + 2*ly*lz*mu6*mx*pz*pow(a6,2)*pow(lambda1,2) - 2*lx*lz*mu6*my*pz*pow(a6,2)*pow(lambda1,2) - a6*ly*mu6*mx*pow(d1,2)*pow(lambda1,2) + a6*lx*mu6*my*pow(d1,2)*pow(lambda1,2) - mu6*my*px*pow(d1,2)*pow(lambda1,2) + mu6*mx*py*pow(d1,2)*pow(lambda1,2) - 3*mu6*my*px*pow(a6,2)*pow(lx,2) + mu6*mx*py*pow(a6,2)*pow(lx,2) - ly*mu6*mx*pow(a6,3)*pow(lx,2) + mu6*my*pow(a6,3)*pow(lx,3) - mu6*my*px*pow(a6,2)*pow(lambda1,2)*pow(ly,2) + 3*mu6*mx*py*pow(a6,2)*pow(lambda1,2)*pow(ly,2) + lx*mu6*my*pow(a6,3)*pow(lambda1,2)*pow(ly,2) - mu6*mx*pow(a6,3)*pow(lambda1,2)*pow(ly,3) - mu6*my*px*pow(a6,2)*pow(lambda1,2)*pow(lz,2) + mu6*mx*py*pow(a6,2)*pow(lambda1,2)*pow(lz,2) - ly*mu6*mx*pow(a6,3)*pow(lambda1,2)*pow(lz,2) + lx*mu6*my*pow(a6,3)*pow(lambda1,2)*pow(lz,2) - 2*a1*a6*d2*lx*mu6*mx*pow(mu1,-1) - 2*a1*a6*d1*lambda1*lx*mu6*mx*pow(mu1,-1) + 2*a1*a6*d4*lx*mu2*mu3*mu6*mx*pow(mu1,-1) + 2*a1*a6*d5*lambda4*lx*mu2*mu3*mu6*mx*pow(mu1,-1) + 4*a1*a6*d6*lambda4*lambda5*lx*mu2*mu3*mu6*mx*pow(mu1,-1) - 2*a1*a6*d5*lambda3*lx*mu2*mu4*mu6*mx*pow(mu1,-1) - 4*a1*a6*d6*lambda3*lambda5*lx*mu2*mu4*mu6*mx*pow(mu1,-1) + 4*a1*a6*d6*lambda3*lambda4*lx*mu2*mu5*mu6*mx*pow(mu1,-1) + 4*a1*a6*d6*lx*mu2*mu3*mu4*mu5*mu6*mx*pow(mu1,-1) - 2*a1*a6*d2*ly*mu6*my*pow(mu1,-1) - 2*a1*a6*d1*lambda1*ly*mu6*my*pow(mu1,-1) + 2*a1*a6*d4*ly*mu2*mu3*mu6*my*pow(mu1,-1) + 2*a1*a6*d5*lambda4*ly*mu2*mu3*mu6*my*pow(mu1,-1) + 4*a1*a6*d6*lambda4*lambda5*ly*mu2*mu3*mu6*my*pow(mu1,-1) - 2*a1*a6*d5*lambda3*ly*mu2*mu4*mu6*my*pow(mu1,-1) - 4*a1*a6*d6*lambda3*lambda5*ly*mu2*mu4*mu6*my*pow(mu1,-1) + 4*a1*a6*d6*lambda3*lambda4*ly*mu2*mu5*mu6*my*pow(mu1,-1) + 4*a1*a6*d6*ly*mu2*mu3*mu4*mu5*mu6*my*pow(mu1,-1) - 4*a1*a6*lambda4*lambda5*lx*mu2*mu3*px*pow(mu1,-1) + 4*a1*a6*lambda3*lambda5*lx*mu2*mu4*px*pow(mu1,-1) - 4*a1*a6*lambda3*lambda4*lx*mu2*mu5*px*pow(mu1,-1) - 4*a1*a6*lx*mu2*mu3*mu4*mu5*px*pow(mu1,-1) + 2*a1*d2*mu6*mx*px*pow(mu1,-1) + 2*a1*d1*lambda1*mu6*mx*px*pow(mu1,-1) + 2*a1*a6*lambda1*lz*mu6*mx*px*pow(mu1,-1) - 2*a1*d4*mu2*mu3*mu6*mx*px*pow(mu1,-1) - 2*a1*d5*lambda4*mu2*mu3*mu6*mx*px*pow(mu1,-1) - 4*a1*d6*lambda4*lambda5*mu2*mu3*mu6*mx*px*pow(mu1,-1) + 2*a1*d5*lambda3*mu2*mu4*mu6*mx*px*pow(mu1,-1) + 4*a1*d6*lambda3*lambda5*mu2*mu4*mu6*mx*px*pow(mu1,-1) - 4*a1*d6*lambda3*lambda4*mu2*mu5*mu6*mx*px*pow(mu1,-1) - 4*a1*d6*mu2*mu3*mu4*mu5*mu6*mx*px*pow(mu1,-1) - 4*a1*a6*lambda1*lx*mu6*mz*px*pow(mu1,-1) - 4*a1*a6*lambda4*lambda5*ly*mu2*mu3*py*pow(mu1,-1) + 4*a1*a6*lambda3*lambda5*ly*mu2*mu4*py*pow(mu1,-1) - 4*a1*a6*lambda3*lambda4*ly*mu2*mu5*py*pow(mu1,-1) - 4*a1*a6*ly*mu2*mu3*mu4*mu5*py*pow(mu1,-1) + 2*a1*d2*mu6*my*py*pow(mu1,-1) + 2*a1*d1*lambda1*mu6*my*py*pow(mu1,-1) + 2*a1*a6*lambda1*lz*mu6*my*py*pow(mu1,-1) - 2*a1*d4*mu2*mu3*mu6*my*py*pow(mu1,-1) - 2*a1*d5*lambda4*mu2*mu3*mu6*my*py*pow(mu1,-1) - 4*a1*d6*lambda4*lambda5*mu2*mu3*mu6*my*py*pow(mu1,-1) + 2*a1*d5*lambda3*mu2*mu4*mu6*my*py*pow(mu1,-1) + 4*a1*d6*lambda3*lambda5*mu2*mu4*mu6*my*py*pow(mu1,-1) - 4*a1*d6*lambda3*lambda4*mu2*mu5*mu6*my*py*pow(mu1,-1) - 4*a1*d6*mu2*mu3*mu4*mu5*mu6*my*py*pow(mu1,-1) - 4*a1*a6*lambda1*ly*mu6*mz*py*pow(mu1,-1) + 2*a1*a6*lambda1*lx*mu6*mx*pz*pow(mu1,-1) + 2*a1*a6*lambda1*ly*mu6*my*pz*pow(mu1,-1) - 2*a1*lambda1*mu6*mx*px*pz*pow(mu1,-1) - 2*a1*lambda1*mu6*my*py*pz*pow(mu1,-1) - 2*a1*lambda1*lx*lz*mu6*mx*pow(a6,2)*pow(mu1,-1) - 2*a1*lambda1*ly*lz*mu6*my*pow(a6,2)*pow(mu1,-1) + 2*a1*lambda4*lambda5*mu2*mu3*pow(a6,2)*pow(lx,2)*pow(mu1,-1) - 2*a1*lambda3*lambda5*mu2*mu4*pow(a6,2)*pow(lx,2)*pow(mu1,-1) + 2*a1*lambda3*lambda4*mu2*mu5*pow(a6,2)*pow(lx,2)*pow(mu1,-1) + 2*a1*mu2*mu3*mu4*mu5*pow(a6,2)*pow(lx,2)*pow(mu1,-1) + 2*a1*lambda1*mu6*mz*pow(a6,2)*pow(lx,2)*pow(mu1,-1) + 2*a1*lambda4*lambda5*mu2*mu3*pow(a6,2)*pow(ly,2)*pow(mu1,-1) - 2*a1*lambda3*lambda5*mu2*mu4*pow(a6,2)*pow(ly,2)*pow(mu1,-1) + 2*a1*lambda3*lambda4*mu2*mu5*pow(a6,2)*pow(ly,2)*pow(mu1,-1) + 2*a1*mu2*mu3*mu4*mu5*pow(a6,2)*pow(ly,2)*pow(mu1,-1) + 2*a1*lambda1*mu6*mz*pow(a6,2)*pow(ly,2)*pow(mu1,-1) - 2*a6*d1*lz*mu6*my*px*pow(mu1,2) + 2*a6*d1*lz*mu6*mx*py*pow(mu1,2) + 2*a6*ly*mu6*my*px*py*pow(mu1,2) + 2*a6*d1*ly*mu6*mx*pz*pow(mu1,2) - 2*a6*d1*lx*mu6*my*pz*pow(mu1,2) + 2*d1*mu6*my*px*pz*pow(mu1,2) + 2*a6*lz*mu6*my*px*pz*pow(mu1,2) - 2*d1*mu6*mx*py*pz*pow(mu1,2) - 2*a6*lz*mu6*mx*py*pz*pow(mu1,2) - 2*d1*ly*lz*mu6*mx*pow(a6,2)*pow(mu1,2) + 2*d1*lx*lz*mu6*my*pow(a6,2)*pow(mu1,2) - 2*lx*ly*mu6*my*py*pow(a6,2)*pow(mu1,2) + 2*ly*lz*mu6*mx*pz*pow(a6,2)*pow(mu1,2) - 2*lx*lz*mu6*my*pz*pow(a6,2)*pow(mu1,2) - a6*ly*mu6*mx*pow(d1,2)*pow(mu1,2) + a6*lx*mu6*my*pow(d1,2)*pow(mu1,2) - mu6*my*px*pow(d1,2)*pow(mu1,2) + mu6*mx*py*pow(d1,2)*pow(mu1,2) - mu6*my*px*pow(a6,2)*pow(ly,2)*pow(mu1,2) + 3*mu6*mx*py*pow(a6,2)*pow(ly,2)*pow(mu1,2) + lx*mu6*my*pow(a6,3)*pow(ly,2)*pow(mu1,2) - mu6*mx*pow(a6,3)*pow(ly,3)*pow(mu1,2) - mu6*my*px*pow(a6,2)*pow(lz,2)*pow(mu1,2) + mu6*mx*py*pow(a6,2)*pow(lz,2)*pow(mu1,2) - ly*mu6*mx*pow(a6,3)*pow(lz,2)*pow(mu1,2) + lx*mu6*my*pow(a6,3)*pow(lz,2)*pow(mu1,2) - 4*a6*d6*lx*mx*my*px*pow(mu6,2) - 2*a6*d1*d6*ly*mx*mz*pow(lambda1,2)*pow(mu6,2) + 2*a6*d1*d6*lx*my*mz*pow(lambda1,2)*pow(mu6,2) - 2*d1*d6*my*mz*px*pow(lambda1,2)*pow(mu6,2) - 2*a6*d6*lz*my*mz*px*pow(lambda1,2)*pow(mu6,2) + 4*a6*d6*ly*mx*my*py*pow(lambda1,2)*pow(mu6,2) + 2*d1*d6*mx*mz*py*pow(lambda1,2)*pow(mu6,2) + 2*a6*d6*lz*mx*mz*py*pow(lambda1,2)*pow(mu6,2) + 2*a6*d6*ly*mx*mz*pz*pow(lambda1,2)*pow(mu6,2) - 2*a6*d6*lx*my*mz*pz*pow(lambda1,2)*pow(mu6,2) + 2*d6*my*mz*px*pz*pow(lambda1,2)*pow(mu6,2) - 2*d6*mx*mz*py*pz*pow(lambda1,2)*pow(mu6,2) - 2*d6*ly*lz*mx*mz*pow(a6,2)*pow(lambda1,2)*pow(mu6,2) + 2*d6*lx*lz*my*mz*pow(a6,2)*pow(lambda1,2)*pow(mu6,2) + 2*d6*mx*my*pow(a6,2)*pow(lx,2)*pow(mu6,2) - 2*d6*mx*my*pow(a6,2)*pow(lambda1,2)*pow(ly,2)*pow(mu6,2) + 2*a1*a6*d6*lambda1*lx*mx*mz*pow(mu1,-1)*pow(mu6,2) + 2*a1*a6*d6*lambda1*ly*my*mz*pow(mu1,-1)*pow(mu6,2) - 2*a1*d6*lambda1*mx*mz*px*pow(mu1,-1)*pow(mu6,2) - 2*a1*d6*lambda1*my*mz*py*pow(mu1,-1)*pow(mu6,2) - 2*a6*d1*d6*ly*mx*mz*pow(mu1,2)*pow(mu6,2) + 2*a6*d1*d6*lx*my*mz*pow(mu1,2)*pow(mu6,2) - 2*d1*d6*my*mz*px*pow(mu1,2)*pow(mu6,2) - 2*a6*d6*lz*my*mz*px*pow(mu1,2)*pow(mu6,2) + 4*a6*d6*ly*mx*my*py*pow(mu1,2)*pow(mu6,2) + 2*d1*d6*mx*mz*py*pow(mu1,2)*pow(mu6,2) + 2*a6*d6*lz*mx*mz*py*pow(mu1,2)*pow(mu6,2) + 2*a6*d6*ly*mx*mz*pz*pow(mu1,2)*pow(mu6,2) - 2*a6*d6*lx*my*mz*pz*pow(mu1,2)*pow(mu6,2) + 2*d6*my*mz*px*pz*pow(mu1,2)*pow(mu6,2) - 2*d6*mx*mz*py*pz*pow(mu1,2)*pow(mu6,2) - 2*d6*ly*lz*mx*mz*pow(a6,2)*pow(mu1,2)*pow(mu6,2) + 2*d6*lx*lz*my*mz*pow(a6,2)*pow(mu1,2)*pow(mu6,2) - 2*d6*mx*my*pow(a6,2)*pow(ly,2)*pow(mu1,2)*pow(mu6,2) + 2*a6*d6*ly*px*pow(mu6,2)*pow(mx,2) + 2*a6*d6*lx*py*pow(mu6,2)*pow(mx,2) - 2*d6*px*py*pow(mu6,2)*pow(mx,2) - 2*d6*lx*ly*pow(a6,2)*pow(mu6,2)*pow(mx,2) - 2*a1*d2*d6*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) - 2*a1*d1*d6*lambda1*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) - 2*a1*a6*d6*lambda1*lz*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) + 2*a1*d4*d6*mu2*mu3*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) + 2*a1*d5*d6*lambda4*mu2*mu3*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) - 2*a1*d5*d6*lambda3*mu2*mu4*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) + 2*a1*d6*lambda1*pz*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) + 2*a1*lambda4*lambda5*mu2*mu3*pow(d6,2)*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) - 2*a1*lambda3*lambda5*mu2*mu4*pow(d6,2)*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) + 2*a1*lambda3*lambda4*mu2*mu5*pow(d6,2)*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) + 2*a1*mu2*mu3*mu4*mu5*pow(d6,2)*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) + a6*lx*my*pow(d6,2)*pow(mu6,3)*pow(mx,2) - my*px*pow(d6,2)*pow(mu6,3)*pow(mx,2) - a6*ly*pow(d6,2)*pow(mu6,3)*pow(mx,3) + py*pow(d6,2)*pow(mu6,3)*pow(mx,3) - 2*a6*d6*ly*px*pow(lambda1,2)*pow(mu6,2)*pow(my,2) - 2*a6*d6*lx*py*pow(lambda1,2)*pow(mu6,2)*pow(my,2) + 2*d6*px*py*pow(lambda1,2)*pow(mu6,2)*pow(my,2) + 2*d6*lx*ly*pow(a6,2)*pow(lambda1,2)*pow(mu6,2)*pow(my,2) - 2*a1*d2*d6*pow(mu1,-1)*pow(mu6,2)*pow(my,2) - 2*a1*d1*d6*lambda1*pow(mu1,-1)*pow(mu6,2)*pow(my,2) - 2*a1*a6*d6*lambda1*lz*pow(mu1,-1)*pow(mu6,2)*pow(my,2) + 2*a1*d4*d6*mu2*mu3*pow(mu1,-1)*pow(mu6,2)*pow(my,2) + 2*a1*d5*d6*lambda4*mu2*mu3*pow(mu1,-1)*pow(mu6,2)*pow(my,2) - 2*a1*d5*d6*lambda3*mu2*mu4*pow(mu1,-1)*pow(mu6,2)*pow(my,2) + 2*a1*d6*lambda1*pz*pow(mu1,-1)*pow(mu6,2)*pow(my,2) + 2*a1*lambda4*lambda5*mu2*mu3*pow(d6,2)*pow(mu1,-1)*pow(mu6,2)*pow(my,2) - 2*a1*lambda3*lambda5*mu2*mu4*pow(d6,2)*pow(mu1,-1)*pow(mu6,2)*pow(my,2) + 2*a1*lambda3*lambda4*mu2*mu5*pow(d6,2)*pow(mu1,-1)*pow(mu6,2)*pow(my,2) + 2*a1*mu2*mu3*mu4*mu5*pow(d6,2)*pow(mu1,-1)*pow(mu6,2)*pow(my,2) - 2*a6*d6*ly*px*pow(mu1,2)*pow(mu6,2)*pow(my,2) - 2*a6*d6*lx*py*pow(mu1,2)*pow(mu6,2)*pow(my,2) + 2*d6*px*py*pow(mu1,2)*pow(mu6,2)*pow(my,2) + 2*d6*lx*ly*pow(a6,2)*pow(mu1,2)*pow(mu6,2)*pow(my,2) - a6*ly*mx*pow(d6,2)*pow(lambda1,2)*pow(mu6,3)*pow(my,2) + mx*py*pow(d6,2)*pow(lambda1,2)*pow(mu6,3)*pow(my,2) - a6*ly*mx*pow(d6,2)*pow(mu1,2)*pow(mu6,3)*pow(my,2) + mx*py*pow(d6,2)*pow(mu1,2)*pow(mu6,3)*pow(my,2) + a6*lx*pow(d6,2)*pow(lambda1,2)*pow(mu6,3)*pow(my,3) - px*pow(d6,2)*pow(lambda1,2)*pow(mu6,3)*pow(my,3) + a6*lx*pow(d6,2)*pow(mu1,2)*pow(mu6,3)*pow(my,3) - px*pow(d6,2)*pow(mu1,2)*pow(mu6,3)*pow(my,3) - a6*ly*mx*pow(d6,2)*pow(lambda1,2)*pow(mu6,3)*pow(mz,2) + a6*lx*my*pow(d6,2)*pow(lambda1,2)*pow(mu6,3)*pow(mz,2) - my*px*pow(d6,2)*pow(lambda1,2)*pow(mu6,3)*pow(mz,2) + mx*py*pow(d6,2)*pow(lambda1,2)*pow(mu6,3)*pow(mz,2) - a6*ly*mx*pow(d6,2)*pow(mu1,2)*pow(mu6,3)*pow(mz,2) + a6*lx*my*pow(d6,2)*pow(mu1,2)*pow(mu6,3)*pow(mz,2) - my*px*pow(d6,2)*pow(mu1,2)*pow(mu6,3)*pow(mz,2) + mx*py*pow(d6,2)*pow(mu1,2)*pow(mu6,3)*pow(mz,2) - a6*ly*mu6*mx*pow(px,2) + 3*a6*lx*mu6*my*pow(px,2) + mu6*mx*py*pow(px,2) + 2*a1*lambda4*lambda5*mu2*mu3*pow(mu1,-1)*pow(px,2) - 2*a1*lambda3*lambda5*mu2*mu4*pow(mu1,-1)*pow(px,2) + 2*a1*lambda3*lambda4*mu2*mu5*pow(mu1,-1)*pow(px,2) + 2*a1*mu2*mu3*mu4*mu5*pow(mu1,-1)*pow(px,2) + 2*a1*lambda1*mu6*mz*pow(mu1,-1)*pow(px,2) + 2*d6*mx*my*pow(mu6,2)*pow(px,2) - mu6*my*pow(px,3) - 3*a6*ly*mu6*mx*pow(lambda1,2)*pow(py,2) + a6*lx*mu6*my*pow(lambda1,2)*pow(py,2) - mu6*my*px*pow(lambda1,2)*pow(py,2) + 2*a1*lambda4*lambda5*mu2*mu3*pow(mu1,-1)*pow(py,2) - 2*a1*lambda3*lambda5*mu2*mu4*pow(mu1,-1)*pow(py,2) + 2*a1*lambda3*lambda4*mu2*mu5*pow(mu1,-1)*pow(py,2) + 2*a1*mu2*mu3*mu4*mu5*pow(mu1,-1)*pow(py,2) + 2*a1*lambda1*mu6*mz*pow(mu1,-1)*pow(py,2) - 3*a6*ly*mu6*mx*pow(mu1,2)*pow(py,2) + a6*lx*mu6*my*pow(mu1,2)*pow(py,2) - mu6*my*px*pow(mu1,2)*pow(py,2) - 2*d6*mx*my*pow(lambda1,2)*pow(mu6,2)*pow(py,2) - 2*d6*mx*my*pow(mu1,2)*pow(mu6,2)*pow(py,2) - lambda2*(pow(mu1,-1)*(2*(d3 + d5*mu3*mu4)*mu6*(-(d2*mu1*(a6*ly*mx - a6*lx*my + my*px - mx*py)) + a1*(a6*lx*mx + a6*ly*my - mx*px - my*py + d6*mu6*(pow(mx,2) + pow(my,2)))) + a1*lambda5*mu3*mu4*(-4*a6*lx*px - 4*a6*ly*py + 4*d6*mu6*(a6*lx*mx + a6*ly*my - mx*px - my*py) + 2*pow(a6,2)*pow(lx,2) + 2*pow(a6,2)*pow(ly,2) + pow(d6,2)*pow(mx,2) + pow(d6,2)*pow(my,2) - cos(2*alpha6)*pow(d6,2)*(pow(mx,2) + pow(my,2)) + 2*pow(px,2) + 2*pow(py,2)) - a1*lambda4*mu3*mu5*(-4*a6*lx*px - 4*a6*ly*py + 4*d6*mu6*(a6*lx*mx + a6*ly*my - mx*px - my*py) + 2*pow(a6,2)*pow(lx,2) + 2*pow(a6,2)*pow(ly,2) + pow(d6,2)*pow(mx,2) + pow(d6,2)*pow(my,2) - cos(2*alpha6)*pow(d6,2)*(pow(mx,2) + pow(my,2)) + 2*pow(px,2) + 2*pow(py,2))) + lambda3*(lambda4*pow(mu1,-1)*(2*d5*mu6*(-(d2*mu1*(a6*ly*mx - a6*lx*my + my*px - mx*py)) + a1*(a6*lx*mx + a6*ly*my - mx*px - my*py + d6*mu6*(pow(mx,2) + pow(my,2)))) + a1*lambda5*(-4*a6*lx*px - 4*a6*ly*py + 4*d6*mu6*(a6*lx*mx + a6*ly*my - mx*px - my*py) + 2*pow(a6,2)*pow(lx,2) + 2*pow(a6,2)*pow(ly,2) + pow(d6,2)*pow(mx,2) + pow(d6,2)*pow(my,2) - cos(2*alpha6)*pow(d6,2)*(pow(mx,2) + pow(my,2)) + 2*pow(px,2) + 2*pow(py,2))) + 2*(d2*d4*mu6*(-(a6*ly*mx) + a6*lx*my - my*px + mx*py) + a1*pow(mu1,-1)*(d4*mu6*(a6*lx*mx + a6*ly*my - mx*px - my*py + d6*mu6*(pow(mx,2) + pow(my,2))) + mu4*mu5*(-2*a6*(lx*px + ly*py) + 2*d6*mu6*(a6*lx*mx + a6*ly*my - mx*px - my*py) + pow(a6,2)*(pow(lx,2) + pow(ly,2)) + pow(d6,2)*pow(mu6,2)*(pow(mx,2) + pow(my,2)) + pow(px,2) + pow(py,2)))))) + mu6*mx*pow(lambda1,2)*pow(py,3) + mu6*mx*pow(mu1,2)*pow(py,3) - a6*ly*mu6*mx*pow(lambda1,2)*pow(pz,2) + a6*lx*mu6*my*pow(lambda1,2)*pow(pz,2) - mu6*my*px*pow(lambda1,2)*pow(pz,2) + mu6*mx*py*pow(lambda1,2)*pow(pz,2) - a6*ly*mu6*mx*pow(mu1,2)*pow(pz,2) + a6*lx*mu6*my*pow(mu1,2)*pow(pz,2) - mu6*my*px*pow(mu1,2)*pow(pz,2) + mu6*mx*py*pow(mu1,2)*pow(pz,2) + a6*lx*mu6*mx*my*nx*pow(d6,2)*sin(2*alpha6) + a6*d1*d6*lx*mz*ny*pow(lambda1,2)*sin(2*alpha6) + a6*d1*d6*lx*my*nz*pow(lambda1,2)*sin(2*alpha6) + d1*d6*mz*nx*py*pow(lambda1,2)*sin(2*alpha6) + a6*d6*lz*mz*nx*py*pow(lambda1,2)*sin(2*alpha6) + d1*d6*mx*nz*py*pow(lambda1,2)*sin(2*alpha6) + a6*d6*lz*mx*nz*py*pow(lambda1,2)*sin(2*alpha6) + a6*d6*ly*mz*nx*pz*pow(lambda1,2)*sin(2*alpha6) + a6*d6*ly*mx*nz*pz*pow(lambda1,2)*sin(2*alpha6) + d6*mz*ny*px*pz*pow(lambda1,2)*sin(2*alpha6) + d6*my*nz*px*pz*pow(lambda1,2)*sin(2*alpha6) + d6*lx*lz*mz*ny*pow(a6,2)*pow(lambda1,2)*sin(2*alpha6) + d6*lx*lz*my*nz*pow(a6,2)*pow(lambda1,2)*sin(2*alpha6) + d6*my*nx*pow(a6,2)*pow(lx,2)*sin(2*alpha6) + d6*mx*ny*pow(a6,2)*pow(lx,2)*sin(2*alpha6) + a1*a6*d6*lambda1*lx*mz*nx*pow(mu1,-1)*sin(2*alpha6) + a1*a6*d6*lambda1*ly*mz*ny*pow(mu1,-1)*sin(2*alpha6) + a1*a6*d6*lambda1*lx*mx*nz*pow(mu1,-1)*sin(2*alpha6) + a1*a6*d6*lambda1*ly*my*nz*pow(mu1,-1)*sin(2*alpha6) + a6*d1*d6*lx*mz*ny*pow(mu1,2)*sin(2*alpha6) + a6*d1*d6*lx*my*nz*pow(mu1,2)*sin(2*alpha6) + d1*d6*mz*nx*py*pow(mu1,2)*sin(2*alpha6) + a6*d6*lz*mz*nx*py*pow(mu1,2)*sin(2*alpha6) + d1*d6*mx*nz*py*pow(mu1,2)*sin(2*alpha6) + a6*d6*lz*mx*nz*py*pow(mu1,2)*sin(2*alpha6) + a6*d6*ly*mz*nx*pz*pow(mu1,2)*sin(2*alpha6) + a6*d6*ly*mx*nz*pz*pow(mu1,2)*sin(2*alpha6) + d6*mz*ny*px*pz*pow(mu1,2)*sin(2*alpha6) + d6*my*nz*px*pz*pow(mu1,2)*sin(2*alpha6) + d6*lx*lz*mz*ny*pow(a6,2)*pow(mu1,2)*sin(2*alpha6) + d6*lx*lz*my*nz*pow(a6,2)*pow(mu1,2)*sin(2*alpha6) + a6*lx*mu6*my*mz*nz*pow(d6,2)*pow(mu1,2)*sin(2*alpha6) + mu6*mx*my*ny*py*pow(d6,2)*pow(mu1,2)*sin(2*alpha6) + mu6*mx*mz*nz*py*pow(d6,2)*pow(mu1,2)*sin(2*alpha6) + d6*my*nx*pow(px,2)*sin(2*alpha6) + d6*mx*ny*pow(px,2)*sin(2*alpha6)) - lambda6*(-16*a1*lambda1*pow(mu1,-1)*(d1*nx*px + d1*ny*py - nx*px*pz - ny*py*pz + a6*(-(d1*(lx*nx + ly*ny)) + lz*nx*px - 2*lx*nz*px + lz*ny*py - 2*ly*nz*py + lx*nx*pz + ly*ny*pz) + pow(a6,2)*(-(lx*lz*nx) + ly*(-(lz*ny) + ly*nz) + nz*pow(lx,2)) + nz*pow(px,2) + nz*pow(py,2)) + 2*(-8*a2*a3*a6*ly*nx + 8*a2*a4*a6*ly*nx + 8*a3*a4*a6*ly*nx - 8*a2*a5*a6*ly*nx - 8*a3*a5*a6*ly*nx + 8*a4*a5*a6*ly*nx + 8*a2*a3*a6*lx*ny - 8*a2*a4*a6*lx*ny - 8*a3*a4*a6*lx*ny + 8*a2*a5*a6*lx*ny + 8*a3*a5*a6*lx*ny - 8*a4*a5*a6*lx*ny - 8*a2*a3*ny*px + 8*a2*a4*ny*px + 8*a3*a4*ny*px - 8*a2*a5*ny*px - 8*a3*a5*ny*px + 8*a4*a5*ny*px + 8*a6*d1*lz*ny*px + 8*a2*a3*nx*py - 8*a2*a4*nx*py - 8*a3*a4*nx*py + 8*a2*a5*nx*py + 8*a3*a5*nx*py - 8*a4*a5*nx*py - 8*a6*d1*lz*nx*py + 8*a6*lx*nx*px*py - 8*a6*ly*ny*px*py - 8*a6*d1*ly*nx*pz + 8*a6*d1*lx*ny*pz - 8*d1*ny*px*pz - 8*a6*lz*ny*px*pz + 8*d1*nx*py*pz + 8*a6*lz*nx*py*pz + 4*a6*ly*nx*pow(a1,2) - 4*a6*lx*ny*pow(a1,2) + 4*ny*px*pow(a1,2) - 4*nx*py*pow(a1,2) - 4*a6*ly*nx*pow(a2,2) + 4*a6*lx*ny*pow(a2,2) - 4*ny*px*pow(a2,2) + 4*nx*py*pow(a2,2) - 4*a6*ly*nx*pow(a3,2) + 4*a6*lx*ny*pow(a3,2) - 4*ny*px*pow(a3,2) + 4*nx*py*pow(a3,2) - 4*a6*ly*nx*pow(a4,2) + 4*a6*lx*ny*pow(a4,2) - 4*ny*px*pow(a4,2) + 4*nx*py*pow(a4,2) - 4*a6*ly*nx*pow(a5,2) + 4*a6*lx*ny*pow(a5,2) - 4*ny*px*pow(a5,2) + 4*nx*py*pow(a5,2) + 8*d1*ly*lz*nx*pow(a6,2) - 8*d1*lx*lz*ny*pow(a6,2) - 8*lx*ly*nx*px*pow(a6,2) + 8*lx*ly*ny*py*pow(a6,2) - 8*ly*lz*nx*pz*pow(a6,2) + 8*lx*lz*ny*pz*pow(a6,2) + 4*d6*mu6*(d1*(mz*ny + my*nz)*px + a6*(d1*ly*(mz*nx + mx*nz) + 2*lx*my*nx*px + 2*lx*mx*ny*px + lz*mz*ny*px + lz*my*nz*px - 2*lx*mx*nx*py + 2*lx*my*ny*py - 2*ly*(mx*nx*px - my*ny*px + my*nx*py + mx*ny*py) + lx*mz*ny*pz + lx*my*nz*pz) + py*(-2*my*ny*px + my*nx*py + mz*nx*pz + mx*(2*nx*px + ny*py + nz*pz)) + ly*(2*lx*mx*nx + ly*my*nx + lz*mz*nx + ly*mx*ny - 2*lx*my*ny + lz*mx*nz)*pow(a6,2)) + 4*a6*ly*nx*pow(d1,2) - 4*a6*lx*ny*pow(d1,2) + 4*ny*px*pow(d1,2) - 4*nx*py*pow(d1,2) - 4*a6*ly*nx*pow(d2,2) + 4*a6*lx*ny*pow(d2,2) - 4*ny*px*pow(d2,2) + 4*nx*py*pow(d2,2) - 4*a6*ly*nx*pow(d3,2) + 4*a6*lx*ny*pow(d3,2) - 4*ny*px*pow(d3,2) + 4*nx*py*pow(d3,2) - 4*a6*ly*nx*pow(d4,2) + 4*a6*lx*ny*pow(d4,2) - 4*ny*px*pow(d4,2) + 4*nx*py*pow(d4,2) - 4*a6*ly*nx*pow(d5,2) + 4*a6*lx*ny*pow(d5,2) - 4*ny*px*pow(d5,2) + 4*nx*py*pow(d5,2) + 4*a6*ly*mx*my*ny*pow(d6,2) + 4*a6*ly*mx*mz*nz*pow(d6,2) - a6*lx*my*mz*nz*pow(d6,2) + 4*mx*my*nx*px*pow(d6,2) + 4*my*mz*nz*px*pow(d6,2) - mx*my*ny*py*pow(d6,2) - mx*mz*nz*py*pow(d6,2) + 12*ny*px*pow(a6,2)*pow(lx,2) - 4*nx*py*pow(a6,2)*pow(lx,2) + 4*ly*nx*pow(a6,3)*pow(lx,2) - 4*ny*pow(a6,3)*pow(lx,3) + 4*ny*px*pow(a6,2)*pow(ly,2) - 12*nx*py*pow(a6,2)*pow(ly,2) - 4*lx*ny*pow(a6,3)*pow(ly,2) + 4*nx*pow(a6,3)*pow(ly,3) + 4*ny*px*pow(a6,2)*pow(lz,2) - 4*nx*py*pow(a6,2)*pow(lz,2) + 4*ly*nx*pow(a6,3)*pow(lz,2) - 4*lx*ny*pow(a6,3)*pow(lz,2) + 6*a6*ly*nx*pow(d6,2)*pow(mx,2) - 2*a6*lx*ny*pow(d6,2)*pow(mx,2) + 2*ny*px*pow(d6,2)*pow(mx,2) - 6*nx*py*pow(d6,2)*pow(mx,2) + 2*a6*ly*nx*pow(d6,2)*pow(my,2) - 6*a6*lx*ny*pow(d6,2)*pow(my,2) + 6*ny*px*pow(d6,2)*pow(my,2) - 2*nx*py*pow(d6,2)*pow(my,2) + 2*a6*ly*nx*pow(d6,2)*pow(mz,2) - 2*a6*lx*ny*pow(d6,2)*pow(mz,2) + 2*ny*px*pow(d6,2)*pow(mz,2) - 2*nx*py*pow(d6,2)*pow(mz,2) + 4*a6*ly*nx*pow(px,2) - 12*a6*lx*ny*pow(px,2) - 4*nx*py*pow(px,2) + 4*ny*pow(px,3) + 12*a6*ly*nx*pow(py,2) - 4*a6*lx*ny*pow(py,2) + 4*ny*px*pow(py,2) - 4*nx*pow(py,3) + 4*a6*ly*nx*pow(pz,2) - 4*a6*lx*ny*pow(pz,2) + 4*ny*px*pow(pz,2) - 4*nx*py*pow(pz,2)) + pow(mu1,-1)*(16*a1*a6*d2*lx*nx + 32*a1*d2*d6*mu6*mx*nx + 16*a1*a6*d2*ly*ny + 32*a1*d2*d6*mu6*my*ny - 16*a1*d2*nx*px - 16*a1*d2*ny*py + 16*a1*d3*lambda2*(a6*lx*nx + a6*ly*ny - nx*px - ny*py) + 16*a1*d4*(a6*lx*nx + a6*ly*ny - nx*px - ny*py)*cos(alpha2 + alpha3) + 16*a1*a6*d5*lx*nx*cos(alpha2 + alpha3 - alpha4) + 16*a1*a6*d5*ly*ny*cos(alpha2 + alpha3 - alpha4) - 16*a1*d5*nx*px*cos(alpha2 + alpha3 - alpha4) - 16*a1*d5*ny*py*cos(alpha2 + alpha3 - alpha4) + 32*a1*a6*d6*lx*nx*cos(alpha2 + alpha3 - alpha4 + alpha5) + 32*a1*a6*d6*ly*ny*cos(alpha2 + alpha3 - alpha4 + alpha5) - 32*a1*d6*nx*px*cos(alpha2 + alpha3 - alpha4 + alpha5) - 32*a1*d6*ny*py*cos(alpha2 + alpha3 - alpha4 + alpha5) - 8*a6*d1*d6*ly*mz*nx*cos(alpha1 + alpha6) - 8*a6*d1*d6*ly*mx*nz*cos(alpha1 + alpha6) + 16*a6*d6*ly*mx*nx*px*cos(alpha1 + alpha6) - 16*a6*d6*lx*my*nx*px*cos(alpha1 + alpha6) - 16*a6*d6*lx*mx*ny*px*cos(alpha1 + alpha6) - 16*a6*d6*ly*my*ny*px*cos(alpha1 + alpha6) - 8*d1*d6*mz*ny*px*cos(alpha1 + alpha6) - 8*a6*d6*lz*mz*ny*px*cos(alpha1 + alpha6) - 8*d1*d6*my*nz*px*cos(alpha1 + alpha6) - 8*a6*d6*lz*my*nz*px*cos(alpha1 + alpha6) + 16*a6*d6*lx*mx*nx*py*cos(alpha1 + alpha6) + 16*a6*d6*ly*my*nx*py*cos(alpha1 + alpha6) + 16*a6*d6*ly*mx*ny*py*cos(alpha1 + alpha6) - 16*a6*d6*lx*my*ny*py*cos(alpha1 + alpha6) - 16*d6*mx*nx*px*py*cos(alpha1 + alpha6) + 16*d6*my*ny*px*py*cos(alpha1 + alpha6) - 8*a6*d6*lx*mz*ny*pz*cos(alpha1 + alpha6) - 8*a6*d6*lx*my*nz*pz*cos(alpha1 + alpha6) - 8*d6*mz*nx*py*pz*cos(alpha1 + alpha6) - 8*d6*mx*nz*py*pz*cos(alpha1 + alpha6) - 16*d6*lx*ly*mx*nx*cos(alpha1 + alpha6)*pow(a6,2) - 8*d6*ly*lz*mz*nx*cos(alpha1 + alpha6)*pow(a6,2) + 16*d6*lx*ly*my*ny*cos(alpha1 + alpha6)*pow(a6,2) - 8*d6*ly*lz*mx*nz*cos(alpha1 + alpha6)*pow(a6,2) - 8*d6*my*nx*cos(alpha1 + alpha6)*pow(a6,2)*pow(ly,2) - 8*d6*mx*ny*cos(alpha1 + alpha6)*pow(a6,2)*pow(ly,2) - 8*d6*my*nx*cos(alpha1 + alpha6)*pow(py,2) - 8*d6*mx*ny*cos(alpha1 + alpha6)*pow(py,2) - 2*a6*lx*my*mz*nz*pow(d6,2)*sin(3*alpha1) - 2*mx*my*ny*py*pow(d6,2)*sin(3*alpha1) - 2*mx*mz*nz*py*pow(d6,2)*sin(3*alpha1) - 8*a6*d2*d3*ly*nx*sin(alpha1 - alpha2) + 8*a6*d2*d3*lx*ny*sin(alpha1 - alpha2) - 8*d2*d3*ny*px*sin(alpha1 - alpha2) + 8*d2*d3*nx*py*sin(alpha1 - alpha2) - 8*a6*d2*d3*ly*nx*sin(alpha1 + alpha2) + 8*a6*d2*d3*lx*ny*sin(alpha1 + alpha2) - 8*d2*d3*ny*px*sin(alpha1 + alpha2) + 8*d2*d3*nx*py*sin(alpha1 + alpha2) - 8*a6*d3*d4*ly*nx*sin(alpha1 - alpha3) + 8*a6*d3*d4*lx*ny*sin(alpha1 - alpha3) - 8*d3*d4*ny*px*sin(alpha1 - alpha3) + 8*d3*d4*nx*py*sin(alpha1 - alpha3) - 8*a6*d2*d4*ly*nx*sin(alpha1 - alpha2 - alpha3) + 8*a6*d2*d4*lx*ny*sin(alpha1 - alpha2 - alpha3) - 8*d2*d4*ny*px*sin(alpha1 - alpha2 - alpha3) + 8*d2*d4*nx*py*sin(alpha1 - alpha2 - alpha3) - 8*a6*d3*d4*ly*nx*sin(alpha1 + alpha3) + 8*a6*d3*d4*lx*ny*sin(alpha1 + alpha3) - 8*d3*d4*ny*px*sin(alpha1 + alpha3) + 8*d3*d4*nx*py*sin(alpha1 + alpha3) - 8*a6*d2*d4*ly*nx*sin(alpha1 + alpha2 + alpha3) + 8*a6*d2*d4*lx*ny*sin(alpha1 + alpha2 + alpha3) - 8*d2*d4*ny*px*sin(alpha1 + alpha2 + alpha3) + 8*d2*d4*nx*py*sin(alpha1 + alpha2 + alpha3) - 8*a6*d4*d5*ly*nx*sin(alpha1 - alpha4) + 8*a6*d4*d5*lx*ny*sin(alpha1 - alpha4) - 8*d4*d5*ny*px*sin(alpha1 - alpha4) + 8*d4*d5*nx*py*sin(alpha1 - alpha4) - 8*a6*d3*d5*ly*nx*sin(alpha1 + alpha3 - alpha4) + 8*a6*d3*d5*lx*ny*sin(alpha1 + alpha3 - alpha4) - 8*d3*d5*ny*px*sin(alpha1 + alpha3 - alpha4) + 8*d3*d5*nx*py*sin(alpha1 + alpha3 - alpha4) - 8*a6*d2*d5*ly*nx*sin(alpha1 + alpha2 + alpha3 - alpha4) + 8*a6*d2*d5*lx*ny*sin(alpha1 + alpha2 + alpha3 - alpha4) - 8*d2*d5*ny*px*sin(alpha1 + alpha2 + alpha3 - alpha4) + 8*d2*d5*nx*py*sin(alpha1 + alpha2 + alpha3 - alpha4) - 8*a6*d4*d5*ly*nx*sin(alpha1 + alpha4) + 8*a6*d4*d5*lx*ny*sin(alpha1 + alpha4) - 8*d4*d5*ny*px*sin(alpha1 + alpha4) + 8*d4*d5*nx*py*sin(alpha1 + alpha4) - 8*a6*d3*d5*ly*nx*sin(alpha1 - alpha3 + alpha4) + 8*a6*d3*d5*lx*ny*sin(alpha1 - alpha3 + alpha4) - 8*d3*d5*ny*px*sin(alpha1 - alpha3 + alpha4) + 8*d3*d5*nx*py*sin(alpha1 - alpha3 + alpha4) - 8*a6*d2*d5*ly*nx*sin(alpha1 - alpha2 - alpha3 + alpha4) + 8*a6*d2*d5*lx*ny*sin(alpha1 - alpha2 - alpha3 + alpha4) - 8*d2*d5*ny*px*sin(alpha1 - alpha2 - alpha3 + alpha4) + 8*d2*d5*nx*py*sin(alpha1 - alpha2 - alpha3 + alpha4) - 4*a6*ly*mx*my*ny*pow(d6,2)*sin(alpha1 - 2*alpha6) - 4*a6*ly*mx*mz*nz*pow(d6,2)*sin(alpha1 - 2*alpha6) + a6*lx*my*mz*nz*pow(d6,2)*sin(alpha1 - 2*alpha6) - 4*mx*my*nx*px*pow(d6,2)*sin(alpha1 - 2*alpha6) - 4*my*mz*nz*px*pow(d6,2)*sin(alpha1 - 2*alpha6) + mx*my*ny*py*pow(d6,2)*sin(alpha1 - 2*alpha6) + mx*mz*nz*py*pow(d6,2)*sin(alpha1 - 2*alpha6) - 6*a6*ly*nx*pow(d6,2)*pow(mx,2)*sin(alpha1 - 2*alpha6) + 2*a6*lx*ny*pow(d6,2)*pow(mx,2)*sin(alpha1 - 2*alpha6) - 2*ny*px*pow(d6,2)*pow(mx,2)*sin(alpha1 - 2*alpha6) + 6*nx*py*pow(d6,2)*pow(mx,2)*sin(alpha1 - 2*alpha6) - 2*a6*ly*nx*pow(d6,2)*pow(my,2)*sin(alpha1 - 2*alpha6) + 6*a6*lx*ny*pow(d6,2)*pow(my,2)*sin(alpha1 - 2*alpha6) - 6*ny*px*pow(d6,2)*pow(my,2)*sin(alpha1 - 2*alpha6) + 2*nx*py*pow(d6,2)*pow(my,2)*sin(alpha1 - 2*alpha6) - 2*a6*ly*nx*pow(d6,2)*pow(mz,2)*sin(alpha1 - 2*alpha6) + 2*a6*lx*ny*pow(d6,2)*pow(mz,2)*sin(alpha1 - 2*alpha6) - 2*ny*px*pow(d6,2)*pow(mz,2)*sin(alpha1 - 2*alpha6) + 2*nx*py*pow(d6,2)*pow(mz,2)*sin(alpha1 - 2*alpha6) + a6*lx*my*mz*nz*pow(d6,2)*sin(3*alpha1 - 2*alpha6) + mx*my*ny*py*pow(d6,2)*sin(3*alpha1 - 2*alpha6) + mx*mz*nz*py*pow(d6,2)*sin(3*alpha1 - 2*alpha6) - 16*a1*d1*d6*mx*nx*sin(alpha1 - alpha6) - 16*a1*a6*d6*lz*mx*nx*sin(alpha1 - alpha6) - 16*a1*d1*d6*my*ny*sin(alpha1 - alpha6) - 16*a1*a6*d6*lz*my*ny*sin(alpha1 - alpha6) - 8*a1*d6*mz*nx*px*sin(alpha1 - alpha6) - 8*a1*d6*mx*nz*px*sin(alpha1 - alpha6) - 8*a1*d6*mz*ny*py*sin(alpha1 - alpha6) - 8*a1*d6*my*nz*py*sin(alpha1 - alpha6) + 16*a1*d6*mx*nx*pz*sin(alpha1 - alpha6) + 16*a1*d6*my*ny*pz*sin(alpha1 - alpha6) - 16*a1*d3*d6*mx*nx*sin(alpha2 - alpha6) - 16*a1*d3*d6*my*ny*sin(alpha2 - alpha6) - 16*a1*d4*d6*mx*nx*sin(alpha2 + alpha3 - alpha6) - 16*a1*d4*d6*my*ny*sin(alpha2 + alpha3 - alpha6) - 16*a1*d5*d6*mx*nx*sin(alpha2 + alpha3 - alpha4 - alpha6) - 16*a1*d5*d6*my*ny*sin(alpha2 + alpha3 - alpha4 - alpha6) - 16*a1*mx*nx*pow(d6,2)*sin(alpha2 + alpha3 - alpha4 + alpha5 - alpha6) - 16*a1*my*ny*pow(d6,2)*sin(alpha2 + alpha3 - alpha4 + alpha5 - alpha6) + 16*a1*d1*d6*mx*nx*sin(alpha1 + alpha6) + 16*a1*a6*d6*lz*mx*nx*sin(alpha1 + alpha6) + 16*a1*d1*d6*my*ny*sin(alpha1 + alpha6) + 16*a1*a6*d6*lz*my*ny*sin(alpha1 + alpha6) + 8*a1*d6*mz*nx*px*sin(alpha1 + alpha6) + 8*a1*d6*mx*nz*px*sin(alpha1 + alpha6) + 8*a1*d6*mz*ny*py*sin(alpha1 + alpha6) + 8*a1*d6*my*nz*py*sin(alpha1 + alpha6) - 16*a1*d6*mx*nx*pz*sin(alpha1 + alpha6) - 16*a1*d6*my*ny*pz*sin(alpha1 + alpha6) + 16*a1*d3*d6*mx*nx*sin(alpha2 + alpha6) + 16*a1*d3*d6*my*ny*sin(alpha2 + alpha6) + 16*a1*d4*d6*mx*nx*sin(alpha2 + alpha3 + alpha6) + 16*a1*d4*d6*my*ny*sin(alpha2 + alpha3 + alpha6) + 16*a1*d5*d6*mx*nx*sin(alpha2 + alpha3 - alpha4 + alpha6) + 16*a1*d5*d6*my*ny*sin(alpha2 + alpha3 - alpha4 + alpha6) + 16*a1*mx*nx*pow(d6,2)*sin(alpha2 + alpha3 - alpha4 + alpha5 + alpha6) + 16*a1*my*ny*pow(d6,2)*sin(alpha2 + alpha3 - alpha4 + alpha5 + alpha6) - 4*a6*ly*mx*my*ny*pow(d6,2)*sin(alpha1 + 2*alpha6) - 4*a6*ly*mx*mz*nz*pow(d6,2)*sin(alpha1 + 2*alpha6) + a6*lx*my*mz*nz*pow(d6,2)*sin(alpha1 + 2*alpha6) - 4*mx*my*nx*px*pow(d6,2)*sin(alpha1 + 2*alpha6) - 4*my*mz*nz*px*pow(d6,2)*sin(alpha1 + 2*alpha6) + mx*my*ny*py*pow(d6,2)*sin(alpha1 + 2*alpha6) + mx*mz*nz*py*pow(d6,2)*sin(alpha1 + 2*alpha6) - 6*a6*ly*nx*pow(d6,2)*pow(mx,2)*sin(alpha1 + 2*alpha6) + 2*a6*lx*ny*pow(d6,2)*pow(mx,2)*sin(alpha1 + 2*alpha6) - 2*ny*px*pow(d6,2)*pow(mx,2)*sin(alpha1 + 2*alpha6) + 6*nx*py*pow(d6,2)*pow(mx,2)*sin(alpha1 + 2*alpha6) - 2*a6*ly*nx*pow(d6,2)*pow(my,2)*sin(alpha1 + 2*alpha6) + 6*a6*lx*ny*pow(d6,2)*pow(my,2)*sin(alpha1 + 2*alpha6) - 6*ny*px*pow(d6,2)*pow(my,2)*sin(alpha1 + 2*alpha6) + 2*nx*py*pow(d6,2)*pow(my,2)*sin(alpha1 + 2*alpha6) - 2*a6*ly*nx*pow(d6,2)*pow(mz,2)*sin(alpha1 + 2*alpha6) + 2*a6*lx*ny*pow(d6,2)*pow(mz,2)*sin(alpha1 + 2*alpha6) - 2*ny*px*pow(d6,2)*pow(mz,2)*sin(alpha1 + 2*alpha6) + 2*nx*py*pow(d6,2)*pow(mz,2)*sin(alpha1 + 2*alpha6) + a6*lx*my*mz*nz*pow(d6,2)*sin(3*alpha1 + 2*alpha6) + mx*my*ny*py*pow(d6,2)*sin(3*alpha1 + 2*alpha6) + mx*mz*nz*py*pow(d6,2)*sin(3*alpha1 + 2*alpha6)))))/8.;

return(dummyexp);

}