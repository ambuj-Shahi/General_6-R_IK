
# include <stdio.h>
# include <math.h>

double cx151(double* coeff)
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



double dummyexp = (mu1*pow(a1,-1)*(8*a1*(a3 - a4 + a5)*lambda1*cos(alpha3 - alpha4 + alpha5) - 4*lambda2*(2*a1*a3*mu6*mz - 2*a1*a4*mu6*mz + 2*a1*a5*mu6*mz + 2*a1*a3*lambda6*nz - 2*a1*a4*lambda6*nz + 2*a1*a5*lambda6*nz - d2*d5*cos(alpha1 - alpha5) - d2*d4*cos(alpha1 + alpha4 - alpha5) - d2*d3*cos(alpha1 - alpha3 + alpha4 - alpha5) + d2*d5*cos(alpha1 + alpha5) + d2*d4*cos(alpha1 - alpha4 + alpha5) + d2*d3*cos(alpha1 + alpha3 - alpha4 + alpha5)) - mu1*(-8*d3*d5*mu5 + 4*d4*d6*(mx*nx + my*ny + mz*nz)*cos(alpha3 - 2*alpha6) + 4*d5*d6*(mx*nx + my*ny + mz*nz)*cos(alpha3 - alpha4 - 2*alpha6) + 4*a6*d4*lx*mx*cos(alpha3 - alpha6) + 4*a6*d4*ly*my*cos(alpha3 - alpha6) + 4*d1*d4*mz*cos(alpha3 - alpha6) + 4*a6*d4*lz*mz*cos(alpha3 - alpha6) - 4*d4*mx*px*cos(alpha3 - alpha6) - 4*d4*my*py*cos(alpha3 - alpha6) - 4*d4*mz*pz*cos(alpha3 - alpha6) + 4*a6*d5*lx*mx*cos(alpha3 - alpha4 - alpha6) + 4*a6*d5*ly*my*cos(alpha3 - alpha4 - alpha6) + 4*d1*d5*mz*cos(alpha3 - alpha4 - alpha6) + 4*a6*d5*lz*mz*cos(alpha3 - alpha4 - alpha6) - 4*d5*mx*px*cos(alpha3 - alpha4 - alpha6) - 4*d5*my*py*cos(alpha3 - alpha4 - alpha6) - 4*d5*mz*pz*cos(alpha3 - alpha4 - alpha6) + 4*a6*d6*lx*mx*cos(alpha3 - alpha4 + alpha5 - alpha6) + 4*a6*d6*ly*my*cos(alpha3 - alpha4 + alpha5 - alpha6) + 4*d1*d6*mz*cos(alpha3 - alpha4 + alpha5 - alpha6) + 4*a6*d6*lz*mz*cos(alpha3 - alpha4 + alpha5 - alpha6) - 4*d6*mx*px*cos(alpha3 - alpha4 + alpha5 - alpha6) - 4*d6*my*py*cos(alpha3 - alpha4 + alpha5 - alpha6) - 4*d6*mz*pz*cos(alpha3 - alpha4 + alpha5 - alpha6) - 4*a6*d4*lx*mx*cos(alpha3 + alpha6) - 4*a6*d4*ly*my*cos(alpha3 + alpha6) - 4*d1*d4*mz*cos(alpha3 + alpha6) - 4*a6*d4*lz*mz*cos(alpha3 + alpha6) + 4*d4*mx*px*cos(alpha3 + alpha6) + 4*d4*my*py*cos(alpha3 + alpha6) + 4*d4*mz*pz*cos(alpha3 + alpha6) - 4*a6*d5*lx*mx*cos(alpha3 - alpha4 + alpha6) - 4*a6*d5*ly*my*cos(alpha3 - alpha4 + alpha6) - 4*d1*d5*mz*cos(alpha3 - alpha4 + alpha6) - 4*a6*d5*lz*mz*cos(alpha3 - alpha4 + alpha6) + 4*d5*mx*px*cos(alpha3 - alpha4 + alpha6) + 4*d5*my*py*cos(alpha3 - alpha4 + alpha6) + 4*d5*mz*pz*cos(alpha3 - alpha4 + alpha6) - 4*a6*d6*lx*mx*cos(alpha3 - alpha4 + alpha5 + alpha6) - 4*a6*d6*ly*my*cos(alpha3 - alpha4 + alpha5 + alpha6) - 4*d1*d6*mz*cos(alpha3 - alpha4 + alpha5 + alpha6) - 4*a6*d6*lz*mz*cos(alpha3 - alpha4 + alpha5 + alpha6) + 4*d6*mx*px*cos(alpha3 - alpha4 + alpha5 + alpha6) + 4*d6*my*py*cos(alpha3 - alpha4 + alpha5 + alpha6) + 4*d6*mz*pz*cos(alpha3 - alpha4 + alpha5 + alpha6) - 4*d4*d6*mx*nx*cos(alpha3 + 2*alpha6) - 4*d4*d6*my*ny*cos(alpha3 + 2*alpha6) - 4*d4*d6*mz*nz*cos(alpha3 + 2*alpha6) - 4*d5*d6*mx*nx*cos(alpha3 - alpha4 + 2*alpha6) - 4*d5*d6*my*ny*cos(alpha3 - alpha4 + 2*alpha6) - 4*d5*d6*mz*nz*cos(alpha3 - alpha4 + 2*alpha6) + 2*mx*nx*cos(alpha3 - alpha4 + alpha5 - 2*alpha6)*pow(d6,2) + 2*my*ny*cos(alpha3 - alpha4 + alpha5 - 2*alpha6)*pow(d6,2) + 2*mz*nz*cos(alpha3 - alpha4 + alpha5 - 2*alpha6)*pow(d6,2) - 2*mx*nx*cos(alpha3 - alpha4 + alpha5 + 2*alpha6)*pow(d6,2) - 2*my*ny*cos(alpha3 - alpha4 + alpha5 + 2*alpha6)*pow(d6,2) - 2*mz*nz*cos(alpha3 - alpha4 + alpha5 + 2*alpha6)*pow(d6,2) + 4*d4*d6*mu3*pow(mx,2) + 4*d4*d6*mu3*pow(my,2) + 4*d4*d6*mu3*pow(mz,2) + 4*d4*d6*mu3*pow(nx,2) + 4*d4*d6*mu3*pow(ny,2) + 4*d4*d6*mu3*pow(nz,2) + 4*d5*d6*pow(mx,2)*sin(alpha3 - alpha4) + 4*d5*d6*pow(my,2)*sin(alpha3 - alpha4) + 4*d5*d6*pow(mz,2)*sin(alpha3 - alpha4) + 4*d5*d6*pow(nx,2)*sin(alpha3 - alpha4) + 4*d5*d6*pow(ny,2)*sin(alpha3 - alpha4) + 4*d5*d6*pow(nz,2)*sin(alpha3 - alpha4) + 8*d4*d5*sin(alpha3 - alpha5) + 4*pow(d5,2)*sin(alpha3 - alpha4 - alpha5) + 8*d3*d4*sin(alpha4 - alpha5) + 4*pow(d4,2)*sin(alpha3 + alpha4 - alpha5) + 8*a3*a4*sin(alpha3 - alpha4 + alpha5) - 8*a3*a5*sin(alpha3 - alpha4 + alpha5) + 8*a4*a5*sin(alpha3 - alpha4 + alpha5) + 8*a6*d1*lz*sin(alpha3 - alpha4 + alpha5) - 8*a6*lx*px*sin(alpha3 - alpha4 + alpha5) - 8*a6*ly*py*sin(alpha3 - alpha4 + alpha5) - 8*d1*pz*sin(alpha3 - alpha4 + alpha5) - 8*a6*lz*pz*sin(alpha3 - alpha4 + alpha5) - 4*pow(a1,2)*sin(alpha3 - alpha4 + alpha5) + 4*pow(a2,2)*sin(alpha3 - alpha4 + alpha5) - 4*pow(a3,2)*sin(alpha3 - alpha4 + alpha5) - 4*pow(a4,2)*sin(alpha3 - alpha4 + alpha5) - 4*pow(a5,2)*sin(alpha3 - alpha4 + alpha5) + 4*pow(d1,2)*sin(alpha3 - alpha4 + alpha5) - 4*pow(d2,2)*sin(alpha3 - alpha4 + alpha5) - 4*pow(d3,2)*sin(alpha3 - alpha4 + alpha5) + 4*pow(a6,2)*pow(lx,2)*sin(alpha3 - alpha4 + alpha5) + 4*pow(a6,2)*pow(ly,2)*sin(alpha3 - alpha4 + alpha5) + 4*pow(a6,2)*pow(lz,2)*sin(alpha3 - alpha4 + alpha5) + 2*pow(d6,2)*pow(mx,2)*sin(alpha3 - alpha4 + alpha5) + 2*pow(d6,2)*pow(my,2)*sin(alpha3 - alpha4 + alpha5) + 2*pow(d6,2)*pow(mz,2)*sin(alpha3 - alpha4 + alpha5) + 2*pow(d6,2)*pow(nx,2)*sin(alpha3 - alpha4 + alpha5) + 2*pow(d6,2)*pow(ny,2)*sin(alpha3 - alpha4 + alpha5) + 2*pow(d6,2)*pow(nz,2)*sin(alpha3 - alpha4 + alpha5) + 4*pow(px,2)*sin(alpha3 - alpha4 + alpha5) + 4*pow(py,2)*sin(alpha3 - alpha4 + alpha5) + 4*pow(pz,2)*sin(alpha3 - alpha4 + alpha5) - 2*d4*d6*pow(mx,2)*sin(alpha3 - 2*alpha6) - 2*d4*d6*pow(my,2)*sin(alpha3 - 2*alpha6) - 2*d4*d6*pow(mz,2)*sin(alpha3 - 2*alpha6) + 2*d4*d6*pow(nx,2)*sin(alpha3 - 2*alpha6) + 2*d4*d6*pow(ny,2)*sin(alpha3 - 2*alpha6) + 2*d4*d6*pow(nz,2)*sin(alpha3 - 2*alpha6) - 2*d5*d6*pow(mx,2)*sin(alpha3 - alpha4 - 2*alpha6) - 2*d5*d6*pow(my,2)*sin(alpha3 - alpha4 - 2*alpha6) - 2*d5*d6*pow(mz,2)*sin(alpha3 - alpha4 - 2*alpha6) + 2*d5*d6*pow(nx,2)*sin(alpha3 - alpha4 - 2*alpha6) + 2*d5*d6*pow(ny,2)*sin(alpha3 - alpha4 - 2*alpha6) + 2*d5*d6*pow(nz,2)*sin(alpha3 - alpha4 - 2*alpha6) - pow(d6,2)*pow(mx,2)*sin(alpha3 - alpha4 + alpha5 - 2*alpha6) - pow(d6,2)*pow(my,2)*sin(alpha3 - alpha4 + alpha5 - 2*alpha6) - pow(d6,2)*pow(mz,2)*sin(alpha3 - alpha4 + alpha5 - 2*alpha6) + pow(d6,2)*pow(nx,2)*sin(alpha3 - alpha4 + alpha5 - 2*alpha6) + pow(d6,2)*pow(ny,2)*sin(alpha3 - alpha4 + alpha5 - 2*alpha6) + pow(d6,2)*pow(nz,2)*sin(alpha3 - alpha4 + alpha5 - 2*alpha6) + 4*a6*d4*lx*nx*sin(alpha3 - alpha6) + 4*a6*d4*ly*ny*sin(alpha3 - alpha6) + 4*d1*d4*nz*sin(alpha3 - alpha6) + 4*a6*d4*lz*nz*sin(alpha3 - alpha6) - 4*d4*nx*px*sin(alpha3 - alpha6) - 4*d4*ny*py*sin(alpha3 - alpha6) - 4*d4*nz*pz*sin(alpha3 - alpha6) + 4*a6*d5*lx*nx*sin(alpha3 - alpha4 - alpha6) + 4*a6*d5*ly*ny*sin(alpha3 - alpha4 - alpha6) + 4*d1*d5*nz*sin(alpha3 - alpha4 - alpha6) + 4*a6*d5*lz*nz*sin(alpha3 - alpha4 - alpha6) - 4*d5*nx*px*sin(alpha3 - alpha4 - alpha6) - 4*d5*ny*py*sin(alpha3 - alpha4 - alpha6) - 4*d5*nz*pz*sin(alpha3 - alpha4 - alpha6) + 4*a6*d6*lx*nx*sin(alpha3 - alpha4 + alpha5 - alpha6) + 4*a6*d6*ly*ny*sin(alpha3 - alpha4 + alpha5 - alpha6) + 4*d1*d6*nz*sin(alpha3 - alpha4 + alpha5 - alpha6) + 4*a6*d6*lz*nz*sin(alpha3 - alpha4 + alpha5 - alpha6) - 4*d6*nx*px*sin(alpha3 - alpha4 + alpha5 - alpha6) - 4*d6*ny*py*sin(alpha3 - alpha4 + alpha5 - alpha6) - 4*d6*nz*pz*sin(alpha3 - alpha4 + alpha5 - alpha6) + 4*a6*d4*lx*nx*sin(alpha3 + alpha6) + 4*a6*d4*ly*ny*sin(alpha3 + alpha6) + 4*d1*d4*nz*sin(alpha3 + alpha6) + 4*a6*d4*lz*nz*sin(alpha3 + alpha6) - 4*d4*nx*px*sin(alpha3 + alpha6) - 4*d4*ny*py*sin(alpha3 + alpha6) - 4*d4*nz*pz*sin(alpha3 + alpha6) + 4*a6*d5*lx*nx*sin(alpha3 - alpha4 + alpha6) + 4*a6*d5*ly*ny*sin(alpha3 - alpha4 + alpha6) + 4*d1*d5*nz*sin(alpha3 - alpha4 + alpha6) + 4*a6*d5*lz*nz*sin(alpha3 - alpha4 + alpha6) - 4*d5*nx*px*sin(alpha3 - alpha4 + alpha6) - 4*d5*ny*py*sin(alpha3 - alpha4 + alpha6) - 4*d5*nz*pz*sin(alpha3 - alpha4 + alpha6) + 4*a6*d6*lx*nx*sin(alpha3 - alpha4 + alpha5 + alpha6) + 4*a6*d6*ly*ny*sin(alpha3 - alpha4 + alpha5 + alpha6) + 4*d1*d6*nz*sin(alpha3 - alpha4 + alpha5 + alpha6) + 4*a6*d6*lz*nz*sin(alpha3 - alpha4 + alpha5 + alpha6) - 4*d6*nx*px*sin(alpha3 - alpha4 + alpha5 + alpha6) - 4*d6*ny*py*sin(alpha3 - alpha4 + alpha5 + alpha6) - 4*d6*nz*pz*sin(alpha3 - alpha4 + alpha5 + alpha6) - 2*d4*d6*pow(mx,2)*sin(alpha3 + 2*alpha6) - 2*d4*d6*pow(my,2)*sin(alpha3 + 2*alpha6) - 2*d4*d6*pow(mz,2)*sin(alpha3 + 2*alpha6) + 2*d4*d6*pow(nx,2)*sin(alpha3 + 2*alpha6) + 2*d4*d6*pow(ny,2)*sin(alpha3 + 2*alpha6) + 2*d4*d6*pow(nz,2)*sin(alpha3 + 2*alpha6) - 2*d5*d6*pow(mx,2)*sin(alpha3 - alpha4 + 2*alpha6) - 2*d5*d6*pow(my,2)*sin(alpha3 - alpha4 + 2*alpha6) - 2*d5*d6*pow(mz,2)*sin(alpha3 - alpha4 + 2*alpha6) + 2*d5*d6*pow(nx,2)*sin(alpha3 - alpha4 + 2*alpha6) + 2*d5*d6*pow(ny,2)*sin(alpha3 - alpha4 + 2*alpha6) + 2*d5*d6*pow(nz,2)*sin(alpha3 - alpha4 + 2*alpha6) - pow(d6,2)*pow(mx,2)*sin(alpha3 - alpha4 + alpha5 + 2*alpha6) - pow(d6,2)*pow(my,2)*sin(alpha3 - alpha4 + alpha5 + 2*alpha6) - pow(d6,2)*pow(mz,2)*sin(alpha3 - alpha4 + alpha5 + 2*alpha6) + pow(d6,2)*pow(nx,2)*sin(alpha3 - alpha4 + alpha5 + 2*alpha6) + pow(d6,2)*pow(ny,2)*sin(alpha3 - alpha4 + alpha5 + 2*alpha6) + pow(d6,2)*pow(nz,2)*sin(alpha3 - alpha4 + alpha5 + 2*alpha6))))/4.;

return(dummyexp);

}