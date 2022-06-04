
# include <stdio.h>
# include <math.h>

double cx057(double* coeff)
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



double dummyexp = -(mu1*pow(a1,-1)*(2*(a2 + a3 + a4 - a5)*d6*(mx*nx + my*ny + mz*nz)*cos(alpha1 - 2*alpha6) + 2*(a2 + a3 + a4 - a5)*(a6*lx*mx + a6*ly*my + d1*mz + a6*lz*mz - mx*px - my*py - mz*pz)*cos(alpha1 - alpha6) - 2*a1*d3*mz*cos(alpha2 - alpha6) - 2*a1*d4*mz*cos(alpha2 + alpha3 - alpha6) - 2*a1*d5*mz*cos(alpha2 + alpha3 + alpha4 - alpha6) - 2*a1*d6*mz*cos(alpha2 + alpha3 + alpha4 - alpha5 - alpha6) - 2*a2*a6*lx*mx*cos(alpha1 + alpha6) - 2*a3*a6*lx*mx*cos(alpha1 + alpha6) - 2*a4*a6*lx*mx*cos(alpha1 + alpha6) + 2*a5*a6*lx*mx*cos(alpha1 + alpha6) - 2*a2*a6*ly*my*cos(alpha1 + alpha6) - 2*a3*a6*ly*my*cos(alpha1 + alpha6) - 2*a4*a6*ly*my*cos(alpha1 + alpha6) + 2*a5*a6*ly*my*cos(alpha1 + alpha6) - 2*a2*d1*mz*cos(alpha1 + alpha6) - 2*a3*d1*mz*cos(alpha1 + alpha6) - 2*a4*d1*mz*cos(alpha1 + alpha6) + 2*a5*d1*mz*cos(alpha1 + alpha6) - 2*a2*a6*lz*mz*cos(alpha1 + alpha6) - 2*a3*a6*lz*mz*cos(alpha1 + alpha6) - 2*a4*a6*lz*mz*cos(alpha1 + alpha6) + 2*a5*a6*lz*mz*cos(alpha1 + alpha6) + 2*a2*mx*px*cos(alpha1 + alpha6) + 2*a3*mx*px*cos(alpha1 + alpha6) + 2*a4*mx*px*cos(alpha1 + alpha6) - 2*a5*mx*px*cos(alpha1 + alpha6) + 2*a2*my*py*cos(alpha1 + alpha6) + 2*a3*my*py*cos(alpha1 + alpha6) + 2*a4*my*py*cos(alpha1 + alpha6) - 2*a5*my*py*cos(alpha1 + alpha6) + 2*a2*mz*pz*cos(alpha1 + alpha6) + 2*a3*mz*pz*cos(alpha1 + alpha6) + 2*a4*mz*pz*cos(alpha1 + alpha6) - 2*a5*mz*pz*cos(alpha1 + alpha6) + 2*a1*d3*mz*cos(alpha2 + alpha6) + 2*a1*d4*mz*cos(alpha2 + alpha3 + alpha6) + 2*a1*d5*mz*cos(alpha2 + alpha3 + alpha4 + alpha6) + 2*a1*d6*mz*cos(alpha2 + alpha3 + alpha4 - alpha5 + alpha6) - 2*a2*d6*mx*nx*cos(alpha1 + 2*alpha6) - 2*a3*d6*mx*nx*cos(alpha1 + 2*alpha6) - 2*a4*d6*mx*nx*cos(alpha1 + 2*alpha6) + 2*a5*d6*mx*nx*cos(alpha1 + 2*alpha6) - 2*a2*d6*my*ny*cos(alpha1 + 2*alpha6) - 2*a3*d6*my*ny*cos(alpha1 + 2*alpha6) - 2*a4*d6*my*ny*cos(alpha1 + 2*alpha6) + 2*a5*d6*my*ny*cos(alpha1 + 2*alpha6) - 2*a2*d6*mz*nz*cos(alpha1 + 2*alpha6) - 2*a3*d6*mz*nz*cos(alpha1 + 2*alpha6) - 2*a4*d6*mz*nz*cos(alpha1 + 2*alpha6) + 2*a5*d6*mz*nz*cos(alpha1 + 2*alpha6) + 2*a2*d6*mu1*pow(mx,2) + 2*a3*d6*mu1*pow(mx,2) + 2*a4*d6*mu1*pow(mx,2) - 2*a5*d6*mu1*pow(mx,2) + 2*a2*d6*mu1*pow(my,2) + 2*a3*d6*mu1*pow(my,2) + 2*a4*d6*mu1*pow(my,2) - 2*a5*d6*mu1*pow(my,2) + 2*a2*d6*mu1*pow(mz,2) + 2*a3*d6*mu1*pow(mz,2) + 2*a4*d6*mu1*pow(mz,2) - 2*a5*d6*mu1*pow(mz,2) + 2*a2*d6*mu1*pow(nx,2) + 2*a3*d6*mu1*pow(nx,2) + 2*a4*d6*mu1*pow(nx,2) - 2*a5*d6*mu1*pow(nx,2) + 2*a2*d6*mu1*pow(ny,2) + 2*a3*d6*mu1*pow(ny,2) + 2*a4*d6*mu1*pow(ny,2) - 2*a5*d6*mu1*pow(ny,2) + 2*a2*d6*mu1*pow(nz,2) + 2*a3*d6*mu1*pow(nz,2) + 2*a4*d6*mu1*pow(nz,2) - 2*a5*d6*mu1*pow(nz,2) - 2*a1*d5*sin(alpha1 - alpha5) + 2*a2*d5*sin(alpha1 - alpha5) + 2*a3*d5*sin(alpha1 - alpha5) + 2*a4*d5*sin(alpha1 - alpha5) - 2*a5*d5*sin(alpha1 - alpha5) - 2*a1*d4*sin(alpha1 + alpha4 - alpha5) + 2*a2*d4*sin(alpha1 + alpha4 - alpha5) + 2*a3*d4*sin(alpha1 + alpha4 - alpha5) + 2*a4*d4*sin(alpha1 + alpha4 - alpha5) - 2*a5*d4*sin(alpha1 + alpha4 - alpha5) - 2*a1*d3*sin(alpha1 + alpha3 + alpha4 - alpha5) + 2*a2*d3*sin(alpha1 + alpha3 + alpha4 - alpha5) + 2*a3*d3*sin(alpha1 + alpha3 + alpha4 - alpha5) + 2*a4*d3*sin(alpha1 + alpha3 + alpha4 - alpha5) - 2*a5*d3*sin(alpha1 + alpha3 + alpha4 - alpha5) - 4*a1*d1*sin(alpha2 + alpha3 + alpha4 - alpha5) - 4*a1*a6*lz*sin(alpha2 + alpha3 + alpha4 - alpha5) + 4*a1*pz*sin(alpha2 + alpha3 + alpha4 - alpha5) - 2*a1*d2*sin(alpha1 + alpha2 + alpha3 + alpha4 - alpha5) + 2*a2*d2*sin(alpha1 + alpha2 + alpha3 + alpha4 - alpha5) + 2*a3*d2*sin(alpha1 + alpha2 + alpha3 + alpha4 - alpha5) + 2*a4*d2*sin(alpha1 + alpha2 + alpha3 + alpha4 - alpha5) - 2*a5*d2*sin(alpha1 + alpha2 + alpha3 + alpha4 - alpha5) + 2*a1*d5*sin(alpha1 + alpha5) + 2*a2*d5*sin(alpha1 + alpha5) + 2*a3*d5*sin(alpha1 + alpha5) + 2*a4*d5*sin(alpha1 + alpha5) - 2*a5*d5*sin(alpha1 + alpha5) + 2*a1*d4*sin(alpha1 - alpha4 + alpha5) + 2*a2*d4*sin(alpha1 - alpha4 + alpha5) + 2*a3*d4*sin(alpha1 - alpha4 + alpha5) + 2*a4*d4*sin(alpha1 - alpha4 + alpha5) - 2*a5*d4*sin(alpha1 - alpha4 + alpha5) + 2*a1*d3*sin(alpha1 - alpha3 - alpha4 + alpha5) + 2*a2*d3*sin(alpha1 - alpha3 - alpha4 + alpha5) + 2*a3*d3*sin(alpha1 - alpha3 - alpha4 + alpha5) + 2*a4*d3*sin(alpha1 - alpha3 - alpha4 + alpha5) - 2*a5*d3*sin(alpha1 - alpha3 - alpha4 + alpha5) + 2*a1*d2*sin(alpha1 - alpha2 - alpha3 - alpha4 + alpha5) + 2*a2*d2*sin(alpha1 - alpha2 - alpha3 - alpha4 + alpha5) + 2*a3*d2*sin(alpha1 - alpha2 - alpha3 - alpha4 + alpha5) + 2*a4*d2*sin(alpha1 - alpha2 - alpha3 - alpha4 + alpha5) - 2*a5*d2*sin(alpha1 - alpha2 - alpha3 - alpha4 + alpha5) - a2*d6*pow(mx,2)*sin(alpha1 - 2*alpha6) - a3*d6*pow(mx,2)*sin(alpha1 - 2*alpha6) - a4*d6*pow(mx,2)*sin(alpha1 - 2*alpha6) + a5*d6*pow(mx,2)*sin(alpha1 - 2*alpha6) - a2*d6*pow(my,2)*sin(alpha1 - 2*alpha6) - a3*d6*pow(my,2)*sin(alpha1 - 2*alpha6) - a4*d6*pow(my,2)*sin(alpha1 - 2*alpha6) + a5*d6*pow(my,2)*sin(alpha1 - 2*alpha6) - a2*d6*pow(mz,2)*sin(alpha1 - 2*alpha6) - a3*d6*pow(mz,2)*sin(alpha1 - 2*alpha6) - a4*d6*pow(mz,2)*sin(alpha1 - 2*alpha6) + a5*d6*pow(mz,2)*sin(alpha1 - 2*alpha6) + a2*d6*pow(nx,2)*sin(alpha1 - 2*alpha6) + a3*d6*pow(nx,2)*sin(alpha1 - 2*alpha6) + a4*d6*pow(nx,2)*sin(alpha1 - 2*alpha6) - a5*d6*pow(nx,2)*sin(alpha1 - 2*alpha6) + a2*d6*pow(ny,2)*sin(alpha1 - 2*alpha6) + a3*d6*pow(ny,2)*sin(alpha1 - 2*alpha6) + a4*d6*pow(ny,2)*sin(alpha1 - 2*alpha6) - a5*d6*pow(ny,2)*sin(alpha1 - 2*alpha6) + a2*d6*pow(nz,2)*sin(alpha1 - 2*alpha6) + a3*d6*pow(nz,2)*sin(alpha1 - 2*alpha6) + a4*d6*pow(nz,2)*sin(alpha1 - 2*alpha6) - a5*d6*pow(nz,2)*sin(alpha1 - 2*alpha6) + 2*a2*a6*lx*nx*sin(alpha1 - alpha6) + 2*a3*a6*lx*nx*sin(alpha1 - alpha6) + 2*a4*a6*lx*nx*sin(alpha1 - alpha6) - 2*a5*a6*lx*nx*sin(alpha1 - alpha6) + 2*a2*a6*ly*ny*sin(alpha1 - alpha6) + 2*a3*a6*ly*ny*sin(alpha1 - alpha6) + 2*a4*a6*ly*ny*sin(alpha1 - alpha6) - 2*a5*a6*ly*ny*sin(alpha1 - alpha6) + 2*a2*d1*nz*sin(alpha1 - alpha6) + 2*a3*d1*nz*sin(alpha1 - alpha6) + 2*a4*d1*nz*sin(alpha1 - alpha6) - 2*a5*d1*nz*sin(alpha1 - alpha6) + 2*a2*a6*lz*nz*sin(alpha1 - alpha6) + 2*a3*a6*lz*nz*sin(alpha1 - alpha6) + 2*a4*a6*lz*nz*sin(alpha1 - alpha6) - 2*a5*a6*lz*nz*sin(alpha1 - alpha6) - 2*a2*nx*px*sin(alpha1 - alpha6) - 2*a3*nx*px*sin(alpha1 - alpha6) - 2*a4*nx*px*sin(alpha1 - alpha6) + 2*a5*nx*px*sin(alpha1 - alpha6) - 2*a2*ny*py*sin(alpha1 - alpha6) - 2*a3*ny*py*sin(alpha1 - alpha6) - 2*a4*ny*py*sin(alpha1 - alpha6) + 2*a5*ny*py*sin(alpha1 - alpha6) - 2*a2*nz*pz*sin(alpha1 - alpha6) - 2*a3*nz*pz*sin(alpha1 - alpha6) - 2*a4*nz*pz*sin(alpha1 - alpha6) + 2*a5*nz*pz*sin(alpha1 - alpha6) - 2*a1*d3*nz*sin(alpha2 - alpha6) - 2*a1*d4*nz*sin(alpha2 + alpha3 - alpha6) - 2*a1*d5*nz*sin(alpha2 + alpha3 + alpha4 - alpha6) - 2*a1*d6*nz*sin(alpha2 + alpha3 + alpha4 - alpha5 - alpha6) + 2*a2*a6*lx*nx*sin(alpha1 + alpha6) + 2*a3*a6*lx*nx*sin(alpha1 + alpha6) + 2*a4*a6*lx*nx*sin(alpha1 + alpha6) - 2*a5*a6*lx*nx*sin(alpha1 + alpha6) + 2*a2*a6*ly*ny*sin(alpha1 + alpha6) + 2*a3*a6*ly*ny*sin(alpha1 + alpha6) + 2*a4*a6*ly*ny*sin(alpha1 + alpha6) - 2*a5*a6*ly*ny*sin(alpha1 + alpha6) + 2*a2*d1*nz*sin(alpha1 + alpha6) + 2*a3*d1*nz*sin(alpha1 + alpha6) + 2*a4*d1*nz*sin(alpha1 + alpha6) - 2*a5*d1*nz*sin(alpha1 + alpha6) + 2*a2*a6*lz*nz*sin(alpha1 + alpha6) + 2*a3*a6*lz*nz*sin(alpha1 + alpha6) + 2*a4*a6*lz*nz*sin(alpha1 + alpha6) - 2*a5*a6*lz*nz*sin(alpha1 + alpha6) - 2*a2*nx*px*sin(alpha1 + alpha6) - 2*a3*nx*px*sin(alpha1 + alpha6) - 2*a4*nx*px*sin(alpha1 + alpha6) + 2*a5*nx*px*sin(alpha1 + alpha6) - 2*a2*ny*py*sin(alpha1 + alpha6) - 2*a3*ny*py*sin(alpha1 + alpha6) - 2*a4*ny*py*sin(alpha1 + alpha6) + 2*a5*ny*py*sin(alpha1 + alpha6) - 2*a2*nz*pz*sin(alpha1 + alpha6) - 2*a3*nz*pz*sin(alpha1 + alpha6) - 2*a4*nz*pz*sin(alpha1 + alpha6) + 2*a5*nz*pz*sin(alpha1 + alpha6) - 2*a1*d3*nz*sin(alpha2 + alpha6) - 2*a1*d4*nz*sin(alpha2 + alpha3 + alpha6) - 2*a1*d5*nz*sin(alpha2 + alpha3 + alpha4 + alpha6) - 2*a1*d6*nz*sin(alpha2 + alpha3 + alpha4 - alpha5 + alpha6) - a2*d6*pow(mx,2)*sin(alpha1 + 2*alpha6) - a3*d6*pow(mx,2)*sin(alpha1 + 2*alpha6) - a4*d6*pow(mx,2)*sin(alpha1 + 2*alpha6) + a5*d6*pow(mx,2)*sin(alpha1 + 2*alpha6) - a2*d6*pow(my,2)*sin(alpha1 + 2*alpha6) - a3*d6*pow(my,2)*sin(alpha1 + 2*alpha6) - a4*d6*pow(my,2)*sin(alpha1 + 2*alpha6) + a5*d6*pow(my,2)*sin(alpha1 + 2*alpha6) - a2*d6*pow(mz,2)*sin(alpha1 + 2*alpha6) - a3*d6*pow(mz,2)*sin(alpha1 + 2*alpha6) - a4*d6*pow(mz,2)*sin(alpha1 + 2*alpha6) + a5*d6*pow(mz,2)*sin(alpha1 + 2*alpha6) + a2*d6*pow(nx,2)*sin(alpha1 + 2*alpha6) + a3*d6*pow(nx,2)*sin(alpha1 + 2*alpha6) + a4*d6*pow(nx,2)*sin(alpha1 + 2*alpha6) - a5*d6*pow(nx,2)*sin(alpha1 + 2*alpha6) + a2*d6*pow(ny,2)*sin(alpha1 + 2*alpha6) + a3*d6*pow(ny,2)*sin(alpha1 + 2*alpha6) + a4*d6*pow(ny,2)*sin(alpha1 + 2*alpha6) - a5*d6*pow(ny,2)*sin(alpha1 + 2*alpha6) + a2*d6*pow(nz,2)*sin(alpha1 + 2*alpha6) + a3*d6*pow(nz,2)*sin(alpha1 + 2*alpha6) + a4*d6*pow(nz,2)*sin(alpha1 + 2*alpha6) - a5*d6*pow(nz,2)*sin(alpha1 + 2*alpha6)))/4.;

return(dummyexp);

}