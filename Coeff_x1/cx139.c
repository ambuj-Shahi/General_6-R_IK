
# include <stdio.h>
# include <math.h>

double cx139(double* coeff)
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



double dummyexp = mu2*pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1)*(-(a3*a6*lx*mx*mz) - a4*a6*lx*mx*mz - a5*a6*lx*mx*mz - a3*a6*ly*my*mz - a4*a6*ly*my*mz - a5*a6*ly*my*mz - a3*a6*lx*nx*nz - a4*a6*lx*nx*nz - a5*a6*lx*nx*nz - a3*a6*ly*ny*nz - a4*a6*ly*ny*nz - a5*a6*ly*ny*nz + a3*mx*mz*px + a4*mx*mz*px + a5*mx*mz*px + a3*nx*nz*px + a4*nx*nz*px + a5*nx*nz*px + a3*my*mz*py + a4*my*mz*py + a5*my*mz*py + a3*ny*nz*py + a4*ny*nz*py + a5*ny*nz*py - d5*(a6*ly*mx - a6*lx*my + my*px - mx*py)*cos(alpha5 - alpha6) - d4*(a6*ly*mx - a6*lx*my + my*px - mx*py)*cos(alpha4 + alpha5 - alpha6) - a6*d3*ly*mx*cos(alpha3 + alpha4 + alpha5 - alpha6) + a6*d3*lx*my*cos(alpha3 + alpha4 + alpha5 - alpha6) - d3*my*px*cos(alpha3 + alpha4 + alpha5 - alpha6) + d3*mx*py*cos(alpha3 + alpha4 + alpha5 - alpha6) + a3*a6*lx*mx*mz*cos(2*alpha6) + a4*a6*lx*mx*mz*cos(2*alpha6) + a5*a6*lx*mx*mz*cos(2*alpha6) + a3*a6*ly*my*mz*cos(2*alpha6) + a4*a6*ly*my*mz*cos(2*alpha6) + a5*a6*ly*my*mz*cos(2*alpha6) - a3*a6*lx*nx*nz*cos(2*alpha6) - a4*a6*lx*nx*nz*cos(2*alpha6) - a5*a6*lx*nx*nz*cos(2*alpha6) - a3*a6*ly*ny*nz*cos(2*alpha6) - a4*a6*ly*ny*nz*cos(2*alpha6) - a5*a6*ly*ny*nz*cos(2*alpha6) - a3*mx*mz*px*cos(2*alpha6) - a4*mx*mz*px*cos(2*alpha6) - a5*mx*mz*px*cos(2*alpha6) + a3*nx*nz*px*cos(2*alpha6) + a4*nx*nz*px*cos(2*alpha6) + a5*nx*nz*px*cos(2*alpha6) - a3*my*mz*py*cos(2*alpha6) - a4*my*mz*py*cos(2*alpha6) - a5*my*mz*py*cos(2*alpha6) + a3*ny*nz*py*cos(2*alpha6) + a4*ny*nz*py*cos(2*alpha6) + a5*ny*nz*py*cos(2*alpha6) + a6*d5*ly*mx*cos(alpha5 + alpha6) - a6*d5*lx*my*cos(alpha5 + alpha6) + d5*my*px*cos(alpha5 + alpha6) - d5*mx*py*cos(alpha5 + alpha6) + a6*d4*ly*mx*cos(alpha4 + alpha5 + alpha6) - a6*d4*lx*my*cos(alpha4 + alpha5 + alpha6) + d4*my*px*cos(alpha4 + alpha5 + alpha6) - d4*mx*py*cos(alpha4 + alpha5 + alpha6) + a6*d3*ly*mx*cos(alpha3 + alpha4 + alpha5 + alpha6) - a6*d3*lx*my*cos(alpha3 + alpha4 + alpha5 + alpha6) + d3*my*px*cos(alpha3 + alpha4 + alpha5 + alpha6) - d3*mx*py*cos(alpha3 + alpha4 + alpha5 + alpha6) + a3*d1*pow(mx,2) + a4*d1*pow(mx,2) + a5*d1*pow(mx,2) + a3*a6*lz*pow(mx,2) + a4*a6*lz*pow(mx,2) + a5*a6*lz*pow(mx,2) - a3*pz*pow(mx,2) - a4*pz*pow(mx,2) - a5*pz*pow(mx,2) - a3*d1*cos(2*alpha6)*pow(mx,2) - a4*d1*cos(2*alpha6)*pow(mx,2) - a5*d1*cos(2*alpha6)*pow(mx,2) - a3*a6*lz*cos(2*alpha6)*pow(mx,2) - a4*a6*lz*cos(2*alpha6)*pow(mx,2) - a5*a6*lz*cos(2*alpha6)*pow(mx,2) + a3*pz*cos(2*alpha6)*pow(mx,2) + a4*pz*cos(2*alpha6)*pow(mx,2) + a5*pz*cos(2*alpha6)*pow(mx,2) + a3*d1*pow(my,2) + a4*d1*pow(my,2) + a5*d1*pow(my,2) + a3*a6*lz*pow(my,2) + a4*a6*lz*pow(my,2) + a5*a6*lz*pow(my,2) - a3*pz*pow(my,2) - a4*pz*pow(my,2) - a5*pz*pow(my,2) - a3*d1*cos(2*alpha6)*pow(my,2) - a4*d1*cos(2*alpha6)*pow(my,2) - a5*d1*cos(2*alpha6)*pow(my,2) - a3*a6*lz*cos(2*alpha6)*pow(my,2) - a4*a6*lz*cos(2*alpha6)*pow(my,2) - a5*a6*lz*cos(2*alpha6)*pow(my,2) + a3*pz*cos(2*alpha6)*pow(my,2) + a4*pz*cos(2*alpha6)*pow(my,2) + a5*pz*cos(2*alpha6)*pow(my,2) + a3*d1*pow(nx,2) + a4*d1*pow(nx,2) + a5*d1*pow(nx,2) + a3*a6*lz*pow(nx,2) + a4*a6*lz*pow(nx,2) + a5*a6*lz*pow(nx,2) - a3*pz*pow(nx,2) - a4*pz*pow(nx,2) - a5*pz*pow(nx,2) + a3*d1*cos(2*alpha6)*pow(nx,2) + a4*d1*cos(2*alpha6)*pow(nx,2) + a5*d1*cos(2*alpha6)*pow(nx,2) + a3*a6*lz*cos(2*alpha6)*pow(nx,2) + a4*a6*lz*cos(2*alpha6)*pow(nx,2) + a5*a6*lz*cos(2*alpha6)*pow(nx,2) - a3*pz*cos(2*alpha6)*pow(nx,2) - a4*pz*cos(2*alpha6)*pow(nx,2) - a5*pz*cos(2*alpha6)*pow(nx,2) + a3*d1*pow(ny,2) + a4*d1*pow(ny,2) + a5*d1*pow(ny,2) + a3*a6*lz*pow(ny,2) + a4*a6*lz*pow(ny,2) + a5*a6*lz*pow(ny,2) - a3*pz*pow(ny,2) - a4*pz*pow(ny,2) - a5*pz*pow(ny,2) + a3*d1*cos(2*alpha6)*pow(ny,2) + a4*d1*cos(2*alpha6)*pow(ny,2) + a5*d1*cos(2*alpha6)*pow(ny,2) + a3*a6*lz*cos(2*alpha6)*pow(ny,2) + a4*a6*lz*cos(2*alpha6)*pow(ny,2) + a5*a6*lz*cos(2*alpha6)*pow(ny,2) - a3*pz*cos(2*alpha6)*pow(ny,2) - a4*pz*cos(2*alpha6)*pow(ny,2) - a5*pz*cos(2*alpha6)*pow(ny,2) - a6*d5*ly*nx*sin(alpha5 - alpha6) + a6*d5*lx*ny*sin(alpha5 - alpha6) - d5*ny*px*sin(alpha5 - alpha6) + d5*nx*py*sin(alpha5 - alpha6) - a6*d4*ly*nx*sin(alpha4 + alpha5 - alpha6) + a6*d4*lx*ny*sin(alpha4 + alpha5 - alpha6) - d4*ny*px*sin(alpha4 + alpha5 - alpha6) + d4*nx*py*sin(alpha4 + alpha5 - alpha6) - a6*d3*ly*nx*sin(alpha3 + alpha4 + alpha5 - alpha6) + a6*d3*lx*ny*sin(alpha3 + alpha4 + alpha5 - alpha6) - d3*ny*px*sin(alpha3 + alpha4 + alpha5 - alpha6) + d3*nx*py*sin(alpha3 + alpha4 + alpha5 - alpha6) + 2*a3*d1*mx*nx*sin(2*alpha6) + 2*a4*d1*mx*nx*sin(2*alpha6) + 2*a5*d1*mx*nx*sin(2*alpha6) + 2*a3*a6*lz*mx*nx*sin(2*alpha6) + 2*a4*a6*lz*mx*nx*sin(2*alpha6) + 2*a5*a6*lz*mx*nx*sin(2*alpha6) - a3*a6*lx*mz*nx*sin(2*alpha6) - a4*a6*lx*mz*nx*sin(2*alpha6) - a5*a6*lx*mz*nx*sin(2*alpha6) + 2*a3*d1*my*ny*sin(2*alpha6) + 2*a4*d1*my*ny*sin(2*alpha6) + 2*a5*d1*my*ny*sin(2*alpha6) + 2*a3*a6*lz*my*ny*sin(2*alpha6) + 2*a4*a6*lz*my*ny*sin(2*alpha6) + 2*a5*a6*lz*my*ny*sin(2*alpha6) - a3*a6*ly*mz*ny*sin(2*alpha6) - a4*a6*ly*mz*ny*sin(2*alpha6) - a5*a6*ly*mz*ny*sin(2*alpha6) - a3*a6*lx*mx*nz*sin(2*alpha6) - a4*a6*lx*mx*nz*sin(2*alpha6) - a5*a6*lx*mx*nz*sin(2*alpha6) - a3*a6*ly*my*nz*sin(2*alpha6) - a4*a6*ly*my*nz*sin(2*alpha6) - a5*a6*ly*my*nz*sin(2*alpha6) + a3*mz*nx*px*sin(2*alpha6) + a4*mz*nx*px*sin(2*alpha6) + a5*mz*nx*px*sin(2*alpha6) + a3*mx*nz*px*sin(2*alpha6) + a4*mx*nz*px*sin(2*alpha6) + a5*mx*nz*px*sin(2*alpha6) + a3*mz*ny*py*sin(2*alpha6) + a4*mz*ny*py*sin(2*alpha6) + a5*mz*ny*py*sin(2*alpha6) + a3*my*nz*py*sin(2*alpha6) + a4*my*nz*py*sin(2*alpha6) + a5*my*nz*py*sin(2*alpha6) - 2*a3*mx*nx*pz*sin(2*alpha6) - 2*a4*mx*nx*pz*sin(2*alpha6) - 2*a5*mx*nx*pz*sin(2*alpha6) - 2*a3*my*ny*pz*sin(2*alpha6) - 2*a4*my*ny*pz*sin(2*alpha6) - 2*a5*my*ny*pz*sin(2*alpha6) - a6*d5*ly*nx*sin(alpha5 + alpha6) + a6*d5*lx*ny*sin(alpha5 + alpha6) - d5*ny*px*sin(alpha5 + alpha6) + d5*nx*py*sin(alpha5 + alpha6) - a6*d4*ly*nx*sin(alpha4 + alpha5 + alpha6) + a6*d4*lx*ny*sin(alpha4 + alpha5 + alpha6) - d4*ny*px*sin(alpha4 + alpha5 + alpha6) + d4*nx*py*sin(alpha4 + alpha5 + alpha6) - a6*d3*ly*nx*sin(alpha3 + alpha4 + alpha5 + alpha6) + a6*d3*lx*ny*sin(alpha3 + alpha4 + alpha5 + alpha6) - d3*ny*px*sin(alpha3 + alpha4 + alpha5 + alpha6) + d3*nx*py*sin(alpha3 + alpha4 + alpha5 + alpha6));

return(dummyexp);

}