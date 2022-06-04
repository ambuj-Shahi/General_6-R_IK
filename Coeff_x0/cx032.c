
# include <stdio.h>
# include <math.h>

double cx032(double* coeff)
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



double dummyexp = -(pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1)*(a5*(2*d1*(mx*nx + my*ny) - a6*(lx*mz*nx + ly*mz*ny - 2*lz*(mx*nx + my*ny) + lx*mx*nz + ly*my*nz) + mz*nx*px + mx*nz*px + mz*ny*py + my*nz*py - 2*mx*nx*pz - 2*my*ny*pz)*cos(alpha2 + alpha3 - alpha4 - 2*alpha6) - a6*d3*ly*nx*cos(alpha2 - alpha5 - alpha6) + a6*d3*lx*ny*cos(alpha2 - alpha5 - alpha6) - d3*ny*px*cos(alpha2 - alpha5 - alpha6) + d3*nx*py*cos(alpha2 - alpha5 - alpha6) - a6*d4*ly*nx*cos(alpha2 + alpha3 - alpha5 - alpha6) + a6*d4*lx*ny*cos(alpha2 + alpha3 - alpha5 - alpha6) - d4*ny*px*cos(alpha2 + alpha3 - alpha5 - alpha6) + d4*nx*py*cos(alpha2 + alpha3 - alpha5 - alpha6) - a6*d5*ly*nx*cos(alpha2 + alpha3 - alpha4 - alpha5 - alpha6) + a6*d5*lx*ny*cos(alpha2 + alpha3 - alpha4 - alpha5 - alpha6) - d5*ny*px*cos(alpha2 + alpha3 - alpha4 - alpha5 - alpha6) + d5*nx*py*cos(alpha2 + alpha3 - alpha4 - alpha5 - alpha6) + a6*d3*ly*nx*cos(alpha2 + alpha5 - alpha6) - a6*d3*lx*ny*cos(alpha2 + alpha5 - alpha6) + d3*ny*px*cos(alpha2 + alpha5 - alpha6) - d3*nx*py*cos(alpha2 + alpha5 - alpha6) + a6*d4*ly*nx*cos(alpha2 + alpha3 + alpha5 - alpha6) - a6*d4*lx*ny*cos(alpha2 + alpha3 + alpha5 - alpha6) + d4*ny*px*cos(alpha2 + alpha3 + alpha5 - alpha6) - d4*nx*py*cos(alpha2 + alpha3 + alpha5 - alpha6) + a6*d5*ly*nx*cos(alpha2 + alpha3 - alpha4 + alpha5 - alpha6) - a6*d5*lx*ny*cos(alpha2 + alpha3 - alpha4 + alpha5 - alpha6) + d5*ny*px*cos(alpha2 + alpha3 - alpha4 + alpha5 - alpha6) - d5*nx*py*cos(alpha2 + alpha3 - alpha4 + alpha5 - alpha6) - a6*d3*ly*nx*cos(alpha2 - alpha5 + alpha6) + a6*d3*lx*ny*cos(alpha2 - alpha5 + alpha6) - d3*ny*px*cos(alpha2 - alpha5 + alpha6) + d3*nx*py*cos(alpha2 - alpha5 + alpha6) - a6*d4*ly*nx*cos(alpha2 + alpha3 - alpha5 + alpha6) + a6*d4*lx*ny*cos(alpha2 + alpha3 - alpha5 + alpha6) - d4*ny*px*cos(alpha2 + alpha3 - alpha5 + alpha6) + d4*nx*py*cos(alpha2 + alpha3 - alpha5 + alpha6) - a6*d5*ly*nx*cos(alpha2 + alpha3 - alpha4 - alpha5 + alpha6) + a6*d5*lx*ny*cos(alpha2 + alpha3 - alpha4 - alpha5 + alpha6) - d5*ny*px*cos(alpha2 + alpha3 - alpha4 - alpha5 + alpha6) + d5*nx*py*cos(alpha2 + alpha3 - alpha4 - alpha5 + alpha6) + a6*d3*ly*nx*cos(alpha2 + alpha5 + alpha6) - a6*d3*lx*ny*cos(alpha2 + alpha5 + alpha6) + d3*ny*px*cos(alpha2 + alpha5 + alpha6) - d3*nx*py*cos(alpha2 + alpha5 + alpha6) + a6*d4*ly*nx*cos(alpha2 + alpha3 + alpha5 + alpha6) - a6*d4*lx*ny*cos(alpha2 + alpha3 + alpha5 + alpha6) + d4*ny*px*cos(alpha2 + alpha3 + alpha5 + alpha6) - d4*nx*py*cos(alpha2 + alpha3 + alpha5 + alpha6) + a6*d5*ly*nx*cos(alpha2 + alpha3 - alpha4 + alpha5 + alpha6) - a6*d5*lx*ny*cos(alpha2 + alpha3 - alpha4 + alpha5 + alpha6) + d5*ny*px*cos(alpha2 + alpha3 - alpha4 + alpha5 + alpha6) - d5*nx*py*cos(alpha2 + alpha3 - alpha4 + alpha5 + alpha6) - 2*a5*d1*mx*nx*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) - 2*a5*a6*lz*mx*nx*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*a6*lx*mz*nx*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) - 2*a5*d1*my*ny*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) - 2*a5*a6*lz*my*ny*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*a6*ly*mz*ny*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*a6*lx*mx*nz*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*a6*ly*my*nz*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*mz*nx*px*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*mx*nz*px*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*mz*ny*py*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*my*nz*py*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) + 2*a5*mx*nx*pz*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) + 2*a5*my*ny*pz*cos(alpha2 + alpha3 - alpha4 + 2*alpha6) - 2*a5*a6*lx*mx*mz*sin(alpha2 + alpha3 - alpha4) - 2*a5*a6*ly*my*mz*sin(alpha2 + alpha3 - alpha4) - 2*a5*a6*lx*nx*nz*sin(alpha2 + alpha3 - alpha4) - 2*a5*a6*ly*ny*nz*sin(alpha2 + alpha3 - alpha4) + 2*a5*mx*mz*px*sin(alpha2 + alpha3 - alpha4) + 2*a5*nx*nz*px*sin(alpha2 + alpha3 - alpha4) + 2*a5*my*mz*py*sin(alpha2 + alpha3 - alpha4) + 2*a5*ny*nz*py*sin(alpha2 + alpha3 - alpha4) + 2*a5*d1*pow(mx,2)*sin(alpha2 + alpha3 - alpha4) + 2*a5*a6*lz*pow(mx,2)*sin(alpha2 + alpha3 - alpha4) - 2*a5*pz*pow(mx,2)*sin(alpha2 + alpha3 - alpha4) + 2*a5*d1*pow(my,2)*sin(alpha2 + alpha3 - alpha4) + 2*a5*a6*lz*pow(my,2)*sin(alpha2 + alpha3 - alpha4) - 2*a5*pz*pow(my,2)*sin(alpha2 + alpha3 - alpha4) + 2*a5*d1*pow(nx,2)*sin(alpha2 + alpha3 - alpha4) + 2*a5*a6*lz*pow(nx,2)*sin(alpha2 + alpha3 - alpha4) - 2*a5*pz*pow(nx,2)*sin(alpha2 + alpha3 - alpha4) + 2*a5*d1*pow(ny,2)*sin(alpha2 + alpha3 - alpha4) + 2*a5*a6*lz*pow(ny,2)*sin(alpha2 + alpha3 - alpha4) - 2*a5*pz*pow(ny,2)*sin(alpha2 + alpha3 - alpha4) + a5*a6*lx*mx*mz*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) + a5*a6*ly*my*mz*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) - a5*a6*lx*nx*nz*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) - a5*a6*ly*ny*nz*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) - a5*mx*mz*px*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) + a5*nx*nz*px*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) - a5*my*mz*py*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) + a5*ny*nz*py*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) - a5*d1*pow(mx,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) - a5*a6*lz*pow(mx,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) + a5*pz*pow(mx,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) - a5*d1*pow(my,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) - a5*a6*lz*pow(my,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) + a5*pz*pow(my,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) + a5*d1*pow(nx,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) + a5*a6*lz*pow(nx,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) - a5*pz*pow(nx,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) + a5*d1*pow(ny,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) + a5*a6*lz*pow(ny,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) - a5*pz*pow(ny,2)*sin(alpha2 + alpha3 - alpha4 - 2*alpha6) + a6*d3*ly*mx*sin(alpha2 - alpha5 - alpha6) - a6*d3*lx*my*sin(alpha2 - alpha5 - alpha6) + d3*my*px*sin(alpha2 - alpha5 - alpha6) - d3*mx*py*sin(alpha2 - alpha5 - alpha6) + a6*d4*ly*mx*sin(alpha2 + alpha3 - alpha5 - alpha6) - a6*d4*lx*my*sin(alpha2 + alpha3 - alpha5 - alpha6) + d4*my*px*sin(alpha2 + alpha3 - alpha5 - alpha6) - d4*mx*py*sin(alpha2 + alpha3 - alpha5 - alpha6) + a6*d5*ly*mx*sin(alpha2 + alpha3 - alpha4 - alpha5 - alpha6) - a6*d5*lx*my*sin(alpha2 + alpha3 - alpha4 - alpha5 - alpha6) + d5*my*px*sin(alpha2 + alpha3 - alpha4 - alpha5 - alpha6) - d5*mx*py*sin(alpha2 + alpha3 - alpha4 - alpha5 - alpha6) - a6*d3*ly*mx*sin(alpha2 + alpha5 - alpha6) + a6*d3*lx*my*sin(alpha2 + alpha5 - alpha6) - d3*my*px*sin(alpha2 + alpha5 - alpha6) + d3*mx*py*sin(alpha2 + alpha5 - alpha6) - a6*d4*ly*mx*sin(alpha2 + alpha3 + alpha5 - alpha6) + a6*d4*lx*my*sin(alpha2 + alpha3 + alpha5 - alpha6) - d4*my*px*sin(alpha2 + alpha3 + alpha5 - alpha6) + d4*mx*py*sin(alpha2 + alpha3 + alpha5 - alpha6) - a6*d5*ly*mx*sin(alpha2 + alpha3 - alpha4 + alpha5 - alpha6) + a6*d5*lx*my*sin(alpha2 + alpha3 - alpha4 + alpha5 - alpha6) - d5*my*px*sin(alpha2 + alpha3 - alpha4 + alpha5 - alpha6) + d5*mx*py*sin(alpha2 + alpha3 - alpha4 + alpha5 - alpha6) - a6*d3*ly*mx*sin(alpha2 - alpha5 + alpha6) + a6*d3*lx*my*sin(alpha2 - alpha5 + alpha6) - d3*my*px*sin(alpha2 - alpha5 + alpha6) + d3*mx*py*sin(alpha2 - alpha5 + alpha6) - a6*d4*ly*mx*sin(alpha2 + alpha3 - alpha5 + alpha6) + a6*d4*lx*my*sin(alpha2 + alpha3 - alpha5 + alpha6) - d4*my*px*sin(alpha2 + alpha3 - alpha5 + alpha6) + d4*mx*py*sin(alpha2 + alpha3 - alpha5 + alpha6) - a6*d5*ly*mx*sin(alpha2 + alpha3 - alpha4 - alpha5 + alpha6) + a6*d5*lx*my*sin(alpha2 + alpha3 - alpha4 - alpha5 + alpha6) - d5*my*px*sin(alpha2 + alpha3 - alpha4 - alpha5 + alpha6) + d5*mx*py*sin(alpha2 + alpha3 - alpha4 - alpha5 + alpha6) + a6*d3*ly*mx*sin(alpha2 + alpha5 + alpha6) - a6*d3*lx*my*sin(alpha2 + alpha5 + alpha6) + d3*my*px*sin(alpha2 + alpha5 + alpha6) - d3*mx*py*sin(alpha2 + alpha5 + alpha6) + a6*d4*ly*mx*sin(alpha2 + alpha3 + alpha5 + alpha6) - a6*d4*lx*my*sin(alpha2 + alpha3 + alpha5 + alpha6) + d4*my*px*sin(alpha2 + alpha3 + alpha5 + alpha6) - d4*mx*py*sin(alpha2 + alpha3 + alpha5 + alpha6) + a6*d5*ly*mx*sin(alpha2 + alpha3 - alpha4 + alpha5 + alpha6) - a6*d5*lx*my*sin(alpha2 + alpha3 - alpha4 + alpha5 + alpha6) + d5*my*px*sin(alpha2 + alpha3 - alpha4 + alpha5 + alpha6) - d5*mx*py*sin(alpha2 + alpha3 - alpha4 + alpha5 + alpha6) + a5*a6*lx*mx*mz*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*a6*ly*my*mz*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*a6*lx*nx*nz*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*a6*ly*ny*nz*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*mx*mz*px*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*nx*nz*px*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*my*mz*py*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*ny*nz*py*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*d1*pow(mx,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*a6*lz*pow(mx,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*pz*pow(mx,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*d1*pow(my,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*a6*lz*pow(my,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*pz*pow(my,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*d1*pow(nx,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*a6*lz*pow(nx,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*pz*pow(nx,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*d1*pow(ny,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) + a5*a6*lz*pow(ny,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6) - a5*pz*pow(ny,2)*sin(alpha2 + alpha3 - alpha4 + 2*alpha6)))/2.;

return(dummyexp);

}