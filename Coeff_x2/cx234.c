
# include <stdio.h>
# include <math.h>

double cx234(double* coeff)
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



double dummyexp = 2*pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1)*(d3*mu2*(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py))*sin(alpha4 - alpha5) + lambda2*mu3*((a4 - a5)*pow(lambda6,2)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))) + mu6*(a6*d5*ly*mu5*mx - a6*d5*lx*mu5*my - a4*a6*lx*mu6*mx*mz + a5*a6*lx*mu6*mx*mz - a4*a6*ly*mu6*my*mz + a5*a6*ly*mu6*my*mz + d5*mu5*my*px + a4*mu6*mx*mz*px - a5*mu6*mx*mz*px - d5*mu5*mx*py + a4*mu6*my*mz*py - a5*mu6*my*mz*py + 2*lambda6*(mx*nx + my*ny)*(a4*(d1 + a6*lz) + a5*pz) + a4*d1*mu6*pow(mx,2) - a5*d1*mu6*pow(mx,2) + a4*a6*lz*mu6*pow(mx,2) - a5*a6*lz*mu6*pow(mx,2) - a4*mu6*pz*pow(mx,2) + a5*mu6*pz*pow(mx,2) + a4*d1*mu6*pow(my,2) - a5*d1*mu6*pow(my,2) + a4*a6*lz*mu6*pow(my,2) - a5*a6*lz*mu6*pow(my,2) - a4*mu6*pz*pow(my,2) + a5*mu6*pz*pow(my,2) - d4*(a6*ly*mx - a6*lx*my + my*px - mx*py)*sin(alpha4 - alpha5)) - lambda6*(-(d5*mu5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)) + mu6*(a5*(2*d1*(mx*nx + my*ny) - a6*(lx*mz*nx + ly*mz*ny - 2*lz*(mx*nx + my*ny) + lx*mx*nz + ly*my*nz) + mz*nx*px + mx*nz*px + mz*ny*py + my*nz*py) + a4*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz)) + d4*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*sin(alpha4 - alpha5))) + lambda3*mu2*(-((a4 - a5)*pow(lambda6,2)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2))))) + mu6*(-(a6*d5*ly*mu5*mx) + a6*d5*lx*mu5*my + a4*a6*lx*mu6*mx*mz - a5*a6*lx*mu6*mx*mz + a4*a6*ly*mu6*my*mz - a5*a6*ly*mu6*my*mz - d5*mu5*my*px - a4*mu6*mx*mz*px + a5*mu6*mx*mz*px + d5*mu5*mx*py - a4*mu6*my*mz*py + a5*mu6*my*mz*py + 2*lambda6*(mx*nx + my*ny)*(a5*(d1 + a6*lz) + a4*pz) - a4*d1*mu6*pow(mx,2) + a5*d1*mu6*pow(mx,2) - a4*a6*lz*mu6*pow(mx,2) + a5*a6*lz*mu6*pow(mx,2) + a4*mu6*pz*pow(mx,2) - a5*mu6*pz*pow(mx,2) - a4*d1*mu6*pow(my,2) + a5*d1*mu6*pow(my,2) - a4*a6*lz*mu6*pow(my,2) + a5*a6*lz*mu6*pow(my,2) + a4*mu6*pz*pow(my,2) - a5*mu6*pz*pow(my,2) + d4*(a6*ly*mx - a6*lx*my + my*px - mx*py)*sin(alpha4 - alpha5)) + lambda6*(-(d5*mu5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)) - mu6*(a4*(2*d1*(mx*nx + my*ny) - a6*(lx*mz*nx + ly*mz*ny - 2*lz*(mx*nx + my*ny) + lx*mx*nz + ly*my*nz) + mz*nx*px + mx*nz*px + mz*ny*py + my*nz*py) + a5*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz)) + d4*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*sin(alpha4 - alpha5))));

return(dummyexp);

}