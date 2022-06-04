
# include <stdio.h>
# include <math.h>

double cx232(double* coeff)
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



double dummyexp = 2*(mu2*(-(lambda6*((d3 + d5*mu3*mu4)*mu5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) + a5*mu3*mu4*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz))) + mu6*(-((d3 + d5*mu3*mu4)*mu5*(a6*ly*mx - a6*lx*my + my*px - mx*py)) + a5*mu3*mu4*(2*lambda6*(d1 + a6*lz)*(mx*nx + my*ny) + mu6*(mx*mz*px + my*mz*py - pz*pow(mx,2) - pz*pow(my,2) + d1*(pow(mx,2) + pow(my,2)) + a6*(-((lx*mx + ly*my)*mz) + lz*(pow(mx,2) + pow(my,2)))))) + a5*mu3*mu4*pow(lambda6,2)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2))))) + lambda3*(d4*mu2*mu5*(mu6*(-(a6*ly*mx) + a6*lx*my - my*px + mx*py) + lambda6*(-(a6*ly*nx) + a6*lx*ny - ny*px + nx*py)) + lambda4*mu2*(-(lambda6*(d5*mu5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) + a5*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz))) + mu6*(2*a5*lambda6*(d1 + a6*lz)*(mx*nx + my*ny) - d5*mu5*(a6*ly*mx - a6*lx*my + my*px - mx*py) + a5*mu6*(mx*mz*px + my*mz*py - pz*pow(mx,2) - pz*pow(my,2) + d1*(pow(mx,2) + pow(my,2)) + a6*(-((lx*mx + ly*my)*mz) + lz*(pow(mx,2) + pow(my,2))))) + a5*pow(lambda6,2)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2))))) + lambda2*mu4*(-(lambda6*(d5*mu5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) + a5*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - mx*nz*px - my*nz*py - mz*(nx*px + ny*py) + 2*mx*nx*pz + 2*my*ny*pz))) + mu6*(2*a5*lambda6*(d1 + a6*lz)*(mx*nx + my*ny) - d5*mu5*(a6*ly*mx - a6*lx*my + my*px - mx*py) + a5*mu6*(mx*mz*px + my*mz*py - pz*pow(mx,2) - pz*pow(my,2) + d1*(pow(mx,2) + pow(my,2)) + a6*(-((lx*mx + ly*my)*mz) + lz*(pow(mx,2) + pow(my,2))))) + a5*pow(lambda6,2)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))))) - lambda2*mu3*(d4*mu5*(mu6*(-(a6*ly*mx) + a6*lx*my - my*px + mx*py) + lambda6*(-(a6*ly*nx) + a6*lx*ny - ny*px + nx*py)) + lambda4*(lambda6*(-(d5*mu5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)) + a5*mu6*(-(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz)) + mz*nx*px + mx*nz*px + mz*ny*py + my*nz*py - 2*my*ny*pz)) + mu6*(-(d5*mu5*(a6*ly*mx - a6*lx*my + my*px - mx*py)) + 2*a5*lambda6*(d1*mx*nx + a6*lz*mx*nx + d1*my*ny + a6*lz*my*ny - mx*nx*pz) + a5*mu6*(mx*mz*px + my*mz*py - pz*pow(mx,2) - pz*pow(my,2) + d1*(pow(mx,2) + pow(my,2)) + a6*(-((lx*mx + ly*my)*mz) + lz*(pow(mx,2) + pow(my,2))))) + a5*pow(lambda6,2)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))))))*pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1);

return(dummyexp);

}