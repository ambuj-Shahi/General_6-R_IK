
# include <stdio.h>
# include <math.h>

double cx164(double* coeff)
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



double dummyexp = mu1*pow(a1,-1)*(4*lambda3*(-(d4*(d2 + d3*lambda2)*lambda5*mu1*mu4) + d2*d5*mu1*mu5 + d3*d5*lambda2*mu1*mu5 + d2*d4*lambda4*mu1*mu5 + d3*d4*lambda2*lambda4*mu1*mu5 + a1*a4*mu6*mz - a1*a5*mu6*mz + a1*(a4 - a5)*lambda6*nz) - 4*a1*lambda1*((a4 - a5)*lambda2*cos(alpha4 - alpha5) - a2*mu2*sin(alpha4 - alpha5)) + mu1*(-4*mu4*(d2*d3*lambda5 + a2*(a4 - a5)*mu2*mu5) + lambda4*(-4*a2*(a4 - a5)*lambda5*mu2 + mu5*(4*d2*d3 + lambda2*(-4*a4*a5 - 4*a6*d1*lz - 4*a6*d6*lx*mu6*mx - 4*a6*d6*ly*mu6*my - 4*d1*d6*mu6*mz - 4*a6*d6*lz*mu6*mz + 4*a6*lx*px + 4*d6*mu6*mx*px + 4*a6*ly*py + 4*d6*mu6*my*py + 4*d1*pz + 4*a6*lz*pz + 4*d6*mu6*mz*pz - 4*d6*lambda6*(d1*nz + a6*(lx*nx + ly*ny + lz*nz) - nx*px - ny*py - nz*pz) + 2*pow(a1,2) + 2*pow(a2,2) - 2*pow(a3,2) + 2*pow(a4,2) + 2*pow(a5,2) - 2*pow(d1,2) + 2*pow(d2,2) + 2*pow(d3,2) + 2*pow(d4,2) + 2*pow(d5,2) - 2*pow(a6,2)*pow(lx,2) - 2*pow(a6,2)*pow(ly,2) - 2*pow(a6,2)*pow(lz,2) - pow(d6,2)*pow(mx,2) - pow(d6,2)*pow(my,2) - pow(d6,2)*pow(mz,2) - pow(d6,2)*pow(nx,2) - pow(d6,2)*pow(ny,2) + cos(2*alpha6)*pow(d6,2)*(pow(mx,2) + pow(my,2) + pow(mz,2) - pow(nx,2) - pow(ny,2) - pow(nz,2)) - pow(d6,2)*pow(nz,2) - 2*pow(px,2) - 2*pow(py,2) - 2*pow(pz,2) - 2*mx*nx*pow(d6,2)*sin(2*alpha6) - 2*my*ny*pow(d6,2)*sin(2*alpha6) - 2*mz*nz*pow(d6,2)*sin(2*alpha6)))) + lambda2*(4*d4*d5*mu5 + mu4*(2*d5*(2*a6*lx*mu6*mx + 2*a6*ly*mu6*my + 2*d1*mu6*mz + 2*a6*lz*mu6*mz - 2*mu6*mx*px - 2*mu6*my*py - 2*mu6*mz*pz + 2*lambda6*(d1*nz + a6*(lx*nx + ly*ny + lz*nz) - nx*px - ny*py - nz*pz) + d6*pow(mx,2) + d6*pow(my,2) + d6*pow(mz,2) + d6*pow(nx,2) + d6*pow(ny,2) - d6*cos(2*alpha6)*(pow(mx,2) + pow(my,2) + pow(mz,2) - pow(nx,2) - pow(ny,2) - pow(nz,2)) + d6*pow(nz,2) + 2*d6*mx*nx*sin(2*alpha6) + 2*d6*my*ny*sin(2*alpha6) + 2*d6*mz*nz*sin(2*alpha6)) + lambda5*(4*a4*a5 + 4*a6*d1*lz + 4*a6*d6*lx*mu6*mx + 4*a6*d6*ly*mu6*my + 4*d1*d6*mu6*mz + 4*a6*d6*lz*mu6*mz - 4*a6*lx*px - 4*d6*mu6*mx*px - 4*a6*ly*py - 4*d6*mu6*my*py - 4*d1*pz - 4*a6*lz*pz - 4*d6*mu6*mz*pz + 4*d6*lambda6*(d1*nz + a6*(lx*nx + ly*ny + lz*nz) - nx*px - ny*py - nz*pz) - 2*pow(a1,2) - 2*pow(a2,2) + 2*pow(a3,2) - 2*pow(a4,2) - 2*pow(a5,2) + 2*pow(d1,2) - 2*pow(d2,2) - 2*pow(d3,2) - 2*pow(d4,2) + 2*pow(d5,2) + 2*pow(a6,2)*pow(lx,2) + 2*pow(a6,2)*pow(ly,2) + 2*pow(a6,2)*pow(lz,2) + pow(d6,2)*pow(mx,2) + pow(d6,2)*pow(my,2) + pow(d6,2)*pow(mz,2) + pow(d6,2)*pow(nx,2) + pow(d6,2)*pow(ny,2) - cos(2*alpha6)*pow(d6,2)*(pow(mx,2) + pow(my,2) + pow(mz,2) - pow(nx,2) - pow(ny,2) - pow(nz,2)) + pow(d6,2)*pow(nz,2) + 2*pow(px,2) + 2*pow(py,2) + 2*pow(pz,2) + 2*mx*nx*pow(d6,2)*sin(2*alpha6) + 2*my*ny*pow(d6,2)*sin(2*alpha6) + 2*mz*nz*pow(d6,2)*sin(2*alpha6))))));

return(dummyexp);

}