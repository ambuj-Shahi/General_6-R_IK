
# include <stdio.h>
# include <math.h>

double cx227(double* coeff)
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



double dummyexp = (pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1)*(2*pow(mu6,2)*((a6*ly*mx - a6*lx*my + my*px - mx*py)*(d1*mz + a6*(lx*mx + ly*my + lz*mz) - mx*px - my*py - mz*pz) + a1*(d2 + d3*lambda2 + d4*cos(alpha2 - alpha3) + d5*cos(alpha2 - alpha3 - alpha4) + d6*cos(alpha2 - alpha3 - alpha4 + alpha5))*pow(mu1,-1)*(pow(mx,2) + pow(my,2)) + a1*lambda1*pow(mu1,-1)*(mx*mz*px + my*mz*py - pz*pow(mx,2) - pz*pow(my,2) + d1*(pow(mx,2) + pow(my,2)) + a6*(-((lx*mx + ly*my)*mz) + lz*(pow(mx,2) + pow(my,2))))) + (mu6*pow(mu1,-1)*(4*a1*(a6*lx*mx + a6*ly*my - mx*px - my*py)*cos(alpha2 - alpha3 - alpha4 + alpha5) - 2*mu1*(a6*ly*mx - a6*lx*my + my*px - mx*py)*(-2*d5*lambda5 - 2*d4*cos(alpha4 - alpha5) - 2*d3*cos(alpha3 + alpha4 - alpha5) - 2*d2*cos(alpha2 - alpha3 - alpha4 + alpha5) - d6*pow(mx,2) + d6*cos(2*alpha6)*pow(mx,2) - d6*pow(my,2) + d6*cos(2*alpha6)*pow(my,2) - d6*pow(mz,2) + d6*cos(2*alpha6)*pow(mz,2) - d6*pow(nx,2) - d6*cos(2*alpha6)*pow(nx,2) - d6*pow(ny,2) - d6*cos(2*alpha6)*pow(ny,2) - d6*pow(nz,2) - d6*cos(2*alpha6)*pow(nz,2))))/2. + 2*pow(lambda6,2)*pow(mu1,-1)*(a6*d1*ly*mu1*nx*nz - a6*d1*lx*mu1*ny*nz + 2*a6*lx*mu1*nx*ny*px + d1*mu1*ny*nz*px + a6*lz*mu1*ny*nz*px - 2*a6*ly*mu1*nx*ny*py - d1*mu1*nx*nz*py - a6*lz*mu1*nx*nz*py - a6*ly*mu1*nx*nz*pz + a6*lx*mu1*ny*nz*pz - mu1*ny*nz*px*pz + mu1*nx*nz*py*pz - a6*d6*lx*mx*nx*ny*cos(alpha1 - alpha6) + a6*d6*ly*my*nx*ny*cos(alpha1 - alpha6) + a6*d6*ly*mz*nx*nz*cos(alpha1 - alpha6) - a6*d6*lx*mz*ny*nz*cos(alpha1 - alpha6) + d6*mx*nx*ny*px*cos(alpha1 - alpha6) + d6*mz*ny*nz*px*cos(alpha1 - alpha6) - d6*my*nx*ny*py*cos(alpha1 - alpha6) - d6*mz*nx*nz*py*cos(alpha1 - alpha6) + a6*d6*lx*mx*nx*ny*cos(alpha1 + alpha6) - a6*d6*ly*my*nx*ny*cos(alpha1 + alpha6) - a6*d6*ly*mz*nx*nz*cos(alpha1 + alpha6) + a6*d6*lx*mz*ny*nz*cos(alpha1 + alpha6) - d6*mx*nx*ny*px*cos(alpha1 + alpha6) - d6*mz*ny*nz*px*cos(alpha1 + alpha6) + d6*my*nx*ny*py*cos(alpha1 + alpha6) + d6*mz*nx*nz*py*cos(alpha1 + alpha6) + ly*lz*mu1*nx*nz*pow(a6,2) - lx*lz*mu1*ny*nz*pow(a6,2) - mu1*nx*ny*pow(a6,2)*pow(lx,2) + mu1*nx*ny*pow(a6,2)*pow(ly,2) + a1*d2*pow(nx,2) - a6*ly*mu1*px*pow(nx,2) - a6*lx*mu1*py*pow(nx,2) + mu1*px*py*pow(nx,2) + a1*d4*cos(alpha2 - alpha3)*pow(nx,2) + a1*d5*cos(alpha2 - alpha3 - alpha4)*pow(nx,2) + a1*d6*cos(alpha2 - alpha3 - alpha4 + alpha5)*pow(nx,2) + a6*d6*ly*mx*cos(alpha1 - alpha6)*pow(nx,2) - d6*mx*py*cos(alpha1 - alpha6)*pow(nx,2) - a6*d6*ly*mx*cos(alpha1 + alpha6)*pow(nx,2) + d6*mx*py*cos(alpha1 + alpha6)*pow(nx,2) + lx*ly*mu1*pow(a6,2)*pow(nx,2) + a1*d2*pow(ny,2) + a6*ly*mu1*px*pow(ny,2) + a6*lx*mu1*py*pow(ny,2) - mu1*px*py*pow(ny,2) + a1*d4*cos(alpha2 - alpha3)*pow(ny,2) + a1*d5*cos(alpha2 - alpha3 - alpha4)*pow(ny,2) + a1*d6*cos(alpha2 - alpha3 - alpha4 + alpha5)*pow(ny,2) - a6*d6*lx*my*cos(alpha1 - alpha6)*pow(ny,2) + d6*my*px*cos(alpha1 - alpha6)*pow(ny,2) + a6*d6*lx*my*cos(alpha1 + alpha6)*pow(ny,2) - d6*my*px*cos(alpha1 + alpha6)*pow(ny,2) - lx*ly*mu1*pow(a6,2)*pow(ny,2) + a1*d3*lambda2*(pow(nx,2) + pow(ny,2)) + a1*lambda1*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))) - mu1*nx*ny*pow(px,2) + mu1*nx*ny*pow(py,2)) + (d1*mz*ny*px + d1*my*nz*px + mz*nx*py*pz + mx*nz*py*pz + a6*(d1*ly*(mz*nx + mx*nz) + (mz*ny + my*nz)*(lz*px + lx*pz)) + ly*(ly*my*nx + lz*mz*nx + ly*mx*ny + lz*mx*nz)*pow(a6,2) + a1*lambda1*(mz*nx*px + mx*nz*px + mz*ny*py + my*nz*py)*pow(mu1,-1) + my*nx*pow(py,2) + mx*ny*pow(py,2))*sin(2*alpha6) + lambda6*(-2*a6*d2*lambda3*lambda4*ly*mu2*mu5*nx + 2*a6*d3*lambda4*ly*mu3*mu5*nx + 2*a6*d2*lambda2*lambda4*ly*mu3*mu5*nx + 2*a6*d4*ly*mu4*mu5*nx + 2*a6*d3*lambda3*ly*mu4*mu5*nx + 2*a6*d2*lambda2*lambda3*ly*mu4*mu5*nx + 2*a6*d2*ly*mu2*mu3*mu4*mu5*nx + 2*a6*d2*lambda3*lambda4*lx*mu2*mu5*ny - 2*a6*d3*lambda4*lx*mu3*mu5*ny - 2*a6*d2*lambda2*lambda4*lx*mu3*mu5*ny - 2*a6*d4*lx*mu4*mu5*ny - 2*a6*d3*lambda3*lx*mu4*mu5*ny - 2*a6*d2*lambda2*lambda3*lx*mu4*mu5*ny - 2*a6*d2*lx*mu2*mu3*mu4*mu5*ny - 2*a6*d1*lx*mu6*mz*ny - 2*a6*d1*lx*mu6*my*nz - 4*a6*ly*mu6*mx*nx*px + 4*a6*lx*mu6*my*nx*px - 2*d2*lambda3*lambda4*mu2*mu5*ny*px + 2*d3*lambda4*mu3*mu5*ny*px + 2*d2*lambda2*lambda4*mu3*mu5*ny*px + 2*d4*mu4*mu5*ny*px + 2*d3*lambda3*mu4*mu5*ny*px + 2*d2*lambda2*lambda3*mu4*mu5*ny*px + 2*d2*mu2*mu3*mu4*mu5*ny*px + 4*a6*lx*mu6*mx*ny*px + 4*a6*ly*mu6*my*ny*px + 2*d2*lambda3*lambda4*mu2*mu5*nx*py - 2*d3*lambda4*mu3*mu5*nx*py - 2*d2*lambda2*lambda4*mu3*mu5*nx*py - 2*d4*mu4*mu5*nx*py - 2*d3*lambda3*mu4*mu5*nx*py - 2*d2*lambda2*lambda3*mu4*mu5*nx*py - 2*d2*mu2*mu3*mu4*mu5*nx*py - 4*a6*lx*mu6*mx*nx*py - 4*a6*ly*mu6*my*nx*py - 2*d1*mu6*mz*nx*py - 2*a6*lz*mu6*mz*nx*py - 4*a6*ly*mu6*mx*ny*py + 4*a6*lx*mu6*my*ny*py - 2*d1*mu6*mx*nz*py - 2*a6*lz*mu6*mx*nz*py + 4*mu6*mx*nx*px*py - 4*mu6*my*ny*px*py - 2*a6*ly*mu6*mz*nx*pz - 2*a6*ly*mu6*mx*nz*pz - 2*mu6*mz*ny*px*pz - 2*mu6*my*nz*px*pz + 4*lx*ly*mu6*mx*nx*pow(a6,2) - 4*lx*ly*mu6*my*ny*pow(a6,2) - 2*lx*lz*mu6*mz*ny*pow(a6,2) - 2*lx*lz*mu6*my*nz*pow(a6,2) - 2*mu6*my*nx*pow(a6,2)*pow(lx,2) - 2*mu6*mx*ny*pow(a6,2)*pow(lx,2) - 2*a1*a6*lambda3*lambda4*lx*mu2*mu5*nx*pow(mu1,-1) + 2*a1*a6*lambda2*lambda4*lx*mu3*mu5*nx*pow(mu1,-1) + 2*a1*a6*lambda2*lambda3*lx*mu4*mu5*nx*pow(mu1,-1) + 2*a1*a6*lx*mu2*mu3*mu4*mu5*nx*pow(mu1,-1) + 4*a1*d2*mu6*mx*nx*pow(mu1,-1) + 4*a1*d1*lambda1*mu6*mx*nx*pow(mu1,-1) + 4*a1*d3*lambda2*mu6*mx*nx*pow(mu1,-1) + 4*a1*d4*lambda2*lambda3*mu6*mx*nx*pow(mu1,-1) + 4*a1*d5*lambda2*lambda3*lambda4*mu6*mx*nx*pow(mu1,-1) + 4*a1*a6*lambda1*lz*mu6*mx*nx*pow(mu1,-1) + 4*a1*d4*mu2*mu3*mu6*mx*nx*pow(mu1,-1) + 4*a1*d5*lambda4*mu2*mu3*mu6*mx*nx*pow(mu1,-1) + 4*a1*d5*lambda3*mu2*mu4*mu6*mx*nx*pow(mu1,-1) - 4*a1*d5*lambda2*mu3*mu4*mu6*mx*nx*pow(mu1,-1) - 4*a1*d6*lambda3*lambda4*mu2*mu5*mu6*mx*nx*pow(mu1,-1) + 4*a1*d6*lambda2*lambda4*mu3*mu5*mu6*mx*nx*pow(mu1,-1) + 4*a1*d6*lambda2*lambda3*mu4*mu5*mu6*mx*nx*pow(mu1,-1) + 4*a1*d6*mu2*mu3*mu4*mu5*mu6*mx*nx*pow(mu1,-1) - 2*a1*a6*lambda1*lx*mu6*mz*nx*pow(mu1,-1) - 2*a1*a6*lambda3*lambda4*ly*mu2*mu5*ny*pow(mu1,-1) + 2*a1*a6*lambda2*lambda4*ly*mu3*mu5*ny*pow(mu1,-1) + 2*a1*a6*lambda2*lambda3*ly*mu4*mu5*ny*pow(mu1,-1) + 2*a1*a6*ly*mu2*mu3*mu4*mu5*ny*pow(mu1,-1) + 4*a1*d2*mu6*my*ny*pow(mu1,-1) + 4*a1*d1*lambda1*mu6*my*ny*pow(mu1,-1) + 4*a1*d3*lambda2*mu6*my*ny*pow(mu1,-1) + 4*a1*d4*lambda2*lambda3*mu6*my*ny*pow(mu1,-1) + 4*a1*d5*lambda2*lambda3*lambda4*mu6*my*ny*pow(mu1,-1) + 4*a1*a6*lambda1*lz*mu6*my*ny*pow(mu1,-1) + 4*a1*d4*mu2*mu3*mu6*my*ny*pow(mu1,-1) + 4*a1*d5*lambda4*mu2*mu3*mu6*my*ny*pow(mu1,-1) + 4*a1*d5*lambda3*mu2*mu4*mu6*my*ny*pow(mu1,-1) - 4*a1*d5*lambda2*mu3*mu4*mu6*my*ny*pow(mu1,-1) - 4*a1*d6*lambda3*lambda4*mu2*mu5*mu6*my*ny*pow(mu1,-1) + 4*a1*d6*lambda2*lambda4*mu3*mu5*mu6*my*ny*pow(mu1,-1) + 4*a1*d6*lambda2*lambda3*mu4*mu5*mu6*my*ny*pow(mu1,-1) + 4*a1*d6*mu2*mu3*mu4*mu5*mu6*my*ny*pow(mu1,-1) - 2*a1*a6*lambda1*ly*mu6*mz*ny*pow(mu1,-1) - 2*a1*a6*lambda1*lx*mu6*mx*nz*pow(mu1,-1) - 2*a1*a6*lambda1*ly*mu6*my*nz*pow(mu1,-1) + 2*a1*lambda3*lambda4*mu2*mu5*nx*px*pow(mu1,-1) - 2*a1*lambda2*lambda4*mu3*mu5*nx*px*pow(mu1,-1) - 2*a1*lambda2*lambda3*mu4*mu5*nx*px*pow(mu1,-1) - 2*a1*mu2*mu3*mu4*mu5*nx*px*pow(mu1,-1) + 2*a1*lambda3*lambda4*mu2*mu5*ny*py*pow(mu1,-1) - 2*a1*lambda2*lambda4*mu3*mu5*ny*py*pow(mu1,-1) - 2*a1*lambda2*lambda3*mu4*mu5*ny*py*pow(mu1,-1) - 2*a1*mu2*mu3*mu4*mu5*ny*py*pow(mu1,-1) - 4*a1*lambda1*mu6*mx*nx*pz*pow(mu1,-1) - 4*a1*lambda1*mu6*my*ny*pz*pow(mu1,-1) - 4*a6*d6*lx*mx*my*nx*pow(mu6,2) + 4*a6*d6*ly*mx*my*ny*pow(mu6,2) + 4*a6*d6*ly*mx*mz*nz*pow(mu6,2) - 4*a6*d6*lx*my*mz*nz*pow(mu6,2) + 4*d6*mx*my*nx*px*pow(mu6,2) + 4*d6*my*mz*nz*px*pow(mu6,2) - 4*d6*mx*my*ny*py*pow(mu6,2) - 4*d6*mx*mz*nz*py*pow(mu6,2) + a6*d6*ly*nx*pow(mx,2) - a6*d6*lx*ny*pow(mx,2) + d6*ny*px*pow(mx,2) - d6*nx*py*pow(mx,2) + 4*a6*d6*ly*nx*pow(mu6,2)*pow(mx,2) - 4*d6*nx*py*pow(mu6,2)*pow(mx,2) + a6*d6*ly*nx*pow(my,2) - a6*d6*lx*ny*pow(my,2) + d6*ny*px*pow(my,2) - d6*nx*py*pow(my,2) - 4*a6*d6*lx*ny*pow(mu6,2)*pow(my,2) + 4*d6*ny*px*pow(mu6,2)*pow(my,2) + a6*d6*ly*nx*pow(mz,2) - a6*d6*lx*ny*pow(mz,2) + d6*ny*px*pow(mz,2) - d6*nx*py*pow(mz,2) - a6*d6*lx*ny*pow(nx,2) + d6*ny*px*pow(nx,2) + a6*d6*ly*pow(nx,3) - d6*py*pow(nx,3) + a6*d6*ly*nx*pow(ny,2) - d6*nx*py*pow(ny,2) - a6*d6*lx*pow(ny,3) + d6*px*pow(ny,3) - d6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*cos(2*alpha6)*(pow(mx,2) + pow(my,2) + pow(mz,2) - pow(nx,2) - pow(ny,2) - pow(nz,2)) + a6*d6*ly*nx*pow(nz,2) - a6*d6*lx*ny*pow(nz,2) + d6*ny*px*pow(nz,2) - d6*nx*py*pow(nz,2) - 2*mu6*my*nx*pow(px,2) - 2*mu6*mx*ny*pow(px,2) + lambda5*pow(mu1,-1)*(2*d5*mu1*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) + 2*a1*(a6*lx*nx + a6*ly*ny - nx*px - ny*py)*cos(alpha2 - alpha3 - alpha4) + a6*d4*ly*nx*sin(alpha1 - alpha4) - a6*d4*lx*ny*sin(alpha1 - alpha4) + d4*ny*px*sin(alpha1 - alpha4) - d4*nx*py*sin(alpha1 - alpha4) + a6*d3*ly*nx*sin(alpha1 - alpha3 - alpha4) - a6*d3*lx*ny*sin(alpha1 - alpha3 - alpha4) + d3*ny*px*sin(alpha1 - alpha3 - alpha4) - d3*nx*py*sin(alpha1 - alpha3 - alpha4) + a6*d2*ly*nx*sin(alpha1 + alpha2 - alpha3 - alpha4) - a6*d2*lx*ny*sin(alpha1 + alpha2 - alpha3 - alpha4) + d2*ny*px*sin(alpha1 + alpha2 - alpha3 - alpha4) - d2*nx*py*sin(alpha1 + alpha2 - alpha3 - alpha4) + a6*d4*ly*nx*sin(alpha1 + alpha4) - a6*d4*lx*ny*sin(alpha1 + alpha4) + d4*ny*px*sin(alpha1 + alpha4) - d4*nx*py*sin(alpha1 + alpha4) + a6*d3*ly*nx*sin(alpha1 + alpha3 + alpha4) - a6*d3*lx*ny*sin(alpha1 + alpha3 + alpha4) + d3*ny*px*sin(alpha1 + alpha3 + alpha4) - d3*nx*py*sin(alpha1 + alpha3 + alpha4) + a6*d2*ly*nx*sin(alpha1 - alpha2 + alpha3 + alpha4) - a6*d2*lx*ny*sin(alpha1 - alpha2 + alpha3 + alpha4) + d2*ny*px*sin(alpha1 - alpha2 + alpha3 + alpha4) - d2*nx*py*sin(alpha1 - alpha2 + alpha3 + alpha4) - 2*a1*d6*mx*nx*sin(alpha2 - alpha3 - alpha4 - alpha6) - 2*a1*d6*my*ny*sin(alpha2 - alpha3 - alpha4 - alpha6) + 2*a1*d6*mx*nx*sin(alpha2 - alpha3 - alpha4 + alpha6) + 2*a1*d6*my*ny*sin(alpha2 - alpha3 - alpha4 + alpha6)))))/2.;

return(dummyexp);

}