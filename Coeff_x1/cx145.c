
# include <stdio.h>
# include <math.h>

double cx145(double* coeff)
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



double dummyexp = 2*pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1)*(8*(a5*(d4 + d3*lambda3)*lambda5*mu2 + (a2*d2*lambda3 + a2*lambda2*(d4 + d3*lambda3) - mu2*(a3*d3*mu3 + a4*d5*mu4))*mu5)*(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)) + lambda4*(8*a2*d5*lambda2*mu5*mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) - 2*a5*mu2*pow(lambda6,2)*(-4*d1*ny*nz*px + 4*d1*nx*nz*py + 4*ny*nz*px*pz - 4*nx*nz*py*pz + 4*(-(ly*nx) + lx*ny)*(lx*nx + ly*ny + lz*nz)*pow(a6,2) - my*mz*nx*nz*pow(d6,2) + mx*mz*ny*nz*pow(d6,2) + nx*ny*pow(d6,2)*pow(mx,2) - nx*ny*pow(d6,2)*pow(my,2) - 4*px*py*pow(nx,2) - mx*my*pow(d6,2)*pow(nx,2) + 4*px*py*pow(ny,2) + mx*my*pow(d6,2)*pow(ny,2) - 4*a6*(d1*(ly*nx - lx*ny)*nz + 2*lx*nx*ny*px + lz*ny*nz*px - lz*nx*nz*py + lx*ny*nz*pz - lx*py*pow(nx,2) + lx*py*pow(ny,2) - ly*(2*nx*ny*py + nx*nz*pz + px*pow(nx,2) - px*pow(ny,2))) + 4*a1*lambda1*pow(mu1,-1)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))) + 4*nx*ny*pow(px,2) - 4*nx*ny*pow(py,2)) + lambda6*(8*a2*d5*lambda2*mu5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) + a5*mu2*(8*a6*d1*ly*mu6*mz*nx - 8*a6*d1*lx*mu6*mz*ny + 8*a6*d1*ly*mu6*mx*nz - 8*a6*d1*lx*mu6*my*nz + 8*d1*mu6*mz*ny*px + 8*a6*lz*mu6*mz*ny*px + 8*d1*mu6*my*nz*px + 8*a6*lz*mu6*my*nz*px - 8*d1*mu6*mz*nx*py - 8*a6*lz*mu6*mz*nx*py - 8*d1*mu6*mx*nz*py - 8*a6*lz*mu6*mx*nz*py + 8*d5*lambda5*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) - 8*a6*ly*mu6*mz*nx*pz + 8*a6*lx*mu6*mz*ny*pz - 8*a6*ly*mu6*mx*nz*pz + 8*a6*lx*mu6*my*nz*pz - 8*mu6*mz*ny*px*pz - 8*mu6*my*nz*px*pz + 8*mu6*mz*nx*py*pz + 8*mu6*mx*nz*py*pz + 8*ly*lz*mu6*mz*nx*pow(a6,2) - 8*lx*lz*mu6*mz*ny*pow(a6,2) + 8*ly*lz*mu6*mx*nz*pow(a6,2) - 8*lx*lz*mu6*my*nz*pow(a6,2) - 2*my*mz*nx*nz*cos(3*alpha6)*pow(d6,2) + 2*mx*mz*ny*nz*cos(3*alpha6)*pow(d6,2) - 8*mu6*my*nx*pow(a6,2)*pow(lx,2) - 8*mu6*mx*ny*pow(a6,2)*pow(lx,2) + 8*mu6*my*nx*pow(a6,2)*pow(ly,2) + 8*mu6*mx*ny*pow(a6,2)*pow(ly,2) + 8*a1*a6*lambda1*lx*mu6*mz*nx*pow(mu1,-1) + 8*a1*a6*lambda1*ly*mu6*mz*ny*pow(mu1,-1) + 8*a1*a6*lambda1*lx*mu6*mx*nz*pow(mu1,-1) + 8*a1*a6*lambda1*ly*mu6*my*nz*pow(mu1,-1) - 8*a1*lambda1*mu6*mz*nx*px*pow(mu1,-1) - 8*a1*lambda1*mu6*mx*nz*px*pow(mu1,-1) - 8*a1*lambda1*mu6*mz*ny*py*pow(mu1,-1) - 8*a1*lambda1*mu6*my*nz*py*pow(mu1,-1) + 4*a6*d6*ly*nx*pow(mx,2) - 4*a6*d6*lx*ny*pow(mx,2) + 4*d6*ny*px*pow(mx,2) - 4*d6*nx*py*pow(mx,2) + 2*nx*ny*cos(3*alpha6)*pow(d6,2)*pow(mx,2) - 2*mu6*ny*pow(d6,2)*pow(mx,3) + 4*a6*d6*ly*nx*pow(my,2) - 4*a6*d6*lx*ny*pow(my,2) + 4*d6*ny*px*pow(my,2) - 4*d6*nx*py*pow(my,2) - 2*mu6*mx*ny*pow(d6,2)*pow(my,2) - 2*nx*ny*cos(3*alpha6)*pow(d6,2)*pow(my,2) + 4*a6*d6*ly*nx*pow(mz,2) - 4*a6*d6*lx*ny*pow(mz,2) + 4*d6*ny*px*pow(mz,2) - 4*d6*nx*py*pow(mz,2) - 2*mu6*mx*ny*pow(d6,2)*pow(mz,2) - 4*a6*d6*lx*ny*pow(nx,2) + 4*d6*ny*px*pow(nx,2) - 2*mx*my*cos(3*alpha6)*pow(d6,2)*pow(nx,2) + 4*a6*d6*ly*pow(nx,3) - 4*d6*py*pow(nx,3) - 2*mu6*my*pow(d6,2)*pow(nx,3) + 4*a6*d6*ly*nx*pow(ny,2) - 4*d6*nx*py*pow(ny,2) - 2*mu6*my*nx*pow(d6,2)*pow(ny,2) + 2*mx*my*cos(3*alpha6)*pow(d6,2)*pow(ny,2) - 4*a6*d6*lx*pow(ny,3) + 4*d6*px*pow(ny,3) + 4*a6*d6*ly*nx*pow(nz,2) - 4*a6*d6*lx*ny*pow(nz,2) + 4*d6*ny*px*pow(nz,2) - 4*d6*nx*py*pow(nz,2) - 2*mu6*my*nx*pow(d6,2)*pow(nz,2) - 4*d6*cos(2*alpha6)*(ny*px*pow(mx,2) + ny*px*pow(my,2) + ny*px*pow(mz,2) + py*pow(nx,3) + nx*py*pow(ny,2) + nx*py*pow(nz,2) + a6*(ly*nx*(pow(mx,2) + pow(my,2) + pow(mz,2)) + lx*ny*(pow(nx,2) + pow(ny,2) + pow(nz,2)))) - 8*mu6*my*nx*pow(px,2) - 8*mu6*mx*ny*pow(px,2) + 8*mu6*my*nx*pow(py,2) + 8*mu6*mx*ny*pow(py,2) - 8*a6*d6*lx*mx*nx*ny*sin(2*alpha6) + 8*a6*d6*ly*my*nx*ny*sin(2*alpha6) + 8*a6*d6*ly*mz*nx*nz*sin(2*alpha6) - 8*a6*d6*lx*mz*ny*nz*sin(2*alpha6) + 8*d6*mx*nx*ny*px*sin(2*alpha6) + 8*d6*mz*ny*nz*px*sin(2*alpha6) - 8*d6*my*nx*ny*py*sin(2*alpha6) - 8*d6*mz*nx*nz*py*sin(2*alpha6) + 8*a6*d6*ly*mx*pow(nx,2)*sin(2*alpha6) - 8*d6*mx*py*pow(nx,2)*sin(2*alpha6) - 8*a6*d6*lx*my*pow(ny,2)*sin(2*alpha6) + 8*d6*my*px*pow(ny,2)*sin(2*alpha6) - my*nx*pow(d6,2)*pow(mx,2)*sin(3*alpha6) + ny*pow(d6,2)*pow(mx,3)*sin(3*alpha6) + mx*ny*pow(d6,2)*pow(my,2)*sin(3*alpha6) - nx*pow(d6,2)*pow(my,3)*sin(3*alpha6) - my*nx*pow(d6,2)*pow(mz,2)*sin(3*alpha6) + mx*ny*pow(d6,2)*pow(mz,2)*sin(3*alpha6) - mx*ny*pow(d6,2)*pow(nx,2)*sin(3*alpha6) + my*pow(d6,2)*pow(nx,3)*sin(3*alpha6) + my*nx*pow(d6,2)*pow(ny,2)*sin(3*alpha6) - mx*pow(d6,2)*pow(ny,3)*sin(3*alpha6) + my*nx*pow(d6,2)*pow(nz,2)*sin(3*alpha6) - mx*ny*pow(d6,2)*pow(nz,2)*sin(3*alpha6))) + a5*mu2*(-2*pow(mu6,2)*(-4*d1*my*mz*px + 4*d1*mx*mz*py + 4*my*mz*px*pz - 4*mx*mz*py*pz + 4*(-(ly*mx) + lx*my)*(lx*mx + ly*my + lz*mz)*pow(a6,2) + my*mz*nx*nz*pow(d6,2) - mx*mz*ny*nz*pow(d6,2) - 4*px*py*pow(mx,2) - nx*ny*pow(d6,2)*pow(mx,2) + 4*px*py*pow(my,2) + nx*ny*pow(d6,2)*pow(my,2) - 4*a6*(d1*(ly*mx - lx*my)*mz + 2*lx*mx*my*px + lz*my*mz*px - lz*mx*mz*py + lx*my*mz*pz - lx*py*pow(mx,2) + lx*py*pow(my,2) - ly*(2*mx*my*py + mx*mz*pz + px*pow(mx,2) - px*pow(my,2))) + 4*a1*lambda1*pow(mu1,-1)*(mx*mz*px + my*mz*py - pz*pow(mx,2) - pz*pow(my,2) + d1*(pow(mx,2) + pow(my,2)) + a6*(-((lx*mx + ly*my)*mz) + lz*(pow(mx,2) + pow(my,2)))) + mx*my*pow(d6,2)*pow(nx,2) - mx*my*pow(d6,2)*pow(ny,2) + 4*mx*my*pow(px,2) - 4*mx*my*pow(py,2)) - 8*a6*ly*mx*nx*px*sin(2*alpha6) + 8*a6*lx*my*nx*px*sin(2*alpha6) + 8*a6*lx*mx*ny*px*sin(2*alpha6) + 8*a6*ly*my*ny*px*sin(2*alpha6) - 8*a6*lx*mx*nx*py*sin(2*alpha6) - 8*a6*ly*my*nx*py*sin(2*alpha6) - 8*a6*ly*mx*ny*py*sin(2*alpha6) + 8*a6*lx*my*ny*py*sin(2*alpha6) + 8*mx*nx*px*py*sin(2*alpha6) - 8*my*ny*px*py*sin(2*alpha6) + 8*lx*ly*mx*nx*pow(a6,2)*sin(2*alpha6) - 8*lx*ly*my*ny*pow(a6,2)*sin(2*alpha6) - 8*a1*d1*lambda1*mx*nx*pow(mu1,-1)*sin(2*alpha6) - 8*a1*a6*lambda1*lz*mx*nx*pow(mu1,-1)*sin(2*alpha6) - 8*a1*d1*lambda1*my*ny*pow(mu1,-1)*sin(2*alpha6) - 8*a1*a6*lambda1*lz*my*ny*pow(mu1,-1)*sin(2*alpha6) + 8*a1*lambda1*mx*nx*pz*pow(mu1,-1)*sin(2*alpha6) + 8*a1*lambda1*my*ny*pz*pow(mu1,-1)*sin(2*alpha6) + my*nx*pow(d6,2)*pow(mx,2)*sin(2*alpha6) + nx*pow(d6,2)*pow(my,3)*sin(2*alpha6) + my*nx*pow(d6,2)*pow(mz,2)*sin(2*alpha6) + mx*ny*pow(d6,2)*pow(nx,2)*sin(2*alpha6) + mx*pow(d6,2)*pow(ny,3)*sin(2*alpha6) + mx*ny*pow(d6,2)*pow(nz,2)*sin(2*alpha6) + mu6*(8*d5*lambda5*(a6*ly*mx - a6*lx*my + my*px - mx*py) + d6*(-4*(a6*ly*mx - a6*lx*my + my*px - mx*py)*cos(2*alpha6)*(pow(mx,2) + pow(my,2) + pow(mz,2) - pow(nx,2) - pow(ny,2) - pow(nz,2)) + d6*(my*nx - mx*ny)*cos(3*alpha6)*(pow(mx,2) + pow(my,2) + pow(mz,2) - pow(nx,2) - pow(ny,2) - pow(nz,2)) + 2*(2*(a6*ly*mx - a6*lx*my + my*px - mx*py)*(pow(mx,2) + pow(my,2) + pow(mz,2) + pow(nx,2) + pow(ny,2) + pow(nz,2)) + 4*(mx*nx + my*ny + mz*nz)*(a6*ly*mx - a6*lx*my + my*px - mx*py)*sin(2*alpha6) + d6*(-(my*nx) + mx*ny)*(mx*nx + my*ny + mz*nz)*sin(3*alpha6)))) + a6*d6*lx*ny*pow(mu6,-1)*pow(mx,2)*sin(4*alpha6) + d6*nx*py*pow(mu6,-1)*pow(mx,2)*sin(4*alpha6) + a6*d6*lx*ny*pow(mu6,-1)*pow(my,2)*sin(4*alpha6) + d6*nx*py*pow(mu6,-1)*pow(my,2)*sin(4*alpha6) + a6*d6*lx*ny*pow(mu6,-1)*pow(mz,2)*sin(4*alpha6) + d6*nx*py*pow(mu6,-1)*pow(mz,2)*sin(4*alpha6) + d6*ny*px*pow(mu6,-1)*pow(nx,2)*sin(4*alpha6) + a6*d6*ly*pow(mu6,-1)*pow(nx,3)*sin(4*alpha6) + a6*d6*ly*nx*pow(mu6,-1)*pow(ny,2)*sin(4*alpha6) + d6*px*pow(mu6,-1)*pow(ny,3)*sin(4*alpha6) + a6*d6*ly*nx*pow(mu6,-1)*pow(nz,2)*sin(4*alpha6) + d6*ny*px*pow(mu6,-1)*pow(nz,2)*sin(4*alpha6))));

return(dummyexp);

}