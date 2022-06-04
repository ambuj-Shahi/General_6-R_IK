
# include <stdio.h>
# include <math.h>

double cx045(double* coeff)
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



double dummyexp = pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1)*(8*(a2 + a3)*a5*lambda5*mu2*mu3*(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)) + mu5*(8*d2*d3*mu3*mu6*(-(a6*ly*mx) + a6*lx*my - my*px + mx*py) + 2*lambda3*mu2*(my*nx - mx*ny)*(mx*nx + my*ny + mz*nz)*pow(d6,3)*pow(lambda6,4) - lambda3*mu2*pow(d6,2)*pow(lambda6,3)*(2*mx*my*nx*px + 2*my*mz*nz*px - 2*mx*my*ny*py - 2*mx*mz*nz*py - 3*ny*px*pow(mx,2) + nx*py*pow(mx,2) - ny*px*pow(my,2) + 3*nx*py*pow(my,2) - 3*ny*px*pow(mz,2) + 3*nx*py*pow(mz,2) + 3*ny*px*pow(nx,2) - 3*py*pow(nx,3) - 3*nx*py*pow(ny,2) + 3*px*pow(ny,3) + 2*d6*mu6*(-(my*nx) + mx*ny)*(pow(mx,2) + pow(my,2) + pow(mz,2) - pow(nx,2) - pow(ny,2) - pow(nz,2)) + 3*ny*px*pow(nz,2) - 3*nx*py*pow(nz,2) + a6*(lx*(-2*mx*my*nx - 2*my*mz*nz + 3*ny*pow(mx,2) + ny*pow(my,2) - 3*ny*(-pow(mz,2) + pow(nx,2) + pow(ny,2) + pow(nz,2))) + ly*(2*mx*(my*ny + mz*nz) - nx*pow(mx,2) + 3*nx*(-pow(my,2) - pow(mz,2) + pow(nx,2) + pow(ny,2) + pow(nz,2))))) + d6*lambda3*mu2*pow(lambda6,2)*(8*a1*lambda1*pow(mu1,-1)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))) + 2*(4*(ny*px - nx*py)*(-(d1*nz) + nx*px + ny*py + nz*pz) + 4*(-(ly*nx) + lx*ny)*(lx*nx + ly*ny + lz*nz)*pow(a6,2) + (-(my*nx) + mx*ny)*(mx*nx + my*ny + mz*nz)*pow(d6,2) - 4*a6*(d1*(ly*nx - lx*ny)*nz + 2*lx*nx*ny*px + lz*ny*nz*px - lz*nx*nz*py + lx*ny*nz*pz - lx*py*pow(nx,2) + lx*py*pow(ny,2) - ly*(2*nx*ny*py + nx*nz*pz + px*pow(nx,2) - px*pow(ny,2)))) + d6*mu6*(-10*mz*ny*nz*px + 10*my*nx*ny*py + 10*mz*nx*nz*py + my*px*pow(mx,2) - py*pow(mx,3) + px*pow(my,3) + my*px*pow(mz,2) - my*px*pow(nx,2) - 11*my*px*pow(ny,2) - my*px*pow(nz,2) + mx*(-10*nx*ny*px + 11*py*pow(nx,2) + py*(-pow(my,2) - pow(mz,2) + pow(ny,2) + pow(nz,2))) + a6*(ly*(-10*nx*(my*ny + mz*nz) + pow(mx,3) + mx*(pow(my,2) + pow(mz,2) - 11*pow(nx,2) - pow(ny,2) - pow(nz,2))) + lx*(10*mx*nx*ny + 10*mz*ny*nz - my*pow(mx,2) - pow(my,3) + my*(-pow(mz,2) + pow(nx,2) + 11*pow(ny,2) + pow(nz,2)))))) - lambda6*(8*d2*d3*mu3*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) + mu2*(-8*d3*(d4 + d5*lambda4)*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) + lambda3*(-8*a2*a3*a6*ly*nx + 8*a6*d1*d6*ly*mu6*mz*nx + 8*a2*a3*a6*lx*ny - 8*a6*d1*d6*lx*mu6*mz*ny + 8*a6*d1*d6*ly*mu6*mx*nz - 8*a6*d1*d6*lx*mu6*my*nz - 8*a2*a3*ny*px + 8*a6*d1*lz*ny*px + 8*d1*d6*mu6*mz*ny*px + 8*a6*d6*lz*mu6*mz*ny*px + 8*d1*d6*mu6*my*nz*px + 8*a6*d6*lz*mu6*my*nz*px + 8*a2*a3*nx*py - 8*a6*d1*lz*nx*py - 8*d1*d6*mu6*mz*nx*py - 8*a6*d6*lz*mu6*mz*nx*py - 8*d1*d6*mu6*mx*nz*py - 8*a6*d6*lz*mu6*mx*nz*py + 8*a6*lx*nx*px*py - 8*a6*ly*ny*px*py - 8*d4*d5*lambda4*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) - 8*a6*d1*ly*nx*pz - 8*a6*d6*ly*mu6*mz*nx*pz + 8*a6*d1*lx*ny*pz + 8*a6*d6*lx*mu6*mz*ny*pz - 8*a6*d6*ly*mu6*mx*nz*pz + 8*a6*d6*lx*mu6*my*nz*pz - 8*d1*ny*px*pz - 8*a6*lz*ny*px*pz - 8*d6*mu6*mz*ny*px*pz - 8*d6*mu6*my*nz*px*pz + 8*d1*nx*py*pz + 8*a6*lz*nx*py*pz + 8*d6*mu6*mz*nx*py*pz + 8*d6*mu6*mx*nz*py*pz + 4*a6*ly*nx*pow(a1,2) - 4*a6*lx*ny*pow(a1,2) + 4*ny*px*pow(a1,2) - 4*nx*py*pow(a1,2) - 4*a6*ly*nx*pow(a2,2) + 4*a6*lx*ny*pow(a2,2) - 4*ny*px*pow(a2,2) + 4*nx*py*pow(a2,2) - 4*a6*ly*nx*pow(a3,2) + 4*a6*lx*ny*pow(a3,2) - 4*ny*px*pow(a3,2) + 4*nx*py*pow(a3,2) + 4*a6*ly*nx*pow(a4,2) - 4*a6*lx*ny*pow(a4,2) + 4*ny*px*pow(a4,2) - 4*nx*py*pow(a4,2) - 4*a6*ly*nx*pow(a5,2) + 4*a6*lx*ny*pow(a5,2) - 4*ny*px*pow(a5,2) + 4*nx*py*pow(a5,2) + 8*d1*ly*lz*nx*pow(a6,2) + 8*d6*ly*lz*mu6*mz*nx*pow(a6,2) - 8*d1*lx*lz*ny*pow(a6,2) - 8*d6*lx*lz*mu6*mz*ny*pow(a6,2) + 8*d6*ly*lz*mu6*mx*nz*pow(a6,2) - 8*d6*lx*lz*mu6*my*nz*pow(a6,2) - 8*lx*ly*nx*px*pow(a6,2) + 8*lx*ly*ny*py*pow(a6,2) - 8*ly*lz*nx*pz*pow(a6,2) + 8*lx*lz*ny*pz*pow(a6,2) + 4*a6*ly*nx*pow(d1,2) - 4*a6*lx*ny*pow(d1,2) + 4*ny*px*pow(d1,2) - 4*nx*py*pow(d1,2) + 4*a6*ly*nx*pow(d2,2) - 4*a6*lx*ny*pow(d2,2) + 4*ny*px*pow(d2,2) - 4*nx*py*pow(d2,2) - 4*a6*ly*nx*pow(d3,2) + 4*a6*lx*ny*pow(d3,2) - 4*ny*px*pow(d3,2) + 4*nx*py*pow(d3,2) - 4*a6*ly*nx*pow(d4,2) + 4*a6*lx*ny*pow(d4,2) - 4*ny*px*pow(d4,2) + 4*nx*py*pow(d4,2) - 4*a6*ly*nx*pow(d5,2) + 4*a6*lx*ny*pow(d5,2) - 4*ny*px*pow(d5,2) + 4*nx*py*pow(d5,2) + 2*a6*lx*mx*my*nx*pow(d6,2) - 2*a6*ly*mx*my*ny*pow(d6,2) - 2*a6*ly*mx*mz*nz*pow(d6,2) + 2*a6*lx*my*mz*nz*pow(d6,2) - 2*mx*my*nx*px*pow(d6,2) - 2*my*mz*nz*px*pow(d6,2) + 2*mx*my*ny*py*pow(d6,2) + 2*mx*mz*nz*py*pow(d6,2) - 8*d6*mu6*my*nx*pow(a6,2)*pow(lx,2) - 8*d6*mu6*mx*ny*pow(a6,2)*pow(lx,2) + 12*ny*px*pow(a6,2)*pow(lx,2) - 4*nx*py*pow(a6,2)*pow(lx,2) + 4*ly*nx*pow(a6,3)*pow(lx,2) - 4*ny*pow(a6,3)*pow(lx,3) + 8*d6*mu6*my*nx*pow(a6,2)*pow(ly,2) + 8*d6*mu6*mx*ny*pow(a6,2)*pow(ly,2) + 4*ny*px*pow(a6,2)*pow(ly,2) - 12*nx*py*pow(a6,2)*pow(ly,2) - 4*lx*ny*pow(a6,3)*pow(ly,2) + 4*nx*pow(a6,3)*pow(ly,3) + 4*ny*px*pow(a6,2)*pow(lz,2) - 4*nx*py*pow(a6,2)*pow(lz,2) + 4*ly*nx*pow(a6,3)*pow(lz,2) - 4*lx*ny*pow(a6,3)*pow(lz,2) - 10*a6*lx*mx*my*nx*pow(d6,2)*pow(mu6,2) + 10*a6*ly*mx*my*ny*pow(d6,2)*pow(mu6,2) + 10*a6*ly*mx*mz*nz*pow(d6,2)*pow(mu6,2) - 10*a6*lx*my*mz*nz*pow(d6,2)*pow(mu6,2) + 10*mx*my*nx*px*pow(d6,2)*pow(mu6,2) + 10*my*mz*nz*px*pow(d6,2)*pow(mu6,2) - 10*mx*my*ny*py*pow(d6,2)*pow(mu6,2) - 10*mx*mz*nz*py*pow(d6,2)*pow(mu6,2) + a6*ly*nx*pow(d6,2)*pow(mx,2) - 3*a6*lx*ny*pow(d6,2)*pow(mx,2) + 3*ny*px*pow(d6,2)*pow(mx,2) - nx*py*pow(d6,2)*pow(mx,2) + 2*mu6*my*nx*pow(d6,3)*pow(mx,2) + 11*a6*ly*nx*pow(d6,2)*pow(mu6,2)*pow(mx,2) - a6*lx*ny*pow(d6,2)*pow(mu6,2)*pow(mx,2) + ny*px*pow(d6,2)*pow(mu6,2)*pow(mx,2) - 11*nx*py*pow(d6,2)*pow(mu6,2)*pow(mx,2) - 2*my*nx*pow(d6,3)*pow(mu6,3)*pow(mx,2) + 2*ny*pow(d6,3)*pow(mu6,3)*pow(mx,3) + 3*a6*ly*nx*pow(d6,2)*pow(my,2) - a6*lx*ny*pow(d6,2)*pow(my,2) + ny*px*pow(d6,2)*pow(my,2) - 3*nx*py*pow(d6,2)*pow(my,2) + a6*ly*nx*pow(d6,2)*pow(mu6,2)*pow(my,2) - 11*a6*lx*ny*pow(d6,2)*pow(mu6,2)*pow(my,2) + 11*ny*px*pow(d6,2)*pow(mu6,2)*pow(my,2) - nx*py*pow(d6,2)*pow(mu6,2)*pow(my,2) + 2*mx*ny*pow(d6,3)*pow(mu6,3)*pow(my,2) + 2*mu6*nx*pow(d6,3)*pow(my,3) - 2*nx*pow(d6,3)*pow(mu6,3)*pow(my,3) + 3*a6*ly*nx*pow(d6,2)*pow(mz,2) - 3*a6*lx*ny*pow(d6,2)*pow(mz,2) + 3*ny*px*pow(d6,2)*pow(mz,2) - 3*nx*py*pow(d6,2)*pow(mz,2) + 2*mu6*my*nx*pow(d6,3)*pow(mz,2) + a6*ly*nx*pow(d6,2)*pow(mu6,2)*pow(mz,2) - a6*lx*ny*pow(d6,2)*pow(mu6,2)*pow(mz,2) + ny*px*pow(d6,2)*pow(mu6,2)*pow(mz,2) - nx*py*pow(d6,2)*pow(mu6,2)*pow(mz,2) - 2*my*nx*pow(d6,3)*pow(mu6,3)*pow(mz,2) + 2*mx*ny*pow(d6,3)*pow(mu6,3)*pow(mz,2) - a6*lx*ny*pow(d6,2)*pow(nx,2) + ny*px*pow(d6,2)*pow(nx,2) + 2*mu6*mx*ny*pow(d6,3)*pow(nx,2) + a6*lx*ny*pow(d6,2)*pow(mu6,2)*pow(nx,2) - ny*px*pow(d6,2)*pow(mu6,2)*pow(nx,2) - 2*mx*ny*pow(d6,3)*pow(mu6,3)*pow(nx,2) + a6*ly*pow(d6,2)*pow(nx,3) - py*pow(d6,2)*pow(nx,3) - a6*ly*pow(d6,2)*pow(mu6,2)*pow(nx,3) + py*pow(d6,2)*pow(mu6,2)*pow(nx,3) + 2*my*pow(d6,3)*pow(mu6,3)*pow(nx,3) + a6*ly*nx*pow(d6,2)*pow(ny,2) - nx*py*pow(d6,2)*pow(ny,2) - a6*ly*nx*pow(d6,2)*pow(mu6,2)*pow(ny,2) + nx*py*pow(d6,2)*pow(mu6,2)*pow(ny,2) + 2*my*nx*pow(d6,3)*pow(mu6,3)*pow(ny,2) - a6*lx*pow(d6,2)*pow(ny,3) + px*pow(d6,2)*pow(ny,3) + 2*mu6*mx*pow(d6,3)*pow(ny,3) + a6*lx*pow(d6,2)*pow(mu6,2)*pow(ny,3) - px*pow(d6,2)*pow(mu6,2)*pow(ny,3) - 2*mx*pow(d6,3)*pow(mu6,3)*pow(ny,3) + a6*ly*nx*pow(d6,2)*pow(nz,2) - a6*lx*ny*pow(d6,2)*pow(nz,2) + ny*px*pow(d6,2)*pow(nz,2) - nx*py*pow(d6,2)*pow(nz,2) + 2*mu6*mx*ny*pow(d6,3)*pow(nz,2) - a6*ly*nx*pow(d6,2)*pow(mu6,2)*pow(nz,2) + a6*lx*ny*pow(d6,2)*pow(mu6,2)*pow(nz,2) - ny*px*pow(d6,2)*pow(mu6,2)*pow(nz,2) + nx*py*pow(d6,2)*pow(mu6,2)*pow(nz,2) + 2*my*nx*pow(d6,3)*pow(mu6,3)*pow(nz,2) - 2*mx*ny*pow(d6,3)*pow(mu6,3)*pow(nz,2) + 4*a6*ly*nx*pow(px,2) - 8*d6*mu6*my*nx*pow(px,2) - 12*a6*lx*ny*pow(px,2) - 8*d6*mu6*mx*ny*pow(px,2) - 4*nx*py*pow(px,2) + 4*ny*pow(px,3) + 12*a6*ly*nx*pow(py,2) + 8*d6*mu6*my*nx*pow(py,2) - 4*a6*lx*ny*pow(py,2) + 8*d6*mu6*mx*ny*pow(py,2) + 4*ny*px*pow(py,2) + 8*a1*lambda1*pow(mu1,-1)*(d1*nx*px + d1*ny*py + d6*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - nz*(mx*px + my*py) - mz*(nx*px + ny*py)) - nx*px*pz - ny*py*pz + a6*(-(d1*(lx*nx + ly*ny)) + lz*nx*px - 2*lx*nz*px + lz*ny*py - 2*ly*nz*py + lx*nx*pz + ly*ny*pz) + pow(a6,2)*(-(lx*lz*nx) + ly*(-(lz*ny) + ly*nz) + nz*pow(lx,2)) + nz*pow(px,2) + nz*pow(py,2)) - 4*nx*pow(py,3) + 4*a6*ly*nx*pow(pz,2) - 4*a6*lx*ny*pow(pz,2) + 4*ny*px*pow(pz,2) - 4*nx*py*pow(pz,2)))) + mu2*(2*lambda3*(-(my*nx) + mx*ny)*(mx*nx + my*ny + mz*nz)*pow(d6,3)*pow(mu6,4) - lambda3*pow(d6,2)*pow(mu6,3)*(2*mz*ny*nz*px - 2*my*nx*ny*py - 2*mz*nx*nz*py + 3*my*px*pow(mx,2) - 3*py*pow(mx,3) + 3*px*pow(my,3) + 3*my*px*pow(mz,2) - 3*my*px*pow(nx,2) - my*px*pow(ny,2) - 3*my*px*pow(nz,2) + mx*(2*nx*ny*px + py*pow(nx,2) + 3*py*(-pow(my,2) - pow(mz,2) + pow(ny,2) + pow(nz,2))) + a6*(ly*(2*nx*(my*ny + mz*nz) + 3*pow(mx,3) + mx*(3*pow(my,2) + 3*pow(mz,2) - pow(nx,2) - 3*pow(ny,2) - 3*pow(nz,2))) + lx*(-2*mx*nx*ny - 2*mz*ny*nz - 3*my*pow(mx,2) - 3*pow(my,3) + my*(-3*pow(mz,2) + 3*pow(nx,2) + pow(ny,2) + 3*pow(nz,2))))) - 2*d6*lambda3*pow(mu6,2)*(4*d1*my*mz*px - 4*d1*mx*mz*py - 4*my*mz*px*pz + 4*mx*mz*py*pz + 4*(ly*mx - lx*my)*(lx*mx + ly*my + lz*mz)*pow(a6,2) - my*mz*nx*nz*pow(d6,2) + mx*mz*ny*nz*pow(d6,2) + 4*px*py*pow(mx,2) + nx*ny*pow(d6,2)*pow(mx,2) - 4*px*py*pow(my,2) - nx*ny*pow(d6,2)*pow(my,2) + 4*a6*(d1*(ly*mx - lx*my)*mz + 2*lx*mx*my*px + lz*my*mz*px - lz*mx*mz*py + lx*my*mz*pz - lx*py*pow(mx,2) + lx*py*pow(my,2) - ly*(2*mx*my*py + mx*mz*pz + px*pow(mx,2) - px*pow(my,2))) - 4*a1*lambda1*pow(mu1,-1)*(mx*mz*px + my*mz*py - pz*pow(mx,2) - pz*pow(my,2) + d1*(pow(mx,2) + pow(my,2)) + a6*(-((lx*mx + ly*my)*mz) + lz*(pow(mx,2) + pow(my,2)))) - mx*my*pow(d6,2)*pow(nx,2) + mx*my*pow(d6,2)*pow(ny,2) - 4*mx*my*pow(px,2) + 4*mx*my*pow(py,2)) - mu6*(-8*d3*(d4 + d5*lambda4)*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda3*(-8*a2*a3*a6*ly*mx + 8*a2*a3*a6*lx*my - 8*a2*a3*my*px + 8*a6*d1*lz*my*px + 8*a2*a3*mx*py - 8*a6*d1*lz*mx*py + 8*a6*lx*mx*px*py - 8*a6*ly*my*px*py - 8*d4*d5*lambda4*(a6*ly*mx - a6*lx*my + my*px - mx*py) - 8*a6*d1*ly*mx*pz + 8*a6*d1*lx*my*pz - 8*d1*my*px*pz - 8*a6*lz*my*px*pz + 8*d1*mx*py*pz + 8*a6*lz*mx*py*pz + 4*a6*ly*mx*pow(a1,2) - 4*a6*lx*my*pow(a1,2) + 4*my*px*pow(a1,2) - 4*mx*py*pow(a1,2) - 4*a6*ly*mx*pow(a2,2) + 4*a6*lx*my*pow(a2,2) - 4*my*px*pow(a2,2) + 4*mx*py*pow(a2,2) - 4*a6*ly*mx*pow(a3,2) + 4*a6*lx*my*pow(a3,2) - 4*my*px*pow(a3,2) + 4*mx*py*pow(a3,2) + 4*a6*ly*mx*pow(a4,2) - 4*a6*lx*my*pow(a4,2) + 4*my*px*pow(a4,2) - 4*mx*py*pow(a4,2) - 4*a6*ly*mx*pow(a5,2) + 4*a6*lx*my*pow(a5,2) - 4*my*px*pow(a5,2) + 4*mx*py*pow(a5,2) + 8*d1*ly*lz*mx*pow(a6,2) - 8*d1*lx*lz*my*pow(a6,2) - 8*lx*ly*mx*px*pow(a6,2) + 8*lx*ly*my*py*pow(a6,2) - 8*ly*lz*mx*pz*pow(a6,2) + 8*lx*lz*my*pz*pow(a6,2) + 4*a6*ly*mx*pow(d1,2) - 4*a6*lx*my*pow(d1,2) + 4*my*px*pow(d1,2) - 4*mx*py*pow(d1,2) + 4*a6*ly*mx*pow(d2,2) - 4*a6*lx*my*pow(d2,2) + 4*my*px*pow(d2,2) - 4*mx*py*pow(d2,2) - 4*a6*ly*mx*pow(d3,2) + 4*a6*lx*my*pow(d3,2) - 4*my*px*pow(d3,2) + 4*mx*py*pow(d3,2) - 4*a6*ly*mx*pow(d4,2) + 4*a6*lx*my*pow(d4,2) - 4*my*px*pow(d4,2) + 4*mx*py*pow(d4,2) - 4*a6*ly*mx*pow(d5,2) + 4*a6*lx*my*pow(d5,2) - 4*my*px*pow(d5,2) + 4*mx*py*pow(d5,2) + 2*a6*lx*mx*nx*ny*pow(d6,2) - 2*a6*ly*my*nx*ny*pow(d6,2) - 2*a6*ly*mz*nx*nz*pow(d6,2) + 2*a6*lx*mz*ny*nz*pow(d6,2) - 2*mx*nx*ny*px*pow(d6,2) - 2*mz*ny*nz*px*pow(d6,2) + 2*my*nx*ny*py*pow(d6,2) + 2*mz*nx*nz*py*pow(d6,2) + 12*my*px*pow(a6,2)*pow(lx,2) - 4*mx*py*pow(a6,2)*pow(lx,2) + 4*ly*mx*pow(a6,3)*pow(lx,2) - 4*my*pow(a6,3)*pow(lx,3) + 4*my*px*pow(a6,2)*pow(ly,2) - 12*mx*py*pow(a6,2)*pow(ly,2) - 4*lx*my*pow(a6,3)*pow(ly,2) + 4*mx*pow(a6,3)*pow(ly,3) + 4*my*px*pow(a6,2)*pow(lz,2) - 4*mx*py*pow(a6,2)*pow(lz,2) + 4*ly*mx*pow(a6,3)*pow(lz,2) - 4*lx*my*pow(a6,3)*pow(lz,2) - a6*lx*my*pow(d6,2)*pow(mx,2) + my*px*pow(d6,2)*pow(mx,2) + a6*ly*pow(d6,2)*pow(mx,3) - py*pow(d6,2)*pow(mx,3) + a6*ly*mx*pow(d6,2)*pow(my,2) - mx*py*pow(d6,2)*pow(my,2) - a6*lx*pow(d6,2)*pow(my,3) + px*pow(d6,2)*pow(my,3) + a6*ly*mx*pow(d6,2)*pow(mz,2) - a6*lx*my*pow(d6,2)*pow(mz,2) + my*px*pow(d6,2)*pow(mz,2) - mx*py*pow(d6,2)*pow(mz,2) + a6*ly*mx*pow(d6,2)*pow(nx,2) - 3*a6*lx*my*pow(d6,2)*pow(nx,2) + 3*my*px*pow(d6,2)*pow(nx,2) - mx*py*pow(d6,2)*pow(nx,2) + 3*a6*ly*mx*pow(d6,2)*pow(ny,2) - a6*lx*my*pow(d6,2)*pow(ny,2) + my*px*pow(d6,2)*pow(ny,2) - 3*mx*py*pow(d6,2)*pow(ny,2) + 3*a6*ly*mx*pow(d6,2)*pow(nz,2) - 3*a6*lx*my*pow(d6,2)*pow(nz,2) + 3*my*px*pow(d6,2)*pow(nz,2) - 3*mx*py*pow(d6,2)*pow(nz,2) + 4*a6*ly*mx*pow(px,2) - 12*a6*lx*my*pow(px,2) - 4*mx*py*pow(px,2) + 4*my*pow(px,3) + 12*a6*ly*mx*pow(py,2) - 4*a6*lx*my*pow(py,2) + 4*my*px*pow(py,2) + 8*a1*lambda1*pow(mu1,-1)*(d1*mx*px + d1*my*py - mx*px*pz - my*py*pz + a6*(-(d1*(lx*mx + ly*my)) + lz*mx*px - 2*lx*mz*px + lz*my*py - 2*ly*mz*py + lx*mx*pz + ly*my*pz) + pow(a6,2)*(-(lx*lz*mx) + ly*(-(lz*my) + ly*mz) + mz*pow(lx,2)) + mz*pow(px,2) + mz*pow(py,2)) - 4*mx*pow(py,3) + 4*a6*ly*mx*pow(pz,2) - 4*a6*lx*my*pow(pz,2) + 4*my*px*pow(pz,2) - 4*mx*py*pow(pz,2))) + d6*lambda3*(-8*mx*nx*px*py + 8*my*ny*px*py + 8*a6*(ly*(mx*nx*px - my*ny*px + my*nx*py + mx*ny*py) - lx*(my*nx*px + mx*ny*px - mx*nx*py + my*ny*py)) + 8*lx*ly*(-(mx*nx) + my*ny)*pow(a6,2) + 8*a1*lambda1*(mx*nx + my*ny)*(d1 + a6*lz - pz)*pow(mu1,-1) + ny*pow(d6,2)*pow(mx,3) + mx*ny*pow(d6,2)*pow(my,2) + mx*ny*pow(d6,2)*pow(mz,2) + my*pow(d6,2)*pow(nx,3) + my*nx*pow(d6,2)*pow(ny,2) + my*nx*pow(d6,2)*pow(nz,2))*sin(2*alpha6))) + lambda2*(-8*(a2 + a3)*a5*lambda3*lambda5*(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)) + mu3*mu5*(2*(my*nx - mx*ny)*(mx*nx + my*ny + mz*nz)*pow(d6,3)*pow(lambda6,4) + 2*(-(my*nx) + mx*ny)*(mx*nx + my*ny + mz*nz)*pow(d6,3)*pow(mu6,4) - pow(d6,2)*pow(lambda6,3)*(2*mx*my*nx*px + 2*my*mz*nz*px - 2*mx*my*ny*py - 2*mx*mz*nz*py - 3*ny*px*pow(mx,2) + nx*py*pow(mx,2) - ny*px*pow(my,2) + 3*nx*py*pow(my,2) - 3*ny*px*pow(mz,2) + 3*nx*py*pow(mz,2) + 3*ny*px*pow(nx,2) - 3*py*pow(nx,3) - 3*nx*py*pow(ny,2) + 3*px*pow(ny,3) + 2*d6*mu6*(-(my*nx) + mx*ny)*(pow(mx,2) + pow(my,2) + pow(mz,2) - pow(nx,2) - pow(ny,2) - pow(nz,2)) + 3*ny*px*pow(nz,2) - 3*nx*py*pow(nz,2) + a6*(lx*(-2*mx*my*nx - 2*my*mz*nz + 3*ny*pow(mx,2) + ny*pow(my,2) - 3*ny*(-pow(mz,2) + pow(nx,2) + pow(ny,2) + pow(nz,2))) + ly*(2*mx*(my*ny + mz*nz) - nx*pow(mx,2) + 3*nx*(-pow(my,2) - pow(mz,2) + pow(nx,2) + pow(ny,2) + pow(nz,2))))) - pow(d6,2)*pow(mu6,3)*(2*mz*ny*nz*px - 2*my*nx*ny*py - 2*mz*nx*nz*py + 3*my*px*pow(mx,2) - 3*py*pow(mx,3) + 3*px*pow(my,3) + 3*my*px*pow(mz,2) - 3*my*px*pow(nx,2) - my*px*pow(ny,2) - 3*my*px*pow(nz,2) + mx*(2*nx*ny*px + py*pow(nx,2) + 3*py*(-pow(my,2) - pow(mz,2) + pow(ny,2) + pow(nz,2))) + a6*(ly*(2*nx*(my*ny + mz*nz) + 3*pow(mx,3) + mx*(3*pow(my,2) + 3*pow(mz,2) - pow(nx,2) - 3*pow(ny,2) - 3*pow(nz,2))) + lx*(-2*mx*nx*ny - 2*mz*ny*nz - 3*my*pow(mx,2) - 3*pow(my,3) + my*(-3*pow(mz,2) + 3*pow(nx,2) + pow(ny,2) + 3*pow(nz,2))))) + d6*pow(lambda6,2)*(8*a1*lambda1*pow(mu1,-1)*(nx*nz*px + ny*nz*py - pz*pow(nx,2) - pz*pow(ny,2) + d1*(pow(nx,2) + pow(ny,2)) + a6*(-((lx*nx + ly*ny)*nz) + lz*(pow(nx,2) + pow(ny,2)))) + 2*(4*(ny*px - nx*py)*(-(d1*nz) + nx*px + ny*py + nz*pz) + 4*(-(ly*nx) + lx*ny)*(lx*nx + ly*ny + lz*nz)*pow(a6,2) + (-(my*nx) + mx*ny)*(mx*nx + my*ny + mz*nz)*pow(d6,2) - 4*a6*(d1*(ly*nx - lx*ny)*nz + 2*lx*nx*ny*px + lz*ny*nz*px - lz*nx*nz*py + lx*ny*nz*pz - lx*py*pow(nx,2) + lx*py*pow(ny,2) - ly*(2*nx*ny*py + nx*nz*pz + px*pow(nx,2) - px*pow(ny,2)))) + d6*mu6*(-10*mz*ny*nz*px + 10*my*nx*ny*py + 10*mz*nx*nz*py + my*px*pow(mx,2) - py*pow(mx,3) + px*pow(my,3) + my*px*pow(mz,2) - my*px*pow(nx,2) - 11*my*px*pow(ny,2) - my*px*pow(nz,2) + mx*(-10*nx*ny*px + 11*py*pow(nx,2) + py*(-pow(my,2) - pow(mz,2) + pow(ny,2) + pow(nz,2))) + a6*(ly*(-10*nx*(my*ny + mz*nz) + pow(mx,3) + mx*(pow(my,2) + pow(mz,2) - 11*pow(nx,2) - pow(ny,2) - pow(nz,2))) + lx*(10*mx*nx*ny + 10*mz*ny*nz - my*pow(mx,2) - pow(my,3) + my*(-pow(mz,2) + pow(nx,2) + 11*pow(ny,2) + pow(nz,2)))))) + 2*d6*pow(mu6,2)*(-4*d1*my*mz*px + 4*d1*mx*mz*py + 4*my*mz*px*pz - 4*mx*mz*py*pz + 4*(-(ly*mx) + lx*my)*(lx*mx + ly*my + lz*mz)*pow(a6,2) + my*mz*nx*nz*pow(d6,2) - mx*mz*ny*nz*pow(d6,2) - 4*px*py*pow(mx,2) - nx*ny*pow(d6,2)*pow(mx,2) + 4*px*py*pow(my,2) + nx*ny*pow(d6,2)*pow(my,2) - 4*a6*(d1*(ly*mx - lx*my)*mz + 2*lx*mx*my*px + lz*my*mz*px - lz*mx*mz*py + lx*my*mz*pz - lx*py*pow(mx,2) + lx*py*pow(my,2) - ly*(2*mx*my*py + mx*mz*pz + px*pow(mx,2) - px*pow(my,2))) + 4*a1*lambda1*pow(mu1,-1)*(mx*mz*px + my*mz*py - pz*pow(mx,2) - pz*pow(my,2) + d1*(pow(mx,2) + pow(my,2)) + a6*(-((lx*mx + ly*my)*mz) + lz*(pow(mx,2) + pow(my,2)))) + mx*my*pow(d6,2)*pow(nx,2) - mx*my*pow(d6,2)*pow(ny,2) + 4*mx*my*pow(px,2) - 4*mx*my*pow(py,2)) - lambda6*(-8*a2*a3*a6*ly*nx + 8*a6*d1*d6*ly*mu6*mz*nx + 8*a2*a3*a6*lx*ny - 8*a6*d1*d6*lx*mu6*mz*ny + 8*a6*d1*d6*ly*mu6*mx*nz - 8*a6*d1*d6*lx*mu6*my*nz - 8*a2*a3*ny*px + 8*a6*d1*lz*ny*px + 8*d1*d6*mu6*mz*ny*px + 8*a6*d6*lz*mu6*mz*ny*px + 8*d1*d6*mu6*my*nz*px + 8*a6*d6*lz*mu6*my*nz*px + 8*a2*a3*nx*py - 8*a6*d1*lz*nx*py - 8*d1*d6*mu6*mz*nx*py - 8*a6*d6*lz*mu6*mz*nx*py - 8*d1*d6*mu6*mx*nz*py - 8*a6*d6*lz*mu6*mx*nz*py + 8*a6*lx*nx*px*py - 8*a6*ly*ny*px*py - 8*d4*d5*lambda4*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) - 8*a6*d1*ly*nx*pz - 8*a6*d6*ly*mu6*mz*nx*pz + 8*a6*d1*lx*ny*pz + 8*a6*d6*lx*mu6*mz*ny*pz - 8*a6*d6*ly*mu6*mx*nz*pz + 8*a6*d6*lx*mu6*my*nz*pz - 8*d1*ny*px*pz - 8*a6*lz*ny*px*pz - 8*d6*mu6*mz*ny*px*pz - 8*d6*mu6*my*nz*px*pz + 8*d1*nx*py*pz + 8*a6*lz*nx*py*pz + 8*d6*mu6*mz*nx*py*pz + 8*d6*mu6*mx*nz*py*pz + 4*a6*ly*nx*pow(a1,2) - 4*a6*lx*ny*pow(a1,2) + 4*ny*px*pow(a1,2) - 4*nx*py*pow(a1,2) - 4*a6*ly*nx*pow(a2,2) + 4*a6*lx*ny*pow(a2,2) - 4*ny*px*pow(a2,2) + 4*nx*py*pow(a2,2) - 4*a6*ly*nx*pow(a3,2) + 4*a6*lx*ny*pow(a3,2) - 4*ny*px*pow(a3,2) + 4*nx*py*pow(a3,2) + 4*a6*ly*nx*pow(a4,2) - 4*a6*lx*ny*pow(a4,2) + 4*ny*px*pow(a4,2) - 4*nx*py*pow(a4,2) - 4*a6*ly*nx*pow(a5,2) + 4*a6*lx*ny*pow(a5,2) - 4*ny*px*pow(a5,2) + 4*nx*py*pow(a5,2) + 8*d1*ly*lz*nx*pow(a6,2) + 8*d6*ly*lz*mu6*mz*nx*pow(a6,2) - 8*d1*lx*lz*ny*pow(a6,2) - 8*d6*lx*lz*mu6*mz*ny*pow(a6,2) + 8*d6*ly*lz*mu6*mx*nz*pow(a6,2) - 8*d6*lx*lz*mu6*my*nz*pow(a6,2) - 8*lx*ly*nx*px*pow(a6,2) + 8*lx*ly*ny*py*pow(a6,2) - 8*ly*lz*nx*pz*pow(a6,2) + 8*lx*lz*ny*pz*pow(a6,2) + 4*a6*ly*nx*pow(d1,2) - 4*a6*lx*ny*pow(d1,2) + 4*ny*px*pow(d1,2) - 4*nx*py*pow(d1,2) + 4*a6*ly*nx*pow(d2,2) - 4*a6*lx*ny*pow(d2,2) + 4*ny*px*pow(d2,2) - 4*nx*py*pow(d2,2) + 4*a6*ly*nx*pow(d3,2) - 4*a6*lx*ny*pow(d3,2) + 4*ny*px*pow(d3,2) - 4*nx*py*pow(d3,2) - 4*a6*ly*nx*pow(d4,2) + 4*a6*lx*ny*pow(d4,2) - 4*ny*px*pow(d4,2) + 4*nx*py*pow(d4,2) - 4*a6*ly*nx*pow(d5,2) + 4*a6*lx*ny*pow(d5,2) - 4*ny*px*pow(d5,2) + 4*nx*py*pow(d5,2) + 2*a6*lx*mx*my*nx*pow(d6,2) - 2*a6*ly*mx*my*ny*pow(d6,2) - 2*a6*ly*mx*mz*nz*pow(d6,2) + 2*a6*lx*my*mz*nz*pow(d6,2) - 2*mx*my*nx*px*pow(d6,2) - 2*my*mz*nz*px*pow(d6,2) + 2*mx*my*ny*py*pow(d6,2) + 2*mx*mz*nz*py*pow(d6,2) - 8*d6*mu6*my*nx*pow(a6,2)*pow(lx,2) - 8*d6*mu6*mx*ny*pow(a6,2)*pow(lx,2) + 12*ny*px*pow(a6,2)*pow(lx,2) - 4*nx*py*pow(a6,2)*pow(lx,2) + 4*ly*nx*pow(a6,3)*pow(lx,2) - 4*ny*pow(a6,3)*pow(lx,3) + 8*d6*mu6*my*nx*pow(a6,2)*pow(ly,2) + 8*d6*mu6*mx*ny*pow(a6,2)*pow(ly,2) + 4*ny*px*pow(a6,2)*pow(ly,2) - 12*nx*py*pow(a6,2)*pow(ly,2) - 4*lx*ny*pow(a6,3)*pow(ly,2) + 4*nx*pow(a6,3)*pow(ly,3) + 4*ny*px*pow(a6,2)*pow(lz,2) - 4*nx*py*pow(a6,2)*pow(lz,2) + 4*ly*nx*pow(a6,3)*pow(lz,2) - 4*lx*ny*pow(a6,3)*pow(lz,2) + a6*ly*nx*pow(d6,2)*pow(mx,2) - 3*a6*lx*ny*pow(d6,2)*pow(mx,2) + 3*ny*px*pow(d6,2)*pow(mx,2) - nx*py*pow(d6,2)*pow(mx,2) + 2*mu6*my*nx*pow(d6,3)*pow(mx,2) + 11*a6*ly*nx*pow(d6,2)*pow(mu6,2)*pow(mx,2) - a6*lx*ny*pow(d6,2)*pow(mu6,2)*pow(mx,2) + ny*px*pow(d6,2)*pow(mu6,2)*pow(mx,2) - 11*nx*py*pow(d6,2)*pow(mu6,2)*pow(mx,2) - 2*my*nx*pow(d6,3)*pow(mu6,3)*pow(mx,2) + 2*ny*pow(d6,3)*pow(mu6,3)*pow(mx,3) + 3*a6*ly*nx*pow(d6,2)*pow(my,2) - a6*lx*ny*pow(d6,2)*pow(my,2) + ny*px*pow(d6,2)*pow(my,2) - 3*nx*py*pow(d6,2)*pow(my,2) + a6*ly*nx*pow(d6,2)*pow(mu6,2)*pow(my,2) - 11*a6*lx*ny*pow(d6,2)*pow(mu6,2)*pow(my,2) + 11*ny*px*pow(d6,2)*pow(mu6,2)*pow(my,2) - nx*py*pow(d6,2)*pow(mu6,2)*pow(my,2) + 2*mx*ny*pow(d6,3)*pow(mu6,3)*pow(my,2) + 2*mu6*nx*pow(d6,3)*pow(my,3) - 2*nx*pow(d6,3)*pow(mu6,3)*pow(my,3) + 3*a6*ly*nx*pow(d6,2)*pow(mz,2) - 3*a6*lx*ny*pow(d6,2)*pow(mz,2) + 3*ny*px*pow(d6,2)*pow(mz,2) - 3*nx*py*pow(d6,2)*pow(mz,2) + 2*mu6*my*nx*pow(d6,3)*pow(mz,2) + a6*ly*nx*pow(d6,2)*pow(mu6,2)*pow(mz,2) - a6*lx*ny*pow(d6,2)*pow(mu6,2)*pow(mz,2) + ny*px*pow(d6,2)*pow(mu6,2)*pow(mz,2) - nx*py*pow(d6,2)*pow(mu6,2)*pow(mz,2) - 2*my*nx*pow(d6,3)*pow(mu6,3)*pow(mz,2) + 2*mx*ny*pow(d6,3)*pow(mu6,3)*pow(mz,2) - a6*lx*ny*pow(d6,2)*pow(nx,2) + ny*px*pow(d6,2)*pow(nx,2) + 2*mu6*mx*ny*pow(d6,3)*pow(nx,2) + a6*lx*ny*pow(d6,2)*pow(mu6,2)*pow(nx,2) - ny*px*pow(d6,2)*pow(mu6,2)*pow(nx,2) - 2*mx*ny*pow(d6,3)*pow(mu6,3)*pow(nx,2) + a6*ly*pow(d6,2)*pow(nx,3) - py*pow(d6,2)*pow(nx,3) - a6*ly*pow(d6,2)*pow(mu6,2)*pow(nx,3) + py*pow(d6,2)*pow(mu6,2)*pow(nx,3) + 2*my*pow(d6,3)*pow(mu6,3)*pow(nx,3) + a6*ly*nx*pow(d6,2)*pow(ny,2) - nx*py*pow(d6,2)*pow(ny,2) - a6*ly*nx*pow(d6,2)*pow(mu6,2)*pow(ny,2) + nx*py*pow(d6,2)*pow(mu6,2)*pow(ny,2) + 2*my*nx*pow(d6,3)*pow(mu6,3)*pow(ny,2) - a6*lx*pow(d6,2)*pow(ny,3) + px*pow(d6,2)*pow(ny,3) + 2*mu6*mx*pow(d6,3)*pow(ny,3) + a6*lx*pow(d6,2)*pow(mu6,2)*pow(ny,3) - px*pow(d6,2)*pow(mu6,2)*pow(ny,3) - 2*mx*pow(d6,3)*pow(mu6,3)*pow(ny,3) + a6*ly*nx*pow(d6,2)*pow(nz,2) - a6*lx*ny*pow(d6,2)*pow(nz,2) + ny*px*pow(d6,2)*pow(nz,2) - nx*py*pow(d6,2)*pow(nz,2) + 2*mu6*mx*ny*pow(d6,3)*pow(nz,2) - a6*ly*nx*pow(d6,2)*pow(mu6,2)*pow(nz,2) + a6*lx*ny*pow(d6,2)*pow(mu6,2)*pow(nz,2) - ny*px*pow(d6,2)*pow(mu6,2)*pow(nz,2) + nx*py*pow(d6,2)*pow(mu6,2)*pow(nz,2) + 2*my*nx*pow(d6,3)*pow(mu6,3)*pow(nz,2) - 2*mx*ny*pow(d6,3)*pow(mu6,3)*pow(nz,2) + 4*a6*ly*nx*pow(px,2) - 8*d6*mu6*my*nx*pow(px,2) - 12*a6*lx*ny*pow(px,2) - 8*d6*mu6*mx*ny*pow(px,2) - 4*nx*py*pow(px,2) + 4*ny*pow(px,3) + 12*a6*ly*nx*pow(py,2) + 8*d6*mu6*my*nx*pow(py,2) - 4*a6*lx*ny*pow(py,2) + 8*d6*mu6*mx*ny*pow(py,2) + 4*ny*px*pow(py,2) + 8*a1*lambda1*pow(mu1,-1)*(d1*nx*px + d1*ny*py + d6*mu6*(a6*(lx*mz*nx + ly*mz*ny + lx*mx*nz + ly*my*nz) - nz*(mx*px + my*py) - mz*(nx*px + ny*py)) - nx*px*pz - ny*py*pz + a6*(-(d1*(lx*nx + ly*ny)) + lz*nx*px - 2*lx*nz*px + lz*ny*py - 2*ly*nz*py + lx*nx*pz + ly*ny*pz) + pow(a6,2)*(-(lx*lz*nx) + ly*(-(lz*ny) + ly*nz) + nz*pow(lx,2)) + nz*pow(px,2) + nz*pow(py,2)) - 4*nx*pow(py,3) + 4*a6*ly*nx*pow(pz,2) - 4*a6*lx*ny*pow(pz,2) + 4*ny*px*pow(pz,2) - 4*nx*py*pow(pz,2)) + d6*(-8*mx*nx*px*py + 8*my*ny*px*py + 8*a6*(ly*(mx*nx*px - my*ny*px + my*nx*py + mx*ny*py) - lx*(my*nx*px + mx*ny*px - mx*nx*py + my*ny*py)) + 8*lx*ly*(-(mx*nx) + my*ny)*pow(a6,2) + 8*a1*lambda1*(mx*nx + my*ny)*(d1 + a6*lz - pz)*pow(mu1,-1) + ny*pow(d6,2)*pow(mx,3) + mx*ny*pow(d6,2)*pow(my,2) + mx*ny*pow(d6,2)*pow(mz,2) + my*pow(d6,2)*pow(nx,3) + my*nx*pow(d6,2)*pow(ny,2) + my*nx*pow(d6,2)*pow(nz,2))*sin(2*alpha6) - mu6*(-8*a2*a3*a6*ly*mx + 8*a2*a3*a6*lx*my - 8*a2*a3*my*px + 8*a6*d1*lz*my*px + 8*a2*a3*mx*py - 8*a6*d1*lz*mx*py + 8*a6*lx*mx*px*py - 8*a6*ly*my*px*py - 8*d4*d5*lambda4*(a6*ly*mx - a6*lx*my + my*px - mx*py) - 8*a6*d1*ly*mx*pz + 8*a6*d1*lx*my*pz - 8*d1*my*px*pz - 8*a6*lz*my*px*pz + 8*d1*mx*py*pz + 8*a6*lz*mx*py*pz + 4*a6*ly*mx*pow(a1,2) - 4*a6*lx*my*pow(a1,2) + 4*my*px*pow(a1,2) - 4*mx*py*pow(a1,2) - 4*a6*ly*mx*pow(a2,2) + 4*a6*lx*my*pow(a2,2) - 4*my*px*pow(a2,2) + 4*mx*py*pow(a2,2) - 4*a6*ly*mx*pow(a3,2) + 4*a6*lx*my*pow(a3,2) - 4*my*px*pow(a3,2) + 4*mx*py*pow(a3,2) + 4*a6*ly*mx*pow(a4,2) - 4*a6*lx*my*pow(a4,2) + 4*my*px*pow(a4,2) - 4*mx*py*pow(a4,2) - 4*a6*ly*mx*pow(a5,2) + 4*a6*lx*my*pow(a5,2) - 4*my*px*pow(a5,2) + 4*mx*py*pow(a5,2) + 8*d1*ly*lz*mx*pow(a6,2) - 8*d1*lx*lz*my*pow(a6,2) - 8*lx*ly*mx*px*pow(a6,2) + 8*lx*ly*my*py*pow(a6,2) - 8*ly*lz*mx*pz*pow(a6,2) + 8*lx*lz*my*pz*pow(a6,2) + 4*a6*ly*mx*pow(d1,2) - 4*a6*lx*my*pow(d1,2) + 4*my*px*pow(d1,2) - 4*mx*py*pow(d1,2) + 4*a6*ly*mx*pow(d2,2) - 4*a6*lx*my*pow(d2,2) + 4*my*px*pow(d2,2) - 4*mx*py*pow(d2,2) + 4*a6*ly*mx*pow(d3,2) - 4*a6*lx*my*pow(d3,2) + 4*my*px*pow(d3,2) - 4*mx*py*pow(d3,2) - 4*a6*ly*mx*pow(d4,2) + 4*a6*lx*my*pow(d4,2) - 4*my*px*pow(d4,2) + 4*mx*py*pow(d4,2) - 4*a6*ly*mx*pow(d5,2) + 4*a6*lx*my*pow(d5,2) - 4*my*px*pow(d5,2) + 4*mx*py*pow(d5,2) + 2*a6*lx*mx*nx*ny*pow(d6,2) - 2*a6*ly*my*nx*ny*pow(d6,2) - 2*a6*ly*mz*nx*nz*pow(d6,2) + 2*a6*lx*mz*ny*nz*pow(d6,2) - 2*mx*nx*ny*px*pow(d6,2) - 2*mz*ny*nz*px*pow(d6,2) + 2*my*nx*ny*py*pow(d6,2) + 2*mz*nx*nz*py*pow(d6,2) + 12*my*px*pow(a6,2)*pow(lx,2) - 4*mx*py*pow(a6,2)*pow(lx,2) + 4*ly*mx*pow(a6,3)*pow(lx,2) - 4*my*pow(a6,3)*pow(lx,3) + 4*my*px*pow(a6,2)*pow(ly,2) - 12*mx*py*pow(a6,2)*pow(ly,2) - 4*lx*my*pow(a6,3)*pow(ly,2) + 4*mx*pow(a6,3)*pow(ly,3) + 4*my*px*pow(a6,2)*pow(lz,2) - 4*mx*py*pow(a6,2)*pow(lz,2) + 4*ly*mx*pow(a6,3)*pow(lz,2) - 4*lx*my*pow(a6,3)*pow(lz,2) - a6*lx*my*pow(d6,2)*pow(mx,2) + my*px*pow(d6,2)*pow(mx,2) + a6*ly*pow(d6,2)*pow(mx,3) - py*pow(d6,2)*pow(mx,3) + a6*ly*mx*pow(d6,2)*pow(my,2) - mx*py*pow(d6,2)*pow(my,2) - a6*lx*pow(d6,2)*pow(my,3) + px*pow(d6,2)*pow(my,3) + a6*ly*mx*pow(d6,2)*pow(mz,2) - a6*lx*my*pow(d6,2)*pow(mz,2) + my*px*pow(d6,2)*pow(mz,2) - mx*py*pow(d6,2)*pow(mz,2) + a6*ly*mx*pow(d6,2)*pow(nx,2) - 3*a6*lx*my*pow(d6,2)*pow(nx,2) + 3*my*px*pow(d6,2)*pow(nx,2) - mx*py*pow(d6,2)*pow(nx,2) + 3*a6*ly*mx*pow(d6,2)*pow(ny,2) - a6*lx*my*pow(d6,2)*pow(ny,2) + my*px*pow(d6,2)*pow(ny,2) - 3*mx*py*pow(d6,2)*pow(ny,2) + 3*a6*ly*mx*pow(d6,2)*pow(nz,2) - 3*a6*lx*my*pow(d6,2)*pow(nz,2) + 3*my*px*pow(d6,2)*pow(nz,2) - 3*mx*py*pow(d6,2)*pow(nz,2) + 4*a6*ly*mx*pow(px,2) - 12*a6*lx*my*pow(px,2) - 4*mx*py*pow(px,2) + 4*my*pow(px,3) + 12*a6*ly*mx*pow(py,2) - 4*a6*lx*my*pow(py,2) + 4*my*px*pow(py,2) + 8*a1*lambda1*pow(mu1,-1)*(d1*mx*px + d1*my*py - mx*px*pz - my*py*pz + a6*(-(d1*(lx*mx + ly*my)) + lz*mx*px - 2*lx*mz*px + lz*my*py - 2*ly*mz*py + lx*mx*pz + ly*my*pz) + pow(a6,2)*(-(lx*lz*mx) + ly*(-(lz*my) + ly*mz) + mz*pow(lx,2)) + mz*pow(px,2) + mz*pow(py,2)) - 4*mx*pow(py,3) + 4*a6*ly*mx*pow(pz,2) - 4*a6*lx*my*pow(pz,2) + 4*my*px*pow(pz,2) - 4*mx*py*pow(pz,2) - 5*a6*lx*mx*my*nx*pow(d6,2)*sin(2*alpha6) + 5*a6*ly*mx*my*ny*pow(d6,2)*sin(2*alpha6) + 5*a6*ly*mx*mz*nz*pow(d6,2)*sin(2*alpha6) - 5*a6*lx*my*mz*nz*pow(d6,2)*sin(2*alpha6) + 5*mx*my*nx*px*pow(d6,2)*sin(2*alpha6) + 5*my*mz*nz*px*pow(d6,2)*sin(2*alpha6) - 5*mx*my*ny*py*pow(d6,2)*sin(2*alpha6) - 5*mx*mz*nz*py*pow(d6,2)*sin(2*alpha6)))));

return(dummyexp);

}