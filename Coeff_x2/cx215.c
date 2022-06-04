
# include <stdio.h>
# include <math.h>

double cx215(double* coeff)
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



double dummyexp = 8*pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1)*((a2 - a3)*a5*lambda4*(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)) + (a1*mu5*pow(mu1,-1)*sin(alpha2 - alpha3)*(4*a6*d6*lx*mu6*mx + 4*a6*d6*ly*mu6*my - 4*a6*lx*px - 4*d6*mu6*mx*px - 4*a6*ly*py - 4*d6*mu6*my*py + 4*d6*lambda6*(a6*lx*nx + a6*ly*ny - nx*px - ny*py) + 2*pow(a6,2)*pow(lx,2) + 2*pow(a6,2)*pow(ly,2) + pow(d6,2)*pow(mx,2) + pow(d6,2)*pow(my,2) + pow(d6,2)*pow(nx,2) - cos(2*alpha6)*pow(d6,2)*(pow(mx,2) + pow(my,2) - pow(nx,2) - pow(ny,2)) + pow(d6,2)*pow(ny,2) + 2*pow(px,2) + 2*pow(py,2) + 2*mx*nx*pow(d6,2)*sin(2*alpha6) + 2*my*ny*pow(d6,2)*sin(2*alpha6)))/2.);

return(dummyexp);

}