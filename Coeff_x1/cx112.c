
# include <stdio.h>
# include <math.h>

double cx112(double* coeff)
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



double dummyexp = 8*pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1)*(a2*a5*a6*ly*mu3*mu4*mu6*mx - a2*a5*a6*lx*mu3*mu4*mu6*my + a2*a5*mu3*mu4*mu6*my*px - a2*a5*mu3*mu4*mu6*mx*py + a2*a5*lambda3*lambda4*(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)) + 2*a1*a6*d6*lx*mu2*mu5*mu6*mx*pow(mu1,-1) + 2*a1*a6*d6*ly*mu2*mu5*mu6*my*pow(mu1,-1) - 2*a1*a6*lx*mu2*mu5*px*pow(mu1,-1) - 2*a1*d6*mu2*mu5*mu6*mx*px*pow(mu1,-1) - 2*a1*a6*ly*mu2*mu5*py*pow(mu1,-1) - 2*a1*d6*mu2*mu5*mu6*my*py*pow(mu1,-1) + a1*mu2*mu5*pow(a6,2)*pow(lx,2)*pow(mu1,-1) + a1*mu2*mu5*pow(a6,2)*pow(ly,2)*pow(mu1,-1) + lambda6*(a2*a5*mu3*mu4*(a6*ly*nx - a6*lx*ny + ny*px - nx*py) + 2*a1*d6*mu2*mu5*(a6*lx*nx + a6*ly*ny - nx*px - ny*py)*pow(mu1,-1)) + a1*mu2*mu5*pow(d6,2)*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) + a1*mu2*mu5*pow(d6,2)*pow(mu1,-1)*pow(mu6,2)*pow(my,2) + a1*mu2*mu5*pow(d6,2)*pow(lambda6,2)*pow(mu1,-1)*(pow(nx,2) + pow(ny,2)) + a1*mu2*mu5*pow(mu1,-1)*pow(px,2) + a1*mu2*mu5*pow(mu1,-1)*pow(py,2) + a1*mu2*mu5*mx*nx*pow(d6,2)*pow(mu1,-1)*sin(2*alpha6) + a1*mu2*mu5*my*ny*pow(d6,2)*pow(mu1,-1)*sin(2*alpha6));

return(dummyexp);

}