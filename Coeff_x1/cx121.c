
# include <stdio.h>
# include <math.h>

double cx121(double* coeff)
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



double dummyexp = 2*pow(mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py) + lambda6*(a6*ly*nx - a6*lx*ny + ny*px - nx*py),-1)*(-(a2*a6*lambda3*lambda5*ly*mu4*mu6*mx) + a2*a6*ly*mu3*mu4*mu5*mu6*mx + a2*a6*lambda3*lambda5*lx*mu4*mu6*my - a2*a6*lx*mu3*mu4*mu5*mu6*my - a2*lambda3*lambda5*mu4*mu6*my*px + a2*mu3*mu4*mu5*mu6*my*px + a2*lambda3*lambda5*mu4*mu6*mx*py - a2*mu3*mu4*mu5*mu6*mx*py + (lambda6*(-4*a1*a4*mu2*mu6*(mx*nx + my*ny) + a2*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*cos(alpha1 - alpha3 + alpha4 - alpha5) - a2*(a6*ly*nx - a6*lx*ny + ny*px - nx*py)*cos(alpha1 + alpha3 - alpha4 + alpha5))*pow(mu1,-1))/2. + a1*a3*mu2*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) - a1*a4*mu2*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) + a1*a5*mu2*pow(mu1,-1)*pow(mu6,2)*pow(mx,2) + a1*a3*mu2*pow(mu1,-1)*pow(mu6,2)*pow(my,2) - a1*a4*mu2*pow(mu1,-1)*pow(mu6,2)*pow(my,2) + a1*a5*mu2*pow(mu1,-1)*pow(mu6,2)*pow(my,2) + a1*(a3 - a4 + a5)*mu2*pow(lambda6,2)*pow(mu1,-1)*(pow(nx,2) + pow(ny,2)) + a2*lambda4*mu6*(a6*ly*mx - a6*lx*my + my*px - mx*py)*sin(alpha3 + alpha5) + a1*a3*mu2*mx*nx*pow(mu1,-1)*sin(2*alpha6) + a1*a5*mu2*mx*nx*pow(mu1,-1)*sin(2*alpha6) + a1*a3*mu2*my*ny*pow(mu1,-1)*sin(2*alpha6) + a1*a5*mu2*my*ny*pow(mu1,-1)*sin(2*alpha6));

return(dummyexp);

}