/** 
 
 @date September 30, 2017

 @file example1.c 
 @brief This file demonstrates the use of the Inverse kinematics function InverseKinematicsGeneral6R, through an example on a general 6-R serial manipulator, taken from
        the 1990 paper by 'Raghavan and Roth'.

 
*/

# include "cx.h"
# include "csvparser.h"
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>              // Standard Library of Complex Numbers 
# include <gsl/gsl_math.h>
# include <gsl/gsl_eigen.h>
# include <gsl/gsl_vector.h>
# include <gsl/gsl_complex.h>
# include <gsl/gsl_complex_math.h>
# include <gsl/gsl_statistics.h>
# include <time.h>


/**

create_csv: Function to create a CSV file. 

@param filename:    Name of the CSV file to be generated.
@param solution:    A 2-D array containing the information to be stored in the CSV file.
@param m:           No. of rows of the above defined 2-D array 'solution'.
@param n:           No. of columns of the above defined 2-D array 'solution'.

              
Variables used:
Variable                 | Description
--------------------     |------------
fp                       | A pointer of type FILE to hold the CSV file named 'filename', which is an user input.


@returns: The CSV file thus created.

*/

void create_csv(char *filename, double _Complex solution[][10], int n, int m)
{
 
    printf("\n Creating %s.csv file",filename);
 
    FILE *fp;
 
    int i, j;
 
    filename=strcat(filename,".csv");
    
    fp=fopen(filename,"w+");

    fprintf(fp,"Solutions, theta1, theta2, theta3, theta4, theta5, theta6");
 
    for(i=0;i<m;i++)
    {
 
        fprintf(fp,"\n%d",i+1);
 
        for(j=0;j<n;j++)
        {
            fprintf(fp,", %f%+fi ", creal(solution[i][j]), cimag(solution[i][j]));
        }
 
    }
 
    fclose(fp);
 
    printf("\n %sfile created",filename);
 
}


/**

main program to demonstrate the usage of the function 'InverseKinematicsGeneral6R'.

For computing the inverse kinematics of any other general 6-R robot, the CSV file containing the D-H parameters needs to be placed at the position of 
'dhmat_Raghavan_Roth.csv' in the program below, and for any other End-effector pose of the robot, the CSV file containing the End-effector pose needs 
to be placed at the position of 'eemat_Raghavan_Roth.csv'. 
The way to input D-H parameters and the End-Effector pose into the corresponding CSV files has been explained in the document 
'IK_6R_Program_Documentation.pdf'.

Variables used:
Variable                   | Description
--------------------       |------------
soln                       | Structure type variable, to store the output of the function 'InverseKinematicsGeneral6R'. Its definition is there in the header file, 'cx.h'.
start, end, cpu_time_used  | To calculate the average computation time of the inverse kinematics solution through the use of the function 'InverseKinematicsGeneral6R'
NN                         | No. of cycles over which the function 'InverseKinematicsGeneral6R'  has been evaluated to get an estimate of its average computation time.





*/



int main() 
{

    struct Solution soln;
    clock_t start, end;
    double cpu_time_used;

    int NN = 200;

    int i = 0, j;

    /*

    Formation of 'dhmat' through the use of a CSV parser with its 1st argument being, the CSV file to be read ('dhmat_Raghavan_Roth.csv' in this case),
    its 2nd argument being the delimiter (comma or semicolon) and its 3rd argument being 0 or 1 (0 -> CSV files without header and 1 -> CSV files with header).
    Here, we are using CSV files without headers with this function.

    */ 

    CsvParser *csvparser1 = CsvParser_new("dhmat_Raghavan_Roth.csv", ",", 0);
    CsvRow *row1;
    double dhmat[10][10];


    while ((row1 = CsvParser_getRow(csvparser1))) 
    {

        const char **rowFields1 = CsvParser_getFields(row1);
        for (j = 0 ; j < CsvParser_getNumFields(row1) ; j++) 
        {
            dhmat[i][j] = atof(rowFields1[j]) ;
        }
    
        CsvParser_destroy_row(row1);
        i++;
    }
    CsvParser_destroy(csvparser1);


    
    /*

    Formation of 'eemat' through the use of the CSV parser in a similar way as the 'dhmat' was formed. 

    */ 

    int ii = 0, jj;
    CsvParser *csvparser2 = CsvParser_new("eemat_Raghavan_Roth.csv", ",", 0);
    CsvRow *row2;
    double eemat[10][10];


    while ((row2 = CsvParser_getRow(csvparser2))) 
    {

        const char **rowFields2 = CsvParser_getFields(row2);
        for (jj = 0 ; jj < CsvParser_getNumFields(row2) ; jj++) 
        {
            eemat[ii][jj] = atof(rowFields2[jj]) ;
        }
    
        CsvParser_destroy_row(row2);
        ii++;
    }
    CsvParser_destroy(csvparser2); 
  

    /*

    Computation of the inverse kinematics solutions using the function 'InverseKinematicsGeneral6R', with the previuosly formed 'dhmat' and 'eemat' as input and storing it
    into structure type variable 'soln'.

    */ 

    start=clock();

    for(i = 0; i < NN; i++)
    {
        soln = InverseKinematicsGeneral6R(dhmat, eemat);
    }
  
    end= clock();

    cpu_time_used = (((double) (end - start)) / CLOCKS_PER_SEC)/((double) (NN));

    printf("\n code took %f seconds to execute \n \n", cpu_time_used);

    int nsolns = 16;   // No. of branches of the inverse kinematics solution. 

    double absresidue[20][10];
    double maxresidue[20];

    for ( i = 0; i < nsolns; i++ ) 
    {
        for ( j = 0; j < 6; j++ ) 
        {
            absresidue[i][j] = cabs(soln.residue[i][j]);
        }
        
    }

    for ( i = 0; i < nsolns; i++ ) 
    {
        double data[6] = {absresidue[i][0], absresidue[i][1], absresidue[i][2], absresidue[i][3], absresidue[i][4], absresidue[i][5]};
        maxresidue[i] = gsl_stats_max(data, 1, 6);
        
    }
  
   // Printing the solutions thus obtained and their corresponding residues:

    for ( i = 0; i < nsolns; i++ ) 
    {
                printf("Solution-%d\n",i+1);
                printf("theta1 = %f%+fi ", creal(soln.thetad[i][0]), cimag(soln.thetad[i][0]));
                printf("theta2 = %f%+fi ", creal(soln.thetad[i][1]), cimag(soln.thetad[i][1]));
                printf("theta3 = %f%+fi ", creal(soln.thetad[i][2]), cimag(soln.thetad[i][2]));
                printf("theta4 = %f%+fi ", creal(soln.thetad[i][3]), cimag(soln.thetad[i][3]));
                printf("theta5 = %f%+fi ", creal(soln.thetad[i][4]), cimag(soln.thetad[i][4]));
                printf("theta6 = %f%+fi\n\n", creal(soln.thetad[i][5]), cimag(soln.thetad[i][5]));
                printf("Residue1 = %f%+fi, Residue2 = %f%+fi, Residue3 = %f%+fi, Residue4 = %f%+fi, Residue5 = %f%+fi, Residue6 = %f%+fi\n\n",\
                        creal(soln.residue[i][0]), cimag(soln.residue[i][0]), creal(soln.residue[i][1]), cimag(soln.residue[i][1]),\
                        creal(soln.residue[i][2]), cimag(soln.residue[i][2]), creal(soln.residue[i][3]), cimag(soln.residue[i][3]),\
                        creal(soln.residue[i][4]), cimag(soln.residue[i][4]), creal(soln.residue[i][5]), cimag(soln.residue[i][5]));

                printf("Maxresidue = %f\n\n", maxresidue[i]);

        printf("\n ********************************************************************************************************** \n");
 
    }
    
    
    char str[100];       // String to store the name of the output CSV file, to be defined by the user.
 
    printf("\n Enter the filename :");
 
    scanf("%[^\n]s",str);
 
    create_csv(str,soln.thetad,6,nsolns);

    return(0);
    
   
}




