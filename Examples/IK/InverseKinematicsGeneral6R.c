/** 
 
 @date September 30, 2017

 @file InverseKinematicsGeneral6R.c
 @brief The program is used to find the Inverse kinematics solutions for a general 6-R serial manipulator.

 
 */

# include "cx.h"
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

/**

	user_atan2_c: gives the arc tangent of z2/z1, taking into account which quadrant the point (z1,z2) is in, where 'z1' and 'z2' are both complex nos.
	

	Variables used:
	Variable | Description
	-------- |------------
	z3       | Complex no. to store the output of the function, user_atan2_c.
	
	@returns z3.
 */

double _Complex user_atan2_c( double _Complex z1,  double _Complex z2)
{
 		
		double _Complex z3 = -I * clog((z1 + I * z2)/csqrt(cpow(z1,2) + cpow(z2,2)));    
		return (z3);
}


/**

	InverseKinematicsGeneral6R: Function to compute the Inverse kinematics solutions for a general 6-R serial manipulator.

	@param dhmat: A 2-D array containing the D-H parameters of the 6-R robot in cosideration.
	@param eemat: A 2-D array containing the End-effector pose of the 6-R robot in consideration. 
	              The End-effector pose being described by the 3 x 3 rotation matix and 3x1 position vector.

				  Further information about dhmat and eemat is given in the 'Input files' section of the document 'IK_6R_Program_Documentation'.


	Variables used:
	Variable                 | Description
	--------------------     |------------
	rtod                     | Constant to convert angle in radians to angle in degrees.
	lambda1 to lambda6       | cos(alpha1) to cos(alpha6),  where alphai is from the D-H table. 
	mu1 to mu6               | sin(alpha1) to sin(alpha6),  where alphai is from the D-H table.
	vars                     | A 2-D array to store the D-H parameters and the End-effector pose of the 6-R robot in cosideration, to be used in other variables.
	cmatx2                   | The 12 x 12 coefficient matrix, 'A'.
	cmatx1                   | The 12 x 12 coefficient matrix, 'B'.
	cmatx0                   | The 12 x 12 coefficient matrix, 'C'. The descriptions of A, B, and C is there in the following function.
	cmatx01 to cmatx06       | The rows of the 6 x 9 submatrix used to form the 12 x 12 matrix, cmatx0.
	cmatx11 to cmatx16       | The rows of the 6 x 9 submatrix used to form the 12 x 12 matrix, cmatx1.
	cmatx21 to cmatx26       | The rows of the 6 x 9 submatrix used to form the 12 x 12 matrix, cmatx2.
	i, j, l, m               | counters used in the function.
	angle                    | User-defined Structure type variable, Solution. Its definition is there in the header file, 'cx.h'



	@returns angle

 */


struct Solution InverseKinematicsGeneral6R(double dhmat[10][10], double eemat[10][10])
{

	
    const double rtod = 180.0/M_PI;          

	double lambda1 = cos(dhmat[1][0]);
	double lambda2 = cos(dhmat[1][1]);
	double lambda3 = cos(dhmat[1][2]);
	double lambda4 = cos(dhmat[1][3]);
	double lambda5 = cos(dhmat[1][4]);
	double lambda6 = cos(dhmat[1][5]);

	double mu1 = sin(dhmat[1][0]);
	double mu2 = sin(dhmat[1][1]);
	double mu3 = sin(dhmat[1][2]);
	double mu4 = sin(dhmat[1][3]);
	double mu5 = sin(dhmat[1][4]);
	double mu6 = sin(dhmat[1][5]);



	double vars[] = { dhmat[0][0], dhmat[0][1], dhmat[0][2], dhmat[0][3], dhmat[0][4], dhmat[0][5],
	 	              dhmat[1][0], dhmat[1][1], dhmat[1][2], dhmat[1][3], dhmat[1][4], dhmat[1][5],
	                  dhmat[2][0], dhmat[2][1], dhmat[2][2], dhmat[2][3], dhmat[2][4], dhmat[2][5],
	                  lambda1, lambda2, lambda3, lambda4, lambda5,lambda6,
	                  eemat[0][0], eemat[0][1], eemat[0][2],
	                  mu1, mu2, mu3, mu4, mu5, mu6,
	                  eemat[1][0], eemat[1][1], eemat[1][2],
	                  eemat[2][0], eemat[2][1], eemat[2][2],
	                  eemat[3][0], eemat[3][1], eemat[3][2] }; 

	
   	/*
   Following the 1990 paper by Raghavan and Roth and 1992 paper by Manocha and Canny on Inverse kinematics of general 6-R robots, 
   one can arrive at the 12 x 12 matrix, whose entries are quadratic polynomial in x3, lets call it 'capsigdd'. When the problem 
   of root finding of the determinant of the matrix, capsigdd, is reduced to an eigenvalue problem, the matrix capsigdd can further 
   be decomposed into (A * x3^2 + B * x3 + C), where A, B and C are 12 x 12 matrices and x3 = tan(theta3 / 2). Here, 'A', 'B', and 'C' 
   is represented by cmatx2, cmatx1, and cmatx0 repectively, each of the entries of which are functions of the D-H parameters and the End-effector pose.  

   	*/

  	/*
   	The 12 x 12 matrix cmatx0 has been formed out of a 6 x 9 submatrix, as can be seen from the paper by Manocha and Canny. Each of the entries of this
   	6 x 9 submatrix has been simplified in the CAS Mathematica and converted into C code, and the expression thus obtained has been used to create a 
   	C function, which is being called here to form the 6 x 9 submatrix. The rows of this submatrix are cmatx01 to cmatx06. Each of the entries of these
   	rows are the aforementioned functions, which have been given a generic nomenclature of cx0ij, where i and j denotes the row position and column position 
   	of the entry in the submatrix, respectively.

   	*/


    // maxexp = highest power of the polynomial entries in the 12 x 12 matrix, capsigdd, which in this case is 2.
	const int maxexp=2; 
    
    // cord = order of the matrix, capsigdd, which is 12.
	const int cord=12;                
    
    int i,j;

	double cmatx01[12]={cx011(vars), cx012(vars), cx013(vars), cx014(vars), cx015(vars), \
	cx016(vars), cx017(vars), cx018(vars), cx019(vars)};

	double cmatx02[12]={cx021(vars), cx022(vars), cx023(vars), cx024(vars), cx025(vars), \
	cx026(vars), cx027(vars), cx028(vars), cx029(vars)};

	double cmatx03[12]={cx031(vars), cx032(vars), cx033(vars), cx034(vars), cx035(vars), \
	cx036(vars), cx037(vars), cx038(vars), cx039(vars)};

	double cmatx04[12]={cx041(vars), cx042(vars), cx043(vars), cx044(vars), cx045(vars), \
	cx046(vars), cx047(vars), cx048(vars), cx049(vars)};

	double cmatx05[12]={cx051(vars), cx052(vars), cx053(vars), cx054(vars), cx055(vars), \
	cx056(vars), cx057(vars), cx058(vars), cx059(vars)};

	double cmatx06[12]={cx061(vars), cx062(vars), cx063(vars), cx064(vars), cx065(vars), \
	cx066(vars), cx067(vars), cx068(vars), cx069(vars)};

    double cmatx0[cord][cord];

	for(i=0;i<cord;i++)
	{
		for(j=0;j<cord;j++)
		{
				cmatx0[i][j]=0;
			
		}		
	}

	for(j=0;j<9;j++)
	{
		cmatx0[0][j]=cmatx01[j];
	}

    for(j=0;j<9;j++)
	{
		cmatx0[1][j]=cmatx02[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx0[2][j]=cmatx03[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx0[3][j]=cmatx04[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx0[4][j]=cmatx05[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx0[5][j]=cmatx06[j];
	}
    
    for(i=6;i<cord;i++)
	{
		for(j=3;j<cord;j++)
		{
				cmatx0[i][j]=cmatx0[i-6][j-3];
			
		}		
	}

	

	/**

    The 12 x 12 matrix cmatx1 has also been formed in the same way as the matrix cmatx0. The rows of the 6 x 9 submatrix used in the formation of cmatx1
    are cmatx11 to cmatx16. Each of the entries of these rows have been given a generic nomenclature of cx1ij, where i and j denotes the row position and 
    column position of the entry in the submatrix, respectively. 

   	*/


	double cmatx11[12]={cx111(vars), cx112(vars), cx113(vars), cx114(vars), cx115(vars), \
	cx116(vars), cx117(vars), cx118(vars), cx119(vars)};

	double cmatx12[12]={cx121(vars), cx122(vars), cx123(vars), cx124(vars), cx125(vars), \
	cx126(vars), cx127(vars), cx128(vars), cx129(vars)};

	double cmatx13[12]={cx131(vars), cx132(vars), cx133(vars), cx134(vars), cx135(vars), \
	cx136(vars), cx137(vars), cx138(vars), cx139(vars)};

	double cmatx14[12]={cx141(vars), cx142(vars), cx143(vars), cx144(vars), cx145(vars), \
	cx146(vars), cx147(vars), cx148(vars), cx149(vars)};

	double cmatx15[12]={cx151(vars), cx152(vars), cx153(vars), cx154(vars), cx155(vars), \
	cx156(vars), cx157(vars), cx158(vars), cx159(vars)};

	double cmatx16[12]={cx161(vars), cx162(vars), cx163(vars), cx164(vars), cx165(vars), \
	cx166(vars), cx167(vars), cx168(vars), cx169(vars)};



	double cmatx1[cord][cord];

	for(i=0;i<cord;i++)
	{
		for(j=0;j<cord;j++)
		{
				cmatx1[i][j]=0;
			
		}		
	}

	for(j=0;j<9;j++)
	{
		cmatx1[0][j]=cmatx11[j];
	}

    for(j=0;j<9;j++)
	{
		cmatx1[1][j]=cmatx12[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx1[2][j]=cmatx13[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx1[3][j]=cmatx14[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx1[4][j]=cmatx15[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx1[5][j]=cmatx16[j];
	}
    
    for(i=6;i<cord;i++)
	{
		for(j=3;j<cord;j++)
		{
				cmatx1[i][j]=cmatx1[i-6][j-3];
			
		}		
	}



	/*

    The 12 x 12 matrix cmatx2 has also been formed in the same way as the matrix cmatx0. The rows of the 6 x 9 submatrix used in the formation of cmatx2
    are cmatx21 to cmatx26. Each of the entries of these rows have been given a generic nomenclature of cx2ij, where i and j denotes the row position and 
    column position of the entry in the submatrix, respectively. 

   	*/


	double cmatx21[12]={cx211(vars), cx212(vars), cx213(vars), cx214(vars), cx215(vars), \
	cx216(vars), cx217(vars), cx218(vars), cx219(vars)};

	double cmatx22[12]={cx221(vars), cx222(vars), cx223(vars), cx224(vars), cx225(vars), \
	cx226(vars), cx227(vars), cx228(vars), cx229(vars)};

	double cmatx23[12]={cx231(vars), cx232(vars), cx233(vars), cx234(vars), cx235(vars), \
	cx236(vars), cx237(vars), cx238(vars), cx239(vars)};

	double cmatx24[12]={cx241(vars), cx242(vars), cx243(vars), cx244(vars), cx245(vars), \
	cx246(vars), cx247(vars), cx248(vars), cx249(vars)};

	double cmatx25[12]={cx251(vars), cx252(vars), cx253(vars), cx254(vars), cx255(vars), \
	cx256(vars), cx257(vars), cx258(vars), cx259(vars)};

	double cmatx26[12]={cx261(vars), cx262(vars), cx263(vars), cx264(vars), cx265(vars), \
	cx266(vars), cx267(vars), cx268(vars), cx269(vars)};


	double cmatx2[cord][cord];

	for(i=0;i<cord;i++)
	{
		for(j=0;j<cord;j++)
		{
				cmatx2[i][j]=0;
			
		}		
	}

	for(j=0;j<9;j++)
	{
		cmatx2[0][j]=cmatx21[j];
	}

    for(j=0;j<9;j++)
	{
		cmatx2[1][j]=cmatx22[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx2[2][j]=cmatx23[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx2[3][j]=cmatx24[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx2[4][j]=cmatx25[j];
	}

	for(j=0;j<9;j++)
	{
		cmatx2[5][j]=cmatx26[j];
	}
    
    for(i=6;i<cord;i++)
	{
		for(j=3;j<cord;j++)
		{
				cmatx2[i][j]=cmatx2[i-6][j-3];
			
		}		
	}


    

    


	// Forming the identity matrix, idmat. 
	double idmat[cord][cord];

	for(i=0;i<cord;i++)
	{
		for(j=0;j<cord;j++)
		{
			if(i==j)
				idmat[i][j]=1;
			else
				idmat[i][j]=0;
		}		
	}

	/* dim of the A and B matrices */
	// mord = order of the matrices M1 and M2, to be used in the generalized eigenvalue problem. For information about the structure of M1 
	// and M2 matrices, look into the aforementioned paper by Manocha and Canny.

	int mord = maxexp*cord;

	double M1[mord][mord];
	double M2[mord][mord];

	
	int l,m;

    // Initialisation of matrices, M1 and M2 as null matrices of order 24 x 24.

	for(i=0;i<mord;i++)
    {
        for(j=0;j<mord;j++)
        {
            M1[i][j]=0;
            M2[i][j]=0;
        }       
    } 

	// Forming the 1st row of the M2 matrix, with 12 x 12 null matrix as its 1st entry and 12 x 12 identity matrix as its 2nd entry.
 

	for(l=0;l<cord;l=l+cord)
	{
			int temp1= l/cord;
			for(i=0;i<cord;i++)
			{
				for(m=0;m<mord;m=m+cord)
				{
					int temp2= m/cord;
					for(j=0;j<cord;j++)
					{
						if(temp2==(temp1+1))
							M2[l+i][m+j]=idmat[i][j];
						else
							M2[l+i][m+j]=0;					
					}
				}
			}		
	}
    
    // Forming the 2nd row of the M2 matrix, with the 12 x 12 matrix, cmatx0, as its 1st entry and 12 x 12 matrix, cmatx1 as its 2nd entry.

	l=cord;
	for(i=0;i<cord;i++)
	{
		for(m=0;m<mord;m=m+cord)
		{
			int temp2= m/cord;
			for(j=0;j<cord;j++)
			{
				if(temp2==0)
					M2[l+i][m+j]=cmatx0[i][j];
				if(temp2==1)
					M2[l+i][m+j]=cmatx1[i][j];			
			}
		}
		
	}	

	// Forming the 1st row of the M1 matrix, with 12 x 12 identity matrix as its 1st entry and 12 x 12 null matrix as its 2nd entry.

	for(l=0;l<cord;l=l+cord)
	{
		int temp1= l/cord;
		for(i=0;i<cord;i++)
		{
			for(m=0;m<mord;m=m+cord)
			{
				int temp2= m/cord;
				for(j=0;j<cord;j++)
				{
					if(temp2==temp1)
						M1[l+i][m+j]=idmat[i][j];
					else
						M1[l+i][m+j]=0;					
				}
			}
		}		
	}

    // Forming the 2nd row of the M1 matrix, with 12 x 12 null matrix as its 1st entry and 12 x 12 matrix, cmatx2 as its 2nd entry.

	l=cord;
	for(i=0;i<cord;i++)
	{
		for(m=0;m<mord;m=m+cord)
		{
			int temp2= m/cord;
			for(j=0;j<cord;j++)
			{
				if(temp2==1)
				M1[l+i][m+j]=-cmatx2[i][j];	
				else
				M1[l+i][m+j]=0;		
			}
		}
	
	}	


	int countA=0;  // Counters to form Vectors A and B.
	int countB=0;
	
    
    /*
    Flattening the 2-D arrays M1 and M2 into 1-D arrays, B and A, respectively to be used by 'gsl_matrix_view_array' function to convert 
    the matrices, M1 and M2, to the form required by the other GSL functions. For more information about these GSL functions and the GSL eigensystem 
    routines, please look into the GSL documentation.

    */

    double A[mord*mord];
	double B[mord*mord];

	for(i=0;i<mord;i++)
	{
		for(j=0;j<mord;j++)
		{
		A[countA++] = M2[i][j];		
		B[countB++] = M1[i][j];
		}
	}
	
 	// Converting the matrices to the  required form 

 	gsl_matrix_view aa = gsl_matrix_view_array (A, mord, mord);

 	gsl_matrix_view bb = gsl_matrix_view_array (B, mord, mord);
	
 	/*  
 	The function 'gsl_eigen_gen' computes the eigenvalues of the real generalized nonsymmetric matrix pair (aa, bb), and
 	stores them as pairs in (alpha, beta), where 'alpha' is complex and 'beta' is real. If 'beta' is non-zero, then (alpha/beta) is an eigenvalue.
 	The eigenvectors are stored in 'evec'. 
 	*/
     

 	// Allocation of memory for alpha and beta 

  	gsl_vector_complex *alpha = gsl_vector_complex_alloc (mord);
  	gsl_vector *beta = gsl_vector_alloc (mord);
    
    // Allocation of memory for evec

  	gsl_matrix_complex *evec = gsl_matrix_complex_alloc (mord, mord);

 	// Allocation of workspace for computing the eigenvalues and eigenvectors of the matrix pair (aa, bb).

  	gsl_eigen_genv_workspace *w =  gsl_eigen_genv_alloc (mord);

 	/* 
 	The function 'gsl_eigen_genv' computes eigenvalues and right eigenvectors of the matrix pair (aa, bb). The eigenvalues are stored in (alpha,beta)
	and eigenvectors are stored in 'evec'.
 	*/
    gsl_eigen_genv(&aa.matrix, &bb.matrix,  alpha, beta, evec, w);

    // This function frees the memory associated with the workspace 'w'.
    gsl_eigen_genv_free (w);


   	int dof = 6;   // dof = No. of joints of the manipulator

	struct Solution angle;

 	/* Getting all the eigenvalues and the corresponding eigenvectors */

    for (i = 0; i < mord; i++)
    {
	
	 	/*
	 	Getting the eigenvectors for each iteration of i from 1 to 24, and storing it in variable 'evec_i'. The structure of the eigenvectors can
	 	be obtained from the paper by Manocha and Canny and it can be seen that the 12th element of that vector is 1 and when the 12th element is 1, 
	 	the 9th and 11th element is 'x4' and 'x5' respectively. 'theta4' and 'theta5' can now be obtained with, theta4 = 2*atan(x4) and theta5 = 2*atan(x5).
        The eigenvalues, x3, of the matrix pair (aa, bb) is (alpha/beta). Thus, theta3 can be obtained with, theta3 = 2*atan(x3).

	 	*/
        
        // The function 'gsl_matrix_complex_column' with argument 'evec' and 'i' gets the ith column vector of 'evec'.

	 	gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);

	    // The eigenvalues exist only if beta is non-zero.

        if(gsl_vector_get (beta, i) != 0)
		{
            
            // Obtaining the eigenvalue, x3

            double x3r = GSL_REAL(gsl_vector_complex_get (alpha, i))/gsl_vector_get (beta, i);   // Real part of x3

			double x3i = GSL_IMAG(gsl_vector_complex_get (alpha, i))/gsl_vector_get (beta, i);   // Imaginary part of x3
            
			double _Complex x3c = x3r + I * x3i;      // x3c = Complex no. representing, x3.


			if( fabsf(x3r) != 0.0f && fabsf(x3i) != 1.0f)    // Discarding the solutions of the equation (1+x^2)^4 = 0. 
			                                                 // i.e. getting the eigenvalues other than I and -I, where I = sqrt(-1). 
			{

				double _Complex th3c = 2 * catan(x3c);     // th3c = Complex no. representing, theta3 in radians.
				double _Complex th3cd = th3c * rtod;       // th3cd = Complex no. representing, theta3 in degrees.
				

	            /*
	            Getting the 12th element of the vector 'evec_i', so that when the 9th and 11th elements are normalized w.r.t. it, we can get the 
	            desired values for x4 and x5.
	            */
				gsl_complex z0 = gsl_vector_complex_get(&evec_i.vector, 11);
				gsl_complex x5val = gsl_vector_complex_get(&evec_i.vector, 10);  // Value of x5 as obtained from the eigenvector
				gsl_complex x4val = gsl_vector_complex_get(&evec_i.vector, 8);   // Value of x4 as obtained from the eigenvector

				gsl_complex x4 = gsl_complex_div(x4val,z0);   // Normalized value of x4 w.r.t. the 12th element of the eigenvector.
				double x4r = GSL_REAL(x4);
				double x4i = GSL_IMAG(x4);
				double _Complex x4c = x4r + I * x4i;          // x4c = Complex no. representing, x4.
				double _Complex th4c = 2 * catan(x4c);        // th4c = Complex no. representing, theta4 in radians.
				double _Complex th4cd = th4c * rtod;		  // th4cd = Complex no. representing, theta4 in degrees.


				gsl_complex x5 = gsl_complex_div(x5val,z0);   // Normalized value of x5 w.r.t. the 12th element of the eigenvector.
				double x5r = GSL_REAL(x5);
				double x5i = GSL_IMAG(x5);
				double _Complex x5c = x5r + I * x5i;          // x5c = Complex no. representing, x5.
				double _Complex th5c = 2 * catan(x5c);        // th5c = Complex no. representing, theta5 in radians.
				double _Complex th5cd = th5c * rtod;		  // th5cd = Complex no. representing, theta5 in degrees.

				/* Calculation of theta1 using the function 'theta1' */

				double _Complex vars2[] = {th3c, th4c, th5c};     
				struct Theta theta1r = theta1(vars, vars2);              
				double _Complex th1c = user_atan2_c(theta1r.costheta, theta1r.sintheta);  // Calculation of theta1 in radians.
				double _Complex th1cd = th1c * rtod;                                      // Calculation of theta1 in degrees.

				/* Calculation of theta2 using the function 'theta2' */

				double _Complex vars3[] = {th1c, th3c, th4c, th5c};     
				struct Theta theta2r = theta2(vars, vars3);              
				double _Complex th2c = user_atan2_c(theta2r.costheta, theta2r.sintheta);  // Calculation of theta2 in radians.
				double _Complex th2cd = th2c * rtod; 									  // Calculation of theta2 in degrees.


				/* Calculation of theta6 using the function 'theta6' */

				double _Complex vars4[] = {th1c, th2c, th3c, th4c, th5c};     
				struct Theta theta6r = theta6(vars, vars4);              
				double _Complex th6c = user_atan2_c(theta6r.costheta, theta6r.sintheta);  // Calculation of theta6 in radians.
				double _Complex th6cd = th6c * rtod;                                      // Calculation of theta6 in degrees.

				     
				struct Residue res16 = residue(vars, vars4);                              // Calculation of residues.
				
				
				double _Complex vars5[] = {th1c, th2c, th3c, th4c, th5c, th6c};           
				double _Complex vars6[] = {th1cd, th2cd, th3cd, th4cd, th5cd, th6cd};
				double _Complex resvec[] = {res16.residue1, res16.residue2, res16.residue3, res16.residue4, res16.residue5, res16.residue6};

				int j;

				for(j=0;j<dof;j++)
				{
					angle.thetar[i][j] = vars5[j];
					angle.thetad[i][j] = vars6[j];
					angle.residue[i][j] = resvec[j];
				}


			}

				
		}

        			         			
	
	}
 
 
 
 	gsl_vector_complex_free(alpha);
 	gsl_vector_free(beta);

 	return(angle);

}

