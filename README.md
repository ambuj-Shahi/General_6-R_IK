# General_6-R_IK

Description
 
 This file documents the use of the function 'InverseKinematicsGeneral6R', to compute the inverse kinematics solutions of general 6-R serial manipulator, 
 along with its dependencies. The function can be located in the 'IK' folder and its usage has been demonstrated through three examples, situated in the 
 folders 'example1', 'example2', and 'example3_Jaco'. In the folder 'example1', the numerical example taken by 'Raghavan and Roth' in their 1990 paper on 
 this topic has been implemented, in the folder 'example2', the numerical example taken by 'Liu and Zhu' in their 2007 paper has been implemented, 
 and in the folder 'example3_Jaco',the inverse kinematics algorithm has been implemented on the general 6-R robot, Jaco. 

Input files and other dependencies:

The inputs to the function 'InverseKinematicsGeneral6R' are the D-H parameters of the robot and its End-effector pose, which have been denoted by 'dhmat' and 'eemat'
in the function. The inputs 'dhmat' and 'eemat' are being taken from separate 'CSV' files. For the input 'dhmat', the corresponding 'CSV' file should have 3 rows and 6 columns, with the '1st' row representing the D-H parameter 'a' of the robot, with the six column entries of that row representing the values of 'a1' to 'a6' for the robot. The values of 'a1' to 'a6' should be in meters. The 2nd row represents the D-H parameter 'alpha' of the robot with the six column entries of that row representing the values of 'alpha1' to 'alpha6' for the robot. The values of 'alpha1' to 'alpha6' should be in radians. The 3rd row represents the D-H parameter 'd' of the robot with the six column entries of that row representing the values of 'd1' to 'd6' for the robot. The values of 'd1' to 'd6' should be in meters. In this algorithm, the robots have been modelled using the classic D-H parameters. So, if the robot has been modelled using modified D-H parameters, then that must be converted into classic D-H parameters, which can then be used as input to the 'CSV' file in the way described above to create the 'dhmat'. In the folder 'example1', one such CSV file can be located named, 'dhmat_Raghavan_Roth.csv'. Similar files are there in the other example folders as well.

For the other input 'eemat' to the function 'InverseKinematicsGeneral6R', the corresponding 'CSV' file should have 4 rows and 3 columns, the 1st row of which represents the 1st column of the rotation matrix, representing the orientation of the end-effector  of the robot w.r.t the base frame. The 2nd row of this 'CSV' file
represents the 2nd column of the aforementioned rotation matrix, the 3rd row of the 'CSV' file represents the 3rd column of the aforementioned rotation matrix 
and the 4th row of the 'CSV' file represents the position vector (Px, Py, Pz) of the end-effector w.r.t. base frame, and values of the entries of the 
position vector should also be in meters. In the folder 'example1', one such 'CSV' file can be located named, 'eemat_Raghavan_Roth.csv'. Similar files are there in the other example folders as well.

Following the 1990 paper by Raghavan and Roth and 1992 paper by Manocha and Canny on Inverse kinematics of general 6-R robots, one can arrive at the 12 x 12 matrix, whose entries are quadratic polynomial in x3, lets call it 'capsigdd'. When the problem of root finding of the determinant of the matrix, capsigdd, is reduced to an eigenvalue problem, the matrix capsigdd can further be decomposed into (A * x3^2 + B * x3 + C), where A, B and C are 12 x 12 matrices and x3 = tan(theta3 / 2). In the function 'InverseKinematicsGeneral6R', 'A', 'B', and 'C' is represented by cmatx2, cmatx1, and cmatx0 repectively, each of the entries of which are functions of the D-H parameters and the End-effector pose.  
The 12 x 12 coefficient matrices 'cmatx0', 'cmatx1', and 'cmatx2' in the function 'InverseKinematicsGeneral6R' has been formed by calling functions from the folders 
Coeff_x0, Coeff_x1, and Coeff_x2 respectively. Each of the files in these folders has been formed out of the expressions generated in CAS Mathematica, after simplifying the entries of the corresponding coefficient matrices.

 Output files:

 The output file is also a 'CSV' file, in which the all the solutions to the inverse kinematics problem for the general 6-R serial manipulator in consideration,
 are stored. In the folder 'example1', one such 'CSV' file can be located, named, 'example1_output.csv'. Similar files are there in the other example folders as well.

 
