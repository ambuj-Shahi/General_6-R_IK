/** 
 
 @date September 30, 2017

 @file cx.h
 @brief This file contains definition of the functions, used to compute the inverse kinematics solution branches of a general 6-R serial manipulator.

 
*/


# ifndef cx_H
# define cx_H



/**

	Solution: Structure type variable, to store the six joint angles in radians, the six joint angles in degrees, and the six residues 
	          for all the branches of the inverse kinematics solution.

	Variables used:
	Variable | Description
	-------- |------------
	thetar   | A 2-D array to store the six joint angles (in radians) for all the branches of the inverse kinematics solution.
	thetad   | A 2-D array to store the six joint angles (in degrees) for all the branches of the inverse kinematics solution.
	residue  | A 2-D array to store the six residues for for all the branches of the inverse kinematics solution.
	

*/

struct Solution
{
	double _Complex thetar[20][10];
	double _Complex thetad[20][10];
	double _Complex residue[20][10];
    
};

/**

	Theta:    Structure type variable, to store the sin(theta) and cos(theta) outputs from the functions, 'theta1', 'theta2' and 'theta6'.

	Variables used:
	Variable | Description
	-------- |------------
	sintheta | Variable to store sin(theta).
	costheta | Variable to store cos(theta).
	

 */

struct Theta
{
	double _Complex sintheta;
	double _Complex costheta;
    
};

/**

	Residue:    Structure type variable, to store residue1 to residue6, calculated from the six corresponding equations mentioned in the paper by Manocha and Canny.

	Variables used:
	Variable | Description
	-------- |------------
	residue1 | Variable to store residue1.
	residue2 | Variable to store residue2.
	residue3 | Variable to store residue3.
	residue4 | Variable to store residue4.
	residue5 | Variable to store residue5.
	residue6 | Variable to store residue6.
	
 */

struct Residue
{
	double _Complex residue1;
	double _Complex residue2;
	double _Complex residue3;
	double _Complex residue4;
	double _Complex residue5;
	double _Complex residue6;
    
};

struct Solution InverseKinematicsGeneral6R(double[10][10], double[10][10]);

struct Theta theta1(double[],double _Complex[]);
struct Theta theta2(double[],double _Complex[]);
struct Theta theta6(double[],double _Complex[]);

struct Residue residue(double[],double _Complex[]);


double cx011(double[]); double cx012(double[]); double \
cx013(double[]); double cx014(double[]); double cx015(double[]); \
double cx016(double[]); double cx017(double[]); double \
cx018(double[]); double cx019(double[]); double cx021(double[]); \
double cx022(double[]); double cx023(double[]); double \
cx024(double[]); double cx025(double[]); double cx026(double[]); \
double cx027(double[]); double cx028(double[]); double \
cx029(double[]); double cx031(double[]); double cx032(double[]); \
double cx033(double[]); double cx034(double[]); double \
cx035(double[]); double cx036(double[]); double cx037(double[]); \
double cx038(double[]); double cx039(double[]); double \
cx041(double[]); double cx042(double[]); double cx043(double[]); \
double cx044(double[]); double cx045(double[]); double \
cx046(double[]); double cx047(double[]); double cx048(double[]); \
double cx049(double[]); double cx051(double[]); double \
cx052(double[]); double cx053(double[]); double cx054(double[]); \
double cx055(double[]); double cx056(double[]); double \
cx057(double[]); double cx058(double[]); double cx059(double[]); \
double cx061(double[]); double cx062(double[]); double \
cx063(double[]); double cx064(double[]); double cx065(double[]); \
double cx066(double[]); double cx067(double[]); double \
cx068(double[]); double cx069(double[]);

double cx111(double[]); double cx112(double[]); double \
cx113(double[]); double cx114(double[]); double cx115(double[]); \
double cx116(double[]); double cx117(double[]); double \
cx118(double[]); double cx119(double[]); double cx121(double[]); \
double cx122(double[]); double cx123(double[]); double \
cx124(double[]); double cx125(double[]); double cx126(double[]); \
double cx127(double[]); double cx128(double[]); double \
cx129(double[]); double cx131(double[]); double cx132(double[]); \
double cx133(double[]); double cx134(double[]); double \
cx135(double[]); double cx136(double[]); double cx137(double[]); \
double cx138(double[]); double cx139(double[]); double \
cx141(double[]); double cx142(double[]); double cx143(double[]); \
double cx144(double[]); double cx145(double[]); double \
cx146(double[]); double cx147(double[]); double cx148(double[]); \
double cx149(double[]); double cx151(double[]); double \
cx152(double[]); double cx153(double[]); double cx154(double[]); \
double cx155(double[]); double cx156(double[]); double \
cx157(double[]); double cx158(double[]); double cx159(double[]); \
double cx161(double[]); double cx162(double[]); double \
cx163(double[]); double cx164(double[]); double cx165(double[]); \
double cx166(double[]); double cx167(double[]); double \
cx168(double[]); double cx169(double[]);

double cx211(double[]); double cx212(double[]); double \
cx213(double[]); double cx214(double[]); double cx215(double[]); \
double cx216(double[]); double cx217(double[]); double \
cx218(double[]); double cx219(double[]); double cx221(double[]); \
double cx222(double[]); double cx223(double[]); double \
cx224(double[]); double cx225(double[]); double cx226(double[]); \
double cx227(double[]); double cx228(double[]); double \
cx229(double[]); double cx231(double[]); double cx232(double[]); \
double cx233(double[]); double cx234(double[]); double \
cx235(double[]); double cx236(double[]); double cx237(double[]); \
double cx238(double[]); double cx239(double[]); double \
cx241(double[]); double cx242(double[]); double cx243(double[]); \
double cx244(double[]); double cx245(double[]); double \
cx246(double[]); double cx247(double[]); double cx248(double[]); \
double cx249(double[]); double cx251(double[]); double \
cx252(double[]); double cx253(double[]); double cx254(double[]); \
double cx255(double[]); double cx256(double[]); double \
cx257(double[]); double cx258(double[]); double cx259(double[]); \
double cx261(double[]); double cx262(double[]); double \
cx263(double[]); double cx264(double[]); double cx265(double[]); \
double cx266(double[]); double cx267(double[]); double \
cx268(double[]); double cx269(double[]);









#endif

