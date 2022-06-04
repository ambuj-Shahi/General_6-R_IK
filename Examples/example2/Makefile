VPATH =  ../../Coeff_x0 ../../Coeff_x1 ../../Coeff_x2 ../Residue ../Theta ../IK
CC=gcc
LD=gcc
RM=rm -f
CFLAGS=-I$(Coeff_x0)/include -I$(Coeff_x1)/include -I$(Coeff_x2)/include -I$(Residue)/include -I$(Theta)/include -I$(IK)/include -Wall -g -O3 


OBJS = example2.o cx011.o cx012.o cx013.o cx014.o cx015.o cx016.o cx017.o cx018.o cx019.o cx021.o cx022.o cx023.o cx024.o cx025.o cx026.o cx027.o cx028.o cx029.o cx031.o cx032.o cx033.o cx034.o cx035.o cx036.o cx037.o cx038.o cx039.o cx041.o cx042.o cx043.o cx044.o cx045.o cx046.o cx047.o cx048.o cx049.o cx051.o cx052.o cx053.o cx054.o cx055.o cx056.o cx057.o cx058.o cx059.o cx061.o cx062.o cx063.o cx064.o cx065.o cx066.o cx067.o cx068.o cx069.o cx111.o cx112.o cx113.o cx114.o cx115.o cx116.o cx117.o cx118.o cx119.o cx121.o cx122.o cx123.o cx124.o cx125.o cx126.o cx127.o cx128.o cx129.o cx131.o cx132.o cx133.o cx134.o cx135.o cx136.o cx137.o cx138.o cx139.o cx141.o cx142.o cx143.o cx144.o cx145.o cx146.o cx147.o cx148.o cx149.o cx151.o cx152.o cx153.o cx154.o cx155.o cx156.o cx157.o cx158.o cx159.o cx161.o cx162.o cx163.o cx164.o cx165.o cx166.o cx167.o cx168.o cx169.o cx211.o cx212.o cx213.o cx214.o cx215.o cx216.o cx217.o cx218.o cx219.o cx221.o cx222.o cx223.o cx224.o cx225.o cx226.o cx227.o cx228.o cx229.o cx231.o cx232.o cx233.o cx234.o cx235.o cx236.o cx237.o cx238.o cx239.o cx241.o cx242.o cx243.o cx244.o cx245.o cx246.o cx247.o cx248.o cx249.o cx251.o cx252.o cx253.o cx254.o cx255.o cx256.o cx257.o cx258.o cx259.o cx261.o cx262.o cx263.o cx264.o cx265.o cx266.o cx267.o cx268.o cx269.o theta1.o theta2.o theta6.o residue.o InverseKinematicsGeneral6R.o csvparser.o 
 
DEPS = cx.h csvparser.h


MAIN = example2

$(MAIN):$(OBJS)
	$(LD) -pg -I/usr/include -L/usr/local/lib $(OBJS) $(CFLAGS) -lgsl -lgslcblas -lm -o $(MAIN) 

%.o: %.c $(DEPS)
	$(LD) $(CFLAGS) -c $< -pg

clean:
	$(RM) *.o stats timings roots testsuite *.bak *~ core *.res
	$(RM)  *.o
