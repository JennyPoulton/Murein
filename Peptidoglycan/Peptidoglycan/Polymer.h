#ifndef POLYMER_HEADER
#define POLYMER_HEADER

#define DIMENSION 100

#include"Monomer.h"

class Polymer
{
private:
	Monomer Murein[DIMENSION][DIMENSION];

	int Numerator[DIMENSION][DIMENSION];
	int Denominator[DIMENSION][DIMENSION];

	double Force_Upwards[DIMENSION][DIMENSION];
	double Force_Downwards[DIMENSION][DIMENSION];

	void Find_Force_Upwards(double Input_Force, int p, int q);
	void Find_Force_Downwards(double Input_Force, int p, int q);

	double Spring_Constant_Verticle;
	double Spring_Constant_Horizontal;

	double Splitter_For_Sorting_G;
	double Splitter_For_Sorting_P;
	double Max_Length_Glycan;
	double Min_Length_Glycan;
	double Max_Length_Peptide;
	double Min_Length_Peptide;

public:
	Polymer();//in which the forces through each spring set to zero and the spring constants in each direction are found, max/min lengths are set

	void Calculate_Spring_Constant_Horizontal();
	void Calculate_Spring_Constant_Vertical();

	void Set_Forces_And_Lengths(double input_force);

	void Calculate_Bond_With_Max_Force();
	void Break_Bond();

	double Return_Spring_Constant_Horizontal();
	double Return_Spring_Constant_Verticle();

	void Sort_Lengths_Into_Groups_For_Histogram(); // ensure if min=max then lengths remain the same
	double Return_Splitter_For_Sorting_G();
	double Return_Splitter_For_Sorting_P();
	double Return_Min_Length_G();
	double Return_Min_Length_P();

	void Set_Length_Glycan(double l, int p, int q);
	double Return_Length_Glycan(int p, int q);
	void Set_Spring_Constant_Glycan(double k, int p, int q);
	double Return_Spring_Constant_Glycan(int p, int q);
	void Set_Number_Bonds_Glycan(int n, int p, int q);
	int Return_Number_Bonds_Glycan(int p, int q);

	void Set_Length_Peptide(double l, int p, int q);
	double Return_Length_Peptide(int p, int q);
	void Set_Spring_Constant_Peptide(double k, int p, int q);
	double Return_Spring_Constant_Peptide(int p, int q);
	void Set_Number_Bonds_Peptide(int n, int p, int q);
	int Return_Number_Bonds_Peptide(int p, int q);
};

#endif