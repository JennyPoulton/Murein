#ifndef MONOMER_HEADER
#define MONOMER_HEADER

#include"Subunit.h"

class Monomer
{
private:
	Subunit Glycan;
	Subunit Peptide;
	
	double Horizontal_Force;
	double Vertical_Force;

public:
	Monomer();//in which the initial lengths, spring constants and number of bonds are set for the 4 subunits, in which the forces are set to zero
	
	void Set_Horizontal_Force(double F);
	double Return_Horizontal_Force();
	
	void Set_Vertical_Force(double F);
	double Return_Vertical_Force();

	void Set_Length_Glycan(double x);
	double Return_Length_Glycan();
	void Set_Spring_Constant_Glycan(double x);
	double Return_Spring_Constant_Glycan();
	void Set_Number_Bonds_Glycan(int x);
	int Return_Number_Bonds_Glycan();

	void Set_Length_Peptide(double x);
	double Return_Length_Peptide();
	void Set_Spring_Constant_Peptide(double x);
	double Return_Spring_Constant_Peptide();
	void Set_Number_Bonds_Peptide(int x);
	int Return_Number_Bonds_Peptide();


	
};

#endif