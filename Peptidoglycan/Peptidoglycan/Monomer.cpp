#include"Monomer.h"
#include<stdlib.h>
#include<time.h>




Monomer::Monomer()
{
	srand(time(NULL));

	Set_Length_Glycan(GLYCAN_LENGTH);
	Set_Length_Peptide(PEPTIDE_LENGTH);
	Set_Spring_Constant_Glycan(GLYCAN_SPRING_CONSTANT);
	Set_Spring_Constant_Peptide(PEPTIDE_SPRING_CONSTANT);

	Set_Number_Bonds_Glycan(0);
	Set_Number_Bonds_Peptide(0);
		
	if (rand()/(double)RAND_MAX <= PEPTIDE_PROB)
	{
		Set_Number_Bonds_Peptide(1);
	}
	if (rand()/(double)RAND_MAX <= GLYCAN_PROB)
	{
		Set_Number_Bonds_Glycan(1);
	}
	
}
//in which the initial lengths, spring constants and number of bonds are set for the 4 subunits, in which the forces are set to zero

void Monomer::Set_Horizontal_Force(double F)
{
	Horizontal_Force = F;
	return;
}

double Monomer::Return_Horizontal_Force()
{
	return Horizontal_Force;
}

void Monomer::Set_Vertical_Force(double F)
{
	Vertical_Force = F;
	return;
}

double Monomer::Return_Vertical_Force()
{
	return Vertical_Force;
}

void Monomer::Set_Length_Glycan(double x)
{
	Glycan.set_length(x);
	return;
}

double Monomer::Return_Length_Glycan()
{
	return Glycan.return_length();
}

void Monomer::Set_Spring_Constant_Glycan(double x)
{
	Glycan.set_spring_constant(x);
	return;
}

double Monomer::Return_Spring_Constant_Glycan()
{
	return Glycan.return_spring_constant();
}

void Monomer::Set_Number_Bonds_Glycan(int x)
{
	Glycan.set_number_of_bonds_created(x);
	return;
}

int Monomer::Return_Number_Bonds_Glycan()
{
	return Glycan.return_number_of_bonds_created();
}

void Monomer::Set_Length_Peptide(double x)
{
	Peptide.set_length(x);
	return;
}

double Monomer::Return_Length_Peptide()
{
	return Peptide.return_length();
}

void Monomer::Set_Spring_Constant_Peptide(double x)
{
	Peptide.set_spring_constant(x);
	return;
}

double Monomer::Return_Spring_Constant_Peptide()
{
	return Peptide.return_spring_constant();
}

void Monomer::Set_Number_Bonds_Peptide(int x)
{
	Peptide.set_number_of_bonds_created(x);
	return;
}

int Monomer::Return_Number_Bonds_Peptide()
{
	return Peptide.return_number_of_bonds_created();
}