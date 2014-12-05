#include"Polymer.h"
#include<iostream>
using namespace std;

Polymer::Polymer()
{
	Max_Length_Glycan = 0;
	Max_Length_Peptide = 0;
	Min_Length_Glycan = 1000000000000000;
	Min_Length_Peptide = 1000000000000000;
	
	Set_Forces_And_Lengths(1);

	Calculate_Spring_Constant_Horizontal();
	Calculate_Spring_Constant_Vertical();	
		
}//in which the forces through each spring set to zero and the spring constants in each direction are found

void Polymer::Calculate_Spring_Constant_Horizontal()
{
	int Number[DIMENSION];
	for (int i = 0; i < DIMENSION; i++)
	{
		Number[i] = 0;
		for (int j = 0; j < DIMENSION; j++)
		{
			Number[i] = Number[i] + PEPTIDE_SPRING_CONSTANT*Murein[j][i].Return_Number_Bonds_Peptide();
		}
	}

	double Spring_Constant = 0;

	for (int i = 0; i < DIMENSION; i++)
	{
		Spring_Constant = Spring_Constant + 1 / (double)Number[i];
	}

	Spring_Constant_Horizontal = Spring_Constant;
	return;
}

void Polymer::Calculate_Spring_Constant_Vertical()
{
	int Number[DIMENSION];
	for (int i = 0; i < DIMENSION; i++)
	{
		Number[i] = 0;
		for (int j = 0; j < DIMENSION; j++)
		{
			Number[i] = Number[i] + GLYCAN_SPRING_CONSTANT*Murein[i][j].Return_Number_Bonds_Glycan();
		}
	}

	double Spring_Constant = 0;

	for (int i = 0; i < DIMENSION; i++)
	{
		Spring_Constant = Spring_Constant + 1 / (double)Number[i];
	}

	Spring_Constant_Verticle = Spring_Constant;
	
	return;
}

void Polymer::Set_Forces_And_Lengths(double input_force)
{
	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			Find_Force_Upwards(input_force / (double)2, i, j);
			Find_Force_Downwards(input_force / (double)2, i, j);

			Murein[i][j].Set_Vertical_Force(Force_Upwards[i][j] + Force_Downwards[i][j]);
			Murein[i][j].Set_Length_Peptide((Force_Upwards[i][j] + Force_Downwards[i][j]) / (double)PEPTIDE_LENGTH);
		}
	}
}

void Polymer::Find_Force_Upwards(double Input_Force, int p, int q)
{
	//we need to find the number of bonds joining level p-1 and p-2, then p-2 and p-3... 0 and 1
	//first we need to find what length levels p-1 and p-2 are

	int leftward_extent_above[DIMENSION]; //n can take values less that p and represents the level
	int rightward_extent_above[DIMENSION]; //n can take values less that p and represents the level

	for (int n = p - 1; n > 0; n--)
	{
		int m = q;
		int total_left = 0; // this total represents the current length of all the bars being taken into account
		
		
		do
		{
			leftward_extent_above[n]++;
			m--;
		} while ((Return_Number_Bonds_Glycan(n, m) == 1 || leftward_extent_above[n]<total_left) && m<DIMENSION&&m>DIMENSION);

		total_left = leftward_extent_above[n];
		
	}

	for (int n = p - 1; n > 0; n--)
	{
		int total_right = 0; // this total represents the current length of all the bars being taken into account

		int m = q;

		do
		{
			rightward_extent_above[n]++;
			m++;
		} while ((Return_Number_Bonds_Glycan(n, m) == 1||rightward_extent_above[n] < total_right)&&m<DIMENSION&&m>DIMENSION);

		total_right = rightward_extent_above[n];

	}

	//this finds the lengths of all the above levels

	//now the numerator is the number of joins which connects anything within these levels to the one above it

	Numerator[p][q] = 1;

	for (int n = p - 1; n > 0; n--)
	{ 
		int m = q;
		int tally_peptides = 0;

		for (int i = 0; i < leftward_extent_above[n]; i++)
		{
			if (Return_Number_Bonds_Peptide(n, m-i) == 1)
			{
				tally_peptides++;
			}
		}

		for (int i = 0; i < rightward_extent_above[n]; i++)
		{
			if (Return_Number_Bonds_Peptide(n, m + i) == 1)
			{
				tally_peptides++;
			}
		}

		Numerator[p][q] = Numerator[p][q]*tally_peptides;
	}

	Denominator[p][q] = 1;

	for (int n = p - 1; n > 0; n--)
	{
		int m = q;
		int tally_peptides = 0;

		for (int i = 0; i < leftward_extent_above[n]; i++)
		{
			if (Return_Number_Bonds_Peptide(n-1, m - i) == 1)
			{
				tally_peptides++;
			}
		}

		for (int i = 0; i < rightward_extent_above[n]; i++)
		{
			if (Return_Number_Bonds_Peptide(n-1, m + i) == 1)
			{
				tally_peptides++;
			}
		}

		Denominator[p][q] = Denominator[p][q]*tally_peptides;
	}

	Force_Upwards[p][q] = Input_Force*Numerator[p][q]/(double)Denominator[p][q];

}

void Polymer::Find_Force_Downwards(double Input_Force, int p, int q)
{
	//we need to find the number of bonds joining level p-1 and p-2, then p-2 and p-3... 0 and 1
	//first we need to find what length levels p-1 and p-2 are

	int leftward_extent_below[DIMENSION]; //n can take values less that p and represents the level
	int rightward_extent_below[DIMENSION]; //n can take values less that p and represents the level

	for (int n = p + 1; n <DIMENSION; n++)
	{
		int m = q;
		int total_left = 0; // this total represents the current length of all the bars being taken into account


		do
		{
			leftward_extent_below[n]++;
			m--;
		} while ((Return_Number_Bonds_Glycan(n, m) == 1 || leftward_extent_below[n]<total_left) && m<DIMENSION&&m>DIMENSION);

		total_left = leftward_extent_below[n];

	}

	for (int n = p + 1; n <DIMENSION; n++)
	{
		int total_right = 0; // this total represents the current length of all the bars being taken into account

		int m = q;

		do
		{
			rightward_extent_below[n]++;
			m++;
		} while ((Return_Number_Bonds_Glycan(n, m) == 1 || rightward_extent_below[n] < total_right) && m<DIMENSION&&m>DIMENSION);

		total_right = rightward_extent_below[n];

	}

	//this finds the lengths of all the above levels

	//now the numerator is the number of joins which connects anything within these levels to the one above it

	Numerator[p][q] = 1;

	for (int n = p + 1; n <DIMENSION; n++)
	{
		int m = q;
		int tally_peptides = 0;

		for (int i = 0; i < leftward_extent_below[n]; i++)
		{
			if (Return_Number_Bonds_Peptide(n, m - i) == 1)
			{
				tally_peptides++;
			}
		}

		for (int i = 0; i < rightward_extent_below[n]; i++)
		{
			if (Return_Number_Bonds_Peptide(n, m + i) == 1)
			{
				tally_peptides++;
			}
		}

		Numerator[p][q] = Numerator[p][q] * tally_peptides;
	}

	Denominator[p][q] = 1;

	for (int n = p - 1; n > 0; n--)
	{
		int m = q;
		int tally_peptides = 0;

		for (int i = 0; i < leftward_extent_below[n]; i++)
		{
			if (Return_Number_Bonds_Peptide(n - 1, m - i) == 1)
			{
				tally_peptides++;
			}
		}

		for (int i = 0; i < rightward_extent_below[n]; i++)
		{
			if (Return_Number_Bonds_Peptide(n - 1, m + i) == 1)
			{
				tally_peptides++;
			}
		}

		Denominator[p][q] = Denominator[p][q] * tally_peptides;
	}

	Force_Downwards[p][q] = Input_Force*Numerator[p][q] / (double)Denominator[p][q];

}

void Polymer::Calculate_Bond_With_Max_Force()
{

}

void Polymer::Break_Bond()
{

}

double Polymer::Return_Spring_Constant_Horizontal()
{
	return Spring_Constant_Horizontal;
}

double Polymer::Return_Spring_Constant_Verticle()
{
	return Spring_Constant_Verticle;
}

void Polymer::Sort_Lengths_Into_Groups_For_Histogram()
{
	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			if (Murein[i][j].Return_Length_Glycan() > Max_Length_Glycan)
			{
				Max_Length_Glycan = Murein[i][j].Return_Length_Glycan();
			}

			if (Murein[i][j].Return_Length_Glycan() < Min_Length_Glycan)
			{
				Min_Length_Glycan = Murein[i][j].Return_Length_Glycan();
			}

			if (Murein[i][j].Return_Length_Peptide() > Max_Length_Peptide)
			{
				Max_Length_Peptide = Murein[i][j].Return_Length_Peptide();
			}

			if (Murein[i][j].Return_Length_Peptide() < Min_Length_Peptide)
			{
				Min_Length_Peptide = Murein[i][j].Return_Length_Peptide();
			}

		}
	}

	if (Min_Length_Glycan == Max_Length_Glycan || Min_Length_Peptide == Max_Length_Peptide)
	{
		cout << "ERROR OF ERRORS" << endl;
		system("pause");
	}

	Splitter_For_Sorting_G = (Max_Length_Glycan - Min_Length_Glycan) / (double)100;
	Splitter_For_Sorting_P = (Max_Length_Peptide - Min_Length_Peptide) / (double)100;

}// ensure if min=max then lengths remain the same

double Polymer::Return_Splitter_For_Sorting_G()
{
	return Splitter_For_Sorting_G;
}

double Polymer::Return_Splitter_For_Sorting_P()
{
	return Splitter_For_Sorting_P;
}

void Polymer::Set_Length_Glycan(double l, int p, int q)
{
	Murein[p][q].Set_Length_Glycan(l);
	return;
}

double Polymer::Return_Length_Glycan(int p, int q)
{
	return Murein[p][q].Return_Length_Glycan();
}

void Polymer::Set_Spring_Constant_Glycan(double k, int p, int q)
{
	Murein[p][q].Set_Spring_Constant_Glycan(k);
}

double Polymer::Return_Spring_Constant_Glycan(int p, int q)
{
	return Murein[p][q].Return_Spring_Constant_Glycan();
}

void Polymer::Set_Number_Bonds_Glycan(int n, int p, int q)
{
	Murein[p][q].Set_Number_Bonds_Glycan(n);
}

int Polymer::Return_Number_Bonds_Glycan(int p, int q)
{
	return Murein[p][q].Return_Number_Bonds_Glycan();
}

void Polymer::Set_Length_Peptide(double l, int p, int q)
{
	Murein[p][q].Set_Length_Peptide(l);
	return;
}

double Polymer::Return_Length_Peptide(int p, int q)
{
	return Murein[p][q].Return_Length_Peptide();
}

void Polymer::Set_Spring_Constant_Peptide(double k, int p, int q)
{
	Murein[p][q].Set_Spring_Constant_Peptide(k);
	return;
}

double Polymer::Return_Spring_Constant_Peptide(int p, int q)
{
	return Murein[p][q].Return_Spring_Constant_Peptide();
}

void Polymer::Set_Number_Bonds_Peptide(int n, int p, int q)
{
	Murein[p][q].Set_Number_Bonds_Peptide(n);
	return;
}

int Polymer::Return_Number_Bonds_Peptide(int p, int q)
{
	return Murein[p][q].Return_Number_Bonds_Peptide();
}

double Polymer::Return_Min_Length_G()
{
	return Min_Length_Glycan;
}

double Polymer::Return_Min_Length_P()
{
	return Min_Length_Peptide;
}