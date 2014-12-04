#include"Polymer.h"
#define PEPTIDEK 6
#define GLYCANK 7
#include<iostream>
using namespace std;

Polymer::Polymer()
{
	Max_Length_Glycan = 0;
	Max_Length_Peptide = 0;
	Min_Length_Glycan = 1000000000000000;
	Min_Length_Peptide = 1000000000000000;
	
	Calculate_All_Forces_Upwards(0);
	Calculate_All_Forces_Downwards(0);
	Calculate_All_Forces_Leftwards(0);
	Calculate_All_Forces_Rightwards(0);

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
			Number[i] = Number[i] + PEPTIDEK*Murein[j][i].Return_Number_Bonds_Peptide();
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
			Number[i] = Number[i] + GLYCANK*Murein[i][j].Return_Number_Bonds_Glycan();
		}
	}

	double Spring_Constant = 0;

	for (int i = 0; i < DIMENSION; i++)
	{
		Spring_Constant = Spring_Constant + 1 / (double)Number[i];
	}

	Spring_Constant_Verticle = Spring_Constant;

	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			int n = 0;
			if (Current_Bar_Peptide[i][j] != 0)
			{
				n = Current_Bar_Peptide[i][j];
			}
			else if (Current_Bar_Peptide[i][j] == 0)
			{
				Current_Bar_Peptide[i][j] = n;
			}
		}
	}

	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			int n = 0;
			if (Bar_Above_Peptide[i][j] != 0)
			{
				n = Bar_Above_Peptide[i][j];
			}
			else if (Bar_Above_Peptide[i][j] == 0)
			{
				Bar_Above_Peptide[i][j] = n;
			}
		}
	}

	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			int n = 0;
			if (Bar_Below_Peptide[i][j] != 0)
			{
				n = Bar_Below_Peptide[i][j];
			}
			else if (Bar_Below_Peptide[i][j] == 0)
			{
				Bar_Below_Peptide[i][j] = n;
			}
		}
	}
	return;
}

void Polymer::Set_Up_Bars_Peptide()
{
	for (int i = 0; i < DIMENSION; i++)
	{
		int n = 0;
		for (int j = 0; j < DIMENSION; j++)
		{
			Current_Bar_Peptide[i][j] = 0;
			Bar_Above_Peptide[i][j] = 0;
			Bar_Below_Peptide[i][j] = 0;

			if (Murein[i][j].Return_Number_Bonds_Peptide() == 1)
			{
				Current_Bar_Peptide[i][n]++;
			}
			else if (Murein[i][j].Return_Number_Bonds_Peptide() == 0)
			{
				n = j + 1;
			}//this gives length of bar at start of bar and zero elseware
		}
	}

	for (int i = 1; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			int n = j;

			if (Current_Bar_Peptide[i][j] != 0)
			{
				do
				{
					n--;
				} while (Current_Bar_Peptide[i - 1][n] != 0);

				int Total = 0;

				do
				{
					Total = Total + Current_Bar_Peptide[i - 1][n];
					n = n + Current_Bar_Peptide[i - 1][n];
				} while (Current_Bar_Peptide[i][n] > Total);

				Bar_Above_Peptide[i][j] = Total;
			}//This gives length of bars above at start of each bar and 0 elseware
		}
	}

	for (int i = 0; i < DIMENSION-1; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			int n = j;

			if (Current_Bar_Peptide[i][j] != 0)
			{
				do
				{
					n--;
				} while (Current_Bar_Peptide[i + 1][n] != 0);

				int Total = 0;
				do
				{
					Total = Total + Current_Bar_Peptide[i + 1][n];
					n = n + Current_Bar_Peptide[i + 1][n];
				} while (Current_Bar_Peptide[i][n] > Total);

				Bar_Below_Peptide[i][j] = Total;
			}//This gives length of bars below at start of each bar and zero elseware
		}
	}
	
	Bar_Above_Peptide[0][0] = DIMENSION;
	Bar_Below_Peptide[DIMENSION - 2][0] = DIMENSION;
}

void Polymer::Calculate_All_Forces_Upwards(double Total_Force)
{
	Set_Up_Bars_Peptide();

	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			int Numerator = 1;
			int Denominator = 1;

			for (int n = i; n > 0; n--)
			{
				Numerator = Numerator*Current_Bar_Peptide[n][j];
				Denominator = Denominator*Bar_Above_Peptide[n][j];
			}

			Murein[i][j].Set_Vertical_Force(Total_Force*(double)Numerator/(double)Denominator);
		}


	}
}

void Polymer::Calculate_All_Forces_Downwards(double Total_Force)
{
	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			int Numerator = 1;
			int Denominator = 1;
			for (int n = i; n < DIMENSION; n++)
			{
				Numerator = Numerator*Current_Bar_Peptide[n][j];
				Denominator = Denominator*Bar_Below_Peptide[n][j];
			}

			Murein[i][j].Set_Vertical_Force(Total_Force*(double)Numerator/(double)Denominator);
		}

		
	}

	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			Murein[i][j].Set_Length_Peptide(Murein[i][j].Return_Vertical_Force()/(double)PEPTIDEK);
		}
	}
}//this will give total force and set lengths

void Polymer::Calculate_All_Forces_Leftwards(double Total_Force)
{

}

void Polymer::Calculate_All_Forces_Rightwards(double Total_Force)
{

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