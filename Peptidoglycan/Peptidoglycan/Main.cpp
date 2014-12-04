#include"Polymer.h"

#include<iostream>
#include<fstream>
#include<ctime>
#include<iomanip>
#include <stdlib.h>

using namespace std;

int main(void)
{
	ofstream Output1("OutputForce.txt");

	if (!Output1.is_open())
	{
		cout << "Error: output file 1 cannot be opened." << endl;
		system("pause");
		return 1;
	}
	else
	{
		cout << "Output 1 file opened successfully." << endl;
	}

	ofstream Output2("OutputSpringConstant.txt");

	if (!Output2.is_open())
	{
		cout << "Error: output file 1 cannot be opened." << endl;
		system("pause");
		return 1;
	}
	else
	{
		cout << "Output 1 file opened successfully." << endl;
	}

	Polymer StaphAureous;

	double Force = 10;

	StaphAureous.Calculate_All_Forces_Upwards(Force/(double)2);
	StaphAureous.Calculate_All_Forces_Downwards(Force/(double)2);
	//StaphAureous.Calculate_All_Forces_Leftwards(Force / (double)2);
	//StaphAureous.Calculate_All_Forces_Rightwards(Force / (double)2);

	StaphAureous.Sort_Lengths_Into_Groups_For_Histogram();

	//double Histogram_Glycan[100];
	double Histogram_Peptide[100];

	for (int i = 0; i < 100; i++)
	{
		//Histogram_Glycan[i]=0;
		Histogram_Peptide[i]=0;
	}

	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			//Histogram_Glycan[(int)((StaphAureous.Return_Length_Glycan(i, j)-StaphAureous.Return_Min_Length_G()) / StaphAureous.Return_Splitter_For_Sorting_G())]++;
			double p = (StaphAureous.Return_Length_Peptide(i, j) - StaphAureous.Return_Min_Length_P());
			double q = StaphAureous.Return_Splitter_For_Sorting_P();

			int r = (int)(p /q);
			Histogram_Peptide[r]++;
		}
	}

	for (int i = 0; i < 100; i++)
	{
		Output1 /*<< i*StaphAureous.Return_Splitter_For_Sorting_G() << "\t" << Histogram_Glycan[i] << "\t"*/ << i*StaphAureous.Return_Splitter_For_Sorting_P() << "\t" << Histogram_Peptide[i] << endl;
	}

	system("pause");
	return 0;


}