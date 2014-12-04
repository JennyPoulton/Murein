#include"Subunit.h"

Subunit::Subunit()
{
	set_length(1);
	
	set_number_of_bonds_created(0);

	set_spring_constant(1);
}

void Subunit::set_length(double l)
{
	length = l;
	return;
}

double Subunit::return_length()
{
	return length;
}

void Subunit::set_spring_constant(double k)
{
	spring_constant = k;
	return;
}

double Subunit::return_spring_constant()
{
	return spring_constant;
}

void Subunit::set_number_of_bonds_created(int n)
{
	number_of_bonds_created=n;
	return;
}

int Subunit::return_number_of_bonds_created()
{
	return number_of_bonds_created;
}