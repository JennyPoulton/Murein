#ifndef SUBUNIT_HEADER
#define SUBUNIT_HEADER

class Subunit
{
	private:
		double spring_constant;
		double length;
		int number_of_bonds_created;

	public:
		Subunit();

		void set_spring_constant(double k);
		double return_spring_constant();

		void set_length(double l);
		double return_length();

		void set_number_of_bonds_created(int n);
		int return_number_of_bonds_created();

};

#endif