#ifndef SUBUNIT_HEADER
#define SUBUNIT_HEADER

#define GLYCAN_LENGTH 1
#define PEPTIDE_LENGTH 2

#define GLYCAN_SPRING_CONSTANT 5
#define PEPTIDE_SPRING_CONSTANT 6

#define PEPTIDE_PROB 1
#define GLYCAN_PROB 1

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