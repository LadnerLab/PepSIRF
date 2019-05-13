#ifndef STUDENT_HH_INCLUDED
#define STUDENT_HH_INCLUDED
#include <string>

class Student
{
	private:
		std::string name;
	public:
		Student( std::string );
		virtual void display();

}

#endif // STUDENT_HH_INCLUDED
