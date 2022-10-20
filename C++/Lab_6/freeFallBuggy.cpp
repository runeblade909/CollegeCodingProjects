// Program for calculating freefall distance as a function of time
// d = 1/2 * g * t^2


#include <iostream>
#include <string>

using namespace std;

int main( )
{
	double d,t;
	const double g = 9.8;
	

	cout << "enter time: " << endl; 

	cin >> t;

	d = 1.0/2 * g * t; 

	cout << "distance = " << d << endl;


	return 0;
}
