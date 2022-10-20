#include <iostream>
#include <cmath>
using namespace std;


int main()
{

	double x, a=2;
	cout << "Please enter a value for x = " ;
	cin >> x; 

	cout << "x to the power of " << a;
	cout << " is equal to " << pow(x,++a) << endl;

	cout << "x to the power of " << a;
	cout << " is equal to " << pow(x,++a) << endl;

	cout << "x to the power of " << a;
	cout << " is equal to " << pow(x,++a) << endl;

	return 0;
}

