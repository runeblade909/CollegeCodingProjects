#include <iostream>
#include <cmath>

using namespace std;

double n, guess, r;


int main() 
{
 cout << "Please enter a number I can find the square root of: ";

 cin >> n;

 guess = n/2;

 r = n/guess;

 guess = (guess +r)/2;

 r = n/guess;

 guess = (guess +r)/2;

 r = n/guess;

 guess = (guess +r)/2;

 r = n/guess;

 guess = (guess +r)/2;

 r = n/guess;

 guess = (guess + r)/2;

 printf("%.2f \n" , guess);
 
 return 0;
}