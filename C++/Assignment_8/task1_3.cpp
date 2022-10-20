#include <iostream>
#include <cmath>

using namespace std;

double weight,METS,Etime,mass,calories;

int main () 
{
 cout << "How much junk in yo trunk? (please enter your weight): " << endl;
 cin  >> weight;
 cout << "How many MET's is your excercise worth? (Because I wanna know, just curious, don't ask me questions):" << endl;
 cin >> METS;
 cout << "How many minutes did ya do it fer?:" << endl;
 cin >> Etime;

 mass = weight/2.2;

 calories = 0.0175*METS*mass*Etime;

 cout << "You burned " << calories << ", you should be proud." << endl;


 return 0;
}