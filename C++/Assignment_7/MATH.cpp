#include <iostream>

using namespace std;

int main()
{

    float number,number2,Add,Sub,Mult,Div;

    cout << "Input 2 numbers for math" << endl;

    cin >> number >> number2;

    Add = number + number2;

    Sub = number - number2;

    Mult = number * number2;

    Div = number / number2;

    cout << "" << number << "+" << number2 <<"="<< Add <<endl;

    cout << "" << number << "-" << number2 <<"="<< Sub <<endl;

    cout << "" << number << "x" << number2 <<"="<< Mult <<endl;

    cout << "" << number << "/" << number2 <<"="<< Div <<endl;

    return 0;




}