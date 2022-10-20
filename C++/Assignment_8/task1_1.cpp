#include <iostream>
#include <cmath>

using namespace std;

int main()
{
 double x,Sout,Hout,Mout,HourSec;
 int Minutes,Hours;
    cout << "Give me a time in seconds, and I will convert it for you!"<<endl;
    cin >> x;

    Minutes = (x)/60;
    
    Hours= Minutes/60;

    HourSec = (Hours * 60 * 60);

    

    Mout = Minutes - (Hours *60)  ;

    Sout = (x - HourSec) -(Mout * 60);

    cout << "So that is " << Hours << " Hours," << Mout << " minutes and " << Sout << " seconds." << endl;

 return 0;
}


