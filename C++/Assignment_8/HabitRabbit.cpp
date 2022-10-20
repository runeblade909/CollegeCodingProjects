#include <iostream>
#include <cmath>

using namespace std;

int main()
{
    int i;
 for(i=1;i<=170;i++){

     if(i%3==0 && i%7==0){

         cout<<"HabitRabbit"<<endl;

        }
        else if(i%3==0){

            cout<<"Habit"<<endl;

        }
        else if (i%7==0){

            cout<<"Rabbit"<<endl;

        }
        else{

            cout<<i<<endl;
            
        }
 }
 return 0;    
}