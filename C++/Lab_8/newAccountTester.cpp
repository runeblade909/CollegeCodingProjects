#include <iostream>
#include <string>
using namespace std;

string accountArray[100];
int N = 100;
string accountName;
int j=0;
string newAccount (string accountArray[],int length);

int main() {   
    string accountArray[100] = {};
    int j=0;
    newAccount(accountArray,N);

    cout << "Your account name is: " << accountArray[j] <<endl;



    return 0;
}
    
    
    //if (menuChoice == 2){}


string newAccount (string accountArray[],int length){
    int j = 0;
    cout << "Please enter an account name: \n";
    cin >> accountArray[j];


    
    return accountArray[j];
    //j=j++;

}

    