#include <iostream>
#include <string>

using namespace std;

int accountName[100] = {0};
int const N=100;
string name;
string password;
int money;

int withdraw(int arr[], int length);

int main(){

    withdraw(accountName,N);



    return 0;
}



int withdraw(int arr[], int length){

    cout << "What is the account name?: " << endl;
    // check if name is in the accountArray matrix
    cin >> name;
    cout << "Enter password: " << endl;
    // check password, if it doesnt match say sorry
    cin >> password;
    cout << "Enter amount you wish to withdraw: " << endl;
    // add deposited ammount to the total
    cin >> money;
    cout <<"Thank you " << name << " you have taken " << money << " dollars from your account." << endl;



}