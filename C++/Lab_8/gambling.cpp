#include <iostream>
#include <string>

using namespace std;

string name;
string password;
string accountArray[100]= {};
int N=100;


string gambling (string arr[],int size);

int main(){



    gambling(accountArray,N);


    return 0;
}


string gambling (string arr[],int size){

cout << "Gambling feature allows you to swap accounts with someone with the same name. Give it a try!" << endl;
// if no one has the same name
//cout <<"sorry, your the only" <<name << endl;
cout << "What is the account name?: " << endl;
    // check if name is in the accountArray matrix
    cin >> name;
    cout << "Enter password: " << endl;
    // check password, if it doesnt match say sorry
    cin >> password; 
    // display balance, if bellow 500 say they are too poor to gamble
    // use rand to show a number 0-9, and that is going to be the index of the gambling account
    // call swap values.
    // show new balance.
    cout << "Now we will switch your balance with someone else!" << endl;
    


}