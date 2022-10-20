#include <iostream>
#include <ctime>
#include <string>


void menuTest(int arr[], int length); // function to open menu
void newAccount (char arr[],int length); // function to open account creator 
void deposit(int arr[], int length); // function to allow deposits
void withdraw (int arr[], int length); // function to allow withdraw
void gamblingldx (int idk); // Gambling???
void swapValues(); // I guess steal from the bank?

int menuChoice ;
int accountName ; 
char accountArray[100];

using namespace std;

    int main(){
    int const N = 6;
    int const M = 100;
    int myArray[] = {1,2,3,4,5,6};
    menuTest(myArray,N);

        if (menuChoice ==1){

            newAccount(accountArray,M);
        }

        if (menuChoice ==2){

            
        }

        if (menuChoice ==3){


        }

        if (menuChoice ==4){


        }

        if (menuChoice ==5){


        }

        if (menuChoice ==6){


        }

    }


    void menuTest(int arr[], int length) {

    cout << "Welcome to your bank, what are we doing today?" << endl;
    cout << "1. Making an account\n2. Deposit\n3. Withdraw\n4. Balance\n5. Gambling\n6. Finish" << endl; 
    cout << "Menu choice: " << endl;
    cin >> menuChoice;

    //Here we will have a loop for each of the 6 options
    cout << "You have chosen to do number " << menuChoice << endl << endl;

}
    void newAccount (char arr[],int length){
    int j = 0;
    cout << "Please enter an account name: \n";
    cin >> accountName;

    accountName = accountArray[j];
    cout << accountArray;
    j=j++;
    

}

    void deposit(int arr[], int length){



    }

    void withdraw (int arr[], int length){


    }

    void gamblingldx (int idk){


    }

    void swapValues(){


        
    }









