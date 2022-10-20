#include<iostream>
#include<map>
#include<cstdlib>
using namespace std;

string accountName[100] = { };
int accountBalance[100] = {0};
int accountPassword[100] = {0};
int gambling[10] = {0, 500, 10, 1000, 20, 5000, 30, 10000, 40, 50000};
int n=0;
map<string, int>accountNo; //stores(accountName-->accountNo)

void newAccount();
void deposit(int a);
void withdraw(int a);
void balance(int a);
int gamblingIdx();
void swapValues(int a, int g);
void gamble(int a);

int main(){
   string s;
   int menuChoice;
   while(true)
   {
       cout<<"\n-- ATM MENU LIST --\n1. Make an account\n2. Deposit\n3. Withdraw"<<endl;
       cout<<"4. Balance\n5. Gambling\n6. Finish\nSelect an option: ";
       cin>>menuChoice;
       if(menuChoice == 6)
           break;
       switch(menuChoice)
       {
           case 1:
               if(n < 100)
                   newAccount();
               else
                   cout<<"The bank account database is full\n";
               break;
          
           case 2:
               cout<<"\nEnter your name: ";
               cin>>ws;
               getline(cin, s);
               if(accountNo.find(s) != accountNo.end())
                   deposit(accountNo[s]);
               else
                   cout<<"No account!\n";
               break;
          
           case 3:
               cout<<"\nEnter your name: ";
               cin>>ws;
               getline(cin, s);
               if(accountNo.find(s) != accountNo.end())
                   withdraw(accountNo[s]);
               else
                   cout<<"No account!\n";
               break;
          
           case 4:
               cout<<"\nEnter your name: ";
               cin>>ws;
               getline(cin, s);
               if(accountNo.find(s) != accountNo.end())
                   balance(accountNo[s]);
               else
                   cout<<"No account!\n";
               break;
          
           case 5:
               cout<<"\nEnter your name: ";
               cin>>ws;
               getline(cin, s);
               if(accountNo.find(s) != accountNo.end())
                   gamble(accountNo[s]);
               else
                   cout<<"No account!\n";
               break;
           default:
               cout<<"Invalid option\n";
       }
   }
   return 0;
}

void newAccount(){
   string s;
   int p;
   cout<<"\nEnter your name: ";
   cin>>ws;
   getline(cin, s);
   if(accountNo.find(s) == accountNo.end())
   {
       accountNo[s] = n;
       accountName[n] = s;
       cout<<"Enter your password: ";
       cin>>p;
       accountPassword[n] = p;
       n++;
   }
   else
       cout<<"Account name already exists\n";
}
void deposit(int a){
   int p;
   cout<<"Enter your password: ";
   cin>>p;
   if(accountPassword[a] == p)
   {
       int deposit;
       cout<<"\nEnter the deposit amount: ";
       cin>>deposit;
       if(deposit >= 0)
           accountBalance[a] += deposit;
       else
           cout<<"Deposit amount cannot be negative\n";
   }
   else
       cout<<"Wrong password!\n";
}

void withdraw(int a){
   int p;
   cout<<"Enter your password: ";
   cin>>p;
   if(accountPassword[a] == p)
   {
       while(true)
       {
           cout<<"\nThe current balance is "<<accountBalance[a]<<".\n";
           int withdraw;
           cout<<"Enter the withdraw amount: ";
           cin>>withdraw;
           if(withdraw>=0 && withdraw<=accountBalance[a])
           {
               accountBalance[a] -= withdraw;
               break;
           }
       }
   }
   else
       cout<<"Wrong password!\n";
}

void balance(int a){
   int p;
   cout<<"Enter your password: ";
   cin>>p;
   if(accountPassword[a] == p)
       cout<<"The balance of the account is "<<accountBalance[a]<<".\n";
   else
       cout<<"Wrong password!\n";
}

int gamblingIdx(){
   return rand()%10;
}
void swapValues(int a, int g){
   int t = accountBalance[a];
   accountBalance[a] = gambling[g];
   gambling[g] = t;
}
void gamble(int a){
   int p;
   cout<<"Enter your password: ";
   cin>>p;
   if(accountPassword[a] == p)
   {
       cout<<"\nThe current balance is "<<accountBalance[a]<<".\n";
       if(accountBalance[a] >= 500)
       {
           int g = gamblingIdx();
           swapValues(a, g);
           cout<<"The balance of the account is "<<accountBalance[a]<<".\n";
       }
       else
           cout<<"The balance is less than the minimum balance for gambling!";
   }
   else
       cout<<"Wrong password!\n";
}