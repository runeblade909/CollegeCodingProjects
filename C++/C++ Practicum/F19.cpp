#include <iostream>
#include <string>
#include <ctime>

using namespace std;
int Array1[29] = {};
int Array2[29] = {};
int const size_A = 30;
int minVal = {};
int maxVal = {};


int fill_array(int arr[], int min_ele, int max_ele, int size);

int fuse_universe(int ying[],int yang[],int size_y);

int main()
{
srand(time(NULL));
cout << "Hey there friend! I need some things from you!\nHow about you go ahead and enter a minumum and maximum value (in that order) and ill make an array for you. Okay? Cool." << endl;
cout << "Enter minumum value: ";
cin >> minVal;
cout << "Enter maximum value: ";
cin >> maxVal;

fill_array(Array1,minVal,maxVal,size_A);
cout << endl;

fill_array(Array2,minVal,maxVal,size_A);
cout << endl << endl;

fuse_universe(Array1,Array2,size_A);




return 0;

}

int fill_array(int arr[], int min_ele, int max_ele, int size)
{



    for (int i=0;i<=size-1;i++)
    {
        int num = rand()%max_ele;
        arr[i] = num;
    }


}

int fuse_universe(int ying[],int yang[],int size_y)
{
    int fusion_success =0;
    cout << "This is our Yang value:" << endl;
    for (int k=0; k<=size_y-1;k++) 
    {
        cout << yang[k] << " ";
    }
    cout << endl;

    cout << "Minus our Ying value:" << endl;
    for (int k=0; k<=size_y-1;k++) 
    {
        cout << ying[k] << " ";
    }
    cout << endl;

    cout << "Which gives us (adding up each number):" << endl;
    for (int i=0;i<=size_y-1;i++)
    {

    fusion_success += (yang[i] -ying[i]);
    cout << fusion_success << " ";
    }

    cout << endl << "Giving us a grand total of: " << fusion_success << endl;

}