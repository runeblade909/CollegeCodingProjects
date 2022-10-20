#include <iostream>
#include <string>
#include <cmath>

using namespace std;

int arr[] = {4,3,1,5,3,3,1,8,3};
int sizeArr = 9;
int arr1[] = {4,8,1,5,8,3,1,8,3};
int arr2[] = {2,2,7,0,1,9,8,5};
int sizeArr2 = 8;
int arr3[] = {4,1,5,4,8};
int sizeArr3 = 5;
int checker( int arr[], int size, int var);

int main()
{

checker(arr,sizeArr,3);
checker(arr1,sizeArr,8);
checker(arr2,sizeArr2,2);
checker(arr3,sizeArr3,4);
cout << endl;


}


int checker( int arr[], int size, int var)
{
    for (int k=0; k<= size-1; k++)
    {
        cout << arr[k] << " ";
    }
    cout << endl << "Now we will find numbers that match our variable number: " << var << endl;

    for (int i = 0; i<= size-1; i++)
    {
        if (arr[i] == var)
        {
            arr[i] = pow(arr[i],3);
        }
        else
        {
            arr[i] = arr[i];
        }
        
        cout << arr[i] << " ";
    }



}