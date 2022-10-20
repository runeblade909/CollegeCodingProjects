#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>

using namespace std;

int array[] = {8,5,4,4,2,4,2,5,9,4,6,81,4,8,12,4,8};
int sizeArr = 17;
int min_ele ={};

void practicuum (int arr[], int size, int &min_ele);



int main(){


    practicuum (array,sizeArr,min_ele);
    cout << "Our minimum element is " << min_ele << endl;

    return 0;
}


void practicuum (int arr[], int size, int &min_ele)
{
    for (int y =0; y<=size-1;y++)
    {
        cout << arr[y] << " ";
    }
    cout << endl;
    int temp;
    for (int k=0; k<=size-1 ; k++)
    {
        for(int i=0; i<=size-2 ;i++)
        {
    
                if (arr[i] > arr[i+1])
            {
               
            temp = arr[i];
            arr[i] = arr[i+1];
            arr[i+1] = temp;
             }
             else
            {
            arr[i] = arr[i];
            }
                
         cout << arr[i] << " ";
        }
        cout << endl;
    }   
    for (int z =0; z<=size-1;z++)
    {
        cout << arr[z] << " ";
    }
    cout << endl;

    cout << "The final array max is: " << arr[size-1] << endl;
    min_ele = arr[0];
}