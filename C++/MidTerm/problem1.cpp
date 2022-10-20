#include <iostream>
using namespace std;





void removeDivByFive(/*fill this with pointer reference to the array, int reference to the length */int*arr, int &size)
{
      int count = 0;
    for(int i = 0;i<size;i++)
    {
        if(arr[i]%5 ==0){
            arr[i] = 0;
            count ++;


        }

    }
    cout << count << endl;
    int *newArr = new int[count];
    int j = 0;
    
    
    for(int i = 0;i<size;i++)
    {
        if(arr[i]!=0){

            newArr[j] = arr[i];
            j = j+1;            
           
        }
        
    }
    
    size = size - count;

    for(int i = 0; i<size; i++)
    {
        //cout<< newArr[i]<<endl;
        arr[i] = newArr[i];

    }


}


int main(int argc, char* argv[])
{
    int size = 5;
    int* arr = new int[size];
    arr[0] = 10;
    arr[1] = 2;
    arr[2] = 3;
    arr[3] = 15;
    arr[4] = 6;
    cout<<"before calling the function"<<endl;
    cout<<"size:"<<size<<endl;
    cout<<"arr:";
    for(int i=0;i<size;i++)
    {
        cout<<arr[i]<<" ";
    }
    cout<<endl;

    removeDivByFive(arr,size);
    cout<<"after calling the function"<<endl;
    cout<<"size:"<<size<<endl;
    cout<<"arr:";
    for(int i=0;i<size;i++)
    {
        cout<<arr[i]<<" ";
    }
    cout<<endl;
}