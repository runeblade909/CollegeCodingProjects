#include <iostream>
#include <string.h>

using namespace std;

bool isNew(int arr[],int k, int n) //this function is just to check if the new element already exists or not, this is just to format the output
{
bool flag = true;
for (int i=0; i<n; i++) // Here we will Pick all elements one by one
{
if(arr[i]==k) // Check if the picked element is already picked
        {
           flag = false;
break;
       }       
}
return flag;
}

bool isPrime(int n) //in this function we will check if a particular number is prime or not, this will come handy when we iterate through the array
{  
   int i=2;
   bool flag = true;
   if(n==1) //since one is neither prime nor composite we will take it as composite
   {
       return false;
   }
   for(i=2;i<n/2;i++) //loop to check if the number is prime
   {
       if(n%i==0)
           flag = false;
   }
   return flag;
     
}
int findMin(int arr[], int arrL) // this function will help us find the smallest prime number
{   int i;
   int min = arr[0];
   for(i=1;i<arrL;i++)
   {
       if(isPrime(arr[i]) && arr[i]<min) //Notice how we are checking if the the number is prime as well if it's less than the current min
       {
           min = arr[i];
       }
   }
   return min;
}

int findMax(int arr[], int arrL) // this function goes same as above except that it checks for max
{
   int i;
   int max = 0;
   for(i=0;i<arrL;i++)
   {
       if(isPrime(arr[i]) && arr[i]>max)
       {
           max = arr[i];
       }
   }
   return max;
}

void multipleList(int min,int max, int arr[],int arrL) //this is the last function in which we will
{   int i;
   int MinMultiples[arrL];
   int MaxMultiples[arrL];
   int MinCount=0; //variable to count the number of elements in minMultiples
   int MaxCount=0; //variable to count the number of elements in maxMultiples
   for(i=0;i<arrL;i++)
   {
       if((arr[i]%min == 0) && isNew(MinMultiples,arr[i],MinCount)) //now this is where we will filter the elements if they are the multiples of min and if they already exist
       {
           MinMultiples[MinCount++] = arr[i];   
       }
       if((arr[i]%max == 0) && isNew(MaxMultiples,arr[i],MaxCount))
       {
           MaxMultiples[MaxCount++] = arr[i];   
       }
   }
  
   cout<<"\nMinMultiples =";
   for(i=0;i<MinCount;i++) //printing both the arrays
   {
       cout<<" "<<MinMultiples[i];
   }
   cout<<"\nMaxMultiples =";
   for(i=0;i<MaxCount;i++)
   {
       cout<<" "<<MaxMultiples[i];
   }
  
}

main()
{
   int arr[]={8,2,6,7,14,6,7,2,1,9,5,6};
   int arrL=12;
   int min = findMin(arr,arrL);
   int max = findMax(arr,arrL);
   cout<<min<<"\n"<<max;
   multipleList(min,max,arr,arrL);
   cout << endl;
}