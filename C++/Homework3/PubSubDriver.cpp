/****************************************************************/
/*                Publisher Subscriber Driver File              */
/****************************************************************/
/* TODO: Implement menu options as described in the writeup     */
/****************************************************************/

#include "PubSub.hpp"
#include "PubSub.cpp"
#include <iostream>
// you may include more libraries as needed

using namespace std;

/*
 * Purpose: displays a menu with options
 * @param none
 * @return none
 */

void menu()
{
	cout << "*----------------------------------------*" << endl;
	cout << "Choose an option:" << endl;
    cout << "1. Publisher (Publish items to the queue)" << endl;
	cout << "2. Subscriber (Retrieve items from the queue)" << endl;
	cout << "3. Return the queue size and exit" << endl;
	cout << "*----------------------------------------*" << endl;
}

// Use getline for reading
int main(int argc, char const *argv[])
{
	// Create an object of the class
   PubSub queue;
   // Declare the required variables
   int choice, num;
   string item;
   // Use switch case to execute the condition
   do {
      menu();
      cin >> choice;

      
     switch (choice) {
	  // Option 1:
      case 1: {

         cout << "Enter the number of items to be published: "<<endl;
         cin >> num;
			for (int i = 0; i<num; i++) {
			cout << "Item" << i + 1 << ":" <<endl;
			cin >> item;
			queue.enqueue(item);
			}

         break;
      }

      // Option 2:
      case 2: {
         cout << "Enter the number of items to be retrieved:" << endl;
         cin >> num;
         
		 if (num > queue.queueSize()){
			
			for (int i = 0; i<num-1; i++) {
			cout << "Retrieved: ";
			cout << queue.peek() << endl;
			//cout << "consumed : \n";
			queue.dequeue();
			}
        	cout << "No more items to retrieve from queue" <<endl;
         }

         else{
			for (int i = 0; i<num; i++) {
			cout << "Retrieved: ";
			cout << queue.peek() << endl;
			//cout << "consumed : \n";
			queue.dequeue();
         	}
         }

         break;
      }

      // Option 3
      case 3: {
         cout << "Number of items in the queue:" << queue.queueSize() << " " << endl;
         break;
      }
    }
	// Ends program when size i asked for (exit)
   } while (choice != 3); 
   return 0;
}






