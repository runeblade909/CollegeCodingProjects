#include <iostream>

using namespace std;

#define SIZE 32

class Stack{
	private:
		int capacity;
		int top;
		char* arr;
	public:
		// Constructor to initialize stack
		Stack() {
			capacity = SIZE;
			top = 0;
			arr = new char[capacity];
		}
		// Destructor to free memory
		~Stack() {
			delete arr;
		}
		void push(char c) {
			if (top < capacity) 
				arr[top++] = c;
		}
		char pop() {
			if (top > 0) 
				return arr[--top];
			else
				return '\0';
		}
};


/*
	ToDo: implement the printMessage method to print hidden messages in our	strings. 
	 - Hidden messages are enclosed in parenthesis 
	 - If no such message can be found, do not print anything
	 - Messages must be output in a readable format (no reverse strings)
*/
void printSecretMessage(char str[]) {
	/* 
		Your code here!
	*/
        int start;
        int end;
        int checker;
        int checker2;
        for(int i = 0; i< SIZE ; i++){

            if(str[i] =='('){
                start=i;
                checker = 50;
            }

            if(str[i] ==')'){
                end=i;
                checker2 = 100;
            }

        }
        int length = 0;

        Stack word;

        length = end-start;

        //cout << length << ' '<<start << ' ' << end;

            if(checker2 == 100 && checker == 50)
            {

                if (length == end){
                    length = length;


                }

                int count = 0;
                for(int i = end-1; count < length-1 ; i-- ){

                    word.push(str[i]);
                    
                    cout << word.pop();
                    count ++;
                    

                }



            }

}

int main() {
	// Test Cases
	int testCases = 7;
	char string[testCases][SIZE] = \
		{
        {"This is a (terces) message"}, \
		{"This is (sdrowkcab)"}, \
		{"I (evol pasta"}, \
		{"This is sdrowrof)"}, \
		{""}, {"()"}, {"($)"}
        };

	for (int i=1; i<=testCases; i++) {
		cout << "Test " << i << endl;
		cout << "\tStack contents: " << string[i-1] << endl;
		cout << "\tSecret message: "; 
		printSecretMessage(string[i-1]);
		cout << endl;
	}

	return 0;
}