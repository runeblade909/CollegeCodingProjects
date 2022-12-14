/****************************************************************/
/*                CountryNetwork Implementation                 */
/****************************************************************/
/* TODO: Implement the member functions of class CountryNetwork */
/*     This class uses a linked-list of Country structs to      */
/*     represent communication paths between nations             */
/****************************************************************/

#include "CountryNetwork.hpp"

using namespace std;

/*
 * Purpose: Constructer for empty linked list
 * @param none
 * @return none
 */
CountryNetwork::CountryNetwork()
{
     head = NULL;
}


/*
 * Purpose: Check if list is empty
 * @return true if empty; else false
 */
bool CountryNetwork::isEmpty()
{
       return (head == NULL);
}


/*
 * Purpose: Add a new Country to the network
 *   between the Country *previous and the Country that follows it in the network.
 * @param previous name of the Country that comes before the new Country
 * @param countryName name of the new Country
 * @return none
 */
void CountryNetwork::insertCountry(Country* previous, string countryName) 
{    
    Country *nn = new Country();
    nn->name = countryName;

    if(previous == NULL) {
        nn->next = head;
        head = nn;
        cout << "adding: " << countryName << " (HEAD) " << endl;
    } else {
        nn->next = previous->next;
        previous->next = nn;
        cout << "adding: " << countryName << " (prev: " << previous->name << ")" << endl;
    }

}


/*
 * Purpose: delete the country in the network with the specified name.
 * @param countryName name of the country to delete in the network
 * @return none
 */
void CountryNetwork::deleteCountry(string countryName) 
{
Country *start = head;
    if(start && start->name.compare(countryName) == 0) {
        //cout << countryName << " deleted." << endl;
        head = head->next;
        return;
    }

    while(start->next != NULL && start->next->name.compare(countryName) != 0) {
        start = start->next;
    }
    if(start->next != NULL) {
        Country *del = start->next;
        //cout << countryName << " deleted." << endl;
        start->next = del->next;
        delete del;
        return;
    } else {
        cout << "Country does not exist." << endl;
    }
}

/*
 * Purpose: populates the network with the predetermined countries
 * @param none
 * @return none
 */

    // Needs to be fixed

void CountryNetwork::loadDefaultSetup() {

    deleteEntireNetwork();
    insertCountry(NULL, "United States");
    insertCountry(NULL, "Canada");
    insertCountry(NULL, "Brazil");
    insertCountry(NULL, "India");
    insertCountry(NULL, "China");
    insertCountry(NULL, "Australia");
    
    
}

/*
 * Purpose: Search the network for the specified country and return a pointer to that node
 * @param countryName name of the country to look for in network
 * @return pointer to node of countryName, or NULL if not found
 * @see insertCountry, deletecountry
 */
Country* CountryNetwork::searchNetwork(string countryName)
{
    Country *start = head;

    while(start != NULL && start->name.compare(countryName) != 0) {
        start = start->next;
    }
    if(start != NULL) {
        return start;
    } else {
        return NULL;
    }
}

/*
 * Purpose: deletes all countries in the network starting at the head country.
 * @param none
 * @return none
 */
void CountryNetwork::deleteEntireNetwork()
{
    Country *start = head;

    while(start != NULL) {
        Country *del = start;
        start = start->next;
        cout << "deleting: " << del->name << endl;
        delete del;
    }

    cout << "Deleted network" << endl;
    head = NULL;
}

/*
 * Purpose: Transmit a message across the network to the
 *   receiver. Msg should be stored in each country it arrives
 *   at, and should increment that country's count.
 * @param receiver name of the country to receive the message
 * @param message the message to send to the receiver
 * @return none
 */
void CountryNetwork::transmitMsg(string receiver, string message)
{
  if(isEmpty()) {
        cout << "Empty list" << endl;
        return;
    }
    if(searchNetwork(receiver) == NULL) {
        cout << endl << "Country not found" << endl;
        return;
    }

    Country *start = head;
    while(start != NULL && start->name.compare(receiver) != 0) {
        start->message = message;
        start->numberMessages += 1;
        cout << start->name << " [# messages received: " << start->numberMessages << "] received: " << start->message << endl;
        start = start->next;
    }
    start->message = message;
    start->numberMessages += 1;
    cout << start->name << " [# messages received: " << start->numberMessages << "] received: " << start->message << endl;
}

/*
 * Purpose: prints the current list nicely
 * @param ptr head of list
 */
void CountryNetwork::printPath() 
{
 if(isEmpty()) {
        cout << "== CURRENT PATH ==" << endl << "nothing in path" << endl << "===" << endl;
        return;
    }

    cout << "== CURRENT PATH ==" << endl;
    Country *start = head;

    while(start != NULL) {
        cout << start->name << " -> ";
        start = start->next;
    }
    cout << "NULL" << endl << "===" << endl;
}
