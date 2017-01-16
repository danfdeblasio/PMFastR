/*
 * Queue driver
 */

#include <stdio.h>

#include <iostream>

#include "simple_queue.h"

using namespace std;

void s_queue::enqueue(pthread_t a){

	cout << "enqueue" << endl;
	if(l==NULL) l = new S_LinkedList(a);
	else l->append(new S_LinkedList(a));
}

pthread_t s_queue::dequeue(){
    
    if(!empty()){
        pthread_t ret = l->a;
        l = l->next;
        return ret;
    }
    return NULL;
    
    /*
	pthread_t ret = new pthread_t[2];
	if(empty()){
		pthread_t[0] = pthread_t[1] = NULL; 
	}else{
		pthread_t[0] = l->a;
		l = l->next;
	}
	return ret[0];
	*/
}

pthread_t s_queue::top(){
	if(empty()){ return NULL;}
	return l->a;
}
