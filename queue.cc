/*
 * Queue driver
 */

#include <stdio.h>

#include <iostream>

#include "queue.h"

using namespace std;

void queue::enqueue(int a,int b){
	if(l==NULL) l = new LinkedList(a,b);
	else l->append(new LinkedList(a,b));
	//cout << "enqueue" << a << " " << b << " " << l->a << " " << l->b << endl;
}

int* queue::dequeue(){
	int* ret = new int[2];
	if(empty()){
		ret[0] = ret[1] = -1; 
	}else{
		ret[0] = l->a;
		ret[1] = l->b;
		l = l->next;
	}
	return ret;
}
