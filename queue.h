/*Queue used to store 2 integers (as an array)

*/


class LinkedList{
	public:
		int a,b;
		LinkedList* next;
		LinkedList(){next=NULL;}
		~LinkedList(){delete next;}
		LinkedList(int ai,int bi){a=ai;b=bi;next=NULL;}
		void append(LinkedList* in){if(next==NULL) next = in;else next->append(in);}
};

class queue{
	public:
		queue(){l = NULL;}
		~queue(){delete l;}
		void enqueue(int a,int b);
		int* dequeue();
		bool empty(){return l==NULL;}
	private:
		LinkedList* l;
};
