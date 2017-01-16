/*Queue used to store 1 integer

*/


class S_LinkedList{
	public:
		pthread_t a;
		S_LinkedList* next;
		S_LinkedList(){next=NULL;}
		~S_LinkedList(){delete next;}
		S_LinkedList(pthread_t ai){a=ai;next=NULL;}
		void append(S_LinkedList* in){if(next==NULL) next = in;else next->append(in);}
};

class s_queue{
	public:
		s_queue(){l = NULL;}
		~s_queue(){delete l;}
		void enqueue(pthread_t a);
		pthread_t dequeue();
		pthread_t top();
		bool empty(){return l==NULL;}
	private:
		S_LinkedList* l;
};
