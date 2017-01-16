all: rna_compare5 rna_compare6 
rna_compare5: fastalignment.5.cc fastalignment.5.h fatest5.cc filebuffer.cc stack.cc stockholmbuffer.cc  weight3.h
	g++ -w  -Wno-deprecated -o rna_compare5 fastalignment.5.cc filebuffer.cc stack.cc stockholmbuffer.cc fatest5.cc 
rna_compare6: fastalignment.6.cc fastalignment.6.h fatest6.cc filebuffer.cc stack.cc stockholmbuffer.cc  weight3.h stk_profile.cc queue.o
	g++ -w  -Wno-deprecated -o rna_compare6 fastalignment.6.cc filebuffer.cc stack.cc stockholmbuffer.cc fatest6.cc stk_profile.cc queue.o
rna_compare8: fastalignment.8.cc fastalignment.8.h fatest6.cc filebuffer.cc stack.cc stockholmbuffer.cc  weight3.h stk_profile.cc queue.o simple_queue.cc simple_queue.h queue.cc
	g++ -w  -Wno-deprecated -lpthread -o rna_compare8 fastalignment.8.cc filebuffer.cc stack.cc stockholmbuffer.cc fatest6.cc stk_profile.cc simple_queue.cc queue.cc
queue.o: queue.cc queue.h
	g++ -w  -Wno-deprecated -c queue.cc
clean:
	rm rna_compare5 rna_compare6 rna_compare8 queue.o 
