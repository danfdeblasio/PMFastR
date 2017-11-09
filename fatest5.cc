#include <stdio.h>
#include <iostream>
#include <string.h>
#include "stockholmbuffer.h"
#include "filebuffer.h"
#include "fastalignment.5.h"
#include "stack.h"
#include <unistd.h>
#include <cstdlib>
using namespace std;

int main(int argc, char* argv[])      {
int total_nuc = 16;
char na_name[17] = {'a', 'c', 'g', 't', 'n', 'r', 'y', 'w', 's', 'm',
	                       'k', 'h', 'b', 'v', 'd','x'};
                                                                                                                                                             
        int inpseq = 0;
        if (argc > 4 || argc < 3)       {
                cout << "Please input like: fastr stk_file fasta_file\n";
                exit(0);
        }
        cout << argv[0] << " " << argv[1] << " " << argv[2] << endl;
        char *genome, *query, *filterfile;
        genome = argv[2];
        query = argv[1];
		//cout << argv[0] << " " << argv[1] <<  " " << argv[2] << endl;
        StockholmBuffer *fb1 = new StockholmBuffer( query,  20000000);
//      fb1->Test_buffer();
//	StockholmBuffer *fb2 = new StockholmBuffer( "rtemp4.stk", 20000000);
//	fb2->Test_buffer();
	FileBuffer *fb2 = new FileBuffer( genome, 2000000);
	//cout << "name2 " << fb2->name;
	//for(int i=0;i<fb1->length;i++) cout << (int)fb1->structure[i];
//	fb2->Test_buffer();
	Alignment *ali;
	if(argc==4) ali = new Alignment(fb1->seq, fb1->structure, fb2->seq, fb1->length, fb2->length-1, fb1->name, fb2->name,argv[3]);
	else ali = new Alignment(fb1->seq, fb1->structure, fb2->seq, fb1->length, fb2->length-1, fb1->name, fb2->name,(char*)"test.stok");
	//cout << getpid() << endl;
	char* c = new char[100];
	sprintf(c,"more /proc/%d/status > %s-out",getpid(),argv[0]);
  system(c);
	delete fb1;
	delete fb2;
	delete ali;
        return 0;
}


