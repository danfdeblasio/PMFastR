#include <stdio.h>
#include <iostream>
#include "string"
#include "stockholmbuffer.h"
#include "filebuffer.h"
#include "fastalignment.6.h"
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
                cout << "Please input like: <executable> stok_file fasta_file [output_stok_file]\n";
                exit(0);
        }
        cout << argv[0] << " " << argv[1] << " " << argv[2] << endl;
        char *genome, *query, *filterfile;
        genome = argv[2];
        query = argv[1];
		//cout << argv[0] << " " << argv[1] <<  " " << argv[2] << endl;
        //StockholmBuffer *fb1 = new StockholmBuffer( query,  20000000);
        
        Stk_profile *stkpf = new Stk_profile(query, 40000, 4000);
        //stkpf->Test_buffer();
        for(int i=0;i<stkpf->length;i++){
		//cout << (int)stkpf->seq[i];
	}	
	//cout << endl;
//      fb1->Test_buffer();
//	StockholmBuffer *fb2 = new StockholmBuffer( "rtemp4.stk", 20000000);
//	fb2->Test_bu1->Test_buffer()ffer();

	FileBuffer *fb2 = new FileBuffer( genome, 6000);
	//cout << "Files In Okay" << endl;
	//for(int i=0;i<fb2->length;i++) cout << na_name[(int)fb2->seq[i]];
//	fb2->Test_buffer();
	Alignment *ali;
	if(argc==4) ali = new Alignment(stkpf->seq, stkpf->structure, fb2->seq, stkpf->length, fb2->length-1, stkpf, fb2->name, argv[3]);
	else ali = new Alignment(stkpf->seq, stkpf->structure, fb2->seq, stkpf->length, fb2->length-1, stkpf, fb2->name,(char*) "test.stok");
	//cout << getpid() << endl;
	//char* c = new char[100];
	//sprintf(c,"more /proc/%d/status > %s-out",getpid(),argv[0]);
	//system(c);
	//delete stkpf;
	//delete fb2;
	//delete ali;
        return 0;
}


