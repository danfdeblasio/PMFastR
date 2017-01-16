//filebuffer.cc 
//to read a sequence from a file.
//later on, need to read a set of sequences from one file.
//

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "stockholmbuffer.h"

using namespace std;

extern int total_nuc;
extern char na_name[17];
StockholmBuffer::StockholmBuffer(char *filename, int maxnum) 	{
	char str[500];
	int n=0;
	int m=0;
	int co,cc;
	co = cc = 0;	
	seq = new char[maxnum];
	structure = new int[maxnum];
	std::ifstream from(filename);
	while (from.getline(str, 450))	{
		if (str[0] == '#') {
			sprintf(name,str);
                        for(int i=1;i<200;i++){
                                name[i-1] = name[i];
                        }	
			continue;
		}
		if (str[0] == '>') 	{
			sprintf(name,str);
			for(int i=1;i<200;i++){
				name[i-1] = name[i];
			}
//			name[5] = '\0';
		} else if (str[0] == 'S'){
			for (int i=1; i<strlen(str); i++)	{
				if (str[i] == '(' || str[i] == '>')	{
					structure[m++] = 1;
					co++;
				} else if (str[i] == '<' || str[i] == ')')	{
					structure[m++] = -1;
					cc++;
				} else if ((str[i] == '.')|| (str[i] >='a' && str[i] <= 'z') || (str[i] >='A' && str[i] <= 'Z') )	{
					structure[m++] = 0;
				}
			}
				
		}
		else {
			for (int i=0; i<strlen(str); i++)	{
				if(str[i]=='U') str[i] = 'T';
				if(str[i]=='u') str[i] = 't';
				if (str[i] >= 'a' && str[i] <= 'z')	{
					seq[n++] = Char2int(str[i]);
				} else if (str[i] >= 'A' && str[i] <= 'Z')  {
					seq[n++] = Char2int(str[i] - 'A' + 'a');
				}
			}
		}
	}
	//cout << co << "<open close>" << cc << endl;
	length = n;
	if (m != m)	{
		cout << "Error: structure  and sequence is not consistant \n";
	}
}

StockholmBuffer::~StockholmBuffer()	{
	delete [] seq;
	delete [] structure;
}

void StockholmBuffer::Test_buffer() 	{
	cout << "file length:" << length << '\n';
	for (int i = 0; i < length; i++)	{
		cout << na_name[seq[i]];
	}
	cout << '\n';
	for (int i = 0; i < length; i++)	{
		cout << structure[i];
	}
	cout << '\n';
}

void StockholmBuffer::Out_put(unsigned long start, int length)	{
	for (int i = 0; i < length; i++)	{
		cout << na_name[seq[i+start]];
		if ( i % 60 == 59)	{
			cout << '\n';
		}
	}
}
	
char StockholmBuffer::Char2int(char c)	{
	int k;
	int i;
	for(i = 0; i < total_nuc; i ++)	{
		if(c == na_name[i])     {
			k = i;
			break;
		}
	}
	if(i == total_nuc)      {
		cout <<"Not found" << c << '\n';
		k = rand() % 4;
	}
       if(k > 3)       {
                k = rand() % 4;
	}

      return(k);
}

