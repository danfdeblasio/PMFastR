//filebuffer.cc 
//to read a sequence from a file.
//later on, need to read a set of sequences from one file.
//

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include "na.h"
#include "filebuffer.h"

using namespace std;
FileBuffer::FileBuffer(char *filename, int muxnum) 	{
	char str[500];
	int n=0;
	
	seq = new char[muxnum];
	ifstream from(filename);
	while (from.getline(str, 450))	{
		if (str[0] == '#') 	{
			continue;
		}
		if (str[0] == '>') 	{
			//cout << str << endl;
			sprintf(name,str);
			for(int i=1;i<200;i++){
				name[i-1] = name[i];
			}
//			name[5] = '\0';
		} else {
			//cout << str << endl;
			for (int i=0; i<strlen(str); i++)	{
				if(str[i]=='U') str[i] = 'T';
				if(str[i]=='u') str[i] = 't';
				if (str[i] >= 'a' && str[i] <= 'z')	{
					//cout << "n=" << n  << " " << str[i] << endl;
					seq[n] = Char2int(str[i]);
					n++;
				} else if (str[i] >= 'A' && str[i] <= 'Z')  {
					//cout << "n=" << n  << " " << str[i] << endl;
					seq[n] = Char2int(str[i] - 'A' + 'a');
					n++;
				}//else cout << "!" << str[i] << n;
				//cout << str[i] << "=(" << n << ")" << na_name[(int)seq[n-1]] << "  ";
			}
		}
	}
	length = n;
	
/*	seq2 = new char[length+1];
	n = 0;
	int j = 0;
	for (int i = length-1; i>=0; i--)	{
		switch (seq[i])		{
			case 0: 
				seq2[n] = 3;
				break;

			case 1:
				j++;
				seq2[n] = 2;
				break;
			case 2:
				j++;
				seq2[n] =1;
				break;
			case 3:
				seq2[n] = 0;
				break;
		}
		n++;
	}
//	cout << "GC " << j << '\n'; 
*/
	//from.close();
	//cout << "Done With FASTA in" <<  endl;

}

FileBuffer::~FileBuffer()	{
	delete [] seq;
	delete [] seq2;
}

void FileBuffer::Test_buffer() 	{
	
	cout << "file length:" << length << '\n';
	for (int i = 0; i < length; i++)	{
		cout << na_name[seq[i]];
	}
	cout << '\n';

/*
	cout << "file length:" << length << '\n';
        for (int i = 0; i < length; i++)        {
                cout << na_name[seq2[i]];
		if ( i % 60 == 59)	{
			cout << '\n';
		}
        }
        cout << '\n';
*/
}

void FileBuffer::Out_put(unsigned long start, int length)	{
	for (int i = 0; i < length; i++)	{
		cout << na_name[seq[i+start]];
		if ( i % 60 == 59)	{
			cout << '\n';
		}
	}
}
	
char FileBuffer::Char2int(char c)	{
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


