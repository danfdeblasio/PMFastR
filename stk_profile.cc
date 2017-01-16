//filebuffer.cc 
//to read a sequence from a file.
//later on, need to read a set of sequences from one file.
//

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "stack.h"
//#include "na.h"
#include "stk_profile.h"
//#include "na.h"

using namespace std;

extern int total_nuc;
extern char na_name[17];

Stk_profile::Stk_profile(char *filename, int maxnum, int maxseq,bool doShrink) 	{
	//char str[700];
	int n=0;
	int m=0;
	m_numseq = maxseq;
	char temp[100], temp2[100];
	char main_str[600];	
	int i,j;
	delq = new queue();
	aln = new char*[maxseq];
	name = new char*[maxseq];
	for (i = 0; i < maxseq ; i++)	{
		aln[i] = new char [maxnum];
		name[i] = new char [100];
	}
	structure = new int[maxnum];
	std::ifstream from(filename);
	while (!from.eof()){
		char str[500];
		from >> str;
		if ((str[0] >= 'a' && str[0] <= 'z') || ( str[0] >= 'A' && str[0] <= 'Z') && !from.eof())	{
			i = 0;
			j = n;
			while ((str[0] >= 'a' && str[0] <= 'z') || ( str[0] >= 'A' && str[0] <= 'Z') && !from.eof())	{
				n = j;
				//sscanf(str, "%s %s", name[i], main_str);
				from >> main_str;
				sprintf(name[i],str);
				//cout << "str in-->" << main_str << endl;
				for (int k=0; k< strlen(main_str); k++)	{
					if(main_str[k] == 'u' || main_str[k] == 'U') main_str[k] = 'T';
					if (main_str[k] >= 'a' && main_str[k] <= 'z')	{
						aln[i][n++] = Char2int(main_str[k]);
					} else if (main_str[k] >= 'A' && main_str[k] <= 'Z')  {
						aln[i][n++] = Char2int(main_str[k] - 'A' + 'a');
						//cout << (char)(main_str[k] - 'A' + 'a') << (int)Char2int(main_str[k] - 'A' + 'a') << (int)aln[i][n-1] << endl;
					} else if (main_str[k] != ' ')	{
						aln[i][n++] = 4;
					}else{ cout << i << " " << n << " " << main_str[k] << endl;}
					//cout << main_str[k];
				}
//				cout << endl;
				i++;
				from >> str;
				//cout << "main-->" <<  main_str << endl;
			}
			numseq = i;
		} 
  
		if (str[0] == '#' && str[1] == '=' && str[2] == 'G'){
                        from >> temp >> main_str;
                        for (int i=0; i<strlen(main_str); i++)  {
                                if (main_str[i] == '<') {
                                        structure[m++] = 1;
                                } else if (main_str[i] == '>')  {
                                        structure[m++] = -1;
                                } else if (main_str[i] == '.')  {
                                        structure[m++] = 0;
                                }
                        }
                        //cout << "main2-->" << main_str << endl;                                                                                                                                     
                }
               //cout << "ERROR!! " << str << endl; 
        sprintf(str,"");
		//cout << "str-->" << str << endl;
		//cout << m << ' ' << n << '\n';
	
	}
	length = n;
	//cout << m << ' ' << n << from.eof() << '\n';
	if (m != n)	{
		cout << "Error: structure  and sequence is not consistant (profile) \n";
		for(int q=0;q<n;q++){
			cout << na_name[aln[0][q]];
		}
		cout << endl;
		for(int q=0;q<m;q++){
			cout << (structure[q]==0)?".":((structure[q]>0)?"<":">");
		}
		cout << endl;
	}
	from.close();
	
	aln_old = new char*[numseq];
	for (i =0; i < numseq; i++)	{
		aln_old[i] = new char[length];
	}
	structure_old = new int[length];
	/*
	for(int i=0;i<length;i++){
		cout << (int)aln[0][i];
	}
	cout << endl;
*/
	//toFile((char*)"before.stok");
	Compact();
	//toFile((char*)"after.stok");
	/*
	for(int i=0;i<length;i++){
		cout << (int)aln[0][i];
	}
	cout << endl;
	
	*/
	gap_cut = .9;
	//cout << "before: " << length;
	if(doShrink) Shrink();
	//cout << " after: " << length << endl;
	profile = new double*[5]; 
	prof_matrix = new int*[5];
	for (i = 0; i < 5; i++)	{
		profile[i] = new double [length];
		prof_matrix[i] = new int [length];
	}
	p_profile = new double*[17];
	prof_pmatrix = new int*[17];
	for (i = 0; i < 17; i++)	{
		p_profile[i] = new double[length];
		prof_pmatrix[i]= new int[length];
	}
	seq = new char[length];
	struc1 = new BPLIST [length];
	Get_structure();
	Build_profile();
	
	int max = 0;
	for(int i=0;i<numseq;i++){
		int k=0;
		for(int j=0;j<100 && name[i][j]!='\0';j++){
			k++;
		}
		if(k>max) max = k;
	}
	max += 4;
	for(int i=0;i<numseq;i++){
		int j;
		for(j=0;j<100 && name[i][j]!='\0';j++){
		}
		for(;j<=max;j++){
			name[i][j] = ' ';	
		}
		name[i][j] = '\0';
	}
	maxName = max;
	
	/*for(int i=0;i<length;i++){
		cout << (int)aln[0][i];
	}
	cout << endl;
*/
	//cout << "got here " << endl;
}


Stk_profile::~Stk_profile()	{
	/*for (int i= 0; i < m_numseq; i++)      {
                delete [] aln[i];
                delete [] name[i];
        }
	
	for (int i=0; i < numseq; i++)	{
		delete [] aln_old[i];
	}
	delete [] aln_old;
	delete [] structure_old;

	for (int i= 0; i < 5; i++)      {
                delete [] profile[i];
		delete [] prof_matrix[i];
        }
        for (int i= 0; i < 17; i++)      {
                delete [] p_profile[i];
                delete [] prof_pmatrix[i];
        }


	delete [] aln;
	delete [] name;
	delete [] profile;
	delete [] prof_matrix;
	delete [] p_profile;
	delete [] prof_pmatrix;
	delete [] structure;
	cout << "c" << endl;
	delete [] seq;
	delete [] struc1;
	cout << "d" << endl;*/
}

void Stk_profile::Test_buffer() 	{
	cout << "file length:" << length << '\n';
	cout << "number of sequence: " <<numseq << '\n';
	for (int i = 0; i < numseq; i++)	{
		for (int j = 0; j < length; j++)	{
			if (aln[i][j] == 4)	{
				cout << '-';
			} else {
				cout << na_name[aln[i][j]];
				}
		}
		cout << '\n';
	}
	
	cout << "consensus:\n";                                                                                                                                                             
        for (int i = 0; i < length; i++)        {
                cout << na_name[seq[i]];
        }
	cout << '\n';
	
	for (int i = 0; i < length; i++)	{
		char sTemp = ' ';
		if( structure[i] == 1 ) sTemp = '<';
		if( structure[i] == -1 ) sTemp = '>';
		if( structure[i] == 0 ) sTemp = ' ';	
		cout << sTemp;
	}
	cout << '\n';

/*	for (int i = 0; i < 5; i++)	{
		for (int j = 0; j < 10; j++)	{
			cout << prof_pmatrix[i][j] <<  ' ';
		}
		cout << '\n';
	}
                                                                                                                                                             
        for (int i = 0; i < 5; i++)     {
                for (int j = 0; j < 10; j++)    {
                        cout << p_profile[i][j] <<  ' ';
                }
                cout << '\n';
        }
*/
}

void Stk_profile::Out_put(unsigned long start, int length)	{
	for (int i = 0; i < length; i++)	{
		cout << na_name[aln[0][i+start]];
		if ( i % 60 == 59)	{
			cout << '\n';
		}
	}
}
	
char Stk_profile::Char2int(char c)	{
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

void Stk_profile::Build_profile()	{
	int i; 
	int j;
	int k;
int weight[5][5] = {     22, -19, -15, -14, -20,
                        -19,  12, -25, -10, -20,
                        -15, -25,  10, -17, -20,
                        -14, -10, -17, 17, -20,
			-20, -20, -20, -20, 0};
                                                                                                                                                             
int  p_w[16][16] = {
                                     
	  /*aa   ac   ag   at   ca   cc   cg   ct   ga   gc   gg   gt   ta   tc   tg   tt*/
/*aa*/ -25, -70, -82, -43, -88,-143, -47,-126, -69, -50, -84, -58, -40, 113, -62, -91,
/*ac*/ -70, -21, -89, -20, -94, -91, -59,-105, -97, -38,-111, -47, -53, -86, -69, -78,
/*ag*/ -82, -89,  -8, -51,-104,-145, -46,-101, -86, -58, -54, -66, -54, -89, -59,-111,
/*at*/ -43, -20, -51,  45, -56, -67,  17, -52, -53,  27, -56,   6,  16, -48,  -5, -30,
/*ca*/ -88, -94,-104, -56, -51,-105, -36, -85, -80, -60,-114, -79, -24, -71, -56, -84,
/*cc*/-143, -91,-145, -67, -105,-36, -57, -58,-124, -37,-126, -79, -69, -74, -84, -54,
/*cg*/ -47, -59, -46,  17, -36, -57,  54, -50, -60,  21, -46,  -3,  28, -49,  13, -37,
/*ct*/-126,-105,-101, -52, -85, -58, -50, -23, -77, -58,-137, -56, -47, -38, -74, -52,
/*ga*/ -69, -97, -86, -53, -80,-124, -60, -77, -11, -49, -87, -61, -59, -66, -76,-115,
/*gc*/ -50, -38, -58,  27, -60, -37,  21, -58, -49,  56, -41,  12,  16, -45,  -8, -39,
/*gg*/ -84,-111, -54, -56,-114,-126, -46,-137, -87, -41, -20, -58, -58,-120, -43,-108,
/*gt*/ -58, -47, -66,   6, -79, -79,  -3, -56, -61,  12, -58,  35,  -6, -53, -21, -45,
/*ta*/ -40, -53, -54,  16, -24, -69,  28, -47, -59,  16, -58,  -6,  50, -30,  11, -34,
/*tc*/ 113, -86, -89, -48, -71, -75, -49, -38, -66, -45,-120, -53, -30, -32, -48, -60,
/*tg*/ -62, -69, -59,  -5, -56, -84,  13, -74, -76,  -8, -43, -21,  11, -48,  34, -43,
/*tt*/ -91, -78,-111, -30, -84, -54, -37, -52,-115, -39,-108, -45, -34, -60, -43,   0
/*                                                                                                                        
 -25, -70, -82, -43, -88,-143, -47,-126, -69, -50, -84, -58, -40, 113, -62, -91,
 -70, -21, -89, -20, -94, -91, -59,-105, -97, -38,-111, -47, -53, -86, -69, -78,
 -82, -89,  -8, -51,-104,-145, -46,-101, -86, -58, -54, -66, -54, -89, -59,-111,
 -43, -20, -51,  45, -56, -67,  17, -52, -53,  27, -56,   6,  16, -48,  -5, -30,
 -88, -94,-104, -56, -51,-105, -36, -85, -80, -60,-114, -79, -24, -71, -56, -84,
-143, -91,-145, -67, -105,-36, -57, -58,-124, -37,-126, -79, -69, -74, -84, -54,
 -47, -59, -46,  17, -36, -57,  54, -50, -60,  21, -46,  -3,  28, -49,  13, -37,
-126,-105,-101, -52, -85, -58, -50, -23, -77, -58,-137, -56, -47, -38, -74, -52,
 -69, -97, -86, -53, -80,-124, -60, -77, -11, -49, -87, -61, -59, -66, -76,-115,
 -50, -38, -58,  27, -60, -37,  21, -58, -49,  56, -41,  12,  16, -45,  -8, -39,
 -84,-111, -54, -56,-114,-126, -46,-137, -87, -41, -20, -58, -58,-120, -43,-108,
 -58, -47, -66,   6, -79, -79,  -3, -56, -61,  12, -58,  35,  -6, -53, -21, -45,
 -40, -53, -54,  16, -24, -69,  28, -47, -59,  16, -58,  -6,  50, -30,  11, -34,
 113, -86, -89, -48, -71, -75, -49, -38, -66, -45,-120, -53, -30, -32, -48, -60,
 -62, -69, -59,  -5, -56, -84,  13, -74, -76,  -8, -43, -21,  11, -48,  34, -43,
 -91, -78,-111, -30, -84, -54, -37, -52,-115, -39,-108, -45, -34, -60, -43,   0*/
};

int pa[16][2] = {0, 0, 0, 1, 0, 2, 0, 3,
                 1, 0, 1, 1, 1, 2, 1, 3,
                 2, 0, 2, 1, 2, 2, 2, 3,
                 3, 0, 3, 1, 3, 2, 3, 3};

	for (i = 0; i < 5 ; i++)	{
		for (j = 0; j < length; j ++)	{
			profile[i][j] = 0.01;
			prof_matrix[i][j] = 0;
		}
	}
        for (i = 0; i < 17 ; i++)        {
                for (j = 0; j < length; j ++)   {
                        p_profile[i][j] = 0.01;
                        prof_pmatrix[i][j] = 0;
                }
        }

 	for (i = 0; i < numseq; i++)	{
		for (j=0; j < length; j++)	{
			profile[aln[i][j]][j]++;
			if (struc1[j].pos >=0)	{
				int temp, t1, t2;
				t1= aln[i][j];
				t2 = aln[i][struc1[j].pos];
				if (t1 < 4 && t2 < 4)	{
					p_profile[t1*4+t2][j]++;
				} else {
					p_profile[16][j]++;
				}
			}
		}
	}
	
	for (i = 0; i < 5 ; i++)        {
                for (j = 0; j < length; j ++)   {
                        profile[i][j] = profile[i][j]/(numseq+0.04);
             }
        }

        for (i = 0; i < 17 ; i++)        {
                for (j = 0; j < length; j ++)   {
                        p_profile[i][j] = p_profile[i][j]/(numseq+0.04);
              }
        }

	for (j = 0; j < length; j++)	{
		double tt = 0;
		int ttp = -1;
		for (i =0; i< 5; i++)	{
			double temp = 0;
			for (k=0; k<5; k++)	{
				temp = temp + profile[k][j]*weight[i][k];
//				prof_matrix[i][j] = prof_matrix[i][j] + (int) profile[k][j]*weight[i][k];
			}
			prof_matrix[i][j] = (int) temp;
//			cout << temp << ' ' << prof_matrix[i][j] << '\n';
		}
		for (i=0; i<4; i++)	{
			if (profile[i][j] >= tt)	{
				ttp = i;
				tt = profile[i][j];
			}
		}
		seq[j] = ttp;
	}

        for (j = 0; j < length; j++)    {
		if (struc1[j].pos >=0 )	{	
			double tt = 0;
			int ttp = -1;
	                for (i =0; i< 16; i++)   {
				double temp = 0;
	                        for (k=0; k<16 ; k++)     {
					  temp = temp + p_profile[k][j]*p_w[i][k];
//	                                  prof_pmatrix[i][j] = prof_pmatrix[i][j] + (int) p_profile[k][j]*p_w[i][k];
          	              }
//		 		prof_pmatrix[i][j] = prof_pmatrix[i][j] + (int) p_profile[16][j]*(-20);
				temp = temp + p_profile[16][j]*(-20);
				prof_pmatrix[i][j] = (int) temp;
				if (p_profile[i][j] >= tt)	{
					ttp = i;
					tt = p_profile[i][j];
				} 
               		}
			prof_pmatrix[16][j] = -10;
			seq[j] = pa[ttp][0];
			seq[struc1[j].pos] = pa[ttp][1];
		}
        }
				

/*	double T = 0;
	double current = 0;
        for (i = 0; i < numseq; i++)    {
                for (j=0; j < length; j++)      {
                        T= profile[aln[i][j]][j]+T;
                }
		if (current ==0 )	{
		 	current = T;
		} else if (T < current)	{
			current = T;
		}
        }
	ave_profile = current;
	cout << ave_profile << " ave_profile T\n";
*/


}

void Stk_profile::Shrink()	{
	int i,j;
	//cout << "iS" << endl;	
	for (i=0; i<numseq; i++)	{
		for (j=0; j<length; j++)	{
			aln_old[i][j] = aln[i][j];
		}
	}
	
	for (i=0; i<length; i++)	{
		structure_old[i] = structure[i];
	}
	length_old = length;

	int c, index;
	index =0;
	for (i=0; i<length_old; i++)	{
		c = 0;
		for (j=0; j < numseq; j++)	{
			if (aln_old[j][i] == 4)	{
				c++;
			}
			
		}
		double count = ((double) c)/numseq;
	//	cout << (count <=  gap_cut) << " " << count << " " << c << " " << numseq << endl;
		if ( count <=  gap_cut || structure_old[i] != 0)	{
			for (j=0; j < numseq; j++)	{
				aln[j][index] = aln_old[j][i];
			}
			structure[index] = structure_old[i];
			index++;
		}else{
			//save deleted info
			//Dan 7/27/08
			delq->enqueue(index,i);
		}
	}
	length = index;
}
	
void Stk_profile::Get_structure() {

       int i,j,temp;
       Stack *stack = new Stack(length);
                                                                                                                                                             
        j = 0;
        for (i = 0; i < length; i++)  {
                if (structure[i] == 0)        {
                        struc1[i].pos = -1;
                        struc1[i].offset = -1;
                } else if (structure[i] == 1) {
                                                                                                                                                             
                        stack->Push(i);
                } else if (structure[i] == -1)        {
                        struc1[i].pos = -1;
                        struc1[i].offset = -1;
                        temp = stack->Pop();
                        struc1[temp].pos = i;
                        struc1[temp].offset = j;
                        j++;
                }
        }
                                                                                                                                                             
        if (!stack->Empty())    {
                cout << "stack is not empty!\n";
        }
                                                                                                                                                             
        bp_num = j;
                                                                                                                                                             
        //delete stack;
                                                                                                                                                             
//      for (i=1; i <=length1; i++)     {
//              cout << stru1[i].pos << '\n';
//      }
//      cout << bp_num << '\n';
                                                                                                                                                             
}

void Stk_profile::Compact(){
	for(int i=0;i<length;i++){
		int ct = 0;
		int last = 0;
		for (int j=0; j < numseq; j++)	{
			if(aln[j][i]!=4){
				ct++;
				last = j;
			}
		}
		
		if(ct==1 && structure[i]==0){
			//cout << ct << " <<>> " << last << " " << i << endl;
			if(aln[last][i-1]==4){
				int k=i-1;
				while(aln[last][k]==4 && structure[k]==0) k--;
				for(int j=k+1;j<i;j++){
					//cout << "in" << j << (int)seq[j] << endl;
					if(/*seq[j]==aln[last][i] && */ structure[j]==0){
						aln[last][j] = aln[last][i];
						RemoveCol(i);
						j = i;
						i--;
					}
				}
			}
			//cout << ct << " >><< " << last << endl;
		}
	}
}

void Stk_profile::RemoveCol(int col){
	length--;
	for(int i=col;i<length;i++){
		for (int j=0; j < numseq; j++)	{
			aln[j][i] = aln[j][i+1];
		}
		structure[i] = structure[i+1];
	}
}

void sToCaps(char* in){
	for(int i=0;i<101;i++){
		if(in[i]>='a' && in[i]<='z') in[i] = in[i]-'a'+'A';
		if(in[i] == 'T') in[i] = 'U';
	}
}

void Stk_profile::toFile(char* fname){
	std::ofstream fout(fname);
	static char aline[101];
	char* glines[numseq];
	for(int t = 0;t<numseq;t++){
		glines[t] = new char[101];
	}
	
	char** name2 = new char*[numseq];
	int j;
	for(int i=0;i<numseq;i++){
		name2[i] = new char[maxName+6];
		sprintf(name2[i],name[i]);
		for(j=0;j<maxName && name2[i][j]!='\0';j++){
		}
		for(;j<maxName+5;j++){
			name2[i][j] = ' ';	
		}
		name2[i][j] = '\0';
		//cout << "bline" << endl;
	}
	char* gName = new char[maxName+6];
	sprintf(gName,"#=GC SS_cons");
	for(j=0;j<maxName+5 && gName[j]!='\0';j++){
	}
	for(;j<maxName+5;j++){
		gName[j] = ' ';	
	}
	gName[j] = '\0';
	
	
	int i=0;
	while(i<length){
		for(j=0;j<50 && i<length;j++){
			for(int k=0;k<numseq;k++){
				glines[k][j] = na_name[aln[k][i]];
				if(glines[k][j] == 'n') glines[k][j] = '-';
			}
			aline[j] = (structure[i]==0)?'.':(structure[i]>0)?'<':'>';
			i++;
		}
		
		for(int j=0;j<numseq;j++){
			sToCaps(glines[j]);
			fout << name2[j] << "  " <<  glines[j] << endl;
			for(int i=0;i<101;i++) glines[j][i] = ' ';
		}
		fout << gName << "  " << aline << endl << endl;
		for(int i=0;i<101;i++) aline[i] = ' ';
	}
}
