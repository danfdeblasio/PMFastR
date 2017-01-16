//Here is the ran alignment class.
//alingmment[1.. length1, 1.. length2]
//the starting point is from 1!!
//

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include "string"
#include "stack.h"
#include "fastalignment.5.h"
#include "weight3.h"

#define debug 0

using namespace std;

extern int total_nuc;
extern char na_name[17];


Alignment::Alignment(char *s1, int *st1, char *s2, int l1, int l2, char* n1, char* n2, char* stokOut)
{

        stokoutfilename = stokOut;
        stokfilenameset = true;

	int i, j,k;
	
	name1 = new char[50];
	sprintf(name1,n1);
	name2 = new char[50];
	sprintf(name2,n2);

	//cout << "names " << name1 << "," << name2 << " ns" << n1 << "," << n2 << endl; 
	length1 = l1;
	length2 = l2;
	
	cout << length1 << ' ' << length2 << " length\n";	
	seq1 = new char [length1+1];
	seq2 = new char [length2+3];
	ss1 = new int [length1+1];
	stru1 = new BPLIST [length1+1];
	
	for (i =0; i < length1; i++)	{
		seq1[i+1] = s1[i];
		ss1[i+1] = st1[i];
		//cout << (int)st1[i];
	}
	//cout << "struc<< " << endl;
	length2++;
	for (i=0; i < length2; i++) 	{
		seq2[i+1] = s2[i];
		//cout << (int)s2[i] << " " << (int)seq2[i+1]<< " " << i << endl;
	}
	Get_structure();
	tlist1 = new TLIST [2*length1+1];
	Build_tree();

	g = -20;
	g2 = -40;
	stemplus = 30;
	loopstart = -20;
	
	//start alignment now.
	
	band1 = band2 = BAND;
//	if(BAND < fabs(length2-length1)){
//		band1 = band2 = (int)fabs(length2-length1) + BAND;
//	}
	//align = new int [length1,length1,length2,length2];
	//align = new int [length1][length1][length2][length2];
	align = new int** [newbp_num+2];
	trace = new int** [newbp_num+2];
	for (i = 0; i < newbp_num+2; i++)	{
		align[i]= new int* [2*band1+2];
		trace[i]= new int* [2*band1+2];
		for (j = 0; j < 2*band1+2; j++)	{
			align[i][j] = new int [2*band1+2];
			trace[i][j] = new int [2*band1+2];
		}
	}
	
	Compute_align();

	mlist = new int [length1+length2];
	Trace_alignment();
	//Get_PSSM(s1, &(s2[local_s-1]), st1, length1, local_e - local_s + 1, mlist, 1, local_s);
	Display(s1, &(s2[local_s-1]), st1, length1, local_e - local_s + 1, mlist, 1, local_s);
	if(local_s>1){
		for(int i=length1+length2-1;i>=local_s-1;i--){
			mlist[i] = mlist[i-local_s];
		}
		for(int i=0;i<local_s-1;i++){
			mlist[i] = 3;
		}
		trace_n += local_s;
	}
	if(stokfilenameset) toStokFile(stokoutfilename,s1, s2, st1, length1, local_e - local_s + 1, mlist, 1, 1,n2);
	else toStokFile((char*)"test.stok",s1, s2, st1, length1, local_e - local_s + 1, mlist, 1, 1,n2);
	
	cout << "End alignment\n";
//	Display(s1, s2, st1, length1, length2, mlist, 1, 1);
	
}



Alignment::~Alignment()
{
        delete [] seq1;
        delete [] seq2;
        delete [] ss1;
        delete [] stru1;
        delete [] tlist1;
	delete [] mlist;
	cout << "Got Here" << endl;
	for (int i= 0; i < newbp_num; i++)	{
		for (int j=0; j < 2*band1+2; j++)	{
			delete [] align[i][j];
			delete [] trace[i][j];
		}
		//cout << "in Here" << endl;
		delete [] align[i];
		delete [] trace[i];
	}
        delete [] align;
	delete [] trace;
	//delete [] stokoutfilename;
	//delete name1;
	//free(name2);
	//cout << name1 << name2 << "last" << endl;
}

int Alignment::getAlignVal(int v, int ti, int tj){
	
	
	if(ti >= tlist1[v].pos1 - band1 && ti <= tlist1[v].pos1 + band1 && tj >= tlist1[v].pos2 - band1 && tj <= tlist1[v].pos2 + band1){
		int i = ti - tlist1[v].pos1 + band1;
		int j = tj - tlist1[v].pos2 + band1;
		if(debug == 3) cout << "getA" << ti << "," << tj << "--" << i << "," << j << endl;
		return align[v][i][j];
	}	
	
	if(ti == tj + 1){
		if(v==0) return 0;
		//if(align[v][ti][tj] != (tlist1[v].pos2 - tlist1[v].pos1 + 1)*g) cout << align[v][ti][tj] << "!=" << (tlist1[v].pos2 - tlist1[v].pos1 + 1)*g << "(" << v << "," << ti << "," << tj << ")" << endl;
		return (tlist1[v].pos2 - tlist1[v].pos1 + 1)*g;
	}
	if(tj>=ti && v==0){
		return (tj-ti + 1) * g;
	}
	//if(align[v][ti][tj] != -6000) cout << align[v][ti][tj] << "!=" << -6000 << "(" << v << "," << ti << "," << tj << ")" << endl;
	return -6000;

}

void Alignment::setAlignVal(int v, int ti, int tj, int val){
	int i = ti - tlist1[v].pos1 + band1;
	int j = tj - tlist1[v].pos2 + band1;
	if(debug == 3) cout << "setA" << ti << "," << tj << "--" << i << "," << j << endl;
	if(i<0 || j<0 || i>2*band1 || j>2*band1){ ;  }
	else align[v][i][j] = val;
	if(debug == 3) cout << "-done-" << endl;
}

int Alignment::getTraceVal(int v, int ti, int tj){
	
	
	if(ti >= tlist1[v].pos1 - band1 && ti <= tlist1[v].pos1 + band1 && tj >= tlist1[v].pos2 - band1 && tj <= tlist1[v].pos2 + band1){
		int i = ti - tlist1[v].pos1 + band1;
		int j = tj - tlist1[v].pos2 + band1;
		if(debug == 3) cout << "getT" << v << "(" << newbp_num+2 << ")" << "--" << ti << "," << tj << "--" << i << "," << j << endl;
		if(debug == 3) cout << trace[v][i][j] << endl;
		return trace[v][i][j];
	}	
	return 0;

}

void Alignment::setTraceVal(int v, int ti, int tj, int val){
	int i = ti - tlist1[v].pos1 + band1;
	int j = tj - tlist1[v].pos2 + band1;
	if(debug == 3) cout << "setT" << ti << "," << tj << "--" << i << "," << j << endl;
	if(i>=0 && j >= 0 && i<=2*band1 && j<=2*band1) trace[v][i][j] = val;
}

void Alignment::Compute_align()
{
	int i, j, k, l, m, n;
	int temp_a;
	int temp_t = 0;
	int temp_b;
	int pair[4][4]={ 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0};
	
	/*cout << "seq1:" ;
	for(int s1 = 0; s1<length1; s1++){
		cout << na_name[(int)seq1[s1]];	
	}
	cout << endl << "seq2: ";
	for(int s2 = 0; s2<length2; s2++) cout << na_name[(int)seq2[s2]];
	cout << endl << endl;
	*/
	//initialization here                                                                                                                                                           
        /*for (i = 1; i<= newbp_num; i++)   {
                for (j=1; j<=2*band1 +1; j++)      {
                        for (k= 1; k<= 2*band1 +1; k++)  {
                        		//cout << i << " " << j << " " << k << endl;
                                align[i][j][k] = -6000;
                        }
                }
        }
        for (i = 1; i <= 2*band1 +1; i++)  {
                for (k = 0; k <= 2*band1 +1 - i;  k++)     {
                        align[0][i][i+k] = (k+1)*g;
                }
        }
        for (i = 0; i <= 2*band1; i++)  {
                for (j = 0; j < newbp_num; j++) {
                        align[j+1][i+1][i] = (tlist1[j+1].pos2 - tlist1[j+1].pos1 + 1)*g;
                }
        }
        for (i = 0; i < 2*band1 +1; i++)   {
                align[0][i+1][i] = 0;
        }
        
        */
         for (i = 1; i<= length2; i++)   {
                for (j=1; j<=length2; j++)      {
                        for (k= 1; k<= newbp_num; k++)  {
                        	
                				if(debug  == 1) cout << k << "<k i>" << i << " " << j << "<j" << endl;
                                setAlignVal(k,i,j,-6000);
                        }
                }
        }
        for (i = 1; i <= length2; i++)  {
                for (k = 0; k <= length2 - i;  k++)     {
                	
                		if(debug == 1) cout << k << "<k i>" << i << endl;
                        setAlignVal(0,i,i+k,(k+1)*g);
                }
        }
        for (i = 0; i <= length2; i++)  {
                for (j = 0; j < newbp_num; j++) {
                        //cout << i << "<i j>" << j << endl;
                        setAlignVal(j+1,i+1,i,(tlist1[j+1].pos2 - tlist1[j+1].pos1 + 1)*g);
                }
        }
        for (i = 0; i < length2; i++)   {
                setAlignVal(0,i+1,i,0);
        }
               
		/*for(int i=0;i<newbp_num+2;i++){
		for(int j=0;j<length2 + 2;j++){ 
		for(int k=0;k<length2 + 2;k++){ 

			if(falses[i][j][k] != -99999999) cout << "WRONG!!!" << i << "," << j << "," << k << " " << getAlignVal(i,j,k) << "!=" << falses[i][j][k] << endl;
		}
		}
		} */                                                                                                                                  
//initialization here
/*	for (i = 1; i <= length2; i++)	{
		for (k = 0; k <= length2 - i;  k++)	{
			align[i][i+k][0] = (k+1)*g;
		}
	}
	for (i = 0; i <= length2; i++)	{
		for (j = 0; j < newbp_num; j++)	{
			align[i+1][i][j+1] = (tlist1[j+1].pos2 - tlist1[j+1].pos1 + 1)*g;
		}
	}
	for (i = 0; i < length2; i++)	{
		align[i+1][i][0] = 0;
	}
*/		
//here is the new alignment:
	
for (j = 1; j <= newbp_num; j++)	{
	int tl = tlist1[j].left;
	int tr = tlist1[j].right;
	int tp1 = tlist1[j].pos1;
	int tp2 = tlist1[j].pos2;
	for (l = 0; l <= length2 - 1; l++)	{
		for (i = 1; i <= length2 - l; i++)	{
				
				temp_a = - 10000;
				

				if ( tp1 >= i - band1 && tp1 <= i + band1 && tp2 <= i+l + band1 && tp2 >= i+l - band1)  {
                                        //banded here.
                                        //if want non-banded alignment, please comment off the above line.

				if (tlist1[j].flag == 1)	{
					//the solid egde
					//1
					if (l != 0) {
						temp_a = getAlignVal(tl,i+1,i+l-1) +/* weight[seq2[i]][seq1[tp1]] 
							+ weight[seq2[i+l]][seq1[tp2]] */
							+ p_w[seq2[i]][seq2[i+l]][seq1[tp1]][seq1[tp2]];
						temp_t = 1;
						if (getTraceVal(tl,i+1,i+l-1) == 1)	{
						       temp_a = temp_a + stemplus;
						}	       
					}
					//2
					temp_b = getAlignVal(tl,i,i+l) + g;
					if (temp_b > temp_a)	{
						temp_a = temp_b;
						temp_t = 2;
					}
					//3
					int ttg = g;
					if (tlist1[j].father != -1)	{
						ttg = g2;
					}
					
					temp_b = getAlignVal(j,i+1,i+l) + ttg;
					if (temp_b > temp_a)	{
						temp_a = temp_b;
						temp_t = 3;
					}
					//4
					temp_b = getAlignVal(j,i,i+l-1) + ttg;
					if (temp_b > temp_a)	{
						temp_a = temp_b;
						temp_t = 4;
					}
					//5
					temp_b = getAlignVal(tl,i+1,i+l) + weight[seq2[i]][seq1[tp1]] + g2;
					if (temp_b > temp_a)	{
						temp_a = temp_b;
						temp_t = 5;
					}

					//6
					temp_b = getAlignVal(tl,i,i+l-1) + weight[seq2[i+l]][seq1[tp2]] + g2;
					if (temp_b > temp_a)	{
						temp_a =  temp_b;
						temp_t = 6;
					}
					
					if (pair[seq2[i]][seq2[i+l]]==1)	{
						temp_b = getAlignVal(j,i+1,i+l-1) + g;
						if (temp_b > temp_a)    {
	                                                temp_a =  temp_b;
        	                                        temp_t = 66;
                		                }
					}

					setAlignVal(j,i,i+l,temp_a);
					setTraceVal(j,i,i+l,temp_t);
					
				} else	if (tlist1[j].numson == 1)	{
					//9
					temp_a = getAlignVal(j,i+1,i+l) + g;
					temp_t = 9;
					//affine gap penalty;
					if (getTraceVal(j,i+1,i+l) != 9)    {
						temp_a = temp_a + loopstart;
					}
											
					//10
					temp_b = getAlignVal(j,i,i+l-1) + g;

					if (getTraceVal(j,i,i+l-1) != 10)   {
						temp_b = temp_b + loopstart;
					}
					
					if (temp_b > temp_a)	{
						temp_a = temp_b;
						temp_t = 10;
					}
					//7
					
					temp_b = getAlignVal(tl,i,i+l) + g;

					if (getTraceVal(tl,i,i+l) != 7)   {
						temp_b = temp_b + loopstart;
					}

					if (temp_b > temp_a)	{
						temp_a = temp_b;
						temp_t = 7;
					}
					//8
					temp_b = getAlignVal(tl,i,i+l-1) + weight[seq2[i+l]][seq1[tp2]];
					if ( temp_b > temp_a)	{
						temp_a = temp_b;
						temp_t = 8;
					}
					
					setAlignVal(j,i,i+l,temp_a);
					setTraceVal(j,i,i+l,temp_t);
					
				} else if (tlist1[j].numson == 2)	{
                                        //
                                        //Here is for banded.
/*                                                                                                                                                             
                                        int bi_t = tlist1[tl].pos2;
                                        int minl = (bi_t > i -1 + band1? bi_t -band1 : i-1);
                                        int maxl = (bi_t < i+l - band1? bi_t + band1: i+l);
                                        for(k = minl; k<= maxl; k++)    {
                                                temp_b = align[i][k][tl] + align[k+1][i+l][tr];
                                                if (temp_b > temp_a)    {
                                                       temp_a = temp_b;
                                                       temp_t = -k;
                                                }
                                        }
*/                                                                                                                                                             
                                      //here for nonbanded.

					//
					for(k = i - 1; k<= i+l; k++)	{
						temp_b = getAlignVal(tl,i,k) + getAlignVal(tr,k+1,i+l);
						if (temp_b > temp_a)	{
						       temp_a = temp_b;
						       temp_t = -k;
					       }
					}

					setAlignVal(j,i,i+l,temp_a);
					setTraceVal(j,i,i+l,temp_t);
				}

				} 		//this for banded.
						       
				if(debug  == 2) cout << i << ' '<< ' '<< i+l << ' '<< j << '\n';
				if(debug == 2) cout << temp_a << '\n';
				if(debug == 2) cout << temp_t << '\n';
			}
		}
	}
	//cout << getAlignVal(newbp_num,1,length2) << " here\n";
	//cout << getTraceVal(newbp_num,1,length2) << " trace\n";

}

//*/

void Alignment::Build_tree()
{
	int i, j;

	tempn = 1;
	
	Binarize(1, length1);
	newbp_num = tempn - 1;

	//here is for the test.
	for (i = 0; i < newbp_num; i++)	{
		if (tlist1[i+1].pos2 -tlist1[i+1].pos1 == 0)	{
			tlist1[i+1].left = 0;
		}
		if (tlist1[i+1].pos2 - tlist1[i+1].pos1 == 1 && tlist1[i+1].flag == 1)	{
			tlist1[i+1].left = 0;
		}
		tlist1[i+1].father = -1;
//		cout << tlist1[i+1].pos1 << "--" << tlist1[i+1].pos2 << "--" <<tlist1[i+1].flag << 
//			"left " << tlist1[i+1].left << " right:" << tlist1[i+1].right << " numson:" << tlist1[i+1].numson <<'\n';
	}
	
	for (i = 1; i < newbp_num; i++)	{
		if (tlist1[i].flag == 1 && tlist1[i+1].flag == 1)	{
			tlist1[i].father = i+1;
		} 
	}
}

int 
Alignment::Get_child (int loc1, int loc2, int *childlist)
{
	int i, j, k;
	j = 0;
	i = loc1;
	while ( i <= loc2)	{
		if (stru1[i].pos >= 0)	{
			childlist[j] = i;
			j++;
			i = stru1[i].pos+1;
		} else {
			i++;
		}
			
	}
	return (j);
}

void Alignment::Binarize(int loc1, int loc2)
{
	int i, j, k;
	int *childlist;
	int fedge;
		
	childlist = new int [bp_num];
	k = Get_child(loc1, loc2, childlist);
	
	int start = loc1;
	if (k > 0)	{
		for (i = 0; i < k; i++)	{
			fedge = 0;
			for (j = start; j < childlist[i]; j++)	{
				fedge = 1;
				tlist1[tempn].pos1 = loc1;
				tlist1[tempn].pos2 = j;
				tlist1[tempn].flag = 0;
				tlist1[tempn].numson = 1;
				tlist1[tempn].left = tempn-1;
				tlist1[tempn].right = -1;
				tempn++;
			}
			int templeft = tempn -1;
			
			Binarize(childlist[i]+1, stru1[childlist[i]].pos-1);
			
			tlist1[tempn].pos1 = childlist[i];
			tlist1[tempn].pos2 = stru1[childlist[i]].pos;
			tlist1[tempn].flag = 1;
			tlist1[tempn].numson = 1;
			tlist1[tempn].left = tempn - 1;
			tlist1[tempn].right = -1;
			tempn++;
			
			if (i != 0 || fedge == 1)	{
				tlist1[tempn].pos1 = loc1;
				tlist1[tempn].pos2 = stru1[childlist[i]].pos;
				tlist1[tempn].flag = 0;
				tlist1[tempn].numson = 2;
				tlist1[tempn].left = templeft;
				tlist1[tempn].right = tempn-1;
				tempn++;
			}
			start = stru1[childlist[i]].pos + 1;
		}
		
		for (j = stru1[childlist[k-1]].pos + 1; j <= loc2; j++)	{
			tlist1[tempn].pos1 = loc1;
			tlist1[tempn].pos2 = j;
			tlist1[tempn].flag = 0;
			tlist1[tempn].numson = 1;
			tlist1[tempn].left = tempn - 1;
			tlist1[tempn].right = -1;
			tempn++;
		}
	} else {
		for (j = loc1; j <= loc2; j++)	{
			tlist1[tempn].pos1 = loc1;
			tlist1[tempn].pos2 = j;
			tlist1[tempn].flag = 0;
			tlist1[tempn].numson = 1;
			tlist1[tempn].left = tempn-1;
			tlist1[tempn].right = -1;
			tempn++;
		}
	}
		
	delete [] childlist;
}

	
void
Alignment::Get_structure()	
{
	int i,j,temp;
	Stack *stack = new Stack(length1);

	int po=0,pu=0;
	j = 0;
	for (i = 1; i <= length1; i++)	{
		if (ss1[i] == 0)	{
			stru1[i].pos = -1;
			stru1[i].offset = -1;
		} else if (ss1[i] == 1)	{
			pu++;
			stack->Push(i);
		} else if (ss1[i] == -1)	{
			stru1[i].pos = -1;
			stru1[i].offset = -1;
			po++;
			temp = stack->Pop();
			stru1[temp].pos = i;
			stru1[temp].offset = j;
			j++;
		}
	}

	if (!stack->Empty())	{
		cout << "stack is not empty! " << po << " " << pu << "\n";
	}
	
	bp_num = j;

	delete stack;

//	for (i=1; i <=length1; i++)	{
//		cout << stru1[i].pos << '\n';
//	}
//	cout << bp_num << '\n';

}

void
Alignment::Find_local(int thre)	{
	
	int i1, j1, i2,j2;
	
	i1 = 2;
	j1 = 1;
	i2 = length2 -1;
	j2 = length2;
	
	int tt = getAlignVal(newbp_num,1,length2);
	for (i1=1; i1 <= thre; i1++)	{
		for (i2 =length2; i2 > length2 -thre; i2--)	{
			if (getAlignVal(newbp_num,i1,i2) > tt)	{
				j1 = i1;
				j2 = i2;
				tt = getAlignVal(newbp_num,i1,i2);
			}
		}
	}
	local_s = j1;
	local_e = j2;

}
void 
Alignment::Trace_alignment()
{
	trace_n = 0;
	cout << length2 << '\n';
	cout << newbp_num << '\n';
	Find_local(band1);	
	//cout << local_s << ' ' << local_e << ' ' << getAlignVal(newbp_num,local_s,local_e) << '\n';
	Trace_back(local_s, local_e, newbp_num);
//	Trace_back(1, length2, newbp_num);
}

void
Alignment::Trace_back(int i, int j, int k)
{

		
	Stack *stack1 = new Stack (length1 + length2);
	
	while ( j >= i && k>= 1)
	{
		int tt = getTraceVal(k,i,j);
		switch (tt)	{
			case 1:
				mlist[trace_n] = 4;
				stack1->Push(4);
				trace_n++;
				i++;
				k = tlist1[k].left;
				j--;
				break;

			case 2:
				mlist[trace_n] = 1;
				stack1->Push(1);
				trace_n++;
				k = tlist1[k].left;
				break;

			case 3: 
				mlist[trace_n] = 3;
				trace_n++;
				i = i + 1;
				break;

			case 4:
				stack1->Push(3);
				j = j - 1;
				break;

			case 5: 
				mlist[trace_n] = 2;
				trace_n++;
				stack1->Push(1);
				i++;
				k = tlist1[k].left;
				break;
			
			case 6:
				mlist[trace_n] = 1;
				trace_n++;
				stack1->Push(2);
				j--;
				k = tlist1[k].left;
				break;

			case 7:
				stack1->Push(1);
				k = tlist1[k].left;
				break;
				
			case 8:
				stack1->Push(2);
				j--;
				k = tlist1[k].left;
				break;
			
			case 9:
				mlist[trace_n] = 3;
				trace_n++;
				i = i + 1;
				break;
				
			case 10:
				stack1->Push(3);
				j = j - 1;
				break;
		
			case 66:
				mlist[trace_n] = 5;
				trace_n++;
				stack1->Push(5);
				i = i + 1;
				j = j - 1;
				break;
				
			default:
				if (tt <= 0)	{
					int k1 = tlist1[k].left;
					int k2 = tlist1[k].right;
					Trace_back(i,-tt, k1);
					Trace_back(-tt+1, j, k2);
				}
				j = i - 1;
				k = 0;
				break;
		}
	}

	if (j >= i)	{
		for (int ii = 0; ii <= j -i; ii++)	{
		      mlist[trace_n] = 3;
		      trace_n++;
		}
	} else if (k >= 1)	{
		int d= tlist1[k].pos2 - tlist1[k].pos1 +1;
		for (int ii = 0; ii < d ; ii++)	{
			mlist[trace_n] = 1;
			trace_n++;
		}
	}
	
	while (!stack1->Empty())    {
		mlist[trace_n] = stack1->Pop();
		trace_n++;
	}
	
	delete stack1;
}
				
/*void
Alignment::Get_PSSM(char *A, char *B, int *SS, int M, int N, int *S, int AP, int BP)
{
	int i; 
	int j;
	int k;
	int length = trace_n;
	
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
	
	for (i = 0; i < 5 ; i++)	{
		for (j = 0; j < trace_n; j ++)	{
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

 	/*for (i = 0; i < 2; i++)	{
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
        }*/
	
//}

void Alignment::toStokFile(char* fname, char *A, char *B, int *SS, int M, int N, int *S, int AP, int BP, char* name2){
	ofstream fout(fname);
	static char aline[101], bline[101], cline[101], dline[101];
	register char *a, *b, *c, *d;
	register int i, j, op;
	int lines, ap, bp;
	int pair[4][4]={ 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0};
	
	for (int i = 0; i < trace_n; i++)	{
		cout << mlist[i];
	}
	cout << '\n';

	ap = AP;
	bp = BP;
	a = aline;
	b = bline;
	c = cline;
	d = dline;
	i = 0;
	j = 0;
	lines = 0;
	int ii = 0;
	int idct = 0;
	int totct = 0;
	int loopct = 0;
	while (( i < M || j < N) && ii < trace_n)	
	{
		//cout << S[ii] << " " << j << " " << na_name[B[j]] << endl;
		int func;
		if (ii >= trace_n){
			if(i<M) func = 1;
			else func = 3;
		}else{
			func = S[ii++];
		}
		//cout << S[ii] << " " << j << " " << na_name[B[j]] << endl;
		totct++;
		switch (func)	{
			case 1:
				*a = na_name[A[i]];
				*b = '-';
				*c = ' ';
				if (SS[i] == 1)	{
				       *d = '<';
				} else if (SS[i] == -1)	{
				       *d = '>';
				} else	{
					*d = '.';
				}

				a++;
				b++;
				c++;
				d++;
				i++;
				break;
			case 2:
				*a = na_name[A[i]];
				*b = na_name[B[j]];
				*c = (*a == *b)? '|' : ' ';
				if(*a==*b) idct++;
				if (SS[i] == 1) {
					*d = '<';
				} else if (SS[i] == -1) {
					*d = '>';
				} else  {
					*d = '.';
					}

				a++;
				b++;
				c++;
				d++;
				i++;
				j++;
				break;
				
			case 3: 
				*a = '-';
				*b = na_name[B[j]];
				*c = ' ';
				*d = '.';
				
				
				a++;
				b++;
				c++;
				d++;
				j++;
				break;
			case 4:
				*a = na_name[A[i]];
				*b = na_name[B[j]];
				*c = (*a == *b)? '*' : '+';
				if(*a==*b) idct++;
				if (SS[i] == 1) {
					*d = '<';
				} else if (SS[i] == -1) {
					*d = '>';
				} else  {
					*d = '.';
					}
					
				a++;
				b++;
				c++;
				d++;
				i++;
				j++;
				break;
			case 5:
			       *a = '-';
			       *b = na_name[B[j]];
			       *c = '+';
			       *d = '.';
			       
			       a++;
			       b++;
			       c++;
			       d++;
			       j++;	
		}
		loopct++;
		if (a >= aline+100 || i >= M && j >= N || ii >= trace_n)	{
			*a = *b = *c = *d = '\0';
			
			
			if(name1[0] == 'A' && name1[1] == 'c' && name1[2] == 'c' && name1[3] == 'e' && name1[4] == 's'){
				for(int j=0;j<81;j++){
					name1[j] = name1[j+19];
				} 
			}
			//cout << name2 << endl;	
			if(name2[0] == 'A' && name2[1] == 'c' && name2[2] == 'c' && name2[3] == 'e' && name2[4] == 's'){
				for(int j=0;j<81;j++){
					name2[j] = name2[j+19];
				} 
			}
			//cout << name2 << endl;

			int j;
			if(!(name1[0]>='a' && name1[0]<='z') && !(name1[0]>='A' && name1[0]<='Z')){
				for(int r=77;r>=0;r--){
					name1[r+4] = name1[r];
				}
				name1[0] = 'S';
				name1[1] = 'e';
				name1[2] = 'q';
				name1[3] = '-';
			}
			for(j=0;j<100 && name1[j]!='\0';j++){
			}
			for(;j<=25;j++){
				name1[j] = ' ';	
			}
			name1[j] = '\0';
			if(strlen(name1)==0) name1[0] = 'F';	
			fout << name1 << '\t' << aline << '\n';
			
			if(!(name2[0]>='a' && name2[0]<='z') && !(name2[0]>='A' && name2[0]<='Z')){
				for(int r=77;r>=0;r--){
					name2[r+4] = name2[r];
				}
				name2[0] = 'S';
				name2[1] = 'e';
				name2[2] = 'q';
				name2[2] = '-';
			}

			for(j=0;j<100 && name2[j]!='\0';j++){
			}
			for(;j<=25;j++){
				name2[j] = ' ';	
			}
			name2[j] = '\0';
			if(strlen(name2)==0) name2[0] = 'F';

			fout << name2 << '\t' << bline << '\n';
			
			
			char* gName = new char[100];
			sprintf(gName,"#=GC SS_cons");
			for(j=0;j<100 && gName[j]!='\0';j++){
			}
			for(;j<=25;j++){
				gName[j] = ' ';	
			}
			gName[j] = '\0';
			fout << gName << '\t' << dline << "\n\n\n";
			ap = AP + i;
			bp = BP + j;
			a = aline;
			b = bline;
			c = cline;
			d = dline;
			loopct = 0;
		}
	}
}


	void
Alignment::Display(char *A, char *B, int *SS, int M, int N, int *S, int AP, int BP)
{
	static char aline[101], bline[101], cline[101], dline[101];
	register char *a, *b, *c, *d;
	register int i, j, op;
	int lines, ap, bp;
	int pair[4][4]={ 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0};
	
	for (int i = 0; i < trace_n; i++)	{
		cout << mlist[i];
	}
	cout << '\n';

	ap = AP;
	bp = BP;
	a = aline;
	b = bline;
	c = cline;
	d = dline;
	i = 0;
	j = 0;
	lines = 0;
	int ii = 0;
	int idct = 0;
	int totct = 0;
	while (( i < M || j < N) && ii < trace_n)	
	{
		//cout << S[ii] << " " << j << " " << na_name[B[j]] << endl;
		totct++;
		switch (S[ii++])	{
			case 1:
				*a = na_name[A[i]];
				*b = '-';
				*c = ' ';
				if (SS[i] == 1)	{
				       *d = '<';
				} else if (SS[i] == -1)	{
				       *d = '>';
				} else	{
					*d = '.';
				}
				a++;
				b++;
				c++;
				d++;
				i++;
				break;
			case 2:
				*a = na_name[A[i]];
				*b = na_name[B[j]];
				*c = (*a == *b)? '|' : ' ';
				if(*a==*b) idct++;
				if (SS[i] == 1) {
					*d = '<';
				} else if (SS[i] == -1) {
					*d = '>';
				} else  {
					*d = '.';
					}
				a++;
				b++;
				c++;
				d++;
				i++;
				j++;
				break;
				
			case 3: 
				*a = '-';
				*b = na_name[B[j]];
				*c = ' ';
				*d = ' ';
				a++;
				b++;
				c++;
				d++;
				j++;
				break;
			case 4:
				*a = na_name[A[i]];
				*b = na_name[B[j]];
				*c = (*a == *b)? '*' : '+';
				if(*a==*b) idct++;
				if (SS[i] == 1) {
					*d = '<';
				} else if (SS[i] == -1) {
					*d = '>';
				} else  {
					*d = '.';
					}
				a++;
				b++;
				c++;
				d++;
				i++;
				j++;
				break;
			case 5:
			       *a = '-';
			       *b = na_name[B[j]];
			       *c = '+';
			       *d = ' ';
			       a++;
			       b++;
			       c++;
			       d++;
			       j++;	
		}

		if (a >= aline+100 || i >= M && j >= N || ii >= trace_n)	{
			*a = *b = *c = *d = '\0';
			cout <<'\n';
			for (b = aline+10; b <=a; b += 10)	{
				cout << "    .    :";
			}
			if (b <= a+5)	
				cout << "    .";
			cout << '\n' << ap << '\n' << dline << '\n' << aline << '\n' << cline << '\n'
				<< bline << '\n' << bp << '\n';
			ap = AP + i;
			bp = BP + j;
			a = aline;
			b = bline;
			c = cline;
			d = dline;
		}
	}
	
	cout << "identiy " << (100*idct)/totct << "% "<< endl;
}
