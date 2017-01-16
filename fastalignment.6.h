//alignment.1.h
#include "stk_profile.h"

class  Alignment{
public:

	Stk_profile *stkpf;
	Alignment(char *s1, int *st1, char *s2, int l1, int l2, Stk_profile* sIn, char *name,char* stokIn);
	~Alignment();

	char* stokoutfilename;
	bool stokfilenameset;

	typedef struct treelist	{
		int pos1;
		int pos2;
		int flag;
		int numson;
		int left;
		int right;
		int father;
	} TLIST;		//base pair list structure.
	
	typedef struct bplist	{
		int pos;
		int offset;
	} BPLIST;

	int length1;	//query with structure sequence length;
	int length2;	//db_subsequnce length;
	char *seq1;
	char *seq2;
	int band1;
	int band2;
	int bp_num;	//base pair unit number;
	int newbp_num;	//binarized list number;
	int *ss1;
	BPLIST *stru1;	//structure of seqquence1;
	TLIST *tlist1;


	int g; 	//gap score;
	int g2;  //gap in stem;
	int stemplus;
	int loopstart;
	int ***align;
	int ***ga, ***gb, ***gc;
	int ***trace;
	int *mlist;
	int trace_n;
//	double **profile;
//    double **p_profile;
//    int **prof_matrix;
//    int **prof_pmatrix;
//    double ave_profile;

	
	int local_s;
	int local_e;
	
	int tempn;
	
	void Build_tree(); 	//binarized tree;	
	void Get_structure(); 	//get the structure from 
	int Get_child(int loc1, int loc2, int *childlist); // get the all china list
	void Binarize(int loc1, int loc2); // binarized tree.
	void Compute_align();	//compute the align[].
	void Compute_align_b();
	int getAlignVal(int, int, int);
	void setAlignVal(int, int, int, int);
	int getTraceVal(int, int, int);
	void setTraceVal(int, int, int, int);

	void Find_local(int thre);
	
	void Trace_back(int i, int j, int k);
	void Trace_alignment();
	void toStokFile(char* fname, char *A, char *B, int *SS, int M, int N, int *S, int AP, int BP);
	void Display(char *A, char *B, int *SS, int M, int N, int *S, int AP, int BP);
//	void Get_PSSM(char *A, char *B, int *SS, int M, int N, int *S, int AP, int BP);
	
	char* name2;
	
	
};

