//compressF: A program that reads a multi-fasta file f, and outputs 2 files
//          f.seq : contains the entire sequence (concatenated), perhaps with 
//          special characters
//          f.ind: Contains the offset, and the header (tab-delimited)
// At some point, when we work with very large databases, we might need to
// split files.
// Currently, It's only reading one file.


/*=================================================================
| Class FileBuffer will provide the basic File Read functionality
| 
==================================================================*/
//#include "rna_compare/na.h"
//

#include "queue.h"
class Stk_profile {
public:
   Stk_profile(char *filename, int maxnum, int maxseq,bool doShrink=true) ; 
   ~Stk_profile();
     	//constructor should load both files in main memory. 
   	//Use buffers to speed up reading
//   bool isEOF();
//   char next_char();
//   bool isEndOfCurrentSequence();
//   string &header();
  	 // Returns the header of the sequence at the current fptr location
    char* seq;  //points to the current location in the memory.
    int* structure;
    typedef struct bplist   {
          int pos;
          int offset;
    } BPLIST;

    
    BPLIST* struc1;		    
    int bp_num;
    char** name;  //the file name;
    int length;
    int m_numseq;
    int numseq;
    char**aln;
    
    int maxName;

    char**aln_old;
    int* structure_old;
    int length_old;
    double gap_cut;

    void Shrink();

    void Test_buffer();
    void Out_put(unsigned long start, int length);
    
    double **profile;
    double **p_profile;
    int **prof_matrix;
    int **prof_pmatrix;
    double ave_profile;
    
    void Build_profile();
    void Get_structure();
    void Compact();
    void RemoveCol(int);
    
    void toFile(char *fname);
    
    queue* delq;


    void Output_keywords(int w, double cutoff);
private:
    char Char2int(char c); //to convert into a number 0, 1, 2, 3...
//  const int total_nuc = 16;
//  char na_name[17] = {'a', 'c', 'g', 't', 'n', 'r', 'y', 'w', 's', 'm',
//	    'k', 'h', 'b', 'v', 'd','x'};
	   
    
};
