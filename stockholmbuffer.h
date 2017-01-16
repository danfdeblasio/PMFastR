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
class StockholmBuffer {
public:
   StockholmBuffer(char *filename, int maxnum) ; 
   ~StockholmBuffer();
     	//constructor should load both files in main memory. 
   	//Use buffers to speed up reading
//   bool isEOF();
//   char next_char();
//   bool isEndOfCurrentSequence();
//   string &header();
  	 // Returns the header of the sequence at the current fptr location
    char* seq;  //points to the current location in the memory.
    int* structure;
    char name[200];  //the file name;
    int length;
    void Test_buffer();
    void Out_put(unsigned long start, int length);
private:
    char Char2int(char c); //to convert into a number 0, 1, 2, 3...
//  const int total_nuc = 16;
//  char na_name[17] = {'a', 'c', 'g', 't', 'n', 'r', 'y', 'w', 's', 'm',
//	    'k', 'h', 'b', 'v', 'd','x'};
	   
    
};
