
#define BAND 300
int bpair[4][4] = {
			-1, -1, -1, +1,
			-1, -1, +1, -1,
			-1, +1, -1, +1,
			+1, -1, +1, -1
		};
		
		
int weight[4][4] = {	 22, -19, -15, -14,
			-19,  12, -25, -10,
			-15, -25,  10, -17,
			-14, -10, -17, 17};  
			
			
			

int  p_w[4][4][4][4] = { 

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
};   





			

			
	
			 

