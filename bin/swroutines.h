#ifndef __SWROUTINES__
#define __SWROUTINES__
#include "ssw.h"
#include <map>
#include <queue>
int Do_SW_Pair(Pair* Pairs,char* Pattern,int Flank_Size,int StringLength,int & Err,int Shift,char Sign,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Current_Score,int & Total_SW_Scans);
//int Do_SW_Pair(Pair* Pairs,char* Pattern,int Flank_Size,int StringLength,int & Err,int Shift,char Sign,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Current_Score);
//int Do_SW(BWT* revfmi,SARange* Head_Hits,char* Pattern,int Flank_Size,int StringLength,int & Err,int Shift,char Sign,int Current_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int & Total_SW_Scans,int & filter);
//void ssw_cigar_print(s_align* a,char* Cigar,Cigar_Info & C,char* Ref,char* Pattern);
//void ssw_cigar_process(s_align* a,Cigar_Info & C,char* Ref,char* Pattern,int StringLength);
//void ssw_cigar_processQ(s_align* a,Cigar_Info & C,char* Ref,int Ref_Off,char* Pattern,int Pat_Off,int StringLength,const char* Qual,char* Cigar,int Clip_H,int Clip_T,bool Hard_Penalty=true);
//void ssw_cigar_print(s_align* a,char* Cigar,Cigar_Info & C,char* Ref,char* Pattern,int Clip_H,int Clip_T,int StringLength);
#endif
