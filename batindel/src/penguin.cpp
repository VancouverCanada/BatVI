//#define NDEBUG
//Routines for pairing...  
//TODO: can we change Get_Basic_MapQ to realign all? 
#define __MAIN_CODE__
#include <algorithm>

//{-----------------------------  INCLUDE FILES  -------------------------------------------------/ 
#include <sstream>
#include <cstdio> 
#include "zlib.h"
#include "math.h"
#include <stdarg.h> 
#include <string> 
#include <getopt.h> 
#include <limits.h> 
#include <string.h> 
#include <stdlib.h> 
#include <xmmintrin.h> 
#include <emmintrin.h> 
#include <ctype.h> 
#include "bfix.h" 
#include "assert.h"
#include "ssw.h"
#include "common.h" 
#include "rqindex.h"
#include "batlib.h"
#include <map>
#include <set>
#include <queue>
#include "global.h"
#include "unistd.h"
#include "swroutines.h"
#include "print.h"
#include "filters.h"
#include <pthread.h>
#include "sched.h"
#include "fastsw.h"
#include "print.h"
extern "C" 
{
	#include "iniparser.h"
	#include <time.h>
	#include "MemManager.h"
	#include "MiscUtilities.h"
	#include "TextConverter.h"
	#include "BWT.h"
}


const unsigned MAPPED         =0;
const unsigned PE_SEQ         =1;
const unsigned PROPER_PAIR    =0x2;
const unsigned NOMAP       	 =0x4;
const unsigned MATE_UNMAPPED  =0x8;
const unsigned QUERY_MINUS    =0x10;
const unsigned MATE_MINUS     =0x20;
const unsigned FIRST_READ     =0x40;
const unsigned SECOND_READ    =0x80;
const unsigned AUX	      =0x800;

//}-----------------------------  INCLUDE FILES  -------------------------------------------------/
int CLIP_SAVE_LENGTH=20;
int MIS_IN_AUX=0;
int TOP_TEN=0;
int SW_SIMILARITY_FOR_RESCUE=60;
int DISC_THRESHOLD=10;
int LENGTH_CUTOFF=0;
int BOOST=0;
int Dummy_Int=0;
int CUT_MAX_SWALIGN=200;
int MODE=INT_MAX;
int SEEDSIZE=INT_MAX;
int INSERT=5;
int DELETE=5;
int INSERTSIZE=500;//INT_MAX;
int STD=50;//INT_MAX;
int SW_THRESHOLD=290;
int READS_TO_ESTIMATE=100000;
int MISINSCAN=0;
int SEG_SIZE=18;
const int ScaleQ[]={0,1.5,1.75,2,3,4};
extern const Alignment Default_Alignment={0};
const int ORGSTRINGLENGTH=2200; 
bool PAIRED=FALSE;
bool Hard_Penalty=true;
//extern const int QUALITYCONVERSIONFACTOR=64;
//extern const int QUALITYSCALEFACTOR=33;

bool DASH_DEL=true;
bool ESTIMATE=true;//false;
bool FASTDECODE=false;
bool DEB=false;
bool FASTSW=true;
bool DEBUG_SEGS=false;
unsigned GENOME_SIZE=0;
int QUALITYCONVERSIONFACTOR=33;
int QUALITYSCALEFACTOR=1;
BATPARAMETERS BP;
int THREAD = 0;
int Top_Penalty;
int JUMP=8;
int INDELGAP=21;//15 good for //17 for 8, 19=g00d for 9
FMFILES FMFiles;

typedef std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> ALIGNMENT_Q;
inline void really_free(std::map<unsigned,Alignment>& to_clear)
{
    std::map<unsigned,Alignment> v;
    v.swap(to_clear);
}

//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*
unsigned uabs(unsigned A,unsigned B);
bool Do_Mismatch_Scan(MEMX & MF,MEMX & MC,LEN & L,BWT* fwfmi,BWT* revfmi,int Start_Mis,int End_Mis,int & Last_Mis,int & Head_Top_Count,Hit_Info & H,int & Quality_Score,READ & R,BATREAD & B,FILE* Mishit_File,unsigned Conversion_Factor,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments);
bool Output_Pair(Alignment A1,Alignment A1_P,Alignment B1,Alignment B1_P,int Read_Length);
void Extend_Right(int Plus_Hits,int Minus_Hits,BWT* revfmi,MEMX & MFL,MEMX & MCL,char* Temp_Current_Tag,int StringLength,int & Err, READ & R,int Mis_In_Anchor,int Current_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int & Tot_SW_Scans,int & Filter );
void Extend_Left(int Plus_Hits,int Minus_Hits,BWT* revfmi,MEMX & MFL,MEMX & MCL,char* Temp_Current_Tag,int StringLength,int & Err, READ & R,int Mis_In_Anchor,int Current_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int & Tot_SW_Scans,int & Filter);
int Do_Indel(MEMX & MFLH,MEMX & MCLH,MEMX & MFLT,MEMX & MCLT,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,int StringLength,Pair* & Pairs,FILE* & Single_File,READ & R,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,const Hit_Info & H);
void Show_Progress(unsigned Percentage);
void Read_INI(char* Config_File,unsigned & MAXCOUNT,FMFILES & F,BATPARAMETERS & BP);
void Get_Bases (unsigned Location,int StringLength,char* Org_String);
void Pair_Reads(BWT* fwfmi,BWT* revfmi,RQINDEX & R,FILE* Data_File, In_File IN,unsigned MAXCOUNT,BATPARAMETERS BP,Pair* Pairs,unsigned Entries);
void Get_Best_Alignment_Pair(Alignment & A,Alignment & B,READ & R,const int StringLength,BATREAD & Read,Hit_Info & H,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool Force_Indel,int & Clip_H,int & Clip_T,char* CIG,bool PRINT,bool DUMMY_FORCED);
void Parse_Command_line(int argc, char* argv[],unsigned & MAXCOUNT,FMFILES & FMFiles,BATPARAMETERS & BP);
int Head_Tail(BWT* revfmi,SARange* Head_Hits,SARange* Tail_Hits,int Insert,int MAXCOUNT,char & In_Large,RQINDEX & R,unsigned Entries,Pair* Pairs,int & Pairs_Index,int & HITS,int & Err,unsigned Conversion_Factor);
void *Map_And_Pair_Solexa(void *T);
void Verbose(char *BWTFILE,char *OCCFILE,char *REVBWTINDEX,char *REVOCCFILE,char *PATTERNFILE,char *HITSFILE,char* LOCATIONFILE,char MAX_MISMATCHES,int Patternfile_Count,char* PATTERNFILE1,char FILETYPE,LEN & L,char FORCESOLID);
void Init(In_File & IN,FMFILES F,RQINDEX R,BATPARAMETERS & BP,char Solid,char Npolicy,LEN & L);
void  Paired_Extension(int Last_MisT,int Last_MisH,char* Fwd_Read,char *Revcomp_Read, RQINDEX & RQ,Pair* & Pairs,SARange* & MFH_Hit_Array,SARange* & MFT_Hit_Array,SARange* & MCH_Hit_Array,SARange* & MCT_Hit_Array,int StringLength,READ & R,int & Err,unsigned Conversion_Factor,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Current_Score,int & Tot_SW_Scans);
void Launch_Threads(int NTHREAD, void* (*Map_t)(void*),Thread_Arg T);
bool Report_SW_Hits(const int Err,READ & R,Final_Hit & Single_File,const int StringLength,BATREAD & Read,Hit_Info & Mismatch_Hit,int Quality_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool Force_Indel,bool PRINT,bool DUMMY_FORCED=false);
bool Get_Info(MEMX & MF,MEMX & MC,int STRINGLENGTH,Hit_Info & H,unsigned Conversion_Factor,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,READ & R);
bool Unique_Hits(int Plus_Hits,int Minus_Hits,SARange & P,SARange & M);
bool Check_Subopt(int & Plus_Hits,int & Minus_Hits,int Top_Mis,int Subopt_Mis,READ & R, int StringLength,MEMX & MF,MEMX & MC,float & Top_Score,float & Top_BQScore,float & Sub_Score,float & Sub_BQScore, int & Quality_Score);
Alignment Realign(Hit_Info &  H,BATREAD & Read,int StringLength,const READ & R,bool Dont_Push_To_Q,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments);
Alignment Realign(Hit_Info &  H,BATREAD & Read,int StringLength,const READ & R,bool Dont_Push_To_Q,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments);
Alignment RealignX(Hit_Info &  H,BATREAD & Read,int StringLength, READ & R,bool Dont_Push_To_Q,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,char* Cigar,int & Clip_H,int & Clip_T,int & Filter=Dummy_Int,bool Do_Filter=false);
bool Remove_Duplicates_And_Filter_Top_Hits(Hit_Info *Hit_List,int & Hit_Count,int StringLength);
float Calc_Top_Score(MEMX & MF,MEMX & MC,float & Top_BQ,int Top_Mis,int StringLength,int Plus_Hits,int Minus_Hits,READ & R);
bool Recover_With_SW(int Plus_Hits,int Minus_Hits,READ & R,BATREAD & BR, int StringLength,MEMX & MF,MEMX & MC,unsigned Conversion_Factor,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,Hit_Info & H);
bool Extend_With_SW(int Plus_Hits,int Minus_Hits,READ & R,BATREAD & BR, int StringLength,MEMX & MF,MEMX & MC,unsigned Conversion_Factor,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,Hit_Info & H,bool Mismatch_Scan=true);
bool Map_One_Seg(READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MF,MEMX & MC,LEN & L,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,Pair* & Pairs,bool PRINT,Hit_Info & H,int & Quality_Score,int Segment_Length,int SEG_SIZE,int SHIFT_SEG_SIZE);
unsigned SA2Loc(SARange S,int Pos,unsigned Conversion_Factor);
void Mode_Parameters(BATPARAMETERS BP);
void Set_Force_Indel(bool & Force_Indel,int Last_Mis,Hit_Info & H,int Avg_Q);
int Calculate_Average_Quality(READ & R);
void Set_Affinity();
Alignment RealignFast(Hit_Info &  H,int StringLength, READ & R,int OFF,int Filter,bool Do_Filter);
Alignment RealignFastMinus(Hit_Info &  H,BATREAD & Read,int StringLength, READ & R,int OFF,int Filter,bool Do_Filter);
void FreeQ(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & t );
bool Report_Single(READ & R,Final_Hit & Single_File,const int StringLength,BATREAD & Read,bool & Print_Status,int Clip_H,int Clip_T,Alignment & A);
void Get_Basic_MapQ(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,Alignment & C1, int & MapQ);
void Adjust_Alignments(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Offset,READ & RTemp, BATREAD & BTemp);
void Pop_And_Realign(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,int Offset,READ & RTemp, BATREAD & BTemp);
bool SW_List(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Offset,READ & RTemp, BATREAD & BTemp,bool & List_Exceeded,int & List_Size);
bool Correct_Orientation(Alignment A,Alignment A_P,int Extra_Bit=0);
bool Find_Paired(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & B,std::map<unsigned,Alignment> & D,std::map<unsigned,Alignment> & D_P,int Extra_Bit=0,int SW_Compare=0);
bool Align_Difference(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,unsigned U);
bool Rescue_Mate(unsigned Loc,char Sign,int StringLength,char* Current_Tag,char* Q,int Flank, int Shift, bool Dont_Push_To_Q,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,char* Cigar,int & Clip_H,int & Clip_T,int & Filter,bool Do_Filter);
void Rescue_One_Side(std::map<unsigned,Alignment> & D,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments_P,READ & RTemp_P,BATREAD & BTemp_P);
void Rescue_One_Side_X(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments_P,READ & RTemp_P,BATREAD & BTemp_P);
bool Full_Rescue(READ & RTemp,READ & RTemp_P,BATREAD & BTemp,BATREAD & BTemp_P,int Read_Length,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments_P,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Good_Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Good_Alignments_P,Hit_Info & H1,Hit_Info & H1_P,FILE* Single_File,int Quality_Score1,int Quality_Score1_P,Alignment & A1,Alignment & A1_P,int MapQ1,int MapQ2,bool Max_Pass,MEMX & MF,MEMX & MC,BWT* fwfmi,BWT* revfmi);
bool Proper_Pair(READ & RTemp,READ & RTemp_P,BATREAD & BTemp,BATREAD & BTemp_P,int Read_Length,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments_P,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Good_Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Good_Alignments_P,Hit_Info & H1,Hit_Info & H1_P,FILE* Single_File,int Quality_Score1,int Quality_Score1_P,Alignment & A1,Alignment & A1_P,MEMX & MF,MEMX & MC,bool Max_Pass);
void Mate_Rescue(READ & RTemp,READ & RTemp_P,BATREAD & BTemp,BATREAD & BTemp_P,int Read_Length,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments_P,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Good_Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Good_Alignments_P,Hit_Info & H1,Hit_Info & H1_P,FILE* Single_File,int Quality_Score1,int Quality_Score1_P,Alignment & A1,int MapQ2,MEMX & MF,MEMX & MC,BWT* fwfmi,BWT* revfmi);
void Remove_Dup_Top(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments,int Gap);
void Fix_Offset(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & T,int Offset,int Neg_Off);
void Print_Pair(FILE* Single_File,Final_Hit & H,Final_Hit & T,READ & R1, READ & R2,MEMX & MF,MEMX & MC,BWT* fwfmi,BWT* revfmi);
int Get_Skip(std::string & S);
bool Check_Proper_Pair(int S1,int S2,unsigned Loc1, unsigned Loc2,int Extra_Bit);
void Estimate_Insert(int & INSERTSIZE,int & STD);
void Rescue_Clip(std::string S,int & P_Clip,int & S_Clip);
bool Map_Clip(MEMX & MF,MEMX & MC,BWT* fwfmi,BWT* revfmi,Final_Hit & H,READ & R1,bool Is_Suffix,Final_Hit & Aux_Hit);
void Print_Aux_Hit(FILE* Single_File,Final_Hit & H,Final_Hit & Aux,READ & R,int Clip1,int Clip2);
pthread_mutex_t Lock_Estimate;
std::vector<int> Estimate;
//}-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*

#undef DEBUG
//bool TESTMODE=false;
const bool EXACT=false;
const int QUAL_UNMAPPED=20000;
int SW_STRING_BUFFER=1400;
const char* Code_To_CharCAPS="ACGT";
const int DISK_BUFFER_SIZE=5000;
const int MAXALN=30000;
const int MAXCOUNT=200;
const int MIN_SCORE_DIFF=10;
bool PRINT_ALL_HITS=false;
int INDEX_RESOLUTION=30000;
int Inter_MM = 1;
int Max_MM_GAP = 2;
int Max_MM_GAP_Adjust;
int Max_MM = 5;
int MAX_SW_HITS = 100;
int MAX_SW_sav;
bool USE_MULTI_OUT=false;//true;
bool MAIN_HEADER_ONLY=false;;
bool STACK_LOWQ=false;
const int DEFAULTFASTAQUAL=35;//40;
LEN L_Main;
In_File IN;
INFILE Head_File;
INFILE Tail_File;
time_t Start_Time,End_Time;
FILE* Main_Out;
FILE* Output_File1;
FILE* Output_File2;
int gap_openP,gap_extensionP;
int MAX_PER_LIST=60;
int SCOREGAP=20;
std::string RGID;
//{---------------------------- GLOBAL VARIABLES -------------------------------------------------


int main(int argc, char* argv[])
{
	FILE* Data_File;
	FILE* Mishit_File;
	MMPool *mmPool;
	unsigned MAXCOUNT=1;
	BP.PAIRING_MODE=NORMALFILEMODE;BP.NTHREADS=0;BP.FORCELENGTH=0;BP.ROLLOVER=FALSE;BP.SCANBOTH=FALSE;BP.ONEFMINDEX =FALSE; BP.MAXHITS=1;BP.CMD_Buffer=new char [5000];


	gap_openP=40,gap_extensionP=6;
	Read_INI(NULL,MAXCOUNT,FMFiles,BP);
	Parse_Command_line(argc,argv,MAXCOUNT,FMFiles,BP);	
	init_SSW();Build_Pow10();
	init_SSW_Clip(1 /*match*/,3 /*mismatch*/,11 /*gap_open*/,4 /*gap_extension*/);
	if(CONFIG_FILE) Read_INI(CONFIG_FILE,MAXCOUNT,FMFiles,BP);
	if (MISC_VERB) fprintf(stderr,"BatINDEL V 1.00\n");
	if (BP.MAX_MISMATCHES != INT_MAX)
	{
		if (BP.MAX_MISMATCHES <=5) Inter_MM=BP.MAX_MISMATCHES;
		else fprintf(stderr,"Error tolerence too high ...\n");
	}
	

	Match_FA_Score= -10*log10(Pr(DEFAULTFASTAQUAL));
	Mis_FA_Score= -10*log10((1-Pr(DEFAULTFASTAQUAL))/3);
	//Head_File.Input_File=File_Open(BP.PATTERNFILE,"r");
	//if(PAIRED) Tail_File.Input_File=File_Open(BP.PATTERNFILE1,"r");
	File_OpenGZ(BP.PATTERNFILE,"r",Head_File);
	if(PAIRED) File_OpenGZ(BP.PATTERNFILE1,"r",Tail_File);
	Analyze_File(Head_File,L_Main);
	IN.HEAD_LENGTH=IN.TAIL_LENGTH=IN.STRINGLENGTH=L_Main.STRINGLENGTH;
	
	Init(IN,FMFiles,RQ,BP,Head_File.SOLID,0,L_Main);

	unsigned Location_Array[80];
	Load_FM_and_Sample(fwfmi,revfmi,mmPool,FMFiles); 
	GENOME_SIZE=revfmi->textLength;assert(GENOME_SIZE==fwfmi->textLength);
	Load_Location(FMFiles.LOCATIONFILE,Annotations,S,E,Location_Array);
	if (NFILE) Load_LocationN(FMFiles.NLOCATIONFILE,Annotations,S,E,Location_Array);
	if(!USE_MULTI_OUT)
	{
		std::string O1=FMFiles.OUTPUTFILE;O1+="_1.fa";
		std::string O2=FMFiles.OUTPUTFILE;O2+="_2.fa";
		Output_File1=File_Open(O1.c_str(),"w");
		Output_File2=File_Open(O2.c_str(),"w");
		std::map <unsigned, Ann_Info> ::iterator S,E;
		S=Annotations.begin();E=Annotations.end();
		unsigned CSize=0;
		while (S!=E)
		{
			Ann_Info T=S->second;
			if(T.Size > CSize) CSize=T.Size;
			//fprintf(Main_Out,"@SQ\tSN:%s\tLN:%d\n",T.Name,T.Size);
			S++;
		}
		//if(NFILE && CSize) fprintf(Main_Out,"@SQ\tSN:%s\tLN:%d\n","NOISE",CSize);
		char Current_Dir[1000];
		if (!getcwd(Current_Dir,990))
		{
			sprintf (Current_Dir,"%d",rand());
		}
		std::ostringstream ostr;
		ostr << (rand()%(INT_MAX));
		RGID=ostr.str();
		//fprintf(Main_Out,"@RG\tID:%s\tSM:%s\tLB:%s\n",RGID.c_str(),Current_Dir,BP.PATTERNFILE);
		//fprintf(Main_Out,"@PG\tID:PEnGuin\tCL:%s",BP.CMD_Buffer);
	}
//********************************************************************************************************
	Thread_Arg T;
	if (THREAD)
	{
		printf("Estimating insert size..\n");
		ESTIMATE=false;
		Launch_Threads(THREAD, Map_And_Pair_Solexa,T);

	}
	else
	{
		Set_Affinity();
		ESTIMATE=false;
		Map_And_Pair_Solexa(NULL);
	}

//********************************************************************************************************

	if (PROGRESSBAR) fprintf(stderr,"\r[++++++++100%%+++++++++]\n");//progress bar....
	if (MISC_VERB)
	{
		fprintf(stderr,"%d / %d Reads / Pairs ...\n",Total_Reads,Missed_Hits);
		//printf("%d Large reads....\n",Large);
		time(&End_Time);fprintf(stderr,"\n Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	}
	if (LOG_SUCCESS_FILE) fprintf(Log_SFile,"DONE\n");
}

//------------------------------- Print /Verify the reads ----------------------------------------------------
void *Map_And_Pair_Solexa(void *T)
{
	Thread_Arg *TA;
	TA=(Thread_Arg*) T;
	int Thread_ID;
	if(T)
	{
		Thread_ID=TA->ThreadID;
	}

	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_Mid;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_End;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_Mid;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_End;

	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_P;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_Mid_P;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_P;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_Mid_P;

	Pair* Pairs;
	FILE* Discordant_File;
	FILE* Single_File;
	FILE* Unmapped_File;
	FILE* Multi_File;
	FILE* Mishit_File;
	MEMLOOK MLook;
	MMPool *mmPool;
	LEN L=L_Main,L_Half,L_Third;L.IGNOREHEAD=0;
	OUTPUT O;
	HEADER Header;
	//time_t Start_Time,End_Time;
	//time(&Start_Time);
	int MF1_Top_End,MC1_Top_End;
	char* Disc_Hit_Buffer;
#ifndef NDEBUG
	//DISC_HIT_BUFFER_ORG=Disc_Hit_Buffer_Org;
#endif
//{--------------------------- INIT STUF ---------------------------------------

	L_Half=L_Third=L;
	LEN mL=L;
	if (BP.FORCELENGTH) 
	{
		Split_Read(BP.FORCELENGTH,L);SEEDSIZE=BP.FORCELENGTH;
	}
	else 
	{
		Split_Read(L.STRINGLENGTH_ORG,L);SEEDSIZE=L.STRINGLENGTH_ORG;
	}
	if (!(Pairs=(Pair*)malloc(sizeof(Pair)*30000))) {if(LOG_SUCCESS_FILE) fprintf(Log_SFile,"Init():malloc error...\n");fprintf(stderr,"Init():malloc error...\n");exit(-1);}
	Split_Read(IN.STRINGLENGTH/2,L_Half);
	Split_Read(IN.STRINGLENGTH/3,L_Third);
	Split_Read(20,mL);
// Initialise indexes 
//
	MLook.Lookupsize=3;//Get_Lookup_Size(BP.MAX_MISMATCHES,L.STRINGLENGTH);
	Build_Tables(fwfmi,revfmi,MLook);


	time(&Start_Time);

// Misc stuff 
	if(THREAD==0 || THREAD==1) {Thread_ID=1;}
	if(Thread_ID==1) Verbose(FMFiles.BWTFILE,FMFiles.OCCFILE,FMFiles.REVBWTINDEX,FMFiles.REVOCCFILE,BP.PATTERNFILE,FMFiles.OUTPUTFILE,FMFiles.LOCATIONFILE,Inter_MM,BP.Patternfile_Count,BP.PATTERNFILE1,Head_File.FILETYPE,L,BP.FORCESOLID);

	if(PRINT_MISHIT) Mishit_File=File_Open(BP.MISFILE1,"w");
	//if (WRITE_SINGLE) Single_File=fopen(FMFiles.SINGLEFILE,"w");else 
	WRITE_DISCORDANT=WRITE_MULTI=WRITE_UNMAPPED=WRITE_SINGLE=TRUE;


	FILE* Log_File=File_Open(LOGFILE,"w");GLog_File=Log_File;
	unsigned Mapped=0,Actual_Tag=0,Large=0,Total_Reads=0;
	int LOOKUPSIZE,MAX_MISMATCHES=BP.MAX_MISMATCHES;
	char ONEFMINDEX=BP.ONEFMINDEX;
	READ R;BATREAD B;
	READ M;
	GUESS G;GUESS G1;
	MEMX MF,MC,MFLH,MCLH,MFLT,MCLT;//MEMX MF1,MC1;
	MEMX MFH,MCH,MFT,MCT;//One third pref/suf..
	MEMX mF,mC;

// calculate read portions 
	LOOKUPSIZE=3;//Get_Lookup_Size(MAX_MISMATCHES,L.STRINGLENGTH);
	Split_Read(SEG_SIZE,L);
	B.IGNOREHEAD=L.IGNOREHEAD;B.StringLength=L.STRINGLENGTH;
	MF.Lookupsize=LOOKUPSIZE;MC.Lookupsize=LOOKUPSIZE;MFLH.Lookupsize=LOOKUPSIZE;MCLH.Lookupsize=LOOKUPSIZE;
	MFLT.Lookupsize=LOOKUPSIZE;MCLT.Lookupsize=LOOKUPSIZE;
	MFH.Lookupsize=LOOKUPSIZE;MCH.Lookupsize=LOOKUPSIZE;MFT.Lookupsize=LOOKUPSIZE;MCT.Lookupsize=LOOKUPSIZE;
	MF.L=MC.L=MFLH.L=MCLH.L=MFLT.L=MCLT.L=MFH.L=MCH.L=MFT.L=MCT.L=L;
	mF.L=mC.L=mL;mF.Lookupsize=LOOKUPSIZE;mC.Lookupsize=LOOKUPSIZE;
// initialise memory structures 
	Copy_MEM(MLook,MF,MC,MAX_MISMATCHES);
	Copy_MEM(MLook,mF,mC,MAX_MISMATCHES);
	int Progress=0;unsigned Number_of_Tags=1000;
	if (PROGRESSBAR) Init_Progress();
	int NCount_For_SW=L.STRINGLENGTH/4;
	unsigned Conversion_Factor;

	bwase_initialize(); 
	INPUT_FILE_TYPE=Head_File.FILETYPE;
	MAX_SW_sav=MAX_SW;
	int Read_Length=Head_File.TAG_COPY_LEN;
	SW_THRESHOLD=80*Read_Length*match/100;
//}--------------------------- INIT STUF ---------------------------------------

	while (Read_Tag(Head_File.Input_File,Tail_File.Input_File,Head_File.GZ_File,Tail_File.GZ_File,Head_File.isGZ,Tail_File.isGZ,Head_File.FILETYPE,R,M))
	{
//Read Head start ..
		RECOVER_N=0;
		R.Tag_Number=FIRST_READ;M.Tag_Number=SECOND_READ;
		if(SEG_SIZE>=Read_Length) SEG_SIZE=Read_Length-1;
		int SHIFT_SEG_SIZE=(2*SEG_SIZE>Read_Length)? Read_Length-SEG_SIZE-1:SEG_SIZE;
		int Hits_Printed=0;
		Progress++;Total_Reads++;
		Actual_Tag++;
		if (READS_TO_PROCESS && READS_TO_PROCESS <= Total_Reads) break;
		if (ESTIMATE && READS_TO_ESTIMATE <= Total_Reads) break;
		if (Thread_ID==1 && Progress>Number_of_Tags && PROGRESSBAR) 
		{
			off64_t Current_Pos=Head_File.isGZ? gztell(Head_File.GZ_File):ftello64(Head_File.Input_File);
			off64_t Average_Length=Current_Pos/Actual_Tag+1;
			Number_of_Tags=(Head_File.File_Size/Average_Length)/20;
			Progress=0;
			Show_Progress(Current_Pos*100/Head_File.File_Size);
		}

		Alignment A1,B1,C1;
		Alignment A1_P,C1_P,C2_P;
		int MapQ1=1,MapQ2=1,MapQ3=1;
		int MapQ1_P=1,MapQ2_P=1,MapQ3_P=1;

		A1.Loc=A1_P.Loc=UINT_MAX;
		Hit_Info H1,H2,H3;int Quality_Score1,Quality_Score2,Quality_Score3;
		Final_Hit Head_Hit,Mate_Hit;

		bool Seg_Mapped=false;
		int Shift=0;
		READ RTemp;BATREAD BTemp;
		READ RTemp_P;BATREAD BTemp_P;
		while(!Seg_Mapped)//Prefix mapped for head..
		{
			RTemp=R;BTemp=B;
			Seg_Mapped=Map_One_Seg(R,B,Conversion_Factor,MF,MC,L,Actual_Tag,Head_Hit,Mishit_File,Alignments,Good_Alignments,Pairs,false,H1,Quality_Score1,Shift,SEG_SIZE,0);
			R=RTemp;B=BTemp;
			if(!Seg_Mapped)//Prefix mapped for Tail..
			{
				RTemp_P=M;BTemp_P=B;
				Seg_Mapped=Map_One_Seg(M,B,Conversion_Factor,MF,MC,L,Actual_Tag,Head_Hit,Mishit_File,Alignments_P,Good_Alignments_P,Pairs,false,H1,Quality_Score1,Shift,SEG_SIZE,0);
				M=RTemp_P;B=BTemp_P;
				if(!Seg_Mapped)//Suffix mapped for head..
				{
					RTemp=R;BTemp=B;
					Seg_Mapped=Map_One_Seg(R,B,Conversion_Factor,MF,MC,L,Actual_Tag,Head_Hit,Mishit_File,Alignments,Good_Alignments,Pairs,false,H1,Quality_Score1,0,SEG_SIZE,Read_Length-SEG_SIZE-Shift-1);
					R=RTemp;B=BTemp;
					if(!Seg_Mapped)//Suffix mapped for Tail..
					{
						RTemp_P=M;BTemp_P=B;
						Seg_Mapped=Map_One_Seg(M,B,Conversion_Factor,MF,MC,L,Actual_Tag,Head_Hit,Mishit_File,Alignments_P,Good_Alignments_P,Pairs,false,H1,Quality_Score1,0,SEG_SIZE,Read_Length-SEG_SIZE-Shift-1);
						M=RTemp_P;B=BTemp_P;
					}
				}
			}

			if(Read_Length/2+5< Shift+SEG_SIZE)
				break;
			else
				Shift+=5;
		}
		if(Seg_Mapped)
		{
			Print_Unmapped(Output_File1,R,Output_File2,M,Read_Length);
		}

	}
}


void Verbose(char *BWTFILE,char *OCCFILE,char *REVBWTINDEX,char *REVOCCFILE,char *PATTERNFILE,char *HITSFILE,char* LOCATIONFILE,char MAX_MISMATCHES,int Patternfile_Count,char* PATTERNFILE1,char FILETYPE,LEN & L,char FORCESOLID)
{
	if (!MISC_VERB) return;
	fprintf(stderr,"BAT ALIGN - Long ..\n");
	fprintf(stderr,"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n");
	fprintf(stderr,"Using the genome files\n %s\t %s\n %s\t %s\n", BWTFILE,OCCFILE,REVBWTINDEX,REVOCCFILE); 
	fprintf(stderr,"Location file: %s\n", LOCATIONFILE); 
	fprintf(stderr,"Query File : %s \t\t Output file: %s\n",PATTERNFILE,HITSFILE);
	if(Patternfile_Count) fprintf(stderr,"Mate File : %s \n ",PATTERNFILE1);
	fprintf(stderr,"Length of Tags: %d\t", L.STRINGLENGTH);
	fprintf(stderr,"Mismatches allowed : %d\n",MAX_MISMATCHES);
	if (SOLID) 
	{
		if (FORCESOLID) fprintf (stderr,"DIBASE-SOLiD reads...\n");
		else fprintf (stderr,"SOLiD reads...\n");
	}
	if (FILETYPE == FQ) fprintf(stderr,"FASTQ file..\n"); else fprintf(stderr,"FASTA file..\n");
	if (UNIQ_HIT) fprintf(stderr,"Unique hit mode..\n");
	if (PRIORITYMODE) fprintf(stderr,"Lowest mismatch mode..\n");
	fprintf(stderr,"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n");
}

void Init(In_File & IN,FMFILES F,RQINDEX R,BATPARAMETERS & BP,char Solid,char Npolicy,LEN & L)
{
	JUMP=BP.INDELSIZE;
	INSERT=DELETE=BP.INDELSIZE;
	INDELGAP=2*JUMP+1;

	SOLID=Solid; 
	srand(0);
	if (SOLID) 
	{
		L.IGNOREHEAD=2;
		PAIRING_SCHEME = 5353;
	}
	//NPOLICY=Npolicy;
	Char_To_Code['A']=0;Char_To_Code['C']=1;Char_To_Code['G']=2;Char_To_Code['T']=3;Char_To_Code['a']=0;Char_To_Code['c']=1;Char_To_Code['g']=2;Char_To_Code['t']=3;
	Char_To_Code['0']=0;Char_To_Code['1']=1;Char_To_Code['2']=2;Char_To_Code['3']=3;Char_To_Code['4']=0;
	Char_To_CodeC['A']=3;Char_To_CodeC['C']=2;Char_To_CodeC['G']=1;Char_To_CodeC['T']=0;Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;Char_To_CodeC['n']=0;Char_To_CodeC['N']=0;//we are using character count to store the fmicode for acgt
	Char_To_CodeC['N']=3;Char_To_CodeC['n']=3;Char_To_CodeC['.']=3;Char_To_CodeC[0]=3;Char_To_CodeC[1]=2;Char_To_CodeC[2]=1;Char_To_CodeC[3]=0;Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;Char_To_CodeC['-']='-';Char_To_CodeC['+']='+';//we are using character count to store the fmicode for acgt
	Char_To_CodeC['0']=3;Char_To_CodeC['1']=2;Char_To_CodeC['2']=1;Char_To_CodeC['3']=0;

	Char_To_CharC['A']='T';Char_To_CharC['C']='G';Char_To_CharC['G']='C';Char_To_CharC['T']='A';Char_To_CharC['a']='t';Char_To_CharC['c']='g';Char_To_CharC['g']='c';Char_To_CharC['t']='a';Char_To_CharC['n']='n';Char_To_CharC['N']='N';//we are using character count to store the fmicode for acgt
	//Char_To_CodeCS['3']=3;Char_To_CodeCS['2']=2;Char_To_CodeCS['1']=1;Char_To_CodeCS['0']=0;Char_To_CodeCS['.']=0;
	Code_To_CodeC[0]=3;Code_To_CodeC[1]=2;Code_To_CodeC[2]=1;Code_To_CodeC[3]=0;
	//if (SOLID) sprintf(Default_Cigar,"%dM",IN.STRINGLENGTH-1);else sprintf(Default_Cigar,"%dM",IN.STRINGLENGTH);
	if (BP.STD)//set a margin of avg+3*std
	{
		BP.FLANKSIZE= 3*BP.STD;//for normal pairing
		BP.SW_FLANKSIZE= 2*IN.STRINGLENGTH+6*BP.STD;
	}
	else if (!BP.SW_FLANKSIZE)
	{
		BP.SW_FLANKSIZE= IN.STRINGLENGTH/2+50;//check 50 bp either side if not specified...
		if (BP.SW_FLANKSIZE > BP.INSERTSIZE) BP.SW_FLANKSIZE=BP.INSERTSIZE;
	}
	if (BP.SW_FLANKSIZE>=SW_STRING_BUFFER) {printf("SW Flank too large.. Defaulting\n");BP.SW_FLANKSIZE= IN.STRINGLENGTH/2+50;}
	if (BP.PLUSSTRAND) BP.PLUSSTRAND=IN.STRINGLENGTH;
	IN.Positive_Head=IN.Tag_Copy;


	//FILE* Original_File=File_Open(F.BINFILE,"rb");
	//Original_Text=(unsigned char*) malloc(Get_File_Size(Original_File));
	//if(!fread(Original_Text,Get_File_Size(Original_File),1,Original_File))fprintf (stderr,"Init(): Error loading genome...\n");

	if (IN.STRINGLENGTH < SW_MIN_MATCH_LEN) SW_MIN_MATCH_LEN=IN.STRINGLENGTH-1;
	if (BP.MISHITS) 
	{
		Mishit_File1=File_Open(BP.MISFILE1,"w"); 
		Mishit_File2=File_Open(BP.MISFILE2,"w"); 
	}

	if(MISC_VERB)
	{
		if (R.COMPRESS) fprintf (stderr,"Using compressed index...\n");
		fprintf(stderr,"Read length %d\n",IN.STRINGLENGTH);
		if (BP.SMITH_WATERMAN)
		{
			fprintf(stderr,"Mate rescue mode : Flank size - %d\n",BP.SW_FLANKSIZE);
		}
		fprintf(stderr,"Max Gap : %d\n",BP.INSERTSIZE+BP.FLANKSIZE);
		fprintf(stderr,"Using Pairing index %s\t %s\n", F.INDFILE,F.BLKFILE); 
		fprintf (stderr,"------------------------------------------------------------------------------------------------------------------\n");
	}

}

void FreeQ(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & t ) 
{
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> tmp; 
	while(!t.empty())
	{
		t.pop();
	}
	std::swap(t,tmp );
}

const int FHSIZE=10;//00;
const int FHSKIP=500;


bool Do_Mismatch_Scan(MEMX & MF,MEMX & MC,LEN & L,BWT* fwfmi,BWT* revfmi,int Start_Mis,int End_Mis,int & Last_Mis,int & Head_Top_Count,Hit_Info & H,int & Quality_Score,READ & R,BATREAD & B,FILE* Mishit_File,unsigned Conversion_Factor,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments)
{

	H.Loc=0;H.Indel=0;H.Score=0;H.QScore=-1;H.Cigar[0]=0;H.Cigar[1]=0;
	Last_Mis= -1;
	Last_Mis=Scan(MF,MC,End_Mis,L,fwfmi,revfmi,Start_Mis,Head_Top_Count,MISINSCAN);
	if(Last_Mis != -1)
	{
		return true;
	}
	else
	{
		return false;
	}

	int Plus_Hits=MF.Hit_Array_Ptr-1,Minus_Hits=MC.Hit_Array_Ptr-1;
	int Sub_OptMF=MF.Hit_Array_Ptr,Sub_OptMC=MC.Hit_Array_Ptr,Sub_OptMIS=0;
	int Sub_Plus_Hits,Sub_Minus_Hits;
	float Top_Score=FLT_MAX,Sub_Opt_Score=FLT_MAX,Sub_BQScore=FLT_MAX;
	float Top_BQScore=FLT_MAX;bool Close_Hits=false;
	Top_Penalty=0;MF.Exact_Match_Forward[L.LH-1].Start=0;MC.Exact_Match_Forward[L.LH-1].Start=0;
	if(Last_Mis==0)
	{
		SARange SA=MF.Exact_Match_Forward[L.LH-1];
		if(SA.Start)//Get stat for later mapQ calculation..
		{
			if (SA.Skip)
				Top_Penalty++;
			else
				Top_Penalty=SA.End-SA.Start;
		}
		SA=MC.Exact_Match_Forward[L.LH-1];
		if(SA.Start)
		{
			if (SA.Skip)
				Top_Penalty++;
			else
				Top_Penalty+=(SA.End-SA.Start);
		}
	}

	if(Plus_Hits+Minus_Hits==1 && Last_Mis!=End_Mis && !HEURISTIC)
	{
		bool Deep_Scan=true;
		if(BOOST>=3)
		{
			if(BOOST>=4)
				Deep_Scan=false;
			else
			{
				if(Last_Mis>=2)
					Deep_Scan=false;

			}

		}

		if(Deep_Scan)
		{
			Sub_OptMIS= Scan(MF,MC,Last_Mis+1,L,fwfmi,revfmi,Last_Mis+1,Head_Top_Count,(BOOST>=2)? 1:INT_MAX);//2);
		}
		else
		{
			Sub_OptMIS=-1;
		}

		if(Sub_OptMIS != -1) 
		{
			Close_Hits=true;
			Sub_Plus_Hits=MF.Hit_Array_Ptr-Plus_Hits-2;Sub_Minus_Hits=MC.Hit_Array_Ptr-Minus_Hits-2;
			if(!Check_Subopt(Plus_Hits,Minus_Hits,Last_Mis,Sub_OptMIS,R,L.STRINGLENGTH,MF,MC,Top_Score,Top_BQScore,Sub_Opt_Score,Sub_BQScore,Quality_Score))
			{
				if(Recover_With_SW(Plus_Hits,Minus_Hits,R,B,L.STRINGLENGTH,MF,MC,Conversion_Factor,Alignments,Good_Alignments,H)) 
				{
					H.Status=SW_RECOVERED;
					Last_Mis=H.Mismatch;Top_Score=H.Score;
					return true;
				}
				else
					H.Status=UNRESOLVED_HIT;
					Last_Mis=H.Mismatch;
					return true;
			}
		}
		MF.Hit_Array_Ptr=Plus_Hits+1,MC.Hit_Array_Ptr=Minus_Hits+1;
	}
	if (Last_Mis == End_Mis || Sub_OptMIS== -1)
	{
		Top_Score= Calc_Top_Score(MF,MC,Top_BQScore,Last_Mis,L.STRINGLENGTH,Plus_Hits,Minus_Hits,R);
		Quality_Score=30;
	}
	Alignment Aln;Aln.Realigned=NO;
	if(Last_Mis!= -1)
	{
		H.Loc=UINT_MAX;
	}	
	else 
	{
	       	H.Status=UNMAPPED;
	}

	char Multi_Hit[MULTIHIT_COUNT],Unique='R';
	H.MH=Multi_Hit+1;Multi_Hit[0]='\t';
	bool RESCUE_MULTI=true;H.Sub_Opt_Hit=0;
	char SW_Status=NO_SW_SCAN;//In multi scan can differentialte hit sufficiantly..

	if(Last_Mis >=0)//If mismatch hits..
	{
		//if(!HEURISTIC) MAX_SW=INT_MAX;
		if(Plus_Hits+Minus_Hits==1)//Unique hits 
		{
			if(Unique_Hits(Plus_Hits,Minus_Hits,MF.Hit_Array[0],MC.Hit_Array[0]))
			{
				if(L.STRINGLENGTH<R.Real_Len)//if long read then extend on the fly..
				{
					if(Extend_With_SW(Plus_Hits,Minus_Hits,R,B,L.STRINGLENGTH,MF,MC,Conversion_Factor,Alignments,Good_Alignments,H))
					{
						Last_Mis=H.Mismatch;Top_Score=H.Score;Top_BQScore=H.BQScore;
						if(Close_Hits)
							H.Status=UNIQUEHIT;
						else
							H.Status=SHARP_UNIQUEHIT;
					}
					else
						assert(false);
				}
				else
				{
				Get_Info(MF,MC,L.STRINGLENGTH,H,Conversion_Factor,Alignments,Good_Alignments,R);//obtain hit details in H;
				Unique='U';Multi_Hit[0]=0;
				if(Close_Hits)
					H.Status=UNIQUEHIT;
				else
					H.Status=SHARP_UNIQUEHIT;
				}
			}
			else
			{
				if(Recover_With_SW(Plus_Hits,Minus_Hits,R,B,L.STRINGLENGTH,MF,MC,Conversion_Factor,Alignments,Good_Alignments,H)) 
				{
					H.Status=SW_RECOVERED;
					Last_Mis=H.Mismatch;Top_Score=H.Score;Top_BQScore=H.BQScore;
					return true;
				}
				Get_Info(MF,MC,L.STRINGLENGTH,H,Conversion_Factor,Alignments,Good_Alignments,R);//obtain hit details in H;
				Last_Mis=H.Mismatch;
				H.Sub_Opt_Hit=H.Loc;H.Status=MULTI_HIT;
				Quality_Score=0;
			}
		}
		else if (RESCUE_MULTI && Phred_Check_Multi(Plus_Hits,Minus_Hits,Last_Mis,R,L.STRINGLENGTH,MF,MC,Top_Score,Quality_Score,SW_Status))
		{
			if(SW_Status==NO_SW_SCAN && Unique_Hits(Plus_Hits,Minus_Hits,MF.Hit_Array[0],MC.Hit_Array[0]))
			{
				//assert(false);
				if(L.STRINGLENGTH<R.Real_Len)//if long read then extend on the fly..
				{
					if(Extend_With_SW(Plus_Hits,Minus_Hits,R,B,L.STRINGLENGTH,MF,MC,Conversion_Factor,Alignments,Good_Alignments,H))
					{
						Last_Mis=H.Mismatch;Top_Score=H.Score;Top_BQScore=H.BQScore;
						H.Status=MULTI_HIT;
					}
					else
						assert(false);
				}
				else
				{
					Get_Info(MF,MC,L.STRINGLENGTH,H,Conversion_Factor,Alignments,Good_Alignments,R);//obtain hit details in H;
					Unique='U';Multi_Hit[0]=0;
					H.Mismatch=Last_Mis;H.Status=MULTI_HIT;
					H.Score= Top_Score;
					H.BQScore=Top_BQScore;
				}
			}
			else
			{
				if(Recover_With_SW(Plus_Hits,Minus_Hits,R,B,L.STRINGLENGTH,MF,MC,Conversion_Factor,Alignments,Good_Alignments,H)) 
				{
					H.Status=SW_RECOVERED;
					Last_Mis=H.Mismatch;Top_Score=H.Score;Top_BQScore=H.BQScore;
				}
				else 
				{
					Get_Info(MF,MC,L.STRINGLENGTH,H,Conversion_Factor,Alignments,Good_Alignments,R);//obtain hit details in H;
					Last_Mis=H.Mismatch;
					H.Sub_Opt_Hit=H.Loc;H.Status=MULTI_HIT;
				}
			}
		}

		H.Mismatch=Last_Mis;
		H.Score= Top_Score;
		H.BQScore=Top_BQScore;
		H.Sub_BQScore=Sub_BQScore;
	}
	else
	{
		MAX_SW=MAX_SW_sav;
	}
	return true;
}

Alignment RealignFast(Hit_Info &  H,int StringLength, READ & R,int OFF,int Filter,bool Do_Filter)
{
	Alignment A;
	A.Clip_H=A.Clip_T=INT_MAX;
	if(!FASTSW)
	{
		A.Score=INT_MAX;
		return A; 
	}
	char Org_String[ORGSTRINGLENGTH];
	char *Quality,Rev_Qual[MAXTAG];
	Cigar_Info Cig_Info;


	char Real_String[R.Real_Len],*Current_Tag;
	if(R.Real_Len>=StringLength)
	{
		R.Tag_Copy[R.Real_Len]=0;
		R.Quality[R.Real_Len]=0;
		if(H.Sign=='+')
		{
			Read2Bin(Real_String,R.Tag_Copy,R.Real_Len);
			Quality=R.Quality;
		}
		else
		{
			assert(false);
		}
		StringLength=R.Real_Len;
	}
	Current_Tag=R.Tag_Copy;



	Get_Bases(H.Org_Loc,StringLength+INDELGAP,Org_String);
	for(int i=0;i<StringLength+INDELGAP;i++)
	{
		*(Org_String+i)="ACGT"[*(Org_String+i)];
	}
	*(Org_String+StringLength+INDELGAP)=0;
	A=Fast_SW(Org_String,Current_Tag,Quality,OFF,'+');
	if(A.SW_Score<match*9*StringLength/10) 
		A.Score=INT_MAX;
	A.Loc=H.Org_Loc+1;
	A.Sign=H.Sign;
	return A;
	
}

Alignment RealignFastMinus(Hit_Info &  H,BATREAD & Read,int StringLength, READ & R,int OFF,int Filter,bool Do_Filter)
{
	Alignment A;
	A.Clip_H=A.Clip_T=INT_MAX;
	if(!FASTSW)
	{
		A.Score=INT_MAX;
		return A; 
	}
	char Org_String[ORGSTRINGLENGTH];
	char *Quality;

	char Real_String[R.Real_Len],*Current_Tag;
	if(R.Real_Len>=StringLength)
	{
		R.Tag_Copy[R.Real_Len]=0;
		R.Quality[R.Real_Len]=0;
		if(H.Sign=='+')
		{
			assert(false);
		}
		else
		{
			Quality=R.Quality;
			H.Org_Loc++;
		}
		StringLength=R.Real_Len;
	}
	Current_Tag=R.Tag_Copy;

	Get_Bases(H.Org_Loc-INDELGAP,StringLength+INDELGAP,Org_String);//get ref bases and convert..
	char Org_Rev[StringLength+INDELGAP+1];
	for(int i=StringLength+INDELGAP-1,j=0;i>=0;i--,j++)
	{
		Org_Rev[j]="ACGT"[3-Org_String[i]];
	}
	Org_Rev[StringLength+INDELGAP]=0;

	A=Fast_SW(Org_Rev,Current_Tag,Quality,OFF,'-');
	if(A.SW_Score<match*9*StringLength/10) 
		A.Score=INT_MAX;
	A.Loc+=H.Org_Loc;
	A.Sign=H.Sign;
	return A;
}

Alignment RealignX(Hit_Info &  H,BATREAD & Read,int StringLength, READ & R,bool Dont_Push_To_Q,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,char* Cigar,int & Clip_H,int & Clip_T,int & Filter,bool Do_Filter)
{
	s_align* Aln;
	Alignment A;
	A.Clip_H=A.Clip_T=INT_MAX;
	char Org_String[ORGSTRINGLENGTH];
	char *Quality,Rev_Qual[MAXTAG];
	Cigar_Info Cig_Info;


	char Real_String[R.Real_Len],*Current_Tag;
	if(R.Real_Len>=StringLength)
	{
		R.Tag_Copy[R.Real_Len]=0;
		R.Quality[R.Real_Len]=0;
		if(H.Sign=='+')
		{
			Read2Bin(Real_String,R.Tag_Copy,R.Real_Len);
			Quality=R.Quality;
		}
		else
		{
			Read2RevCBin(Real_String,R.Tag_Copy,R.Real_Len);
			Reverse_Quality(Rev_Qual,R,R.Real_Len);
			H.Org_Loc-=(R.Real_Len-StringLength)+INDELGAP-1;
			//Offset=R.Real_Len-StringLength;
			Quality=Rev_Qual;
		}
		StringLength=R.Real_Len;
	}
	Current_Tag=Real_String;



	int Jump=0;if(H.Sign=='-') Jump=0 +JUMP;
	s_profile* p = ssw_init((int8_t*)Current_Tag, StringLength, mata, n, 1);
	Get_Bases(H.Org_Loc+Jump,StringLength+INDELGAP,Org_String);
	Aln=mengyao_ssw_core(Org_String,StringLength, Current_Tag,StringLength+INDELGAP,Filter,0/*DP*/, p);
	//if(Aln->score1 >= ACC_SCORE)
	if(Aln->score1 >= Filter)
	{
		A.Clip_H=Aln->read_begin1;A.Clip_T=0;
		if(Aln->read_end1!=StringLength-1) A.Clip_T=StringLength-1-Aln->read_end1;

		A.SW_Score=Aln->score1;
		ssw_cigar_processQ(Aln,Cig_Info,Org_String,Aln->ref_begin1,Current_Tag,Aln->read_begin1,StringLength,Quality,A.Cigar,A.Clip_H,A.Clip_T);
		A.Loc=H.Org_Loc+Jump+Aln->ref_begin1;//+Offset;
		A.Score= -Cig_Info.Score;
		A.QScore=Cig_Info.QScore;
		A.BQScore=Cig_Info.BQScore;
		A.Mismatch=Cig_Info.Mis;
		A.Indel=Cig_Info.Indel_Count;
		A.Sign=H.Sign;
		A.Realigned=1;
		if(Do_Filter && BOOST && Aln->score1 >Filter) 
		{
			Filter=Aln->score1;
		}
		if(!Dont_Push_To_Q)
		{
			Good_Alignments.push(A);
		}
	}
	else
	{
		A.Loc=INT_MAX;
		if(BOOST && Aln->score1>=0)
		{
			A.Score= -INT_MIN;
			A.Realigned=1;
			if(!Dont_Push_To_Q)
			{
				Good_Alignments.push(A);
			}
		}
	}
	align_destroy(Aln);
	init_destroy(p); 
	return A;
	
}

inline void Copy_A_to_H(Alignment & A,Hit_Info & H)
{
	H.Mismatch=A.Mismatch;
	H.Indel=A.Indel;
	H.Score= -A.Score;
	H.QScore= A.QScore;
	H.BQScore= A.BQScore;
	H.SW_Score= A.SW_Score;
	H.Loc=H.Org_Loc= A.Loc-1;
	H.Clip_H=A.Clip_H;H.Clip_T=A.Clip_T;
	if(A.Cigar[0])
		strcpy(H.Cigar,A.Cigar);
	else
	{
		H.Cigar[0]=0;		
		H.Cigar[1]='Z';
		A.Cigar[1]='Z';
	}

}

bool Extend_With_SW(int Plus_Hits,int Minus_Hits,READ & R,BATREAD & BR, int StringLength,MEMX & MF,MEMX & MC,unsigned Conversion_Factor,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,Hit_Info & H,bool Mismatch_Scan)
{
	int Hits=0,Clip_T,Clip_H;
	int Filter=0;
	if(Plus_Hits)
	{
		if(Mismatch_Scan) assert(!Minus_Hits);
		for(int i=0;MF.Hit_Array[i].Start;i++)
		{
			SARange SA=MF.Hit_Array[i];
			assert(SA.End-SA.Start>=0); Hits+=(SA.End-SA.Start+1);
			Alignment A;
			A.Clip_H=A.Clip_T=INT_MAX;
			for (int j=0;j<=(SA.End-SA.Start);j++)
			{
				unsigned Loc=SA2Loc(SA,j,Conversion_Factor);
				H.Loc =H.Org_Loc= Loc;
				H.Sign='+';
				Hit_Info H2=H;
				A=RealignFast(H2,StringLength,R,INT_MAX,Filter,true);
				if(A.Score==INT_MAX)
				{
					RealignX(H,BR,StringLength,R,false,Good_Alignments,NULL,Clip_H,Clip_T,Filter,true);
				}
				else
				{
					assert(A.Score<=0);
					A.Realigned=1;A.Clip_T=A.Clip_H=0;
					Good_Alignments.push(A);
				}
			}
		}
	}
	if(Minus_Hits)
	{
		if(Mismatch_Scan) assert(!Plus_Hits);
		for(int i=0;MC.Hit_Array[i].Start;i++)
		{
			SARange SA=MC.Hit_Array[i];
			Hits+=(SA.End-SA.Start+1);
			assert(SA.End-SA.Start>=0);
			for (int j=0;j<=(SA.End-SA.Start);j++)
			{

				unsigned Loc=SA2Loc(SA,j,Conversion_Factor);
				H.Loc =H.Org_Loc= Loc;H.Sign='-';
				Hit_Info H2=H;
				H2.Loc =H2.Org_Loc= Loc-(R.Real_Len-SEEDSIZE);
				Alignment B=RealignFastMinus(H2,BR,StringLength,R,INT_MAX,Filter,true);
				if(B.Score==INT_MAX)
				{
					Alignment A=RealignX(H,BR,StringLength,R,false,Good_Alignments,NULL,Clip_H,Clip_T,Filter,true);
				}
				else
				{
					assert(B.Score<=0);
					B.Realigned=1;B.Clip_T=B.Clip_H=0;
					Good_Alignments.push(B);
				}
			}
		}
	}
	if(!Good_Alignments.empty())
	{
		Alignment A=Good_Alignments.top();
		assert(!Mismatch_Scan || Good_Alignments.size()==1);//there is a unique best alignment
		Copy_A_to_H(A,H);
		Alignments.push(A);
		return true;
	}
	else
	{
		return false;
	}
}
//Tries to recover sw hits from multi hits. returns false if unmappable..
bool Recover_With_SW(int Plus_Hits,int Minus_Hits,READ & R,BATREAD & BR, int StringLength,MEMX & MF,MEMX & MC,unsigned Conversion_Factor,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,Hit_Info & H)
{
	int Hits=0;
	int Filter=0;
	if(Plus_Hits+Minus_Hits==1)
	{
		if(Plus_Hits)//Skip past sentinel..
		{
			if(MF.Hit_Array_Ptr>1)
				for(int i=1;MF.Hit_Array[i+1].Start;i++)
				{
					MF.Hit_Array[i]=MF.Hit_Array[i+1];
				}
		}
		else
		{
			if(MC.Hit_Array_Ptr>1)
				for(int i=1;MC.Hit_Array[i+1].Start;i++)
				{
					MC.Hit_Array[i]=MC.Hit_Array[i+1];
				}
		}
	}

	//might be suboptimal only hits..
	if(!Plus_Hits)
	{
		if(MF.Hit_Array_Ptr>1)
			for(int i=0;MF.Hit_Array[i+1].Start;i++)
			{
				MF.Hit_Array[i]=MF.Hit_Array[i+1];
			}
		Plus_Hits++;
	}
	if(!Minus_Hits)
	{
		if(MC.Hit_Array_Ptr>1)
			for(int i=0;MC.Hit_Array[i+1].Start;i++)
			{
				MC.Hit_Array[i]=MC.Hit_Array[i+1];
			}
		Minus_Hits++;
	}

	int Clip_T,Clip_H;
	if(Plus_Hits)
	{
		for(int i=0;MF.Hit_Array[i].Start;i++)
		{
			SARange SA=MF.Hit_Array[i];
			assert(SA.End-SA.Start>=0);
			Hits+=(SA.End-SA.Start+1);
			for (int j=0;j<=(SA.End-SA.Start);j++)
			{
				unsigned Loc=SA2Loc(SA,j,Conversion_Factor);
				H.Loc =H.Org_Loc= Loc;
				H.Sign='+';
				Hit_Info H2=H;
				Alignment A;
				A=RealignFast(H2,StringLength,R,INT_MAX,Filter,true);
				if(A.Score==INT_MAX)
				{
					RealignX(H,BR,StringLength,R,false,Good_Alignments,NULL,Clip_H,Clip_T,Filter,true);
				}
				else
				{
					assert(A.Score<=0);
					A.Realigned=1;
					Good_Alignments.push(A);
				}
			}
		}
	}
	if(Minus_Hits)
	{
		for(int i=0;MC.Hit_Array[i].Start;i++)
		{
			SARange SA=MC.Hit_Array[i];
			Hits+=(SA.End-SA.Start+1);
			assert(SA.End-SA.Start>=0);
			for (int j=0;j<=(SA.End-SA.Start);j++)
			{
				unsigned Loc=SA2Loc(SA,j,Conversion_Factor);
				H.Loc =H.Org_Loc= Loc;
				H.Sign='-';
				Hit_Info H2=H;
				H2.Loc =H2.Org_Loc= Loc-(R.Real_Len-SEEDSIZE);
				Alignment B=RealignFastMinus(H2,BR,StringLength,R,INT_MAX,Filter,true);
				if(B.Score==INT_MAX)
				{
					Alignment A=RealignX(H,BR,StringLength,R,false,Good_Alignments,NULL,Clip_H,Clip_T,Filter,true);
					//RealignX(H,BR,StringLength,R,false,Alignments,Good_Alignments,NULL,Clip_H,Clip_T,Filter,true);
				}
				else
				{
					assert(B.Score<=0);
					B.Realigned=1;B.Clip_T=B.Clip_H=0;
					Good_Alignments.push(B);
				}
			}
		}
	}
	if(!Good_Alignments.empty())
	{
		Alignment A=Good_Alignments.top();
		if(Good_Alignments.size()==1)//there is a unique best alignment
		{
			A.Sub_Opt_Score=INT_MAX;
		}
		else
		{
			if(PAIRED)
			{
				Good_Alignments.pop();
				Alignments=Good_Alignments;
				Alignment B=Good_Alignments.top();
				//Alignments.push(B);
				Good_Alignments.push(A);//
				A.Sub_Opt_Score= -B.Score;
			}
			else
			{
				Good_Alignments.pop();
				Alignment B=Good_Alignments.top();
				Alignments.push(B);
				Good_Alignments.push(A);//
				A.Sub_Opt_Score= -B.Score;
			}
		}
		Copy_A_to_H(A,H);
		Alignments.push(A);
		return true;
	}
	else
	{
		return false;
	}
}

unsigned SA2Loc(SARange S,int Pos,unsigned Conversion_Factor)
{
	if (S.Start==S.End) 
	{
		return S.Start;
	}
	else
	{
		assert (S.Start);
		return Conversion_Factor-BWTSaValue(revfmi,S.Start+Pos);
	}
	
}

void Launch_Threads(int NTHREAD, void* (*Map_t)(void*),Thread_Arg T)
{
	Threading* Thread_Info=(Threading*) malloc(sizeof(Threading)*NTHREAD);
	int Thread_Num=0;
	pthread_attr_t Attrib;
	pthread_attr_init(&Attrib);
	pthread_attr_setdetachstate(&Attrib, PTHREAD_CREATE_JOINABLE);

	for (int i=0;i<NTHREAD;i++)
	{
		T.ThreadID=i;
		Thread_Info[i].Arg=T;
		//if(!(Thread_Info[i].r=pthread_create(&Thread_Info[i].Thread,NULL,Map_t,(void*) &Thread_Info[i].Arg))) Thread_Num++;
		Thread_Info[i].r=pthread_create(&Thread_Info[i].Thread,&Attrib,Map_t,(void*) &Thread_Info[i].Arg);
		if(Thread_Info[i].r) {printf("Launch_Threads():Cannot create thread..\n");exit(-1);} else Thread_Num++;
	}
	printf("%d Threads runnning ...\n",Thread_Num);
	pthread_attr_destroy(&Attrib);

	for (int i=0;i<NTHREAD;i++)
	{
		pthread_join(Thread_Info[i].Thread,NULL);
	}
}


void Set_Affinity()
{
	cpu_set_t Set;
	CPU_ZERO(&Set);

	if(sched_getaffinity(0,sizeof(cpu_set_t),&Set)<0)
	{
		printf("Affinity could not be get..\n");
	}
	else
	{
		for (int i=0;i<CPU_SETSIZE;i++)
		{
			if(CPU_ISSET(i,&Set))
			{
				printf("Bound to %d\n",i);
				CPU_ZERO(&Set);
				CPU_SET(i, &Set);
				if(sched_setaffinity(0, sizeof(Set), &Set)<0)
				{
					printf("Affinity could not be set..\n");
				}
				return;
			}
		}
	}

}


bool Map_One_Seg(READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MF,MEMX & MC,LEN & L,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,Pair* & Pairs,bool PRINT,Hit_Info & H,int & Quality_Score,int Segment_Length,int SEG_SIZE,int SHIFT_SEG_SIZE)
{
		READ RTemp=R;BATREAD BTemp=B;
		char T_EOS=R.Tag_Copy[SEG_SIZE+1],T_EOQ=R.Quality[SEG_SIZE+1];
		if(SEG_SIZE)
		{
			//for(int i=0;i<=SHIFT_SEG_SIZE;i++)
			for(int i=0;i<=SEG_SIZE;i++)
			{
				assert((R.Tag_Copy[i+SHIFT_SEG_SIZE]>='A' && R.Tag_Copy[i+SHIFT_SEG_SIZE]<='t'));
				R.Tag_Copy[i]=R.Tag_Copy[i+SHIFT_SEG_SIZE];R.Quality[i]=R.Quality[i+SHIFT_SEG_SIZE];
			}

			R.Tag_Copy[SEG_SIZE+1]=0;
		}

		R.Real_Len=0;
		for(;R.Tag_Copy[R.Real_Len]!=0 && R.Tag_Copy[R.Real_Len]!='\n';R.Real_Len++);
		if(Segment_Length)
			R.Real_Len=Segment_Length;
		Conversion_Factor=revfmi->textLength-L.STRINGLENGTH;
		IN.Positive_Head=R.Tag_Copy;
		R.Read_Number=Actual_Tag;

		B.StringLength=SEG_SIZE;
		Process_Read(R,B,MF,MC);
		MF.Strand='+';MF.Larger_Than_Ten=0;MC.Strand='-';MC.Larger_Than_Ten=0;
		MF.Extend=MC.Extend=false;

		FreeQ(Alignments);
		FreeQ(Good_Alignments);
		if (R.NCount > NCOUNT) 
		{
			R.Tag_Copy[SEG_SIZE+1]=T_EOS;R.Quality[SEG_SIZE+1]=T_EOQ;
			R=RTemp;B=BTemp;
			return false;
		}

		int Last_Mis;
		int Head_Top_Count,Tail_Top_Count;//# of top hits in H/T
		int Head_Subopt_Count,Tail_Subopt_Count;//# of top hits in H/T
		Quality_Score=QUAL_UNMAPPED;H.Status=UNMAPPED;


		bool Hit_Found=Do_Mismatch_Scan(MF,MC,L,fwfmi,revfmi,0,Inter_MM,Last_Mis,Head_Top_Count,H,Quality_Score,R,B,Mishit_File,Conversion_Factor,Alignments,Good_Alignments);
		R.Tag_Copy[SEG_SIZE+1]=T_EOS;R.Quality[SEG_SIZE+1]=T_EOQ;
		R=RTemp;B=BTemp;
		return Hit_Found;
}

