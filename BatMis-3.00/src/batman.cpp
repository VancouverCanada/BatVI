//{-----------------------------  INCLUDE FILES  -------------------------------------------------/
#define  NDEBUG
#include "config.h"
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <signal.h> 
#include "assert.h"
//#include <dvec.h>
#include <getopt.h>
#include "zlib.h"
#include <queue>
#include <ctype.h>
#include <map>

extern "C" 
{
	#include "iniparser.h"
	#include <time.h>
	#include "MemManager.h"
	#include "MiscUtilities.h"
	#include "TextConverter.h"
	#include "BWT.h"
}
//}-----------------------------  INCLUDE FILES  -------------------------------------------------/
using namespace std;

//{-----------------------------  DEFINES  -------------------------------------------------/
//#define SOLID_DBG
#define ARRAY_BOUND_BD 0
#define SACUT INT_MAX//20000
//#define EXTCUT  INT_MAX//250
#define TAB	1
#define FQ	2
#define FA	3
#define TWOFILE	4
#define GISMODE 100
#define DEFAULT 0
#define BTS_PER_LOC 8
#define BITMASK 255 //2^BTS_PER_LOC -1
#define DEEP 1
#define PAIREND 2
#define PAIR_END_SEPERATOR '\t'
#define MAXDES 500 
#define MAXTAG 256
#define INDELMARK (MAXTAG-1)
#define INSERTMARK 63
#define DELETEMARK 70
#define EXTRA 0//Stuff like cr/lf at string seperators in the input file
#define INTEGERSIZE 8 //integer size
#define PACKEDBITSIZE 2
#define FALSE 0
#define TRUE 1
#define NOMISMATCHES 100 
#define MAX_MISMATCHES_BOUND 16 //upper boud for mismatch number....
#define REVERSE 1
#define FORWARD 0
#define LEFTPASS 1
#define LEFTPASS1 11
#define LEFTPASS1Y 110
#define LEFTPASSX 111
#define LEFTPASSY 114
#define LEFTPASSXY 113
#define RIGHTPASS 2
#define RIGHTPASS1 21
#define RIGHTPASS1Y 31
#define RIGHTPASSX 121
#define RIGHTPASSXY 101
#define SAINTERVAL 8
#define USER 0
#define SENDMAIL 1
#define MAIL 2


#define START_OF_MARK 9//18  //Start of scanning the tag for the pruning of one mismatch...
//#define BRANCHTHRESHOLD 80 //30 //Threshold at which to check the BWT instead of branching
#define BRANCHTHRESHOLD 0 //30 //Threshold at which to check the BWT instead of branching
//}-----------------------------  DEFINES  -------------------------------------------------/
//{-----------------------------  STRUCTS  -------------------------------------------------/
struct FastaSeq 
{
      int   len; /* the actual string length of seq */
      char* seq; /* the sequence itself */
};

typedef struct gz_stream {
    z_stream stream;
    int      z_err;   /* error code for last stream operation */
    int      z_eof;   /* set if end of input file */
    FILE     *file;   /* .gz file */
    Byte     *inbuf;  /* input buffer */
    Byte     *outbuf; /* output buffer */
    uLong    crc;     /* crc32 of uncompressed data */
    char     *msg;    /* error message */
    char     *path;   /* path name for debugging only */
    int      transparent; /* 1 if input file is not a .gz file */
    char     mode;    /* 'w' or 'r' */
    z_off_t  start;   /* start of compressed data in file (header skipped) */
    z_off_t  in;      /* bytes into deflate or inflate */
    z_off_t  out;     /* bytes out of deflate or inflate */
    int      back;    /* one character push-back */
    int      last;    /* true if push-back is last character */
} gz_stream;


struct Header
{
	char ID[3] ;
	unsigned MAXHITS;
	char FILETYPE;
	char HITMODE;
	char IGNOREHEAD;
	char Index_Count;
	int Tag_Length;
	char Print_Desc;//print the tag desc?
}__attribute__((__packed__));

struct SARange
{

	unsigned long Start;
	unsigned long End;
	int Level;//at what relative node are we?
	char Mismatches;//number of mismatches at this level
	unsigned char Skip;
	int Tag;//Tag number
	unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];//BTS_PER_LOC|BTS_PER_LOC|...
	unsigned Mismatch_Char;//2|2|...
	
};

struct Range
{
	unsigned Start;
	unsigned End;
	int Label;//Final Label of the range...
};

struct Mismatches_Record
{

	int 	Gap;
	unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];//BTS_PER_LOC|BTS_PER_LOC|...
	unsigned Mismatch_Char;//2|2|...


}__attribute__((__packed__));

struct Mismatches_Record_GIS
{

	unsigned Mismatch_Char;//2|2|...
	unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];//BTS_PER_LOC|BTS_PER_LOC|...

}__attribute__((__packed__));


struct Output_Record
{
	unsigned Tag;
	unsigned  Start;
	char Index;
	unsigned char Skip;
	char Mismatches;
	int Gap;
}__attribute__((__packed__));

struct Branches
{
	long Is_Branch [4];
};
//}-----------------------------  STRUCTS  -------------------------------------------------/*

//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*

void Convert_To_Reverse(SARange &Tag);
void Build_Tables();
void Allocate_Memory();
void Init_Variables();
void Load_Indexes();
void Verbose(FILE* out);
void Open_Files();
void File_OpenZ(const char* File_Name,const char* Mode,gzFile & Handle);
void Build_Preindex_Forward(Range Range, int Level, int Bit);
void Build_Preindex_Backward(Range Range, int Level, int Bit);
void Parse_Command_line(int argc, char* argv[]);
void Read_INI();
void Zig_Zag();
void Seed_Scan();
void Extend_Left_Half();
void Print_Hits();
void High_Mismatch_Scan();
void Quality_Calc();

void Left_To_Right();
void Right_To_Left();
void Left_To_RightX();
void Right_To_LeftX();
void Left_To_RightY();
void Right_To_LeftY();

void One_Branch(struct SARange Tag,BWT *fmi);
void Branch_Detect (const struct SARange Tag,BWT *fmi,int Start);
void Branch_Detect_Backwards (const struct SARange Tag,BWT *fmi,int Start);

void Print_LocationX (struct SARange & Tag);
void Show_Progress(unsigned Percentage);
void Get_SARange_Fast( long New_Char,struct SARange & Range,BWT *fmi);

void Reverse(struct SARange & Tag,int Start,int StringLength);
void Backwards(struct SARange & Tag,int Start,int StringLength);

void Search_Exact(struct SARange & Tag,int Start,int StringLength,BWT *fmi);
void Search_Forwards_Exact(struct SARange & Tag,int Start,int StringLength,BWT *fmi);
char Search_Forwards_ExactX(SARange & Tag, int Start,int StringLength);
void Search_Forwards(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_ForwardsX(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Forwards_Indel(SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Forwards_0X(struct SARange & Tag,int Start,int StringLength,BWT *fmi);
void Search_Forwards_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Forwards_OneSAX(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Forwards_OneSA_Indel(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_X01(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_X01_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_XL01(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_XL01_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_01X(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_01X_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_01LX(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_01LX_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_11X(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Half_Tag_11X(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);

void Search_Brute(int Count1,int Count2,int Start,int StringLength,BWT *fmi);
void Emailer(int argc,char* argv[]);
void Search_X11(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Half_Tag_X11(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Backwards_Exact(struct SARange & Tag,int Start,int StringLength,BWT *fmi);
void Search_Backwards_Exact_X0(struct SARange & Tag,int Start,int StringLength,BWT *fmi);
void Search_Backwards(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_BackwardsX(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Backwards_Indel(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Backwards_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Backwards_OneSAX(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Backwards_X10(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Backwards_X10_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Backwards_XL10(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_Backwards_XL10_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_10X(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_10X_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_10LX(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void Search_10LX_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi);
void fprintfX(void* Handle,char* Format, char* String);
void Sig_Handle(int signo);

//----------DUST stuff ------------------------
int dust(FastaSeq* fa);
int wo(int len, char* s, int* beg, int* end);
void wo1(int len, char* s, int ivv); 
inline void addDustRgn(FastaSeq* fa, int f, int t); 
//----------DUST stuff ------------------------

SARange Get_SARange( long New_Char,struct SARange Range,BWT *fmi);
char Guess_Orientation();
FILE* File_Open(const char* File_Name,const char* Mode);
BWT* initFMI(const char* BWTCodeFileName,const char* BWTOccValueFileName);
size_t fwriteX ( const void * ptr, size_t size, size_t count, void * stream );
//}-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*


//{---------------------------- GLOBAL VARIABLES -------------------------------------------------
char debug;
char Patternfile_Count;
char Head_Printed=TRUE;
char Hit_In;
char AMBI=FALSE;
char Tag_Number;
char HIGHSCAN;//inside high mismatch scan..
char Random_ArrayS[]="tatacgataggacaatgtcttcgaagcccacgcggtaagccggtcattgcggttgtgcgaacactatcagcctcgctgcatggttaccctgggtggataggacgtttgcccgacattttgacacgcataaaaggtctgtagtgggggtggcacaccataaaccctggggcggctccacgatcgtaaaatcctgcgatctg";
char Random_ArrayC[]="30301203022010032313312002111012122300211223103321223323212001013031021131213210322330111322232203022012333211120103333201012103000022313230232222232210101103000111322221221311012031230000311321203232";
int Random_Pointer=0;
int Last_Mismatch_Written;
int Argc;
char **Argv;
char Mismatches_InN[15];
char Dustiness[MAXTAG+1];
char Mis_Bound[]={0,0,0,0,0,0,0,2,2,2,3,3,3,3,3,3};
char Mismatch_Check[MAXTAG+1];
char Mismatches_Reduced=FALSE;

FastaSeq fa;
FILE* Input_FileO;
//FILE* Output_File;
FILE* Output_FileF;
gzFile Output_FileG;
void* Output_File;
FILE* Unique_File;
FILE* Mishit_File;
FILE* Log_File;
FILE* Ambiguous_File;
FILE* BLAST_File;
//gzFile Input_File;
//gzFile Mate_File;
FILE* Input_File;
FILE* Mate_File;

Output_Record Record;
Header Header;
MMPool *mmPool;
Mismatches_Record Mismatches;
Mismatches_Record_GIS MismatchesGIS;
BWT *fwfmi,*revfmi;
off64_t File_Size;
time_t Start_Time,End_Time;

SARange Cache_SF[MAXTAG+1];
SARange* BMHStack;//[2*STRINGLENGTHMAX];
SARange* FSHStack;
SARange* FSHStackX0X;
SARange* FSSStack;
SARange* FSSStackX;
SARange* BMStack;
SARange* BMStackX;
SARange* BMStack_X11;
SARange* BMStack_X11H;
SARange* PSBStack;

SARange* Exact_Match_Left;
SARange* Exact_Match_Right;
SARange* Exact_Match_Forward;
SARange* Exact_Match_ForwardF;
SARange* Exact_Match_ForwardC;
SARange* Exact_Match_Backward;
SARange* Left_Mishits;//stores mismatches in first half
SARange* Left_MishitsF;//stores mismatches in first half
SARange* Left_MishitsC;//stores mismatches in first half
SARange* Right_Mishits;//stores mismatches in first half
SARange* Right_MishitsF;//stores mismatches in first half
SARange* Right_MishitsC;//stores mismatches in first half
SARange* Mismatches_Forward;//stores possible 2-mismatches
SARange* Mismatches_ForwardF;
SARange* Mismatches_ForwardC;
SARange* Mismatches_Backward;//stores possible 2-mismatches
SARange* Mismatches_BackwardF;//stores possible 2-mismatches
SARange* Mismatches_BackwardC;//stores possible 2-mismatches
SARange* Two_Mismatches_At_End_Forward;//stores 2-mismatches whose last mismatch occurs at the last nuce..
SARange* Two_Mismatches_At_End_ForwardF;//stores 2-mismatches whose last mismatch occurs at the last nuce..
SARange* Two_Mismatches_At_End_ForwardC;//stores 2-mismatches whose last mismatch occurs at the last nuce..
SARange* Two_Mismatches_At_End;//stores 2-mismatches whose last mismatch occurs at the last nuce..
SARange* Two_Mismatches_At_EndF;//stores 2-mismatches whose last mismatch occurs at the last nuce..
SARange* Two_Mismatches_At_EndC;//stores 2-mismatches whose last mismatch occurs at the last nuce..
SARange* Possible_20;//stores possible 2-mismatches in right half... 
SARange* Possible_20F;//stores possible 2-mismatches in right half... 
SARange* Possible_20C;//stores possible 2-mismatches in right half... 
SARange* Possible_02;//stores possible 2-mismatches in right half... 
SARange* Possible_02F;//stores possible 2-mismatches in right half... 
SARange* Possible_02C;//stores possible 2-mismatches in right half... 

SARange Branch_Ranges[4];
SARange Temp_Branch_Ranges[4];
SARange Temp_Branch_RangesFX[4];
SARange Temp_Branch_RangesBX[4];
SARange Temp_Branch_Ranges2[4];


int Left_Mishits_Pointer;
int Right_Mishits_Pointer;
int Mismatches_Forward_Pointer;
int Mismatches_Backward_Pointer;
int Two_Mismatches_At_End_Pointer; 
int Two_Mismatches_At_End_Forward_Pointer=0;
int Possible_20_Pointer=0;
int Possible_02_Pointer=0;
int Best_Pos;//start of best quality position..
int FMIndex;//keeps track of which fm-index to use...
int One_Mismatch_Already;
int EXTCUT=INT_MAX;//251
int c;//to hold a character temporarily
int word = 3;//dust 
int window = 64; 
int window2 = 32;
int Actual_Tag;
int Stat_Size;
int NCount;
int LOOKUPSIZE=6;
int PLOOKUPSIZE,PRLOOKUPSIZE;
int HITMODE = DEFAULT;
int STRINGLENGTH=36;// 36//36//6+EXTRA//6+EXTRA//36
int TAG_COPY_LEN;
int STRINGLENGTHl=36;// 36//36//6+EXTRA//6+EXTRA//36
int STRINGLENGTHr=36;// 36//36//6+EXTRA//6+EXTRA//36
int STRINGLENGTHt=36;// 36//36//6+EXTRA//6+EXTRA//36
int STRINGLENGTHt1=36;// 36//36//6+EXTRA//6+EXTRA//36
int STRINGLENGTHt2=36;// 36//36//6+EXTRA//6+EXTRA//36
int PAIR_LENGTH_RIGHT,PAIR_LENGTH_LEFT;
int HALFSTRINGLENGTH=18;// 18
int QUARTERSTRINGLENGTH=9;// 9
int LH,RH,LHQL,LHQR,RHQL,RHQR;
int LHt,RHt,LHQLt,LHQRt,RHQLt,RHQRt;
int LHt1,RHt1,LHQLt1,LHQRt1,RHQLt1,RHQRt1;
int LHt2,RHt2,LHQLt2,LHQRt2,RHQLt2,RHQRt2;
int LHl,RHl,LHQLl,LHQRl,RHQLl,RHQRl;
int LHr,RHr,LHQLr,LHQRr,RHQLr,RHQRr;
int PLH,PRH,PLHQL,PLHQR,PRHQL,PRHQR;
int PRLH,PRRH,PRLHQL,PRLHQR,PRRHQL,PRRHQR;
int MAX_HIGH_QUALITY_LENGTH=0;
int MAX_MISMATCHES_H=0;
int IGNOREHEAD=0;
int FORCELENGTH=0;
int THRESHOLD012=50;
int THRESHOLD01=40;
int SLIDER=0;
int SLIDERBLAST=0;
int SLIDERHIT=0;
int Possible_04_Pointer,Possible_40_Pointer,Possible_50_Pointer;
int ARRAY_BOUND,END_BOUND;
int NBOUND = 0;//200;//numbar of N's allowed...
int SEEDSIZE;
int CHOPWHITE=0;
int Delta;
int Max_Seeked;
int Ext_Scan=0;
int mv, iv, jv;//DUST variables...

unsigned char N[500];
unsigned Write_Buf_Size;
unsigned FWDInverseSA0;
unsigned Branch_Characters[4],Temp_BCR[4],Temp_BC[4],Temp_BCB[4],Temp_BC1[4],Temp_BC2[4],Temp_BCX[4],Temp_BCX1[4];//counts the characters in a SARange.
unsigned Hits,Hits1=0,Total_Hits=0,First_Pass_Hits=0,Tags_From_Head=0;
unsigned Total_Tags=0;
unsigned Tags_Read;
unsigned Number_of_Tags,Average_Length;
unsigned Forward_Start_Lookup[4],Forward_End_Lookup[4],Backward_Start_Lookup[4],Backward_End_Lookup[4];
unsigned* Forward_Start_LookupX;
unsigned* Forward_End_LookupX;
unsigned* Backward_Start_LookupX;
unsigned* Backward_End_LookupX;
unsigned SOURCELENGTH;
unsigned MAXHITS =1;
unsigned DISKBUFFERSIZE =1000;
unsigned MAX_TAGS_TO_PROCESS=FALSE; 

char Description[MAXDES+1];
char Temp_Description[MAXDES+1];
char Plus[MAXDES+1];
char Original[MAXTAG+1];
char Tag_Copy[MAXTAG+1];
char Quality_Copy[MAXTAG+1];
char Complement[MAXTAG+1];
char Mishits_Tag[MAXTAG+1];
char Quality[MAXTAG];
char Fasta_String[MAXTAG+1];
char Mismatch_Init[MAX_MISMATCHES_BOUND];
char Quality_Bound;
char Maximum_Quality;
char Sum;
char Print_Header;//For single end tags. Determines if the header is printed.
char New_Record;
char ZigZag;
char In_Mismatch;
char In_MismatchX;
unsigned Extentions;
unsigned SATot,S;
char Last_In_Mis;
char Last_Mismatch;
unsigned short Stats1[7];
unsigned short Stats2[7];
unsigned short* Stats=Stats1;
char DUST=FALSE;
char GIS=TRUE;
char STATMODE=FALSE;//TRUE;
char* Random_Array=Random_ArrayS;
char* Write_Buffer_Ptr;
char* Last_Write_Buffer;
char* Write_Buffer;

char Minimum_Quality,Minimum_Quality_Location;
char Translated_String[MAXTAG +1];//STRINGLENGTHMAX+1];
char Char_To_Code[256];
//char Char_Count[256];
char Char_To_CodeC[256];
char Quality_Count[256];
char All_Zero[256];
char Temp_Char_Array[MAX_MISMATCHES_BOUND];
char Direction=0;//Direction =0,1 scan org, then comple. Direction =1,2 scan compl, then org
char First_Pass=TRUE;
char First_Pass_Printed=FALSE;
char Rollover_Step=FALSE;
char Tag_Stat_Bad=FALSE;
char TagProcessed=FALSE;
char Tag_Printed=FALSE;
char Extending_Tag=FALSE;
char Larger_Than_Ten=FALSE;
unsigned Text=0,T=0;
//char TAG_GUESS=TRUE;
char TAG_GUESS=FALSE;
char SCANCOMPONLY=FALSE;
char ROLLOVER=FALSE;
char FILTER_AMBIGUOUS=FALSE;//decides whether to filter the ambiguous entries to ambiguous.fq
char INPUT_ZIPPED;
char OUTPUT_ZIPPED=FALSE;
char PRINT_MISHITS=FALSE;
char PRINT_DESC=TRUE;//FALSE;
char INDELSCAN=FALSE;
char MAX_MISMATCHES=4;
char MAX_ORG;
char ORG_MISMATCHES=4;
char NISMISMATCH=FALSE;//TRUE;
char NPOLICY=FALSE;
char MAX_MISMATCHES_IN_EXTENSION=5;
char USEQUALITYH=TRUE;
char NOQFILT=FALSE;
char RANDOMLETTER=FALSE;
char LEAST_MISMATCH=FALSE;

char DIRECTIONS = 2;//DIRECTION=1 scan one leg,=2 scan original then comp, =3 scan comple, then original.
char SEED=FALSE;
char SOLID=FALSE;
char FORCESOLID=FALSE;
char BWTFILE_DEFAULT[] = "genome.bwt"; 
char OCCFILE_DEFAULT[] ="genome.fmv";
char REVBWTINDEX_DEFAULT[] ="revgenome.bwt";
char REVOCCFILE_DEFAULT[] ="revgenome.fmv";
char GENOMEFILE_DEFAULT[]="genome";
char PATTERNFILE_DEFAULT[]="tags.fq";
char HITSFILE_DEFAULT[]="hits.txt";
char UNIQUEFILE_DEFAULT[]="unique.txt";
char AMBIGUOUSFILE_DEFAULT[]="ambiguous.txt";
char BLASTFILE_DEFAULT[]="BLAST.txt";
char MISHITFILE_DEFAULT[]="mishits.fq";
char LOGFILE_DEFAULT[]="batman.log";
char MAILCLIENT = MAIL;
char EMAIL=FALSE;
char EMAIL_CFG_DEFAULT[]="email.cfg";
char *EMAIL_CFG=EMAIL_CFG_DEFAULT;
char EMAIL_CLIENT_DEFAULT[]="./email\0";
char *EMAIL_CLIENT=EMAIL_CLIENT_DEFAULT;
#ifdef MMX
char* Version = "BatMis version(MMX) V3.00 [ILLUMINA/SOLiD]\n";
#else
char* Version = "BatMis Aligner version V3.00 [ILLUMINA/SOLiD]\n";
#endif
char UNIQUEHITS = FALSE;//do we need to separate uniquehits to a file?
char USEQUALITY=FALSE;
char HEURISTIC=FALSE;
char ALLHITS=FALSE;
char SUPERACCURATE=FALSE;
char COUNT_ALLHITS=TRUE;
char ONEFMINDEX =FALSE;
char MISMATCHES_TO_TRY_FIRST_LEFT=4;//4;
char MISMATCHES_TO_TRY_FIRST_RIGHT=4;//4;
char NORMAL_TAGS=TRUE;
char SOLIDMARK=0;
char FILETYPE;
char LOG=TRUE;
char PAIRING_TYPE;
char PRINTBLANKS=FALSE;
//char FILTERUNIQUEHITS=FALSE;
char BESTHIT=FALSE;
char SELECTBEST=TRUE;
char SCANBOTH=FALSE;
//char ZIGZAGMODE= FALSE;
char ZIGZAGMODE= TRUE;
char FILTERUNIQUEHITS=TRUE;
char MAXSPECIFIED=FALSE;//maxhits need to be found...

char* Source;
char* Guessed;
char* Guess_Complement;
char* Guessed_NLocation;
char* Guessed_NLocation_Complement;
char* Do_Branch;//decides on quality whether to branch or not...
char* Low_Quality;
char* Low_QualityF;
char* Low_QualityC;
char* Guessed_Quality;
char* Guessed_Complement_Quality;
char* Do_All;
char* Read_Buffer;
char* BWTFILE ; 
char* OCCFILE ;
char* REVBWTINDEX;
char* REVOCCFILE;
char* GENOMEFILE;
const char* Code_To_Char="acgt";
char* Current_Tag=Original;
char* NLocations;
//char* NLocationsF;
//char* NLocationsC;
char* PATTERNFILE;
char* PATTERNFILE1;
char* HITSFILE;
char* AMBIGUOUSFILE=AMBIGUOUSFILE_DEFAULT;
char* BLASTFILE=BLASTFILE_DEFAULT;
char* MISHITFILE=MISHITFILE_DEFAULT;
char* LOGFILE=LOGFILE_DEFAULT;
char* GENFILE;

map<unsigned,char> Multi_Hit;
 
//}---------------------------- GLOBAL VARIABLES -------------------------------------------------

//{---------------------------- Command Line  -------------------------------------------------
option Long_Options[]=
{
{"help",0,NULL,'h'},
{"query",1,NULL,'q'},
{"output",1,NULL,'o'},
{"usequality",0,NULL,'p'},
{"genome",1,NULL,'g'},
{"buffersize",1,NULL,'b'},
{"maxhits",1,NULL,'m'},
{"singleindex",0,NULL,'s'},
{"mishits",optional_argument,NULL,'x'},
{"maxmismatches",1,NULL,'n'},
{"scancomplement",0,NULL,'c'},
{"noquality",0,NULL,'p'},
{"indelscan",0,NULL,'i'},
{"ignorehead",1,NULL,'I'},
{"forcelength",1,NULL,'F'},
{"rollover",0,NULL,'r'},
{"printblanks",0,NULL,'B'},
{"scanboth",0,NULL,'M'},
{"maxtags",1,NULL,'t'}, 
{"zip",0,NULL,'z'}, 
{"statmode",0,NULL,'S'},
{"email",1,NULL,'e'},
{"zigzag",0,NULL,'Z'}, 
{"threshold1",1,NULL,'w'},
{"threshold2",1,NULL,'W'},
{"filteruniquehits",0,NULL,'U'}, 
{"nbound",1,NULL,'C'}, 
{"dust",0,NULL,'d'}, 
{"randomletter",0,NULL,'R'},
{"nismismatch",0,NULL,'X'}, 
{"ignoren",0,NULL,'y'},
{"forcesolid",0,NULL,'L'}, 
{"slider",1,NULL,'a'}, 
{"maxext",1,NULL,'E'},
{"logfile",1,NULL,'l'},
{"leastmismatch",0,NULL,'k'},
{"chopwhitespace",0,NULL,'u'},
{"guess",0,NULL,'G'},


{0,0,0,0}
};

char *myopts[] = { 
#define HEU_LENGTH 0 
"length",
 #define HEU_MISMATCHES 1 
"mm", 
 NULL}; 

//}---------------------------- Command Line -------------------------------------------------

int main(int argc, char* argv[])
{
	Argc=argc;Argv=argv;
//{-----------------------------  INITIALIZE ----------------------------------------------

	signal(SIGSEGV, Sig_Handle); 
	printf("%s",Version);
	Read_INI();
	Parse_Command_line(argc,argv);	
	Open_Files();
	Init_Variables();
	Allocate_Memory();
	Load_Indexes();	
	Build_Tables();
	Verbose(stdout);Verbose(Log_File);
	
//}-----------------------------  INITIALIZE ----------------------------------------------
	time(&Start_Time);
	struct SARange Range,TRange;
	int Start;
	int Progress=0,String_Length_O;
	Actual_Tag= -1;
	for (int i=0;i<MAX_MISMATCHES_BOUND;i++) Mismatch_Init;
	//if (100 == MAX_MISMATCHES_IN_EXTENSION) MAX_MISMATCHES_IN_EXTENSION=3;
	if (MAX_MISMATCHES<10 && MAX_MISMATCHES >5) MAX_MISMATCHES_IN_EXTENSION=MAX_MISMATCHES/2;
//---------------------------------------------------------------------

	if (!NORMAL_TAGS){ HITMODE=PAIREND;}
	if(SLIDER)
	{
		String_Length_O=STRINGLENGTH;
		STRINGLENGTH=SLIDER;
	}
	Header.ID[0]='B';//"BAT";
	Header.ID[1]='A';//"BAT";
	Header.ID[2]='T';//"BAT";
	Header.MAXHITS=MAXHITS;
	Header.FILETYPE=FILETYPE;
	Header.HITMODE = HITMODE;
	Header.IGNOREHEAD = IGNOREHEAD;
	Header.Tag_Length=STRINGLENGTH;
	if (SOLID) SOLIDMARK=100; else SOLIDMARK=0;
	Header.Index_Count=ONEFMINDEX+SOLIDMARK;
	if (GIS) PRINT_DESC=GISMODE;
	Header.Print_Desc=PRINT_DESC;
	fwriteX(&Header,sizeof(Header),1,Output_File);
	fwriteX(&MAX_MISMATCHES,sizeof(MAX_MISMATCHES),1,Output_File);
	fwriteX(&TAG_COPY_LEN,sizeof(int),1,Output_File);
	fwriteX(&ROLLOVER,sizeof(char),1,Output_File);
	fwriteX(&SCANBOTH,sizeof(char),1,Output_File);

	if(!NORMAL_TAGS)
	{
		fwriteX(&PAIR_LENGTH_LEFT,sizeof(int),1,Output_File);
		fwriteX(&PAIR_LENGTH_RIGHT,sizeof(int),1,Output_File);
	}
	if(UNIQUEHITS)
	{
		Header.HITMODE = DEFAULT;
		fwriteX(&Header,sizeof(Header),1,Unique_File);
	}
//---------------------------------------------------------------------

	for (int i=1;i<256;i++)
	{
		All_Zero[i]=0;
	}
	printf("======================]\r[");//progress bar....
	fflush(stdout);

	for (int i=0;i<STRINGLENGTH;i++) Do_All[i]=TRUE;
	Tag_Printed=TRUE;//avoid first printblank...
	Hits=0;char Init=TRUE; 

	bool In_Sliding_Win=false; int Slide_Ptr=0; int Slide_Win_Length=36;
	char* Current_Tag_O=Current_Tag;char Old_Plus;

	for(;;)//Tag Processing loop.....
	{
		if(!Tag_Stat_Bad && !FILTERUNIQUEHITS) Total_Hits=Total_Hits+Hits;

		if(!GIS)
		{
			if (PRINTBLANKS && !Tag_Printed && First_Pass)// && !First_Pass_Printed)
			{
				if(!NORMAL_TAGS && First_Pass_Printed){}
				else
				{
					New_Record='@';
					fwriteX(&New_Record,1,1,Output_File);//write new record marker... later make gap global

					if(PRINT_DESC) fprintfX(Output_File,"%s",Description);
					if (!NORMAL_TAGS)
					{
						fprintfX(Output_File,"%s",Tag_Copy);
					}
					else
					{
						for( int i=0;i<STRINGLENGTH; i++) Translated_String[i]=Code_To_Char[Current_Tag[i]];
						Translated_String[STRINGLENGTH]='\n';//Current_Tag[STRINGLENGTH];
						fwriteX(Translated_String,1,STRINGLENGTH+1,Output_File);//write tag...
					}
				}
			}
		}

		Current_Tag=Original+IGNOREHEAD;
		Mismatches_Reduced=FALSE;
		Direction=0;//this is the default for brute force scan...
		HIGHSCAN=FALSE;
		Actual_Tag++;
		Progress++;
		char QUAL_CALCULATED=FALSE;
		if (Progress==Number_of_Tags) 
		{
			if (MAX_TAGS_TO_PROCESS)
			{
				Number_of_Tags=(MAX_TAGS_TO_PROCESS)/20;
				Progress=0;
				Show_Progress(Actual_Tag*100/MAX_TAGS_TO_PROCESS);
			}
			else
			{
			//off64_t Current_Pos=ftello64(Input_FileO);
			off64_t Current_Pos=ftello64(Input_File);
			Average_Length=Current_Pos/Actual_Tag+1;//+1 avoids divide by zero..
			Number_of_Tags=(File_Size/Average_Length)/20;
			Progress=0;
			Show_Progress(Current_Pos*100/File_Size);
			}
		}
		Print_Header=FALSE;//Header not printed yet...

		if(NORMAL_TAGS)
		{
			if (GIS && !Init) Print_Hits();
			if (MAX_TAGS_TO_PROCESS && Actual_Tag == MAX_TAGS_TO_PROCESS) break;
			Tag_Stat_Bad=FALSE;Tags_From_Head=0;
			Hits=0;Hits1=0;First_Pass_Hits=0;
			Tag_Printed=FALSE;Stats=Stats1;
			Stats[0]=Stats[1]=Stats[2]=Stats[3]=Stats[4]=Stats[5]=Stats[6]=0;
			Head_Printed=FALSE;
			Tag_Number=1;
			//if (gzgets(Input_File,Description,MAXDES)!=0)// read a tag...
			if(SLIDER)
			{
				if (SLIDERHIT) {In_Sliding_Win=0;SLIDERHIT=0;}
				if (In_Sliding_Win)
				{
					if (Slide_Ptr + STRINGLENGTH >= String_Length_O) {Slide_Ptr=0;In_Sliding_Win=false;Current_Tag=Current_Tag_O;}
					else 
					{
						//Current_Tag++;
						Slide_Ptr++;
						for(int i=1;i<=String_Length_O;i++) Original[i-1]=Current_Tag[i];
					}
				}
			}

			if (In_Sliding_Win || fgets(Description,MAXDES,Input_File)!=0)// read a tag... (assume lazy eval)
			{
				if (CHOPWHITE) 
				{
					char* p;
					for(p=Description;*p!=' '&&*p!='\t'&&*p!='\n'&&*p;p++); 
					*p++='\n';*p=0;
				}
				Init=FALSE;//not initial pass
				//gzgets(Input_File,Current_Tag-IGNOREHEAD,MAXDES);//tag
				if(!In_Sliding_Win) 
				{
					fgets(Current_Tag-IGNOREHEAD,MAXDES,Input_File);//tag
					if (SLIDER) 
					{
						strcpy(Temp_Description,Description);
					}
					strcpy(Tag_Copy,Current_Tag-IGNOREHEAD);
				}
				if (FILETYPE == FQ)
				{
					if(!In_Sliding_Win)
					{
						//gzgets(Input_File,Plus,MAXTAG);//plus
						fgets(Plus,MAXTAG,Input_File);//plus
						//gzgets(Input_File,Quality,MAXTAG);//phred
						fgets(Quality,MAXTAG,Input_File);//phred
					}
				}
				NCount=0;int j=0;
				if (SLIDER) 
				{
					sprintf(Description,"%c%d%s",Temp_Description[0],Slide_Ptr,Temp_Description);
					if (Slide_Ptr)
					{
						assert(Current_Tag[STRINGLENGTH-1]!='+');
						Current_Tag[STRINGLENGTH-1]=Old_Plus;
					}
				}

				In_Sliding_Win=SLIDER;
				for (unsigned i=0;i<=STRINGLENGTH-1;i++)
				{
					if (Current_Tag[i] == 'n' || Current_Tag[i]=='N'||Current_Tag[i] == '.')
					{
						N[j++]=i;NLocations[i]=TRUE;NCount++;
						//Current_Tag[i]=Random_Array[Random_Pointer++];N[j++]=Current_Tag[i];
						if (RANDOMLETTER) Current_Tag[i]=Random_Array[Random_Pointer++];N[j++]=Current_Tag[i];
						if (Random_Pointer==sizeof(Random_Array)-1) Random_Pointer=0; 
					}
					else NLocations[i]=FALSE;
					Current_Tag[i]=Char_To_Code[Current_Tag[i]];
					if (SOLID) Complement[STRINGLENGTH-1-i]=Current_Tag[i];
					else Complement[STRINGLENGTH-1-i]=Char_To_CodeC[Current_Tag[i]];
				}
				if(NCount>NBOUND) {if(PRINT_MISHITS) fprintf(Mishit_File,"%s%s%s%s", Description,Tag_Copy,Plus,Quality);continue;}
				if (NISMISMATCH && MAX_MISMATCHES <= NCount ) continue;//Too many N's than mismatches...
				Old_Plus=Current_Tag[STRINGLENGTH];
				Current_Tag[STRINGLENGTH]='+';
				Complement[STRINGLENGTH]='-';

			}
			else {break;}
		}
		else //pair end tags...
		{
			if (First_Pass )//Process head...
			{
				if (GIS) Print_Hits();
				if (MAX_TAGS_TO_PROCESS && Actual_Tag == MAX_TAGS_TO_PROCESS) break;
				Head_Printed=FALSE;
				Tag_Stat_Bad=FALSE;Tags_From_Head=0;
				Tag_Number=1;Stats=Stats1;
				Hits1=0;Hits=0;First_Pass_Hits=0;
				Stats[0]=Stats[1]=Stats[2]=Stats[3]=Stats[4]=Stats[5]=Stats[6]=0;
				Tag_Printed=FALSE;
				First_Pass_Printed=FALSE;
				STRINGLENGTH=PAIR_LENGTH_LEFT;
				//if (gzgets(Input_File,Description,MAXDES)!=0)//des .... read a tag
				if (fgets(Description,MAXDES,Input_File)!=0)//des .... read a tag
				{
					NCount=0;	
					//gzgets(Input_File,Current_Tag-IGNOREHEAD,MAXDES);//tag
					fgets(Current_Tag-IGNOREHEAD,MAXDES,Input_File);//tag
					if(PAIRING_TYPE==TWOFILE)
					{
						//gzgets(Mate_File,Description,MAXDES);
						fgets(Description,MAXDES,Mate_File);
						//gzgets(Mate_File,Current_Tag+STRINGLENGTH+1,MAXDES);//tag
						fgets(Current_Tag+STRINGLENGTH+1,MAXDES,Mate_File);//tag
						Current_Tag[STRINGLENGTH]='\t';
					}
					strcpy(Tag_Copy,Current_Tag-IGNOREHEAD);int j=0;
					for (unsigned i=0;i<=STRINGLENGTH-1;i++)
					{
						if (Current_Tag[i] == 'n' || Current_Tag[i]=='N'|| Current_Tag[i]=='.')
						{
							N[j++]=i;NLocations[i]=TRUE;NCount++;
							//Current_Tag[i]=Random_Array[Random_Pointer++];N[j++]=Current_Tag[i];
							if(RANDOMLETTER) Current_Tag[i]=Random_Array[Random_Pointer++];N[j++]=Current_Tag[i];
							if (Random_Pointer==sizeof(Random_Array)-1) Random_Pointer=0; 
						}
						else NLocations[i]=FALSE;
						Current_Tag[i]=Char_To_Code[Current_Tag[i]];
						if (SOLID) Complement[STRINGLENGTH-1-i]=Current_Tag[i];
						else Complement[STRINGLENGTH-1-i]=Char_To_CodeC[Current_Tag[i]];//change later...
					}
					if(NCount>NBOUND) continue;
					Current_Tag[STRINGLENGTH]='+';
					Complement[STRINGLENGTH]='-';
				}
				else break;

				if (FILETYPE == FQ)
				{
					//gzgets(Input_File,Plus,MAXTAG);//plus
					fgets(Plus,MAXTAG,Input_File);//plus
					//gzgets(Input_File,Quality_Copy,MAXTAG);//phred
					fgets(Quality_Copy,MAXTAG,Input_File);//phred
					strcpy(Quality,Quality_Copy+IGNOREHEAD);
					if(PAIRING_TYPE==TWOFILE)
					{
						//gzgets(Mate_File,Plus,MAXTAG);//plus
						fgets(Plus,MAXTAG,Mate_File);//plus
						//gzgets(Mate_File,Quality_Copy+STRINGLENGTH+1+IGNOREHEAD,MAXTAG);//phred
						fgets(Quality_Copy+STRINGLENGTH+1+IGNOREHEAD,MAXTAG,Mate_File);//phred
						Quality_Copy[STRINGLENGTH+IGNOREHEAD]='\t';
					}
				}
				LH=PLH;LHQL=PLHQL;LHQR=PLHQR;RH=PRH;RHQL=PRHQL;RHQR=PRHQR;//LOOKUPSIZE=PLOOKUPSIZE;
				First_Pass=FALSE;
			}

			else//Process tail
			{
				if (MAX_TAGS_TO_PROCESS && Actual_Tag == MAX_TAGS_TO_PROCESS) break;
				Hits1+=Hits;Hits=0;First_Pass_Hits=0;
				Tag_Number=2;Stats=Stats2;
				Stats[0]=Stats[1]=Stats[2]=Stats[3]=Stats[4]=Stats[5]=Stats[6]=0;
				if (Tag_Stat_Bad) {First_Pass=TRUE;continue;}; 

				STRINGLENGTH=PAIR_LENGTH_RIGHT;
				if (FILETYPE == FQ) strcpy(Quality,Quality_Copy+PAIR_LENGTH_LEFT+2*IGNOREHEAD+1);
				NCount=0;int j=0;
				for (unsigned i=0;i<=STRINGLENGTH-1;i++)
				{
					Current_Tag[i]=Char_To_Code[Current_Tag[i+PAIR_LENGTH_LEFT+1+IGNOREHEAD]];
					if (Current_Tag[i] == 'n' || Current_Tag[i]=='N'|| Current_Tag[i]=='.')
					{
						N[j++]=i;NLocations[i]=TRUE;NCount++;
						if (RANDOMLETTER) Current_Tag[i]=Random_Array[Random_Pointer++];N[j++]=Current_Tag[i];
						if (Random_Pointer==sizeof(Random_Array)-1) Random_Pointer=0; 
					}
					else NLocations[i]=FALSE;
					if (SOLID) Complement[STRINGLENGTH-1-i]=Current_Tag[i];//change later...
					else Complement[STRINGLENGTH-1-i]=Char_To_CodeC[Current_Tag[i]];//change later...
				}
				if(NCount>NBOUND) {if(PRINT_MISHITS && !Hits1) fprintf(Mishit_File,"%s%s%s%s", Description,Tag_Copy,Plus,Quality);continue;}
				Current_Tag[STRINGLENGTH]='+';
				Complement[STRINGLENGTH]='-';
				LH=PRLH;LHQL=PRLHQL;LHQR=PRLHQR;RH=PRRH;RHQL=PRRHQL;RHQR=PRRHQR;//LOOKUPSIZE=PRLOOKUPSIZE;
				First_Pass=TRUE;
			}
		}

		if (DUST)
		{
			fa.len=STRINGLENGTH;
			for (int i=0;i<STRINGLENGTH;i++) Translated_String[i]=Code_To_Char[Current_Tag[i]];
			if (dust(&fa) >5) continue;
		}
//----------------------------------------------------------------------------------------------------------------------------------------
//{------------------------------------- QUALITY CHECK --------------------------------------------------------------------------------------------------
		//if (USEQUALITY || (MAX_MISMATCHES>5 && FILETYPE==FQ))
		if (USEQUALITY || (USEQUALITYH))
		{
			Quality_Calc();
		}
//}------------------------------------- QUALITY CHECK --------------------------------------------------------------------------------------------------
		char Reverse_Exact_Not_Scanned=TRUE;//in tag guess, we did not find an exact match...
		In_Mismatch=0;
		In_MismatchX=0;
		Last_In_Mis=100;
//{------------------------------------------TAG GUESS ---------------------------------------------------------------------------------------------

		//Rollover_Step=FALSE;
		//TagProcessed=FALSE;

		if (TAG_GUESS)
		{
			if (!Guess_Orientation()) //Exact match found...
			{
				if(AMBI){AMBI=FALSE;continue;}
				if(MAX_MISMATCHES == 0) continue;
				if (Direction) Reverse_Exact_Not_Scanned=FALSE;
			}//else tag direction guessed only...
			else
			{
				Print_Header=FALSE;//Do not print description...
			}
		}
		else if (ZIGZAGMODE)
		{
			Zig_Zag();

			if(PRINT_MISHITS && MAX_MISMATCHES <=5)
			{
				if (Tag_Stat_Bad) fprintf(Mishit_File,"%s%sBAD%s%s", Description,Tag_Copy,Plus,Quality);
				else if(!Hits) fprintf(Mishit_File,"%s%s%s%s", Description,Tag_Copy,Plus,Quality);
			}
			/*if(PRINT_MISHITS)
			{
#ifdef SOLID_DBG
				char *i;
				for(i=Description;*i !='\n' && *i !='\r';i++);*i=0;
				if (Hits) fprintf(Mishit_File,"M\t%s\t%d\n", Description,Actual_Tag);
				else fprintf(Mishit_File,"U\t%s\t%d\n", Description,Actual_Tag);
#else
				if(!Hits) fprintf(Mishit_File,"%s%s%s%s", Description,Tag_Copy,Plus,Quality);
#endif
			}*/
			continue;
		}
		if(Tag_Stat_Bad) 
		{
			continue;
		}

		if (SEED)
		{
			Seed_Scan();
			if(PRINT_MISHITS && !Hits) fprintf(Mishit_File,"%s%s%s%s", Description,Tag_Copy,Plus,Quality);
			continue;
		}
//}------------------------------------------TAG GUESS ---------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------
		Rollover_Step=FALSE;
		TagProcessed=FALSE;

		while(Hits<MAXHITS || SCANBOTH) 
		{
			if (TagProcessed)
			{
				if (!ROLLOVER)
				{
					if(Hits) break;
				}
				else
				{ // rolling over ...
					if(Rollover_Step) break;//already in a rollover...
					else
					{
						if (SCANBOTH) {Last_Mismatch=In_Mismatch;Last_Write_Buffer=Write_Buffer_Ptr;if(!FILTERUNIQUEHITS) Total_Hits=Total_Hits+Hits;Hits1+=Hits;First_Pass_Hits=Hits;Hits=0;}
						Rollover_Step=TRUE;
					}
				}
			}
			TagProcessed=TRUE;
			//In_Mismatch=0;
			if (BESTHIT && Rollover_Step && Last_Mismatch == 0) {In_Mismatch=6;continue;}
			if (Direction == DIRECTIONS)//have tried both directions.. 
			{
				if(PRINT_MISHITS && MAX_MISMATCHES <=5) fprintf(Mishit_File,"%s%s%s%s", Description,Tag_Copy,Plus,Quality);
				break;
			}
//----------------------------------------------------------------------------------------------------------------------------------------
			if (Direction == 1)//complementary direction
			{
				Current_Tag=Complement;
				Low_Quality=IGNOREHEAD+Low_QualityC;
				//NLocations=NLocationsC;
				Exact_Match_Forward=Exact_Match_ForwardC;
				Current_Tag[STRINGLENGTH]='-';
				if(Hit_In !=2) Print_Header=TRUE;//we did not start with reverse complement first...
				//Guess_Complement=Current_Tag;
				if(Reverse_Exact_Not_Scanned) //Forward scan found an exact match, and so reverse complement was not scanned for exact match...
				{
					In_Mismatch=0;
					FMIndex=REVERSE;
					if(LOOKUPSIZE ==3)
					{
						c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
					}
					else
					{
						c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
					}

					Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Tag=Actual_Tag;
					Range.Mismatches=0;Range.Level=LOOKUPSIZE+1; Range.Skip=0;
					memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
					Search_Forwards_Exact(Range,-1,STRINGLENGTH,revfmi);//Find exact matches and report... if not found get the range for 0|?
					if (Hits ==MAXHITS) continue;
				}
			}
			else if (Direction ==2)//happens if the first dirction tried is comp, and then we have to try original dirn.
			{
				Current_Tag=Original+IGNOREHEAD;
				Low_Quality=Low_QualityF+IGNOREHEAD;
				Exact_Match_Forward=Exact_Match_ForwardF;
				Current_Tag[STRINGLENGTH]='+';
				//Guess_Complement=Current_Tag;
			}
			


			Direction++;
			Left_Mishits_Pointer=0;
			Right_Mishits_Pointer=0;
			Possible_20_Pointer=0;
			Possible_02_Pointer=0;
			Mismatches_Forward_Pointer=0;//first node where SA range was not found, all other nodes will not have matches..
			Mismatches_Backward_Pointer=0;
			Two_Mismatches_At_End_Pointer=0;
			Two_Mismatches_At_End_Forward_Pointer=0;
			Start=1-2;//Adjust for offsets...
			FMIndex=REVERSE;

			if (SCANBOTH && Hits >= MAXHITS) continue;
			if( !TAG_GUESS)
			{
				In_Mismatch=0;
				if(LOOKUPSIZE ==3)
				{
					c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
				}
				else
				{
					c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
				}
				Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Tag=Actual_Tag;
				Range.Level=LOOKUPSIZE+1;Range.Mismatches=0;Range.Skip=0;
				memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
				Search_Forwards_Exact(Range,Start,STRINGLENGTH,revfmi);//Find exact matches and report... if not found get the range for 0|?
				if(MAXHITS==Hits) continue;  
				if(!NCount && MAX_MISMATCHES == 0) continue;
			}
			if (FILTERUNIQUEHITS && Hits) continue;

//{------------------------------------------- ONE MISMATCH ---------------------------------------------------------------------------------------------
			//One mismatches...
			if (BESTHIT && Rollover_Step && Last_Mismatch == 1) {In_Mismatch=6;continue;}
			In_Mismatch=1;
			Range=Exact_Match_Forward[Start+LH];
			if(Range.Start && Range.Tag == Actual_Tag)//if there are hits of the form 0|?
			{
				Range.Level=1;
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					Search_Forwards(Range,1,LH+1,RH,revfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				Search_Forwards(Range,1,LH+1,RH,revfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
				if(MAXHITS==Hits) continue;


			}		
			FMIndex=FORWARD;
			if(LOOKUPSIZE ==3)
			{
				c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
			}
			else
			{
				c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
			}
			Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];Range.Level=LOOKUPSIZE+1;
			Range.Mismatches=0;Range.Tag=Actual_Tag;Range.Skip=0;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
			Search_Backwards_Exact( Range,STRINGLENGTH,RH,fwfmi);//Backward scan for ?|0
			if(Range.Start)//if there are possible hits of the form ?|0
			{
				Range.Level=1;
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					Search_Backwards(Range,1,LH,LH,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				Search_Backwards(Range,1,LH,LH,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
				if(MAXHITS==Hits) continue;
			}
			if (NPOLICY && NCount){if (!NISMISMATCH && MAX_MISMATCHES + NCount == 1) continue;else if (MAX_MISMATCHES==NCount) continue;}
			else if (MAX_MISMATCHES == 1 ) continue;
			if (FILTERUNIQUEHITS && Hits) continue;
//}------------------------------------------- ONE MISMATCH ---------------------------------------------------------------------------------------------

//{------------------------------------------- TWO MISMATCH ---------------------------------------------------------------------------------------------
			if (BESTHIT && Rollover_Step && Last_Mismatch == 2) {In_Mismatch=6;continue;}
			In_Mismatch=2;
			FMIndex=REVERSE;
			if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
			{
				for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
				{
					Two_Mismatches_At_End_Forward[i].Mismatch_Pos[1]= (STRINGLENGTH-1);//mismatches of the form 0|2, with last mismatch at the end...
					Print_LocationX(Two_Mismatches_At_End_Forward[i]);
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}
			Two_Mismatches_At_End_Forward_Pointer=0;

			FMIndex=FORWARD;
			if(Two_Mismatches_At_End_Pointer)
			{
				Do_Branch=Do_All;
				for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
				{
					Print_LocationX(Two_Mismatches_At_End[i]);//Mismatches of the form 2|0, with one mismatch at the first position
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;

			}

			Two_Mismatches_At_End_Pointer=0;
			FMIndex=REVERSE;
			int Possible_03_Pointer=Mismatches_Forward_Pointer;
			if(Mismatches_Forward_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=Possible_03_Pointer-1;i>=0;i--)
					{
						Search_Forwards(Mismatches_Forward[i],2,LH+1,RH,revfmi);//scan for possible two mismatches of the form 0|2, and store candidates for 0|3
						if(MAXHITS==Hits) break;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=Possible_03_Pointer-1;i>=0;i--)
				{
					Search_Forwards(Mismatches_Forward[i],2,LH+1,RH,revfmi);//scan for possible two mismatches of the form 0|2, and store candidates for 0|3
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}

			FMIndex=FORWARD;
			int Possible_30_Pointer=Mismatches_Backward_Pointer;
			if(Mismatches_Backward_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=Possible_30_Pointer-1;i>=0;i--)
					{
						Search_Backwards(Mismatches_Backward[i],2,LH,LH,fwfmi);//scan for possible two mismatches of the form 2|0, and stores the candidates for 3|0
						if(MAXHITS==Hits) break;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=Possible_30_Pointer-1;i>=0;i--)
				{
					Search_Backwards(Mismatches_Backward[i],2,LH,LH,fwfmi);//scan for possible two mismatches of the form 2|0, and stores the candidates for 3|0
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}

//----------------------------------------------------------------------------------------------------------------------------------------
			if(LOOKUPSIZE==3)
			{
				c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
			}
			else
			{
				c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
			}
			Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];
			Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;//Range.Tag=Actual_Tag;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;

			Search_Backwards_Exact_X0( Range,STRINGLENGTH,RHQR,fwfmi);// ?|?|0
			Range.Level=1;

			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				Search_Backwards_X10(Range,1,LH + RHQL, RHQL,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
				if(MAXHITS==Hits) continue;
			}

			Do_Branch=Do_All;
			Search_Backwards_X10(Range,1,LH + RHQL, RHQL,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
			if(MAXHITS==Hits) continue;
//----------------------------------------------------------------------------------------------------------------------------------------
			if(LOOKUPSIZE==3)
			{
				c=Current_Tag[LH+0] | (Current_Tag[LH+1]<<2) | (Current_Tag[LH+2]<<4);// | (Current_Tag[LH+3]<<6) | Current_Tag[LH+4]<<8 | (Current_Tag[LH+5]<<10);//Use lookup table...
			}
			else
			{
				c=Current_Tag[LH+0] | (Current_Tag[LH+1]<<2) | (Current_Tag[LH+2]<<4) | (Current_Tag[LH+3]<<6) | Current_Tag[LH+4]<<8 | (Current_Tag[LH+5]<<10);//Use lookup table...
			}
			Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];
			Range.Mismatches=0;Range.Level=LOOKUPSIZE+1; Range.Skip=0;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
			Search_Forwards_0X(Range,LH+1,RHQL,revfmi);
			Range.Level=1;TRange=Range;
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;

				Search_X01(Range,1,LH + RHQL +1,RHQR,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
				if(MAXHITS==Hits) continue;
			}

			Do_Branch=Do_All;
			Search_X01(TRange,1,LH + RHQL +1,RHQR,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
			if(MAXHITS==Hits) continue;
			if (NPOLICY && NCount){if (!NISMISMATCH && MAX_MISMATCHES + NCount == 2) continue;else if (MAX_MISMATCHES==NCount) continue;}
			else if( MAX_MISMATCHES ==2) continue;
			if (FILTERUNIQUEHITS && Hits) continue;
//}------------------------------------------- TWO MISMATCH ---------------------------------------------------------------------------------------------

//{------------------------------------------- THREE MISMATCH ---------------------------------------------------------------------------------------------
			//Find three mismatches....
			if (BESTHIT && Rollover_Step && Last_Mismatch == 3) {In_Mismatch=6;continue;}
			In_Mismatch=3;
			FMIndex=REVERSE;
			if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
			{
				for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
				{
					Print_LocationX(Two_Mismatches_At_End_Forward[i]);//mismatches of the form 0|3, with last mismatch at the end...
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}
			Two_Mismatches_At_End_Forward_Pointer=0;

			FMIndex=FORWARD;
			if(Two_Mismatches_At_End_Pointer)
			{
				for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
				{
					Print_LocationX(Two_Mismatches_At_End[i]);//Mismatches of the form 3|0, with one mismatch at the first position
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}
			Two_Mismatches_At_End_Pointer=0;

			FMIndex=REVERSE;
			Possible_04_Pointer=Mismatches_Forward_Pointer;
			if(Mismatches_Forward_Pointer!=Possible_03_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=Possible_04_Pointer-1;i>=Possible_03_Pointer;i--)
					{
						Search_Forwards(Mismatches_Forward[i],3,LH+1,RH,revfmi);//scan for possible three mismatches of the form 0|3, and finds mismatches of the form 1|2, stores possibles in the form 1|3
						if(MAXHITS==Hits) break;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=Possible_04_Pointer-1;i>=Possible_03_Pointer;i--)
				{
					Search_Forwards(Mismatches_Forward[i],3,LH+1,RH,revfmi);//scan for possible three mismatches of the form 0|3, and finds mismatches of the form 1|2, stores possibles in the form 1|3
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}

			FMIndex=FORWARD;
			Possible_40_Pointer=Mismatches_Backward_Pointer;
			if(Mismatches_Backward_Pointer!=Possible_30_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=Possible_40_Pointer-1;i>=Possible_30_Pointer;i--)
					{
						Search_Backwards(Mismatches_Backward[i],3,LH,LH,fwfmi);//scan for possible mismatches of the form 3|0, 2|1 and sotres the candidates for 4|0, 3|1
						if(MAXHITS==Hits) break;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=Possible_40_Pointer-1;i>=Possible_30_Pointer;i--)
				{
					Search_Backwards(Mismatches_Backward[i],3,LH,LH,fwfmi);//scan for possible mismatches of the form 3|0, 2|1 and sotres the candidates for 4|0, 3|1
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;

			}

//----------------------------------------------------------------------------------------------------------------------------------------
			FMIndex=REVERSE;
			if(LOOKUPSIZE==3)
			{
				c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
			}
			else
			{
				c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
			}
			Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0; Range.Skip=0;
			Search_Forwards_0X(Range,1,LHQL,revfmi);
			Range.Level=1;
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				Search_01X(Range,1,LHQL +1,LHQR,revfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
				if(MAXHITS==Hits) continue;
			}

			Do_Branch=Do_All;
			Search_01X(Range,1,LHQL +1,LHQR,revfmi);
			if(MAXHITS==Hits) continue;
//----------------------------------------------------------------------------------------------------------------------------------------
			if(LOOKUPSIZE==3)
			{
				c=Current_Tag[LH-1-0] | (Current_Tag[LH-1-1]<<2) | (Current_Tag[LH-1-2]<<4);// | (Current_Tag[LH-1-3]<<6) | Current_Tag[LH-1-4]<<8 | (Current_Tag[LH-1-5]<<10);//Use lookup table...
			}
			else
			{
				c=Current_Tag[LH-1-0] | (Current_Tag[LH-1-1]<<2) | (Current_Tag[LH-1-2]<<4) | (Current_Tag[LH-1-3]<<6) | Current_Tag[LH-1-4]<<8 | (Current_Tag[LH-1-5]<<10);//Use lookup table...
			}

			Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;//Range.Tag=Actual_Tag;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
			Search_Backwards_Exact_X0( Range,LH,LHQR,fwfmi);// ?|0|?
			Range.Level=1;
			if(USEQUALITY)
			{
				TRange=Range;
				Do_Branch=Low_Quality;
				Search_10X(TRange,1,LHQL, LHQL,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3 
				if(MAXHITS==Hits) continue;
			}
			Do_Branch=Do_All;
			Search_10X(Range,1,LHQL, LHQL,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
			if(MAXHITS==Hits) continue;
			if (NPOLICY && NCount){if (!NISMISMATCH && MAX_MISMATCHES + NCount == 3) continue;else if (MAX_MISMATCHES==NCount) continue;}
			else if( MAX_MISMATCHES ==3 && !INDELSCAN) continue;
			if (FILTERUNIQUEHITS && Hits) continue;
//}------------------------------------------- THREE MISMATCH ---------------------------------------------------------------------------------------------
//{-------------------------------------------  INDEL  ---------------------------------------------------------------------------------------------
			if (INDELSCAN)
			{
				FMIndex=REVERSE;
				Range=Exact_Match_Forward[Start+LH];
				memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
				if(Range.Start && Range.Tag == Actual_Tag)//if there are hits of the form 0|?
				{
					Range.Level=1;
					Do_Branch=Do_All;
					Search_Forwards_Indel(Range,1,LH+1,RH,revfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
					if(MAXHITS==Hits) continue;
				}

				FMIndex=FORWARD;
				if(LOOKUPSIZE ==3)
				{
					c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
				}
				else
				{
					c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
				}
				Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;Range.Tag=Actual_Tag;
				memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
				Search_Backwards_Exact( Range,STRINGLENGTH,RH,fwfmi);//Backward scan for ?|0
				if(Range.Start)//if there are possible hits of the form ?|0
				{
					Range.Level=1;

					Do_Branch=Do_All;
					Search_Backwards_Indel(Range,0,LH,LH,fwfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
					if(MAXHITS==Hits) continue;
				}
				if (NPOLICY && NCount){if (!NISMISMATCH && MAX_MISMATCHES + NCount == 3) continue;else if (MAX_MISMATCHES==NCount) continue;}
				else if( MAX_MISMATCHES ==3 ) continue;
			}
			if (FILTERUNIQUEHITS && Hits) continue;

//}-------------------------------------------  INDEL  ---------------------------------------------------------------------------------------------

//{------------------------------------------- FOUR MISMATCH ---------------------------------------------------------------------------------------------
			if (BESTHIT && Rollover_Step && Last_Mismatch == 4) {In_Mismatch=6;continue;}
			In_Mismatch=4;
			FMIndex=REVERSE;
			if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
			{
				for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
				{
					Print_LocationX(Two_Mismatches_At_End_Forward[i]);//mismatches of the form 0|4, with last mismatch at the end...
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}

			Two_Mismatches_At_End_Forward_Pointer=0;
			FMIndex=FORWARD;
			if(Two_Mismatches_At_End_Pointer)
			{
				for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
				{
					Print_LocationX(Two_Mismatches_At_End[i]);//mismatches of the form 0|4, with one mismatch at the start...
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}
			Two_Mismatches_At_End_Pointer=0;

			FMIndex=REVERSE;
			int Possible_05_Pointer=Mismatches_Forward_Pointer;
			int Wrap=FALSE;
			if(Mismatches_Forward_Pointer)
			{
				if (Possible_04_Pointer > 46000) {Mismatches_Forward_Pointer=0;Wrap=TRUE;}
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=Possible_04_Pointer;i<Possible_05_Pointer;i++)//Mismatches_Forward_Pointer;i++)
					{
						Search_Forwards(Mismatches_Forward[i],4,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|4, and finds mismatches of the form 1|3, stores possibles in the form 1|4
						if(MAXHITS==Hits) break;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=Possible_04_Pointer;i<Possible_05_Pointer;i++)//Mismatches_Forward_Pointer;i++)
				{
					Search_Forwards(Mismatches_Forward[i],4,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|4, and finds mismatches of the form 1|3, stores possibles in the form 1|4
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}
			if(Wrap) Possible_05_Pointer=0;
			int Mismatches_Forward_Pointer_Last4=Mismatches_Forward_Pointer;

			FMIndex=FORWARD;
			Possible_50_Pointer=Mismatches_Backward_Pointer;
			if(Mismatches_Backward_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=Possible_50_Pointer-1;i>=Possible_40_Pointer;i--)//Mismatches_Backward_Pointer-1;i>=0;i--)
					{
						Search_Backwards(Mismatches_Backward[i],4,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
						if(MAXHITS==Hits) break;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=Possible_50_Pointer-1;i>=Possible_40_Pointer;i--)//Mismatches_Backward_Pointer-1;i>=0;i--)
				{
					Search_Backwards(Mismatches_Backward[i],4,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;

			}

			FMIndex=REVERSE;
			int Left_Mishits_Pointer_1=Left_Mishits_Pointer;//Polish some more....
			if(Left_Mishits_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=0;i<Left_Mishits_Pointer;i++)
					{
						TRange=Left_Mishits[i];
						if (Left_Mishits[i].Level != LH+1)
						{
							Search_Exact(Left_Mishits[i],-1,LH,revfmi);
							Left_Mishits[i].Level=LH+1;//search exact does not update level..
						}
						if (Left_Mishits[i].Start)
						{
							Search_Forwards(Left_Mishits[i],4,1,STRINGLENGTH,revfmi);//find mismatches of the form 022 form, stores possibles of the form 023
						}
						Left_Mishits[i]=TRange;
						if(MAXHITS==Hits) break;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=0;i<Left_Mishits_Pointer;i++)
				{
					if (Left_Mishits[i].Level != LH+1)
					{
						Search_Exact(Left_Mishits[i],-1,LH,revfmi);
						Left_Mishits[i].Level=LH+1;//search exact does not update level..
					}
					if (Left_Mishits[i].Start)
					{
						Search_Forwards(Left_Mishits[i],4,1,STRINGLENGTH,revfmi);//find mismatches of the form 022 form, stores possibles of the form 023
					}
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;

			}


			int Mismatches_Forward_Pointer_Last5=Mismatches_Forward_Pointer;
			if( Right_Mishits_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=0;i<Right_Mishits_Pointer;i++)
					{
						TRange=Right_Mishits[i];
						if(Right_Mishits[i].Level!=LHQL) 
						{
							Search_Backwards_Exact( Right_Mishits[i],LHQL,LHQL,fwfmi);//finds mismatches of the form 202, stores possibles of the form 203
						}
						if(Right_Mishits[i].Start)
						{	
							//Right_Mishits[i].Skip=0;
							Backwards(Right_Mishits[i],1,LH);
							if(Right_Mishits[i].Start)
							{
								Right_Mishits[i].Level=1;
								Search_Forwards(Right_Mishits[i],4,LH+1,RH,revfmi);
								if(MAXHITS==Hits) break;
							}
						}
						Right_Mishits[i]=TRange;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=0;i<Right_Mishits_Pointer;i++)
				{

					if(Right_Mishits[i].Level!=LHQL) 
					{
						Search_Backwards_Exact( Right_Mishits[i],LHQL,LHQL,fwfmi);//finds mismatches of the form 202, stores possibles of the form 203
					}
					if(Right_Mishits[i].Start)
					{	
						Backwards(Right_Mishits[i],1,LH);
						if(Right_Mishits[i].Start)
						{
							Right_Mishits[i].Level=1;
							Search_Forwards(Right_Mishits[i],4,LH+1,RH,revfmi);
							if(MAXHITS==Hits) break;
						}
					}
				}
				if(MAXHITS==Hits) continue;
			}

			FMIndex=REVERSE;
			int LHQLrx=LHQL/2;
			if (LHQL % 2) LHQLrx++; int LHQRrx=LHQL-LHQLrx;

			if (LOOKUPSIZE >= LHQLrx)
			{
				Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
				memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
			}
			else
			{
				if(LOOKUPSIZE==3)
				{
					c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
				}
				else
				{
					c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
				}
				Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;
				memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0; Range.Skip=0;
			}
			Search_Forwards_0X(Range,1,LHQLrx,revfmi);
			Range.Level=1;
			if(USEQUALITY)
			{
				TRange=Range;
				Do_Branch=Low_Quality;
				Search_01LX(Range,1,LHQLrx +1,LHQRrx,revfmi);
				if(MAXHITS==Hits) continue;
				Range=TRange;
			}

			Do_Branch=Do_All;
			Search_01LX(Range,1,LHQLrx +1,LHQRrx,revfmi);
			if(MAXHITS==Hits) continue;
//--------------------------------
			if (LOOKUPSIZE >= LHQRrx)
			{
				Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
				memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
			}
			else
			{
				if(LOOKUPSIZE==3)
				{
					c=Current_Tag[LHQL-1-0] | (Current_Tag[LHQL-1-1]<<2) | (Current_Tag[LHQL-1-2]<<4);// | (Current_Tag[LH-1-3]<<6) | Current_Tag[LH-1-4]<<8 | (Current_Tag[LH-1-5]<<10);//Use lookup table...
				}
				else
				{
					c=Current_Tag[LHQL-1-0] | (Current_Tag[LHQL-1-1]<<2) | (Current_Tag[LHQL-1-2]<<4) | (Current_Tag[LHQL-1-3]<<6) | Current_Tag[LHQL-1-4]<<8 | (Current_Tag[LHQL-1-5]<<10);//Use lookup table...
				}

				Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;//Range.Tag=Actual_Tag;
				memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
			}
			Search_Backwards_Exact_X0( Range,LHQL,LHQRrx,fwfmi);// ?|0|?
			Range.Level=1;
			if(USEQUALITY)
			{
				TRange=Range;
				Do_Branch=Low_Quality;
				Search_10LX(TRange,1,LHQLrx, LHQLrx,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
				if(MAXHITS==Hits) continue;
			}
			Do_Branch=Do_All;
			Search_10LX(Range,1,LHQLrx, LHQLrx,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
			if(MAXHITS==Hits) continue;

			if (NPOLICY && NCount){if (!NISMISMATCH && MAX_MISMATCHES + NCount == 4) continue;else if (MAX_MISMATCHES==NCount) continue;}
			else if( MAX_MISMATCHES ==4) continue;
			if (FILTERUNIQUEHITS && Hits) continue;
//}------------------------------------------- FOUR MISMATCH ---------------------------------------------------------------------------------------------
//{------------------------------------------- FIVE MISMATCH ---------------------------------------------------------------------------------------------

			if (BESTHIT && Rollover_Step && Last_Mismatch == 5) {In_Mismatch=6;continue;}
			In_Mismatch=5;
			FMIndex=REVERSE;
			if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
			{
				for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
				{
					Print_LocationX(Two_Mismatches_At_End_Forward[i]);//mismatches of the form 0|4, with last mismatch at the end...
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}
			Two_Mismatches_At_End_Forward_Pointer=0;


			FMIndex=FORWARD;
			if(Two_Mismatches_At_End_Pointer)
			{
				for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
				{
					Print_LocationX(Two_Mismatches_At_End[i]);//mismatches of the form 0|5, with one mismatch at the start...
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}
			Two_Mismatches_At_End_Pointer=0;

			FMIndex=REVERSE;
			if(Mismatches_Forward_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=Possible_05_Pointer;i<Mismatches_Forward_Pointer_Last4;i++)//Mismatches_Forward_Pointer;i++)
					{
						Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
						if(MAXHITS==Hits) break;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=Possible_05_Pointer;i<Mismatches_Forward_Pointer_Last4;i++)//Mismatches_Forward_Pointer;i++)
				{
					Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;

			}


			FMIndex=REVERSE;
			if(Mismatches_Forward_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=Mismatches_Forward_Pointer_Last4;i<Mismatches_Forward_Pointer_Last5;i++)//Mismatches_Forward_Pointer;i++)
					{
						Search_Forwards(Mismatches_Forward[i],5,1,STRINGLENGTH,revfmi);//scan for possible five mismatches of the form 0|5, and finds mismatches of the form 1|4,2|3 
						if(MAXHITS==Hits) break;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=Mismatches_Forward_Pointer_Last4;i<Mismatches_Forward_Pointer_Last5;i++)//Mismatches_Forward_Pointer;i++)
				{
					Search_Forwards(Mismatches_Forward[i],5,1,STRINGLENGTH,revfmi);//scan for possible five mismatches of the form 0|5, and finds mismatches of the form 1|4,2|3 
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}

			FMIndex=REVERSE;
			if(Mismatches_Forward_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=Mismatches_Forward_Pointer_Last5;i<Mismatches_Forward_Pointer;i++)//Mismatches_Forward_Pointer;i++)
					{
						Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
						if(MAXHITS==Hits) break;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=Mismatches_Forward_Pointer_Last5;i<Mismatches_Forward_Pointer;i++)//Mismatches_Forward_Pointer;i++)
				{
					Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;
			}

			FMIndex=FORWARD;
			if(Mismatches_Backward_Pointer!=Possible_50_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=Mismatches_Backward_Pointer-1;i>=Possible_50_Pointer;i--)
					{
						Search_Backwards(Mismatches_Backward[i],5,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
						if(MAXHITS==Hits) break;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=Mismatches_Backward_Pointer-1;i>=Possible_50_Pointer;i--)
				{
					Search_Backwards(Mismatches_Backward[i],5,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) continue;

			}
			FMIndex=FORWARD;
			if (Possible_20_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=Possible_20_Pointer-1;i>=0;i--)
					{
						if(Possible_20[i].Level!=RHQL) 
						{
							Possible_20[i].Level++;
							Search_Backwards_Exact( Possible_20[i],LH+RHQL,RHQL,fwfmi);
						}

						if(Possible_20[i].Start)
						{	
							Possible_20[i].Level=1;
							Search_Backwards(Possible_20[i],5,LH,LH,fwfmi);
							if(MAXHITS==Hits) break;
						}
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=Possible_20_Pointer-1;i>=0;i--)
				{
					if(Possible_20[i].Level!=RHQL) 
					{
						Possible_20[i].Level++; 
						Search_Backwards_Exact( Possible_20[i],LH+RHQL,RHQL,fwfmi);
					}
					
					if(Possible_20[i].Start)
					{	
						Possible_20[i].Level=1;
						Search_Backwards(Possible_20[i],5,LH,LH,fwfmi);
						if(MAXHITS==Hits) break;
					}
				}
				if(MAXHITS==Hits) continue;

			}

			if(Possible_02_Pointer)
			{
				if(USEQUALITY)
				{
					Do_Branch=Low_Quality;
					for(int i=0;i<Possible_02_Pointer;i++)
					{
						TRange=Possible_02[i];
						if(Possible_02[i].Level!=RHQR) 
						{
							Search_Exact( Possible_02[i],LH + RHQL-1 ,RHQR,revfmi);//finds mismatches of the form 202, stores possibles of the form 203
						}
						if(Possible_02[i].Start) 
						{
							Reverse(Possible_02[i],STRINGLENGTH,RH);
							if(Possible_02[i].Start)
							{
								Possible_02[i].Level=1;
								Search_Backwards(Possible_02[i],5,LH,LH,fwfmi);
							}
						}
						if(MAXHITS==Hits) break;
						Possible_02[i]=TRange;
					}
					if(MAXHITS==Hits) continue;
				}

				Do_Branch=Do_All;
				for(int i=0;i<Possible_02_Pointer;i++)
				{

					if(Possible_02[i].Level!=RHQR) 
					{
						Search_Exact( Possible_02[i],LH + RHQL-1 ,RHQR,revfmi);//finds mismatches of the form 202, stores possibles of the form 203
					}
					if(Possible_02[i].Start) 
					{
						Reverse(Possible_02[i],STRINGLENGTH,RH);
						if(Possible_02[i].Start)
						{
							Possible_02[i].Level=1;
							Search_Backwards(Possible_02[i],5,LH,LH,fwfmi);
						}
					}
					if(MAXHITS==Hits) break;
				}

				if(MAXHITS==Hits) continue;
			}

			int RHQLlx=RHQR/2;
			if (RHQR % 2) RHQLlx++; int RHQRlx=RHQR-RHQLlx;
			if (LOOKUPSIZE >= RHQRlx)
			{
				Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
				memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
			}
			else
			{
				if(LOOKUPSIZE==3)
				{
					c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
				}
				else
				{
					c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
				}
				Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];
				Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;//Range.Tag=Actual_Tag;
				memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
			}

			Search_Backwards_Exact_X0( Range,STRINGLENGTH,RHQRlx,fwfmi);// ?|?|0
			Range.Level=1;

			if(USEQUALITY)
			{
				TRange=Range;
				Do_Branch=Low_Quality;
				Search_Backwards_XL10(TRange,1,LH + RHQL+RHQLlx, RHQLlx,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
				if(MAXHITS==Hits) continue;
			}

			Do_Branch=Do_All;
			Search_Backwards_XL10(Range,1,LH + RHQL+RHQLlx, RHQLlx,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
			if(MAXHITS==Hits) continue;



			if (LOOKUPSIZE >= RHQLlx)
			{
				Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
				memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
			}
			else
			{
				if(LOOKUPSIZE==3)
				{
					c=Current_Tag[LH+RHQL+0] | (Current_Tag[LH+RHQL+1]<<2) | (Current_Tag[LH+RHQL+2]<<4);// | (Current_Tag[LH+3]<<6) | Current_Tag[LH+4]<<8 | (Current_Tag[LH+5]<<10);//Use lookup table...
				}
				else
				{
					c=Current_Tag[LH+RHQL+0] | (Current_Tag[LH+RHQL+1]<<2) | (Current_Tag[LH+RHQL+2]<<4) | (Current_Tag[LH+RHQL+3]<<6) | Current_Tag[LH+RHQL+4]<<8 | (Current_Tag[LH+RHQL+5]<<10);//Use lookup table...
				}
				Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];
				Range.Mismatches=0;Range.Level=LOOKUPSIZE+1; Range.Skip=0;
				memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
			}

			Search_Forwards_0X(Range,LH+RHQL+1,RHQLlx,revfmi);
			Range.Level=1;
			TRange=Range;
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;

				Search_XL01(Range,1,LH + RHQL+RHQLlx +1,RHQRlx,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
				if(MAXHITS==Hits) continue;
			}

			Do_Branch=Do_All;
			Search_XL01(TRange,1,LH + RHQL+RHQLlx +1,RHQRlx,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
			if(MAXHITS==Hits) continue;

			if (FILTERUNIQUEHITS && Hits) continue;
		}

//}------------------------------------------- FIVE MISMATCH ---------------------------------------------------------------------------------------------

//{------------------------------------------- MULTI MISMATCH ---------------------------------------------------------------------------------------------
		if (NPOLICY && NCount){if (!NISMISMATCH && MAX_MISMATCHES + NCount == 5) continue;else if (MAX_MISMATCHES==NCount) continue;}
		if (MAX_MISMATCHES >5 && !Hits)
		{
			if (SCANCOMPONLY)
			{
				Guessed=Original+IGNOREHEAD;
				Guessed_Quality=Low_QualityF+IGNOREHEAD;
				Guess_Complement=Complement;
			}
			MAX_ORG=MAX_MISMATCHES;
			int AB=ARRAY_BOUND;int EB=END_BOUND;
			if (ARRAY_BOUND_BD) {ARRAY_BOUND=ARRAY_BOUND_BD;END_BOUND=ARRAY_BOUND_BD;}
			unsigned TTotal_Hits=Total_Hits,TTotal_Tags=Total_Tags;
			S=SATot=0;
			In_Mismatch=6;
			Ext_Scan=0;
			Current_Tag=Guessed;
			Low_Quality=Guessed_Quality;
			Extentions=0;

			HIGHSCAN=TRUE;
			Last_Mismatch_Written=0;
			int H=MAXHITS;if(!SUPERACCURATE) MAXHITS=INT_MAX;//10000;
			Current_Tag=Guessed;
			Low_Quality=Guessed_Quality;
			char L;
			if (Low_Quality[STRINGLENGTH]) L=TRUE; else L=FALSE;
			while(Hits < MAXHITS && Ext_Scan<2) 
			{

				Ext_Scan++;
				if (Ext_Scan==1)
				{
					if (L) 
					{
						Current_Tag=Guessed;
						Low_Quality=Guessed_Quality;
						Left_To_Right();
						if(Hits == MAXHITS) break;
						if(SCANCOMPONLY) continue; 
						Current_Tag=Guess_Complement;
						Low_Quality=Guessed_Complement_Quality;
						Right_To_Left();
					}
					else
					{
						Current_Tag=Guessed;
						Low_Quality=Guessed_Quality;
						Right_To_Left();
						if(Hits==MAXHITS) break;
						if(SCANCOMPONLY) continue; 
						Current_Tag=Guess_Complement;
						Low_Quality=Guessed_Complement_Quality;
						Left_To_Right();
					}
				}
				else	
				{
					if (L) 
					{
						Current_Tag=Guessed;
						Low_Quality=Guessed_Quality;
						Right_To_Left();
						if(Hits == MAXHITS) break;
						if(SCANCOMPONLY) continue; 
						Current_Tag=Guess_Complement;
						Low_Quality=Guessed_Complement_Quality;
						Left_To_Right();
					}
					else
					{
						Current_Tag=Guessed;
						Low_Quality=Guessed_Quality;
						Left_To_Right();
						if(Hits==MAXHITS) break;
						if(SCANCOMPONLY) continue; 
						Current_Tag=Guess_Complement;
						Low_Quality=Guessed_Complement_Quality;
						Right_To_Left();
					}
				}


				//Low_Quality=Guessed_Complement_Quality;
			}

			if (Tag_Stat_Bad)
			{
				Total_Hits=TTotal_Hits;Total_Tags=TTotal_Tags;
			}

			Ext_Scan=0;
			//int Ext_Scan=0;
			Current_Tag=Guessed;
			Low_Quality=Guessed_Quality;
			Larger_Than_Ten = TRUE;
			if (!Hits || (Last_Mismatch_Written >12))
			{
				if (!Tag_Stat_Bad && MAX_MISMATCHES >10)
				{
					while(Hits < MAXHITS && Ext_Scan<2) 
					{
						Ext_Scan++;

						Right_To_LeftX();
						if(Hits==MAXHITS) break;
						Left_To_RightX();
						if(Hits==MAXHITS) break;
						Left_To_RightY();
						if(Hits==MAXHITS) break;
						Right_To_LeftY();

						Current_Tag=Guess_Complement;
						Low_Quality=Guessed_Complement_Quality;
					}
				}
			}
			Larger_Than_Ten = FALSE;
			if(ARRAY_BOUND_BD) {ARRAY_BOUND=AB;END_BOUND=EB;}
			MAXHITS=H;MAX_MISMATCHES=MAX_ORG;
			if(PRINT_MISHITS )
			{
				if (Tag_Stat_Bad) fprintf(Mishit_File,"%s%sBAD%s%s", Description,Tag_Copy,Plus,Quality);
				else if(!Hits) fprintf(Mishit_File,"%s%s%s%s", Description,Tag_Copy,Plus,Quality);
			}
		}

//}------------------------------------------- MULTI MISMATCH ---------------------------------------------------------------------------------------------
	}
	if(HITMODE == DEFAULT && !GIS)
	{
		Record.Start=0;
		fwriteX(&Record,sizeof(Record),1,Output_File);//write sentinel..
	}
	else
	{
		char End_Mark='&';
		fwriteX(&End_Mark,1,1,Output_File);//write sentinel..
	}	

	printf("\r[++++++++100%%+++++++++]\n");fprintf(Log_File,"100%%\n");//progress bar....
	if(FILTERUNIQUEHITS) Total_Tags++;
	if(NORMAL_TAGS) 
	{
		printf("Total Tags/Hits : %d/%d\t Tags parsed : %d\n",Total_Tags,Total_Hits,Actual_Tag); 
		fprintf(Log_File,"Total Tags/Hits : %d/%d\t Tags parsed : %d\n",Total_Tags,Total_Hits,Actual_Tag); 
	}
	else 
	{
		printf("Total Tags/Hits : %d/%d\t Tags parsed : %d [%d]\n",Total_Tags,Total_Hits,Actual_Tag,Actual_Tag/2);
		fprintf(Log_File,"Total Tags/Hits : %d/%d\t Tags parsed : %d [%d]\n",Total_Tags,Total_Hits,Actual_Tag,Actual_Tag/2);
	}

	time(&End_Time);printf("\n Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));fprintf(Log_File,"\n Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	if(OUTPUT_ZIPPED) gzclose(Output_File);
	if (EMAIL) Emailer(argc,argv);



}

//------------------multi mismatch extensions -----------
void Right_To_Left()
{
	Current_Tag += LH;
	STRINGLENGTHt=STRINGLENGTH;LHt=LH;RHt=RH;LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
	STRINGLENGTH=STRINGLENGTHl;LH=LHl;RH=RHl;LHQL=LHQLl;LHQR=LHQRl;RHQL=RHQLl;RHQR=RHQRl;
	Extending_Tag=RIGHTPASS;
	Extend_Left_Half();
	Extending_Tag=FALSE;
	STRINGLENGTH=STRINGLENGTHt;LH=LHt;RH=RHt;LHQL=LHQLt;LHQR=LHQRt;RHQL=RHQLt;RHQR=RHQRt;
	Current_Tag -= LH; 
}

void Left_To_Right()
{
	STRINGLENGTHt=STRINGLENGTH;LHt=LH;RHt=RH;LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
	STRINGLENGTH=STRINGLENGTHr;LH=LHr;RH=RHr;LHQL=LHQLr;LHQR=LHQRr;RHQL=RHQLr;RHQR=RHQRr;
	Extending_Tag=LEFTPASS;
	Extend_Left_Half();
	Extending_Tag=FALSE;
	STRINGLENGTH=STRINGLENGTHt;LH=LHt;RH=RHt;LHQL=LHQLt;LHQR=LHQRt;RHQL=RHQLt;RHQR=RHQRt;
}

void Left_To_RightX()//super high mis on left..
{
	STRINGLENGTHt=STRINGLENGTH;LHt=LH;RHt=RH;LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
	STRINGLENGTHt1=STRINGLENGTHr;LHt1=LHr;RHt1=RHr;LHQLt1=LHQLr;LHQRt1=LHQRr;RHQLt1=RHQLr;RHQRt1=RHQRr;

	STRINGLENGTH=LHr;
	LH=STRINGLENGTH/2;//calculate tag portions...
	LHQL=LH/2;
	if ((STRINGLENGTH % 2)) {LH++;LHQL++;}	
	LHQR=LH-LHQL;
	RH=STRINGLENGTH-LH;
	RHQL=RH/2;RHQR=RH-RHQL;
	Extending_Tag=LEFTPASS1;

	Extend_Left_Half();
	Extending_Tag=FALSE;
	STRINGLENGTH=STRINGLENGTHt;LH=LHt;RH=RHt;LHQL=LHQLt;LHQR=LHQRt;RHQL=RHQLt;RHQR=RHQRt;
}

void Right_To_LeftX()//search right half of left half for >5 and extend
{
	STRINGLENGTHt=STRINGLENGTH;LHt=LH;RHt=RH;LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
	STRINGLENGTHt1=STRINGLENGTHr;LHt1=LHr;RHt1=RHr;LHQLt1=LHQLr;LHQRt1=LHQRr;RHQLt1=RHQLr;RHQRt1=RHQRr;

	STRINGLENGTH=RHr;
	LH=STRINGLENGTH/2;//calculate tag portions...
	LHQL=LH/2;
	if ((STRINGLENGTH % 2)) {LH++;LHQL++;}	
	LHQR=LH-LHQL;
	RH=STRINGLENGTH-LH;
	RHQL=RH/2;RHQR=RH-RHQL;
	STRINGLENGTHt2=STRINGLENGTH;LHt2=LH;RHt2=RH;LHQLt2=LHQL;LHQRt2=LHQR;RHQLt2=RHQL;RHQRt2=RHQR;
	Current_Tag += LHr;

	Extending_Tag=RIGHTPASS1;
	Extend_Left_Half();
	Extending_Tag=FALSE;
	STRINGLENGTH=STRINGLENGTHt;LH=LHt;RH=RHt;LHQL=LHQLt;LHQR=LHQRt;RHQL=RHQLt;RHQR=RHQRt;
	Current_Tag -= LHr; 
}

//-------------- ---------------
void Right_To_LeftY()
{
	STRINGLENGTHt=STRINGLENGTH;LHt=LH;RHt=RH;LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
	STRINGLENGTHt1=STRINGLENGTHl;LHt1=LHl;RHt1=RHl;LHQLt1=LHQLl;LHQRt1=LHQRl;RHQLt1=RHQLl;RHQRt1=RHQRl;

	Current_Tag += LH;
	STRINGLENGTH=RHl;
	LH=STRINGLENGTH/2;//calculate tag portions...
	LHQL=LH/2;
	if ((STRINGLENGTH % 2)) {LH++;LHQL++;}	
	LHQR=LH-LHQL;
	RH=STRINGLENGTH-LH;
	RHQL=RH/2;RHQR=RH-RHQL;
	STRINGLENGTHt2=STRINGLENGTH;LHt2=LH;RHt2=RH;LHQLt2=LHQL;LHQRt2=LHQR;RHQLt2=RHQL;RHQRt2=RHQR;
	Current_Tag += LHl;

	Extending_Tag=RIGHTPASS1Y;
	Extend_Left_Half();
	Extending_Tag=FALSE;
	STRINGLENGTH=STRINGLENGTHt;LH=LHt;RH=RHt;LHQL=LHQLt;LHQR=LHQRt;RHQL=RHQLt;RHQR=RHQRt;
	Current_Tag -= LHl; 
	Current_Tag -= LH;
}
//-------------- ---------------
void Left_To_RightY()
{
	STRINGLENGTHt=STRINGLENGTH;LHt=LH;RHt=RH;LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
	STRINGLENGTHt1=STRINGLENGTHl;LHt1=LHl;RHt1=RHl;LHQLt1=LHQLl;LHQRt1=LHQRl;RHQLt1=RHQLl;RHQRt1=RHQRl;

	Current_Tag += LH;
	STRINGLENGTH=LHl;
	LH=STRINGLENGTH/2;//calculate tag portions...
	LHQL=LH/2;
	if ((STRINGLENGTH % 2)) {LH++;LHQL++;}	
	LHQR=LH-LHQL;
	RH=STRINGLENGTH-LH;
	RHQL=RH/2;RHQR=RH-RHQL;
	Extending_Tag=LEFTPASS1Y;

	Extend_Left_Half();
	Extending_Tag=FALSE;
	STRINGLENGTH=STRINGLENGTHt;LH=LHt;RH=RHt;LHQL=LHQLt;LHQR=LHQRt;RHQL=RHQLt;RHQR=RHQRt;
	Current_Tag -= LH;
}

//------------------multi mismatch extensions -----------

void Emailer(int argc,char* argv[])
{
	gzFile Email;
	char* Command_Line,*Command_Ptr;
	if(MAILCLIENT == USER)
	{
		File_OpenZ(EMAIL_CFG,"r",Email);
		Command_Line=Command_Ptr=(char*)calloc(1,5000);
		gzgets(Email,Command_Ptr,4000);//Read client
		while(*Command_Ptr!= '\r' && *Command_Ptr!='\n'&&*Command_Ptr!=0) Command_Ptr++;*(Command_Ptr++)=' ';*(Command_Ptr++)='"';
		gzgets(Email,Command_Ptr,4000);//To:
		while(*Command_Ptr!= '\r' && *Command_Ptr!='\n'&&*Command_Ptr!=0) Command_Ptr++;*(Command_Ptr++)='"';*Command_Ptr=' ';Command_Ptr++;*(Command_Ptr++)='"';
		gzgets(Email,Command_Ptr,4000);//Subject:
		while(*Command_Ptr!= '\r' && *Command_Ptr!='\n'&&*Command_Ptr!=0) Command_Ptr++;*(Command_Ptr++)='"';*Command_Ptr=' ';Command_Ptr++;*(Command_Ptr++)='"';
		Command_Ptr +=gzread(Email,Command_Ptr,4000);//Body
		while(*Command_Ptr!= '\r' && *Command_Ptr!='\n'&&*Command_Ptr!=0) Command_Ptr++;*(Command_Ptr++)='"';*Command_Ptr=' ';Command_Ptr++;

		*(Command_Ptr++)='"';
		for (int i = 1; i < argc; ) 
		{
			Command_Ptr += sprintf(Command_Ptr,"%s", argv[i]);
			if (++i < argc) *(Command_Ptr++)=' ';
		}
		*(Command_Ptr++)='"';

		while(*Command_Ptr!= '\r' && *Command_Ptr!='\n'&&*Command_Ptr!=0) Command_Ptr++;
	}
	else if (MAILCLIENT == MAIL)
	{
		File_OpenZ(EMAIL_CFG,"r",Email);
		Command_Line=Command_Ptr=(char*)calloc(1,5000);
		char* Subj_Ptr=(char*)calloc(1,5000);
		char* To_Ptr=(char*)calloc(1,5000);
		//strcpy(Command_Line,"mail -s");
		strcpy(Command_Line,"echo -e \"");
		while(*Command_Ptr!= '"') Command_Ptr++;Command_Ptr++;
		gzgets(Email,To_Ptr,5000);//To:
		gzgets(Email,Subj_Ptr,5000);//Subject:
		Command_Ptr +=gzread(Email,Command_Ptr,4000);//Body

		for (int i = 1; i < argc; ) 
		{
			Command_Ptr += sprintf(Command_Ptr,"%s", argv[i]);
			if (++i < argc) *(Command_Ptr++)=' ';
		}

		strcpy(Command_Ptr,"\"|mail -s \"");Command_Ptr +=11;
		strcpy(Command_Ptr,Subj_Ptr);
		while(*Command_Ptr!= '\r' && *Command_Ptr!='\n'&&*Command_Ptr!=0) Command_Ptr++;*(Command_Ptr++)='"';*(Command_Ptr++)=' ';*(Command_Ptr++)='"';
		strcpy(Command_Ptr,To_Ptr);
		while(*Command_Ptr!= '\r' && *Command_Ptr!='\n'&&*Command_Ptr!=0) Command_Ptr++;*(Command_Ptr++)='"';
		
	}
	else 
	{
		File_OpenZ(EMAIL_CFG,"r",Email);
		Command_Line=Command_Ptr=(char*)calloc(1,5000);
		char* Subj_Ptr=(char*)calloc(1,5000);
		char* To_Ptr=(char*)calloc(1,5000);
		strcpy(Command_Line,"echo -e \"");
		while(*Command_Ptr!= '"') Command_Ptr++;Command_Ptr++;
		To_Ptr[0]='T';To_Ptr[1]='o';To_Ptr[2]=':';gzgets(Email,To_Ptr+3,5000);//To:
		Subj_Ptr[0]='S';Subj_Ptr[1]='u';Subj_Ptr[2]='b';Subj_Ptr[3]='j';Subj_Ptr[4]='e';Subj_Ptr[5]='c';Subj_Ptr[6]='t';Subj_Ptr[7]=':';gzgets(Email,Subj_Ptr+8,5000);//Subject:
		strcpy(Command_Ptr,To_Ptr);
		while(*Command_Ptr!= '\r' && *Command_Ptr!='\n'&&*Command_Ptr!=0) Command_Ptr++;*(Command_Ptr++)='\n';
		strcpy(Command_Ptr,Subj_Ptr);
		while(*Command_Ptr!= '\r' && *Command_Ptr!='\n'&&*Command_Ptr!=0) Command_Ptr++;*(Command_Ptr++)='\n';*(Command_Ptr++)='\n';
		Command_Ptr +=gzread(Email,Command_Ptr,4000);//Body

		for (int i = 1; i < argc; ) 
		{
			Command_Ptr += sprintf(Command_Ptr,"%s", argv[i]);
			if (++i < argc) *(Command_Ptr++)=' ';
		}

		strcpy(Command_Ptr,"\n.\n\"|sendmail -t");
		
	}
	system (Command_Line);

}

void Print_Hits()
{
	char Print_Head=FALSE;
	char* Filter_Buffer=Write_Buffer;

	if (FILTERUNIQUEHITS) 
	{
		if (Write_Buffer_Ptr == Write_Buffer) {Tag_Stat_Bad=TRUE;Total_Tags--;}
		else if (BESTHIT)
		{
			Tag_Printed=TRUE;
			if (Last_Mismatch > In_Mismatch)
			{
				Filter_Buffer=Last_Write_Buffer;
			}
			else
			{
				Write_Buffer_Ptr=Last_Write_Buffer;
			}
			Total_Hits++;
		}
		else
		{
			Total_Hits++;
			if(Last_Write_Buffer!=Write_Buffer && Write_Buffer_Ptr!= Last_Write_Buffer) Total_Hits++; 
		}
	}
	if (Tag_Stat_Bad) {Tag_Printed=FALSE;if (PRINTBLANKS) Print_Head=TRUE;}
	else
	{
		if(!Head_Printed)
		{
			if (Tag_Printed) Print_Head=TRUE;
			else if (NORMAL_TAGS)
			{
				if (PRINTBLANKS) Print_Head=TRUE;
			}
			else
			{
				if((2==Tag_Number && PRINTBLANKS))  Print_Head=TRUE;
			}
		}
	}
	
	if(Print_Head)
	{
		if (SLIDERBLAST) 
		{
			SLIDERHIT++;
			fprintfX(BLAST_File,"%s",Description);
			fprintfX(BLAST_File,"%s",Tag_Copy);
			if(FILETYPE == FQ) fprintfX(BLAST_File,"+\n%s",Quality); 
		}
		fprintfX(Output_File,"@%s",Description);
		fwriteX(Stats1,sizeof(unsigned short),Stat_Size,Output_File);
		if(!NORMAL_TAGS)
		{
			fwriteX(Stats2,sizeof(unsigned short),Stat_Size,Output_File);
		}
		fwriteX(Tag_Copy,1,TAG_COPY_LEN,Output_File);//Write_Buffer_Ptr += TAG_COPY_LEN;//write tag...
		if (FILETYPE==FQ) 
		{
			if(NORMAL_TAGS) fwriteX(Quality,1,TAG_COPY_LEN,Output_File);
			else fwriteX(Quality_Copy,1,TAG_COPY_LEN,Output_File);
		}
		Head_Printed=TRUE;
		
	}

	if(Tag_Printed)
	{
		fwriteX(Filter_Buffer,Write_Buffer_Ptr-Filter_Buffer,1,Output_File);
		//else fwriteX(Write_Buffer,Write_Buffer_Ptr-Write_Buffer,1,Output_File);
	}
        Write_Buffer_Ptr=Write_Buffer;
	Last_Write_Buffer=Write_Buffer_Ptr;
}

char Guess_Orientation()
{
	int Complement_Length;
	int Org_Length;
	Hit_In=0;
	SARange Range;
	FMIndex=REVERSE;

//-------------------------------------------------- NORMAL SCAN --------------------------------------------------------------------
	Exact_Match_Forward=Exact_Match_ForwardF;
	Current_Tag=Original+IGNOREHEAD;
	Direction = 0;//forward scn first...
	DIRECTIONS=2;

	if(LOOKUPSIZE ==3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}

	Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Tag=Actual_Tag;
	Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;Range.Skip=0;
	memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
	Search_Forwards_Exact(Range,-1,STRINGLENGTH,revfmi);//Find exact matches and report... if not found get the range for 0|?
	if(Hits) 
	{
		Delta=STRINGLENGTH;
		Guessed=Original+IGNOREHEAD;
		Guessed_Quality=Low_QualityF+IGNOREHEAD;
		Guess_Complement=Complement;
		Guessed_Complement_Quality=Low_QualityC+IGNOREHEAD;
		if (MAX_MISMATCHES != 0) {Hit_In=1;return FALSE;}
		if (!ROLLOVER ) return FALSE;
	}
	Org_Length=Range.Level;
//--------------------------------------------------REVERSE COMPLEMENT SCAN --------------------------------------------------------------------
	Exact_Match_Forward=Exact_Match_ForwardC;
	Current_Tag=Complement;
	//NLocations=NLocationsC;
	Direction =1;//reverse complement scan first...
	DIRECTIONS=3;

	if(LOOKUPSIZE ==3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}

	Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;Range.Tag=Actual_Tag;
	memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
	Search_Forwards_Exact(Range,-1,STRINGLENGTH,revfmi);//Find exact matches and report... if not found get the range for 0|?
	if(Hits || MAX_MISMATCHES == 0) {Hit_In=2;return FALSE;}
	Complement_Length=Range.Level;
//-------------------------------------------------- EXACT SCAN END --------------------------------------------------------------------
	if (Complement_Length > Org_Length)//first try the complement
	{
		Max_Seeked=Complement_Length;
		Delta=Complement_Length-Org_Length;
		Guessed=Complement;
		//Guessed_NLocation=NLocationsC;
		Guessed_Quality=Low_QualityC+IGNOREHEAD;
		Guess_Complement=Original+IGNOREHEAD;
		Guessed_Complement_Quality=Low_QualityF+IGNOREHEAD;
		Direction =1;
		DIRECTIONS=3;
	}
	else
	{

		Delta=Org_Length-Complement_Length;
		Guessed=Original+IGNOREHEAD;
		Guessed_Quality=Low_QualityF+IGNOREHEAD;
		Guess_Complement=Complement;
		//Guessed_NLocation_Complement=NLocationsC;
		Guessed_Complement_Quality=Low_QualityC+IGNOREHEAD;
		if (FILTER_AMBIGUOUS && Complement_Length == Org_Length )
		{
			Current_Tag=Original+IGNOREHEAD;
			for (int i=0;i<=STRINGLENGTH-1;i++)//convert and write fasta
			{
				Current_Tag[i]=Code_To_Char[Current_Tag[i]];
			}
			Current_Tag[STRINGLENGTH]=0;
			fprintf(Ambiguous_File,"%s%s\n+\n%s", Description,Current_Tag,Quality);
			AMBI=TRUE;
			return FALSE;
		}
		Current_Tag=Original+IGNOREHEAD;
		Exact_Match_Forward=Exact_Match_ForwardF;
		Direction = 0;
		DIRECTIONS=2;
	}
	return TRUE;

}


void Build_Preindex_Backward(Range Range, int Level, int Bit)
{

	if (LOOKUPSIZE==Level) 
	{
		assert(Range.Start<=Range.End);
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label
		if(Range.Start)
		{
			Range.Start= fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.Start , Bit) + 1;
			Range.End= fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.End+1, Bit);
			if (Range.Start>Range.End) Range.Start=0;
		}
		Backward_Start_LookupX[Range.Label]=Range.Start;
		Backward_End_LookupX[Range.Label]=Range.End;
		assert(Backward_Start_LookupX[Range.Label]<=Backward_End_LookupX[Range.Label]);
	}
	else
	{


		assert(Range.Start<=Range.End);
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label 
		if(Range.Start)
		{
			Range.Start = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.Start , Bit) + 1;
			Range.End = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.End+1, Bit);
			if (Range.Start>Range.End) Range.Start=0;
		}
		Level ++;
		for ( int i=0;i<4;i++)
		{
			Build_Preindex_Backward( Range, Level,i);
		}

	}

}

void Build_Preindex_Forward(Range Range, int Level, int Bit)
{

	if (LOOKUPSIZE==Level) 
	{
		assert(Range.Start<=Range.End);
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label
		if(Range.Start)
		{
			Range.Start= revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.Start , Bit) + 1;
			Range.End= revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.End+1, Bit);
			if (Range.Start>Range.End) Range.Start=0;
		}
		Forward_Start_LookupX[Range.Label]=Range.Start;
		Forward_End_LookupX[Range.Label]=Range.End;
		assert(Forward_Start_LookupX[Range.Label] <=Forward_End_LookupX[Range.Label]);
	}
	else
	{


		assert(Range.Start<=Range.End);
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label 
		if (Range.Start)
		{
			Range.Start = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.Start , Bit) + 1;
			Range.End = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.End+1, Bit);
			if (Range.Start>Range.End) Range.Start=0;
		}
		Level ++;
		for ( int i=0;i<4;i++)
		{
			Build_Preindex_Forward( Range, Level,i);
		}

	}

}

void Seed_Scan()
{
	int Ext_Scan=0;
	Current_Tag=Guessed;
	Low_Quality=Guessed_Quality;
	Do_Branch=Do_All;

	STRINGLENGTHl=SEEDSIZE;//calculate tag portions...
	LHl=STRINGLENGTHl/2;//calculate tag portions...
	LHQLl=LHl/2;
	if ((STRINGLENGTHl % 2)) {LHl++;LHQLl++;}	
	LHQRl=LHl-LHQLl;
	RHl=STRINGLENGTHl-LHl;
	RHQLl=RHl/2;RHQRl=RHl-RHQLl;

	STRINGLENGTHr=SEEDSIZE;//calculate tag portions...
	LHr=STRINGLENGTHr/2;//calculate tag portions...
	LHQLr=LHr/2;
	if ((STRINGLENGTHr % 2)) {LHr++;LHQLr++;}	
	LHQRr=LHr-LHQLr;
	RHr=STRINGLENGTHr-LHr;
	RHQLr=RHr/2;RHQRr=RHr-RHQLr;

	while(Hits<MAXHITS && Ext_Scan<2) 
	{
		Ext_Scan++;
		if (Current_Tag[STRINGLENGTHt] == '+')//plus has more chance of errors at the end..
		{
			//STRINGLENGTHt=STRINGLENGTH;LHt=LH;RHt=RH;LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
			STRINGLENGTHt=STRINGLENGTH;RHt=SEEDSIZE;//LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
			STRINGLENGTH=STRINGLENGTHr;LH=LHr;RH=RHr;LHQL=LHQLr;LHQR=LHQRr;RHQL=RHQLr;RHQR=RHQRr;
			Extending_Tag=LEFTPASS;
			Extend_Left_Half();
			Extending_Tag=FALSE;
			STRINGLENGTH=STRINGLENGTHt;//LH=LHt;RH=RHt;LHQL=LHQLt;LHQR=LHQRt;RHQL=RHQLt;RHQR=RHQRt;
			if(MAXHITS==Hits) break;
		}
		else
		{

			Current_Tag += STRINGLENGTH-SEEDSIZE;
			//STRINGLENGTHt=STRINGLENGTH;LHt=LH;RHt=RH;LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
			STRINGLENGTHt=STRINGLENGTH;LHt=SEEDSIZE;//STRINGLENGTH-30;RHt=30;//LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
			STRINGLENGTH=STRINGLENGTHr;LH=LHr;RH=RHr;LHQL=LHQLr;LHQR=LHQRr;RHQL=RHQLr;RHQR=RHQRr;
			Extending_Tag=RIGHTPASS;
			Extend_Left_Half();
			Extending_Tag=FALSE;
			STRINGLENGTH=STRINGLENGTHt;//LH=LHt;RH=RHt;LHQL=LHQLt;LHQR=LHQRt;RHQL=RHQLt;RHQR=RHQRt;
			Current_Tag -= LH; 
			if(MAXHITS==Hits) break;
		}
		//if(Delta < 5)
		{
			Current_Tag=Guess_Complement;
		}
		//else Ext_Scan=2;

	}

}

void High_Mismatch_Scan()
{
	In_Mismatch=6;
	int Ext_Scan=0;
	Current_Tag=Guessed;
	Low_Quality=Guessed_Quality;
	while(Hits<MAXHITS && Ext_Scan<2) 
	{
		Ext_Scan++;
		STRINGLENGTHt=STRINGLENGTH;LHt=LH;RHt=RH;LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
		STRINGLENGTH=STRINGLENGTHr;LH=LHr;RH=RHr;LHQL=LHQLr;LHQR=LHQRr;RHQL=RHQLr;RHQR=RHQRr;
		Extending_Tag=LEFTPASS;
		Extend_Left_Half();
		Extending_Tag=FALSE;
		STRINGLENGTH=STRINGLENGTHt;LH=LHt;RH=RHt;LHQL=LHQLt;LHQR=LHQRt;RHQL=RHQLt;RHQR=RHQRt;
		if(MAXHITS==Hits) break;

		if(SUPERACCURATE)
		{
			Current_Tag += LH;
			STRINGLENGTHt=STRINGLENGTH;LHt=LH;RHt=RH;LHQLt=LHQL;LHQRt=LHQR;RHQLt=RHQL;RHQRt=RHQR;
			STRINGLENGTH=STRINGLENGTHl;LH=LHl;RH=RHl;LHQL=LHQLl;LHQR=LHQRl;RHQL=RHQLl;RHQR=RHQRl;
			Extending_Tag=RIGHTPASS;
			Extend_Left_Half();
			Extending_Tag=FALSE;
			STRINGLENGTH=STRINGLENGTHt;LH=LHt;RH=RHt;LHQL=LHQLt;LHQR=LHQRt;RHQL=RHQLt;RHQR=RHQRt;
			Current_Tag -= LH; 
			if(MAXHITS==Hits) break;

			if (MAX_MISMATCHES >10) 
			{
				int Brute_Large = MAX_MISMATCHES/2;
				Brute_Large += MAX_MISMATCHES % 2;
				Do_Branch=Do_All;
				Search_Brute(5,Brute_Large,STRINGLENGTH,RH,fwfmi);//finds mismatches of the form 112, stores possibles 113 and in the left half tag,
			}
		}

		Current_Tag=Guess_Complement;
	}

	if(PRINT_MISHITS )
	{
		if (Tag_Stat_Bad) fprintf(Mishit_File,"%s%sBAD%s%s", Description,Tag_Copy,Plus,Quality);
		else if(!Hits) fprintf(Mishit_File,"%s%s%s%s", Description,Tag_Copy,Plus,Quality);
	}
	//if(PRINT_MISHITS && !Hits) fprintf(Mishit_File,"%s%s%s%s", Description,Tag_Copy,Plus,Quality);
}

//{-----------------------------  BACKWARD SEARCH ROUTINE  -------------------------------------------------/

void Search_Backwards_Exact_X0(struct SARange & Tag,int Start,int StringLength,BWT *fmi)
{
////can optimise..	
	if (!Tag.Start) return;
	for(;;)
	{
		Get_SARange_Fast(Current_Tag[Start-Tag.Level],Tag,fmi);
		if (Tag.Start!=0)
		{
			if(Tag.Level== StringLength)
			{
				return;
			}
			else Tag.Level++;
		} 
		else
		{
			return;//No hit
		}
	}
}

void Search_10X(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	
	if (!Tag.Start) return;
	struct SARange Range,Temp_Range=Tag;
	int BMHStack_Top=0;
	BMHStack[0]=Tag;
	while(BMHStack_Top!=-1)//While Stack non-empty....
	{
		Range=BMHStack[BMHStack_Top];
		BMHStack_Top--;	//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_10X_OneSA(Range,Count,Start,StringLength,fmi);
			if(MAXHITS==Hits) return;
		}
		else
		{
			Branch_Detect_Backwards(Range,fmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start , Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)//only one mismatch allowed here...
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)//a tag of the form ?|1|0
							{
								//Temp_Range.Skip=0;
								SATot += Temp_Range.End-Temp_Range.Start+1;
								S++;
								Backwards(Temp_Range,1,LH);
								if (Temp_Range.Start)
								{
									Temp_Range.Level=1;
									Temp_BC[0]=Branch_Characters[0];Temp_BC[1]=Branch_Characters[1];Temp_BC[2]=Branch_Characters[2];Temp_BC[3]=Branch_Characters[3];
									memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
									Search_Forwards(Temp_Range,3,LH+1,RH,revfmi);
									if(MAXHITS==Hits) return;
									Branch_Characters[0]=Temp_BC[0];Branch_Characters[1]=Temp_BC[1];Branch_Characters[2]=Temp_BC[2];Branch_Characters[3]=Temp_BC[3];
									memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
								}

							}
							else continue;
						}
						else
						{
							BMHStack_Top++;//Push range
							Temp_Range.Level++;
							BMHStack[BMHStack_Top]=Temp_Range;
						}
					}
					else //store mismatches for later use...
					{
						if(Right_Mishits_Pointer < ARRAY_BOUND)
						{
							if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							Right_Mishits[Right_Mishits_Pointer]=Temp_Range;
							Right_Mishits_Pointer++;
						}
						if (HIGHSCAN && !Hits && Mismatches_Forward_Pointer+Mismatches_Backward_Pointer+Left_Mishits_Pointer+Right_Mishits_Pointer>EXTCUT) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
					}
					
				} 
			}
		}
	}
	return;
}

void Search_10X_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{

	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= FWDInverseSA0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if ( !Do_Branch[Start-Tag.Level] && Now != Current_Tag[Start-Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)//a tag of the form ?|1|0 , remove zero mismatch
				{
					SATot++;S++;
					//Tag.Skip=0;
					Backwards(Tag,1,LH);
					if(Tag.Start)
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						Tag.Level=1;
						Search_Forwards(Tag,3,LH+1,RH,revfmi);
					}
					return;

				}
				else return;
			}
			else { Tag.Level++;continue; }
		} 
		else //store mismatches for later use...
		{
			if(Right_Mishits_Pointer < ARRAY_BOUND)
			{
				Tag.End=Tag.Start;
				if (Tag.Level!=StringLength) Tag.Level++; 
				Right_Mishits[Right_Mishits_Pointer]=Tag;
				Right_Mishits_Pointer++;
			}
			if (HIGHSCAN && !Hits && Mismatches_Forward_Pointer+Mismatches_Backward_Pointer+Left_Mishits_Pointer+Right_Mishits_Pointer>EXTCUT) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
			return;
		}
	}
}

void Search_Backwards_X10(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	int BMHStack_Top=0;
	BMHStack[0]=Tag;
	struct SARange Range,Temp_Range;
	while(BMHStack_Top!=-1)//While Stack non-empty....
	{
		Range=BMHStack[BMHStack_Top];
		BMHStack_Top--;	//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Backwards_X10_OneSA(Range,Count,Start,StringLength,fmi);
			if(MAXHITS==Hits) return;
		}
		else
		{
			if (SEED && Range.End-Range.Start >20) continue;//@@@@
			Branch_Detect_Backwards(Range,fmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start , Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)//only one mismatch allowed here...
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)//a tag of the form ?|1|0
							{
								SATot += Temp_Range.End-Temp_Range.Start+1;
								S++;
								Temp_Range.Level=1;
								Temp_BC[0]=Branch_Characters[0];Temp_BC[1]=Branch_Characters[1];Temp_BC[2]=Branch_Characters[2];Temp_BC[3]=Branch_Characters[3];
								memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
								Search_Backwards(Temp_Range,2,LH,LH,fwfmi);
								if(MAXHITS==Hits) return;
								Branch_Characters[0]=Temp_BC[0];Branch_Characters[1]=Temp_BC[1];Branch_Characters[2]=Temp_BC[2];Branch_Characters[3]=Temp_BC[3];
								memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
							}
							else continue;
						}
						else
						{
							BMHStack_Top++;//Push range
							Temp_Range.Level++;
							BMHStack[BMHStack_Top]=Temp_Range;
						}
					}
					else
					{
						//if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
						//Temp_Range.Level++; 
						if (Possible_20_Pointer < END_BOUND)
						{
							Possible_20[Possible_20_Pointer]=Temp_Range;
							Possible_20_Pointer++;
						}
					}
				} 
			}
		}
	}
	return;
}

void Search_Backwards_X10_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{

	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}


	for(;;)
	{
		Index=Tag.Start;
		if (Index >= FWDInverseSA0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if ( !Do_Branch[Start-Tag.Level] && Now != Current_Tag[Start-Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)//only one mismatch allowed here...
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)//a tag of the form ?|1|0 , remove zero mismatch
				{
					SATot ++;
					S++;
					if(!Tag.Skip) Tag.End=Tag.Start;
					Tag.Level=1;
					Search_Backwards(Tag,2,LH,LH,fwfmi);
					return;
				}
				else return;
			}
			else { Tag.Level++;continue; }
		} 
		else
		{
			if(Possible_20_Pointer < END_BOUND)
			{
				if (!Tag.Skip) Tag.End=Tag.Start;
				//if (Tag.Level!=StringLength) Tag.Level++; 
				Possible_20[Possible_20_Pointer]=Tag;
				Possible_20_Pointer++;
			}
			return;
		}
	}
}
//}-----------------------------  BACKWARD SEARCH ROUTINE -------------------------------------------------/

//{-----------------------------  FORWARD SEARCH ROUTINES  -------------------------------------------------/
void Search_Forwards_0X(struct SARange & Tag,int Start,int StringLength,BWT *fmi)
{

	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offsets...
	for(;;)	
	{
		Get_SARange_Fast(Current_Tag[Start+Tag.Level],Tag,fmi);
		if (Tag.Start!=0)
		{
			if(Tag.Level== StringLength)
			{
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
			return;
		}
	}
}

void Search_X11(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	int BMStack_Top=0;
	BMStack_X11[0]=Tag;
	struct SARange Range,Temp_Range;
	while(BMStack_Top!=-1)//While Stack non-empty....
	{
		Range=BMStack_X11[BMStack_Top];
		BMStack_Top--;	//Pop the range
		Branch_Detect_Backwards(Range,fmi,Start);

		for(int Branch=0;Branch<4;Branch++)
		{
			if (Branch_Characters[Branch])//This character actually branches
			{
				Temp_Range=Range;//adjust
				Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
				Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

				if (Current_Tag[Start-Temp_Range.Level] != Branch)
				{
					Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
					Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
					Temp_Range.Mismatches++;
				}

				if (Temp_Range.Mismatches<=1)//we are guaranteed a valid SA range, check only for mismatches
				{
					if(Temp_Range.Level== StringLength)
					{
						if(Temp_Range.Mismatches)
						{
							Temp_Range.Level++;
							Temp_BC2[0]=Branch_Characters[0];Temp_BC2[1]=Branch_Characters[1];Temp_BC2[2]=Branch_Characters[2];Temp_BC2[3]=Branch_Characters[3];
							memcpy(Temp_Branch_Ranges2,Branch_Ranges,4*sizeof(SARange));//X|8
							Search_Half_Tag_X11(Temp_Range,2,STRINGLENGTH,RH,fwfmi);
							if (!FILTERUNIQUEHITS && !Larger_Than_Ten && Last_Mismatch_Written>5 && Last_Mismatch_Written <=10 ) return;
							if(MAXHITS==Hits) return;
							memcpy(Branch_Ranges,Temp_Branch_Ranges2,4*sizeof(SARange));
							Branch_Characters[0]=Temp_BC2[0];Branch_Characters[1]=Temp_BC2[1];Branch_Characters[2]=Temp_BC2[2];Branch_Characters[3]=Temp_BC2[3];

						}
						else continue;
					}
					else
					{
						BMStack_Top++;//Push range
						Temp_Range.Level++;
						BMStack_X11[BMStack_Top]=Temp_Range;
					}
				}
			} 
		}
	}
	return;
}

void Search_Half_Tag_X11_OneSA(struct SARange & Tag,char* Current_Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}


	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if (!Do_Branch[Start-Tag.Level] && Current_Tag[Start-Tag.Level]!=Now) return;  
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength && Tag.Mismatches==2)
			{
				SATot++;S++;
				if (S > SACUT && HIGHSCAN) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
				if (!Tag.Skip) Tag.End=Tag.Start;
				Tag.Level=1;
				Search_Backwards(Tag,5,LH,LH,fwfmi);//LH,fwfmi);
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else return;
	}
}

void Search_Half_Tag_X11(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;

	int BMStack_Top=0;
	BMStack_X11H[0]=Tag;
	struct SARange Range,Temp_Range;
	while(BMStack_Top!=-1)//While Stack non-empty....
	{
		Range=BMStack_X11H[BMStack_Top];
		BMStack_Top--;	//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Half_Tag_X11_OneSA(Range,Current_Tag,Count,Start,StringLength,fmi);
			if (!FILTERUNIQUEHITS && !Larger_Than_Ten && Last_Mismatch_Written>5 && Last_Mismatch_Written <=10 ) return;
			if(MAXHITS==Hits) return;
		}
		else
		{
			Branch_Detect_Backwards(Range,fmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=2)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength && Temp_Range.Mismatches == Count)
						{
							SATot += Temp_Range.End-Temp_Range.Start+1;
							S++;
							if (S > SACUT && HIGHSCAN) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
							Temp_Range.Level=1;
							Temp_BC1[0]=Branch_Characters[0];Temp_BC1[1]=Branch_Characters[1];Temp_BC1[2]=Branch_Characters[2];Temp_BC1[3]=Branch_Characters[3];
							memcpy(Temp_Branch_Ranges2,Branch_Ranges,4*sizeof(SARange));//X|16
							Search_Backwards(Temp_Range,5,LH,LH,fwfmi);//LH,fwfmi);
							if (!FILTERUNIQUEHITS && !Larger_Than_Ten && Last_Mismatch_Written>5 && Last_Mismatch_Written <=10 ) return;
							if(MAXHITS==Hits) return;
							memcpy(Branch_Ranges,Temp_Branch_Ranges2,4*sizeof(SARange));
							Branch_Characters[0]=Temp_BC1[0];Branch_Characters[1]=Temp_BC1[1];Branch_Characters[2]=Temp_BC1[2];Branch_Characters[3]=Temp_BC1[3];
						}
						else
						{
							BMStack_Top++;//Push range
							Temp_Range.Level++;
							BMStack_X11H[BMStack_Top]=Temp_Range;
						}
					}
				} 
			}
		}
	}
	return;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Search_Forwards_ExactX
 *  Description:  forward seach for exact occurence for string in Current_Tag in revfmi
 *  		  Start is the 1-indexed location of string position in Current_Tag
 *  		  String length is the length of substring to search...
 *  		  Start+Tag.Length-2 last position...
 *  		  return TRUE if full StringLength scanned.
 * =====================================================================================
 */
char Search_Forwards_ExactX(SARange & Tag, int Start,int StringLength)
{
	unsigned Index,First,Last;
	SARange Temp;//temporary tag to save last tag details before failure...

	char Now;
	Tag.Level=1;
	Start -= 2;//accomodate 1 based offset.. Start is 1 based and so is tag Level
	for(;;)	
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
			{
				Tag.Skip++;Tag.End=Tag.Start;
			}

			for(;;)
			{
				Index=Tag.Start;
				if (Index >= revfmi->inverseSa0) Index--;//adjust for missing $
				Now=revfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
				if (Current_Tag[Start+Tag.Level] == Now)
				{
					Tag.Start = revfmi->cumulativeFreq[Now] + BWTOccValue(revfmi, Tag.Start, Now) + 1;
					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0) 
					{
						Tag.Skip++;Tag.End=Tag.Start;
					}
					else {Tag.End=Tag.Start;}

					Cache_SF[Tag.Level]=Tag;
					if(Tag.Level== StringLength)
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						return TRUE; 	
					}
					else {Tag.Level++;continue;}
				} 
				else//mismatch...
				{
					if (!Tag.Skip) Tag.End=Tag.Start;
					return FALSE;	
				}
			}
		}
		else//SA range has sevaral possible hits... 
		{
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= revfmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
				for (unsigned Pos=First;Pos<=Last;Pos++)
				{
					Now=revfmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;	
				}

				Now=Current_Tag[Tag.Level+Start];
				if (Branch_Characters[Now])//we have a match... 
				{
					Tag.Start = revfmi->cumulativeFreq[Now] + BWTOccValue(revfmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges
				}
				else//mismatch..
				{
					return FALSE;
				}
			} 
			else
			{
				Temp=Tag;
				Get_SARange_Fast(Current_Tag[Start+Tag.Level],Tag,revfmi);
				if (!Tag.Start) {Tag=Temp;return FALSE;}
			}

			Cache_SF[Tag.Level]=Tag;
			if(Tag.Level== StringLength)
			{
				return TRUE;
			}
			else {Tag.Level++;continue;}

		}
	}
}

/*void Search_11X(const struct SARange & TagF,int Count,int Start,int StringLength,BWT *fmi)////////////////
{
	SARange Tag;
	Tag.Start=1;Tag.End=SOURCELENGTH;Tag.Skip=0;Tag.Mismatches=0;Tag.Mismatch_Char=0;Tag.Level=1;
	if(Search_Forwards_ExactX(Tag,1,LHQL/2))//search + for exact..
	{
		Tag.Level++;
		Start=Start-2;//Adjust for offset difference
		int FSHStack_Top=0;
		FSHStackX0X[0]=Tag;
		struct SARange Range,Temp_Range;
		while(FSHStack_Top!=-1)//While Stack non-empty....
		{
			Range=FSHStackX0X[FSHStack_Top];
			FSHStack_Top--;		//Pop the range

			Branch_Detect(Range,revfmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=1)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							Temp_Range.Level=1;
							Temp_BC2[0]=Branch_Characters[0];Temp_BC2[1]=Branch_Characters[1];Temp_BC2[2]=Branch_Characters[2];Temp_BC2[3]=Branch_Characters[3];
							memcpy(Temp_Branch_Ranges2,Branch_Ranges,4*sizeof(SARange));
							Search_Half_Tag_11X(Temp_Range,2,LHQL +1,LHQR,revfmi);
							if(MAXHITS==Hits) return;
							memcpy(Branch_Ranges,Temp_Branch_Ranges2,4*sizeof(SARange));
							Branch_Characters[0]=Temp_BC2[0];Branch_Characters[1]=Temp_BC2[1];Branch_Characters[2]=Temp_BC2[2];Branch_Characters[3]=Temp_BC2[3];
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							FSHStackX0X[FSHStack_Top]=Temp_Range;
						}
					}
				} 
			}

		}
	}
	//else//search 1|X
	{
		Tag.Start=1;Tag.End=SOURCELENGTH;Tag.Skip=0;Tag.Mismatches=0;Tag.Mismatch_Char=0;Tag.Level=1;
		Search_Backwards_Exact_X0( Tag,LHQL,LHQL/2,fwfmi);// ?|?|0
		Tag.Level=1;
		Start=LHQL-(LHQL/2);
		StringLength=LHQL-(LHQL/2);

		if (!Tag.Start) return;
		int BMHStack_Top=0;
		BMHStack[0]=Tag;
		struct SARange Range,Temp_Range;
		while(BMHStack_Top!=-1)//While Stack non-empty....
		{
			Range=BMHStack[BMHStack_Top];
			BMHStack_Top--;	//Pop the range

			//if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
			//{
			//	Search_Backwards_X10_OneSA(Range,Count,Start,StringLength,fmi);
			//	if(MAXHITS==Hits) return;
			//}
			//else
			{
				Branch_Detect_Backwards(Range,fwfmi,Start);
				for(int Branch=0;Branch<4;Branch++)
				{
					if (Branch_Characters[Branch])
					{
						Temp_Range=Range;//adjust
						Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start , Branch) + 1;
						Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

						if (Current_Tag[Start-Temp_Range.Level] != Branch)//only one mismatch allowed here...
						{
							Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
							Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
							Temp_Range.Mismatches++;
						}

						if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
						{
							if(Temp_Range.Level== StringLength)
							{
								if(Temp_Range.Mismatches)//a tag of the form ?|1|0
								{
									Tag.Start=1;Tag.End=SOURCELENGTH;Tag.Skip=0;Tag.Mismatches=0;Tag.Mismatch_Char=0;Tag.Level=1;
									Backwards(Temp_Range,1,LHQL);

									Temp_Range.Level=1;
									Temp_BC2[0]=Branch_Characters[0];Temp_BC2[1]=Branch_Characters[1];Temp_BC2[2]=Branch_Characters[2];Temp_BC2[3]=Branch_Characters[3];
									memcpy(Temp_Branch_Ranges2,Branch_Ranges,4*sizeof(SARange));
									Search_Half_Tag_11X(Temp_Range,2,LHQL +1,LHQR,revfmi);
									if(MAXHITS==Hits) return;
									memcpy(Branch_Ranges,Temp_Branch_Ranges2,4*sizeof(SARange));
									Branch_Characters[0]=Temp_BC2[0];Branch_Characters[1]=Temp_BC2[1];Branch_Characters[2]=Temp_BC2[2];Branch_Characters[3]=Temp_BC2[3];

								}
								else continue;
							}
							else
							{
								BMHStack_Top++;//Push range
								Temp_Range.Level++;
								BMHStack[BMHStack_Top]=Temp_Range;
							}
						}
						else
						{
							//if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							//Temp_Range.Level++; 
							if (Possible_20_Pointer < END_BOUND)
							{
								Possible_20[Possible_20_Pointer]=Temp_Range;
								Possible_20_Pointer++;
							}
						}
					} 
				}
			}
		}
	}
	return;
}*/

void Search_11X(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)////////////////
{

	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSHStack_Top=0;
	FSHStackX0X[0]=Tag;
	struct SARange Range,Temp_Range;
	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=FSHStackX0X[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		Branch_Detect(Range,revfmi,Start);
		for(int Branch=0;Branch<4;Branch++)
		{
			if (Branch_Characters[Branch])//This character actually branches
			{
				Temp_Range=Range;
				Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
				Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

				if (Current_Tag[Temp_Range.Level+Start]!=Branch)
				{
					Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
					Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
					Temp_Range.Mismatches++;
				}

				if (Temp_Range.Mismatches<=1)//we are guaranteed a valid SA range, check only for mismatches
				{
					if(Temp_Range.Level== StringLength)
					{
						if (Temp_Range.Mismatches == 1)
						{
							Temp_Range.Level=1;
							Temp_BC2[0]=Branch_Characters[0];Temp_BC2[1]=Branch_Characters[1];Temp_BC2[2]=Branch_Characters[2];Temp_BC2[3]=Branch_Characters[3];
							memcpy(Temp_Branch_Ranges2,Branch_Ranges,4*sizeof(SARange));
							Search_Half_Tag_11X(Temp_Range,2,LHQL +1,LHQR,revfmi);
							if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
							if(MAXHITS==Hits) return;
							memcpy(Branch_Ranges,Temp_Branch_Ranges2,4*sizeof(SARange));
							Branch_Characters[0]=Temp_BC2[0];Branch_Characters[1]=Temp_BC2[1];Branch_Characters[2]=Temp_BC2[2];Branch_Characters[3]=Temp_BC2[3];
						}
					}
					else
					{
						FSHStack_Top++;//Push range
						Temp_Range.Level++;
						FSHStackX0X[FSHStack_Top]=Temp_Range;
					}
				}
			} 
		}

	}
	return;
}

void Search_Half_Tag_11X_OneSA(struct SARange & Tag,char* Current_Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if ( !Do_Branch[Start+Tag.Level] && Now != Current_Tag[Start+Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start+Tag.Level] != Now) 
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
		}
		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)// && Tag.Mismatches==Count)
			{
				if (Tag.Mismatches != Count) return;
				SATot++;S++;
				if (S > SACUT && HIGHSCAN) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
				if(!Tag.Skip) Tag.End=Tag.Start;
				Tag.Level=1;
				Search_Forwards(Tag,4,LH+1,RH,revfmi);
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
			/*if (Left_Mishits_Pointer < ARRAY_BOUND)//Bounds Check...
			{
				if(!Tag.Skip) Tag.End=Tag.Start;
				//Tag.Level=LHQR+Tag.Level;
				//if (Tag.Level!=StringLength) Tag.Level++; 
				//Tag.Level=LHQR+Tag.Level;//+1;
				Tag.Level += Start+2;
				Left_Mishits[Left_Mishits_Pointer]=Tag;
				Left_Mishits_Pointer++;
			}*/
			return;
		}
	}
}

void Search_Half_Tag_11X(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSHStack_Top=0;
	FSHStack[0]=Tag;
	struct SARange Range,Temp_Range;
	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=FSHStack[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Half_Tag_11X_OneSA(Range,Current_Tag,Count,Start,StringLength,revfmi);
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
			if(MAXHITS==Hits) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;
					Temp_Range.End = Branch_Ranges[Branch].End;

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength )//&& Temp_Range.Mismatches==Count)
						{
							if (Temp_Range.Mismatches != Count) continue;
							SATot += Temp_Range.End-Temp_Range.Start+1;
							S++;
							if (S > SACUT && HIGHSCAN) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
							Temp_Range.Level=1;
							Temp_BC1[0]=Branch_Characters[0];Temp_BC1[1]=Branch_Characters[1];Temp_BC1[2]=Branch_Characters[2];Temp_BC1[3]=Branch_Characters[3];
							memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
							Search_Forwards(Temp_Range,4,LH+1,RH,revfmi);
							if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
							if(MAXHITS==Hits) return;
							memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
							Branch_Characters[0]=Temp_BC1[0];Branch_Characters[1]=Temp_BC1[1];Branch_Characters[2]=Temp_BC1[2];Branch_Characters[3]=Temp_BC1[3];
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							FSHStack[FSHStack_Top]=Temp_Range;
						}
					}
					else //store mismatches for later use...
					{
						/*if(Left_Mishits_Pointer < ARRAY_BOUND) 
						{
							//Temp_Range.Level=LHQR+Temp_Range.Level;
							//if (Temp_Range.Level!=StringLength) Temp_Range.Level++;
							//Temp_Range.Level=LHQR+Temp_Range.Level;//+1;
							Temp_Range.Level += Start+2;
							Left_Mishits[Left_Mishits_Pointer]=Temp_Range;
							Left_Mishits_Pointer++;
						}*/
					}
				} 
			}

		}
	}
	return;
}
//---new...

void Search_XL01(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSHStack_Top=0;
	FSHStack[0]=Tag;
	struct SARange Range,Temp_Range;
	SARange Temp_Branch_Ranges[4];
	unsigned Temp_BC1[4];

	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=FSHStack[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_XL01_OneSA(Range,Count,Start,StringLength,revfmi);
			if (!Larger_Than_Ten && !FILTERUNIQUEHITS && Last_Mismatch_Written>5 && Last_Mismatch_Written <=10 ) return;
			if(MAXHITS==Hits) return;
		}
		else
		{
			if (SEED && Range.End-Range.Start >20) continue;//@@@@
			Branch_Detect(Range,revfmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)
							{
								int T=RH;RH=RHQR;
								Reverse(Temp_Range,STRINGLENGTH,RH);
								RH=T;
								if(Temp_Range.Start)
								{
									Temp_BC1[0]=Branch_Characters[0];Temp_BC1[1]=Branch_Characters[1];Temp_BC1[2]=Branch_Characters[2];Temp_BC1[3]=Branch_Characters[3];
									Temp_Range.Level=RHQR+1;//Temp_Range.Level=1;
									memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
									Search_Half_Tag_X11(Temp_Range,2,STRINGLENGTH,RH,fwfmi);
									if (!Larger_Than_Ten && !FILTERUNIQUEHITS && Last_Mismatch_Written>5 && Last_Mismatch_Written <=10 ) return;
									memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
									if(MAXHITS==Hits) return;
									Branch_Characters[0]=Temp_BC1[0];Branch_Characters[1]=Temp_BC1[1];Branch_Characters[2]=Temp_BC1[2];Branch_Characters[3]=Temp_BC1[3];

								}
							}
							else continue;
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							FSHStack[FSHStack_Top]=Temp_Range;
						}
					}
					else
					{
						if(Possible_02_Pointer < END_BOUND)
						{
							if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							Possible_02[Possible_02_Pointer]=Temp_Range;
							Possible_02_Pointer++;
						}
					}
				} 
			}

		}
	}
	return;
}

void Search_XL01_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if ( !Do_Branch[Start+Tag.Level] && Now != Current_Tag[Start+Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start+Tag.Level] != Now) 
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
		}
		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)
				{
					//Tag.End=Tag.Start;
					//Tag.Skip=0;
					int T=RH;RH=RHQR;
					Reverse(Tag,STRINGLENGTH,RH);
					RH=T;
					if(Tag.Start)
					{
						if(!Tag.Skip) Tag.End=Tag.Start;
						Tag.Level=RHQR+1;//Temp_Range.Level=1;
						Search_Half_Tag_X11(Tag,2,STRINGLENGTH,RH,fwfmi);
					}
					return;
				}
				else return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
			if(Possible_02_Pointer < END_BOUND)
			{
				if(!Tag.Skip)Tag.End=Tag.Start;
				if (Tag.Level!=StringLength) Tag.Level++; 
				Possible_02[Possible_02_Pointer]=Tag;
				Possible_02_Pointer++;
			}
			return;
		}
	}
}
void Search_Backwards_XL10_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{

	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}


	for(;;)
	{
		Index=Tag.Start;
		if (Index >= FWDInverseSA0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if ( !Do_Branch[Start-Tag.Level] && Now != Current_Tag[Start-Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)//only one mismatch allowed here...
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)//a tag of the form ?|1|0 , remove zero mismatch
				{
					if(!Tag.Skip) Tag.End=Tag.Start;
					Tag.Level=RHQR+1;//Temp_Range.Level=1;
					Search_Half_Tag_X11(Tag,2,STRINGLENGTH,RH,fwfmi);
					return;
				}
				else return;
			}
			else { Tag.Level++;continue; }
		} 
		else
		{
			if(Possible_20_Pointer < END_BOUND)
			{
				if (!Tag.Skip) Tag.End=Tag.Start;
				//if (Tag.Level!=StringLength) Tag.Level++; 
				Possible_20[Possible_20_Pointer]=Tag;
				Possible_20_Pointer++;
			}
			return;
		}
	}
}

void Search_Backwards_XL10(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	int BMHStack_Top=0;
	BMHStack[0]=Tag;
	struct SARange Range,Temp_Range;
	SARange  Temp_Branch_Ranges2[4];
	unsigned Temp_BC[4];

	while(BMHStack_Top!=-1)//While Stack non-empty....
	{
		Range=BMHStack[BMHStack_Top];
		BMHStack_Top--;	//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Backwards_XL10_OneSA(Range,Count,Start,StringLength,fmi);
			if (!Larger_Than_Ten && !FILTERUNIQUEHITS && Last_Mismatch_Written>5 && Last_Mismatch_Written <=10 ) return;
			if(MAXHITS==Hits) return;
		}
		else
		{
			if (SEED && Range.End-Range.Start >20) continue;//@@@@
			Branch_Detect_Backwards(Range,fmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start , Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)//only one mismatch allowed here...
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)//a tag of the form ?|1|0
							{
								Temp_Range.Level=RHQR+1;//Temp_Range.Level=1;
								Temp_BC[0]=Branch_Characters[0];Temp_BC[1]=Branch_Characters[1];Temp_BC[2]=Branch_Characters[2];Temp_BC[3]=Branch_Characters[3];
								memcpy(Temp_Branch_Ranges2,Branch_Ranges,4*sizeof(SARange));
								Search_Half_Tag_X11(Temp_Range,2,STRINGLENGTH,RH,fwfmi);
								if (!Larger_Than_Ten && !FILTERUNIQUEHITS && Last_Mismatch_Written>5 && Last_Mismatch_Written <=10 ) return;
								if(MAXHITS==Hits) return;
								Branch_Characters[0]=Temp_BC[0];Branch_Characters[1]=Temp_BC[1];Branch_Characters[2]=Temp_BC[2];Branch_Characters[3]=Temp_BC[3];
								memcpy(Branch_Ranges,Temp_Branch_Ranges2,4*sizeof(SARange));
							}
							else continue;

						}
						else
						{
							BMHStack_Top++;//Push range
							Temp_Range.Level++;
							BMHStack[BMHStack_Top]=Temp_Range;
						}
					}
					else
					{
						//if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
						//Temp_Range.Level++; 
						if (Possible_20_Pointer < END_BOUND)
						{
							Possible_20[Possible_20_Pointer]=Temp_Range;
							Possible_20_Pointer++;
						}
					}
				} 
			}
		}
	}
	return;
}

void Search_10LX(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	
	if (!Tag.Start) return;
	struct SARange Range,Temp_Range=Tag;
	int BMHStack_Top=0;
	BMHStack[0]=Tag;
	SARange Temp_Branch_Ranges[4];
	unsigned Temp_BC[4];
	while(BMHStack_Top!=-1)//While Stack non-empty....
	{
		Range=BMHStack[BMHStack_Top];
		BMHStack_Top--;	//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_10LX_OneSA(Range,Count,Start,StringLength,fmi);
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
			if(MAXHITS==Hits) return;
		}
		else
		{
			Branch_Detect_Backwards(Range,fmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start , Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)//only one mismatch allowed here...
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)//a tag of the form ?|1|0
							{
								//Temp_Range.Skip=0;
								Backwards(Temp_Range,1,LHQL);
								if (Temp_Range.Start)
								{
									Temp_Range.Level=1;
									Temp_BC[0]=Branch_Characters[0];Temp_BC[1]=Branch_Characters[1];Temp_BC[2]=Branch_Characters[2];Temp_BC[3]=Branch_Characters[3];
									memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
									Search_Half_Tag_11X(Temp_Range,2,LHQL +1,LHQR,revfmi);
									if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
									if(MAXHITS==Hits) return;
									Branch_Characters[0]=Temp_BC[0];Branch_Characters[1]=Temp_BC[1];Branch_Characters[2]=Temp_BC[2];Branch_Characters[3]=Temp_BC[3];
									memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
								}

							}
							else continue;
						}
						else
						{
							BMHStack_Top++;//Push range
							Temp_Range.Level++;
							BMHStack[BMHStack_Top]=Temp_Range;
						}
					}
					else //store mismatches for later use...
					{
						/*if(Right_Mishits_Pointer < ARRAY_BOUND)
						{
							if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							Right_Mishits[Right_Mishits_Pointer]=Temp_Range;
							Right_Mishits_Pointer++;
						}
						if (HIGHSCAN && Mismatches_Forward_Pointer+Mismatches_Backward_Pointer+Left_Mishits_Pointer+Right_Mishits_Pointer>EXTCUT) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}*/
					}
					
				} 
			}
		}
	}
	return;
}

void Search_10LX_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{

	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= FWDInverseSA0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if ( !Do_Branch[Start-Tag.Level] && Now != Current_Tag[Start-Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)//a tag of the form ?|1|0 , remove zero mismatch
				{
					//Tag.Skip=0;
					Backwards(Tag,1,LHQL);
					if(Tag.Start)
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						Tag.Level=1;
						Search_Half_Tag_11X(Tag,2,LHQL +1,LHQR,revfmi);
					}
					return;

				}
				else return;
			}
			else { Tag.Level++;continue; }
		} 
		else //store mismatches for later use...
		{
			/*if(Right_Mishits_Pointer < ARRAY_BOUND)
			{
				Tag.End=Tag.Start;
				if (Tag.Level!=StringLength) Tag.Level++; 
				Right_Mishits[Right_Mishits_Pointer]=Tag;
				Right_Mishits_Pointer++;
			}*/
			return;
		}
	}
}
void Search_01LX_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if ( !Do_Branch[Start+Tag.Level] && Now != Current_Tag[Start+Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start+Tag.Level] != Now) 
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
		}
		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)
				{
					if (!Tag.Skip) Tag.End=Tag.Start;
					Tag.Level=1;
					Search_Half_Tag_11X(Tag,2,LHQL +1,LHQR,revfmi);
					return;
				}
				else return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
//new...
			if(Left_Mishits_Pointer < ARRAY_BOUND)
			{
				/*if (Tag.Level!=StringLength) Tag.Level++; 
				Left_Mishits[Left_Mishits_Pointer]=Tag;
				Left_Mishits_Pointer++;*/

				//if(!Tag.Skip) Tag.End=Tag.Start;
				//Tag.Level=LHQR+Tag.Level;
				//if (Tag.Level!=StringLength) Tag.Level++; 
				//Tag.Level=LHQR+Tag.Level;//+1;
				//Tag.Level += Start+2;
				//Left_Mishits[Left_Mishits_Pointer]=Tag;
				//Left_Mishits_Pointer++;
			}
//new...
			return;
		}
	}
}
void Search_01LX(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSHStack_Top=0;
	FSHStackX0X[0]=Tag;
	struct SARange Range,Temp_Range;
	SARange Temp_Branch_Ranges2[4];
	unsigned Temp_BC2[4]; 

	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=FSHStackX0X[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_01LX_OneSA(Range,Count,Start,StringLength,revfmi);
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
			if(MAXHITS==Hits) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;
					Temp_Range.End = Branch_Ranges[Branch].End;

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)
							{
								Temp_Range.Level=1;
								Temp_BC2[0]=Branch_Characters[0];Temp_BC2[1]=Branch_Characters[1];Temp_BC2[2]=Branch_Characters[2];Temp_BC2[3]=Branch_Characters[3];
								memcpy(Temp_Branch_Ranges2,Branch_Ranges,4*sizeof(SARange));
								Search_Half_Tag_11X(Temp_Range,2,LHQL +1,LHQR,revfmi);
								if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
								if(MAXHITS==Hits) return;
								memcpy(Branch_Ranges,Temp_Branch_Ranges2,4*sizeof(SARange));
								Branch_Characters[0]=Temp_BC2[0];Branch_Characters[1]=Temp_BC2[1];Branch_Characters[2]=Temp_BC2[2];Branch_Characters[3]=Temp_BC2[3];
							}
							else continue;
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							FSHStackX0X[FSHStack_Top]=Temp_Range;
						}
					}
////new...
					else //store mismatches for later use...
					{
						if(Left_Mishits_Pointer < ARRAY_BOUND)
						{
							//Temp_Range.Level += Start+2;
							//Left_Mishits[Left_Mishits_Pointer]=Temp_Range;
							//Left_Mishits_Pointer++;
							/*if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							Left_Mishits[Left_Mishits_Pointer]=Temp_Range;
							Left_Mishits_Pointer++;*/
						}
					}
////new...

				} 
			}

		}
	}
	return;
}


//--new ...
void Search_01X(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSHStack_Top=0;
	FSHStack[0]=Tag;
	struct SARange Range,Temp_Range;
	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=FSHStack[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_01X_OneSA(Range,Count,Start,StringLength,revfmi);
			if(MAXHITS==Hits) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;
					Temp_Range.End = Branch_Ranges[Branch].End;

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)
							{
								SATot += Temp_Range.End-Temp_Range.Start+1;
								S++;
								Temp_Range.Level=1;
								Temp_BC1[0]=Branch_Characters[0];Temp_BC1[1]=Branch_Characters[1];Temp_BC1[2]=Branch_Characters[2];Temp_BC1[3]=Branch_Characters[3];
								memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
								Search_Forwards(Temp_Range,3,LH+1,RH,revfmi);
								memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
								if(MAXHITS==Hits) return;
								Branch_Characters[0]=Temp_BC1[0];Branch_Characters[1]=Temp_BC1[1];Branch_Characters[2]=Temp_BC1[2];Branch_Characters[3]=Temp_BC1[3];
							}
							else continue;
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							FSHStack[FSHStack_Top]=Temp_Range;
						}
					}
////new...
					else //store mismatches for later use...
					{
						if(Left_Mishits_Pointer < ARRAY_BOUND)
						{
							Temp_Range.Level += Start+2;
							Left_Mishits[Left_Mishits_Pointer]=Temp_Range;
							Left_Mishits_Pointer++;
							/*if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							Left_Mishits[Left_Mishits_Pointer]=Temp_Range;
							Left_Mishits_Pointer++;*/
						}
					}
////new...

				} 
			}

		}
	}
	return;
}

void Search_01X_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if ( !Do_Branch[Start+Tag.Level] && Now != Current_Tag[Start+Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start+Tag.Level] != Now) 
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
		}
		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)
				{
					SATot++;S++;
					if (!Tag.Skip) Tag.End=Tag.Start;
					Tag.Level=1;
					Search_Forwards(Tag,3,LH+1,RH,revfmi);
					return;
				}
				else return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
//new...
			if(Left_Mishits_Pointer < ARRAY_BOUND)
			{
				/*if (Tag.Level!=StringLength) Tag.Level++; 
				Left_Mishits[Left_Mishits_Pointer]=Tag;
				Left_Mishits_Pointer++;*/

				if(!Tag.Skip) Tag.End=Tag.Start;
				//Tag.Level=LHQR+Tag.Level;
				//if (Tag.Level!=StringLength) Tag.Level++; 
				//Tag.Level=LHQR+Tag.Level;//+1;
				Tag.Level += Start+2;
				Left_Mishits[Left_Mishits_Pointer]=Tag;
				Left_Mishits_Pointer++;
			}
//new...
			return;
		}
	}
}

void Search_X01(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSHStack_Top=0;
	FSHStack[0]=Tag;
	struct SARange Range,Temp_Range;
	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=FSHStack[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_X01_OneSA(Range,Count,Start,StringLength,revfmi);
			if(MAXHITS==Hits) return;
		}
		else
		{
			if (SEED && Range.End-Range.Start >20) continue;//@@@@
			Branch_Detect(Range,revfmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)
							{
								Reverse(Temp_Range,STRINGLENGTH,RH);
								if(Temp_Range.Start)
								{
									SATot += Temp_Range.End-Temp_Range.Start+1;
									S++;
									Temp_Range.Level=1;
									Temp_BC1[0]=Branch_Characters[0];Temp_BC1[1]=Branch_Characters[1];Temp_BC1[2]=Branch_Characters[2];Temp_BC1[3]=Branch_Characters[3];
									Search_Backwards(Temp_Range,2,LH,LH,fwfmi);
									if(MAXHITS==Hits) return;
									Branch_Characters[0]=Temp_BC1[0];Branch_Characters[1]=Temp_BC1[1];Branch_Characters[2]=Temp_BC1[2];Branch_Characters[3]=Temp_BC1[3];
								}
							}
							else continue;
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							FSHStack[FSHStack_Top]=Temp_Range;
						}
					}
					else
					{
						if(Possible_02_Pointer < END_BOUND)
						{
							if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							Possible_02[Possible_02_Pointer]=Temp_Range;
							Possible_02_Pointer++;
						}
					}
				} 
			}

		}
	}
	return;
}

void Search_X01_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if ( !Do_Branch[Start+Tag.Level] && Now != Current_Tag[Start+Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start+Tag.Level] != Now) 
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
		}
		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)
				{
					SATot++;
					S++;
					//Tag.End=Tag.Start;
					//Tag.Skip=0;
					Reverse(Tag,STRINGLENGTH,RH);
					if(Tag.Start)
					{
						Tag.Level=1;
						Search_Backwards(Tag,2,LH,LH,fwfmi);
					}
					return;
				}
				else return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
			if(Possible_02_Pointer < END_BOUND)
			{
				if(!Tag.Skip)Tag.End=Tag.Start;
				if (Tag.Level!=StringLength) Tag.Level++; 
				Possible_02[Possible_02_Pointer]=Tag;
				Possible_02_Pointer++;
			}
			return;
		}
	}
}

void Backwards_Small(struct SARange & Tag,int Start,int StringLength)
{	

	unsigned Temp=0;
	char Mismatch_Count=Tag.Mismatches;
	int Temp_Pos=0;//New_Char;
	unsigned pos;
	Start=Start-2;
	
	for( int i=0;i<Mismatch_Count;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Temp=Temp | (Current_Tag[pos]<<i*2);
		Current_Tag[pos]=Tag.Mismatch_Char>>(2*i) & 3;
	}

	Temp_BCB[0]=Branch_Characters[0];Temp_BCB[1]=Branch_Characters[1];Temp_BCB[2]=Branch_Characters[2];Temp_BCB[3]=Branch_Characters[3];
	memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
	{
		/*if(LOOKUPSIZE==3)
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		Tag.Start=Forward_Start_LookupX[c];Tag.End=Forward_End_LookupX[c];
		Tag.Level=LOOKUPSIZE + 1; Tag.Skip=0;*/
		Tag.Start=1;Tag.End=SOURCELENGTH;Tag.Skip=0;Tag.Mismatches=0;Tag.Mismatch_Char=0;Tag.Level=1;
		Search_Exact(Tag,Start,StringLength,revfmi);
	}
	Branch_Characters[0]=Temp_BCB[0];Branch_Characters[1]=Temp_BCB[1];Branch_Characters[2]=Temp_BCB[2];Branch_Characters[3]=Temp_BCB[3];
	memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
	for( int i=0;i<Tag.Mismatches;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Current_Tag[pos]=(Temp>>(2*i)) & 3;
	}
	return;
}
void Backwards(struct SARange & Tag,int Start,int StringLength)
{	

	unsigned Temp=0;
	char Mismatch_Count=Tag.Mismatches;
	int Temp_Pos=0;//New_Char;
	unsigned pos;
	Start=Start-2;
	SARange Temp_Branch_Ranges[4];
	unsigned Temp_BCB[4];
	
	for( int i=0;i<Mismatch_Count;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Temp=Temp | (Current_Tag[pos]<<i*2);
		Current_Tag[pos]=Tag.Mismatch_Char>>(2*i) & 3;
	}

	Temp_BCB[0]=Branch_Characters[0];Temp_BCB[1]=Branch_Characters[1];Temp_BCB[2]=Branch_Characters[2];Temp_BCB[3]=Branch_Characters[3];
	memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
	{
		if(LOOKUPSIZE==3)
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		Tag.Start=Forward_Start_LookupX[c];Tag.End=Forward_End_LookupX[c];
		Tag.Level=LOOKUPSIZE + 1; Tag.Skip=0;
		Search_Exact(Tag,Start,StringLength,revfmi);
	}
	Branch_Characters[0]=Temp_BCB[0];Branch_Characters[1]=Temp_BCB[1];Branch_Characters[2]=Temp_BCB[2];Branch_Characters[3]=Temp_BCB[3];
	memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
	for( int i=0;i<Tag.Mismatches;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Current_Tag[pos]=(Temp>>(2*i)) & 3;
	}
	return;
}



void Reverse(struct SARange & Tag,int Start,int StringLength)
{	
	unsigned Temp=0;
	char New_Char;
	char Mismatch_Count=Tag.Mismatches;
	unsigned pos;
	unsigned Temp_BCR[4];

	for( int i=0;i<Mismatch_Count;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Temp=Temp | (Current_Tag[pos]<<i*2);
		Current_Tag[pos]=Tag.Mismatch_Char>>(2*i) & 3;
	}
	Temp_BCR[0]=Branch_Characters[0];Temp_BCR[1]=Branch_Characters[1];Temp_BCR[2]=Branch_Characters[2];Temp_BCR[3]=Branch_Characters[3];
	{
		if(LOOKUPSIZE==3)
		{
			c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		Tag.Start=Backward_Start_LookupX[c];Tag.End=Backward_End_LookupX[c];
		Tag.Level=LOOKUPSIZE + 1;Tag.Skip=0;
		Search_Backwards_Exact( Tag,STRINGLENGTH,RH,fwfmi);//Backward scan for ?|0
	}
	Branch_Characters[0]=Temp_BCR[0];Branch_Characters[1]=Temp_BCR[1];Branch_Characters[2]=Temp_BCR[2];Branch_Characters[3]=Temp_BCR[3];
	for( int i=0;i<Tag.Mismatches;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Current_Tag[pos]=(Temp>>(2*i)) & 3;
	}
	return;
}

void Search_Exact(struct SARange & Tag,int Start,int StringLength,BWT *fmi)
{
	int Level;
	unsigned long Index,Now,First,Last;
	if (!Tag.Start) return;

	for(;;)	
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
			{
				Tag.Skip++;Tag.End=Tag.Start;
			}
			Level=Tag.Level;
			for(;;)
			{
				Index=Tag.Start;
				if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
				Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
				if (Current_Tag[Start+Level] == Now)
				{
					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;
					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0) 
					{
						Tag.Skip++;Tag.End=Tag.Start;
					}
					//Tag.End=Tag.Start;
					if(Level== StringLength)
					{
						if(!Tag.Skip) Tag.End=Tag.Start; 
						return;	
					}
					else {Level++;continue;}
				} 
				else//mismatch...
				{
					Tag.Start=0;//Tag.End=0;
					return;	
				}
			}
		}
		else//SA range has sevaral possible hits... 
		{
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
				for (unsigned long Pos=First;Pos<=Last;Pos++)
				{
					Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;	
				}

				Now=Current_Tag[Tag.Level+Start];
				if (Branch_Characters[Now])//we have a match... 
				{
					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges
				}
				else//mismatch..
				{
					Tag.Start=0;
				}
			} 
			else
			{
				Get_SARange_Fast(Current_Tag[Start+Tag.Level],Tag,fmi);
			}

			if (Tag.Start!=0)
			{
				if(Tag.Level== StringLength)
				{
					return;
				}
				else {Tag.Level++;continue;}
			} 
			else//Mismatch
			{
				return;
			}

		}
	}
}

void Search_Forwards_Exact(struct SARange & Tag,int Start,int StringLength,BWT *fmi)
{
	int Level;
	Exact_Match_Forward[Start+LH].Start=0;
	unsigned long Index,Now,First,Last;
	if (!Tag.Start) return;
	for(;;)	
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
			{
				Tag.Skip++;Tag.End=Tag.Start;
			}

			for(;;)
			{
				
				Index=Tag.Start;
				if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
				Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
				if (Current_Tag[Start+Tag.Level] == Now)
				{
					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;
					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0) 
					{
						Tag.Skip++;Tag.End=Tag.Start;
					}
					//Tag.End=Tag.Start;
					Exact_Match_Forward[Start+Tag.Level]=Tag;
					if (!Tag.Skip) Exact_Match_Forward[Start+Tag.Level].End=Tag.Start;

					if(Tag.Level== StringLength)
					{
						//if (Tag.Skip) Tag.Start= Tag.End; else Tag.End=Tag.Start;
						if (!Tag.Skip) Tag.End=Tag.Start;
						Print_LocationX(Tag);
						return;	
					}
					else {Tag.Level++;continue;}
				} 
				else//mismatch...
				{
					Exact_Match_Forward[0]=Tag;//save old location for heuristics...
					Tag.Start=0;//Tag.End=0;
					Exact_Match_Forward[Start+Tag.Level]=Tag;
					return;	
				}
			}
		}
		else//SA range has sevaral possible hits... 
		{
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
				for (unsigned long Pos=First;Pos<=Last;Pos++)
				{
					Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;	
				}

				Now=Current_Tag[Tag.Level+Start];
				if (Branch_Characters[Now])//we have a match... 
				{
					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges
				}
				else//mismatch..
				{
					Exact_Match_Forward[0]=Tag;//save old location for heuristics...
					Tag.Start=0;
				}
			} 
			else
			{
				Exact_Match_Forward[0]=Tag;//save old location for heuristics...
				Get_SARange_Fast(Current_Tag[Start+Tag.Level],Tag,fmi);
			}

			Exact_Match_Forward[Tag.Level+Start]=Tag;
			if (Tag.Start!=0)
			{
				if(Tag.Level== StringLength)
				{
					Print_LocationX(Tag);
					return;
				}
				else {Tag.Level++;continue;}
			} 
			else//Mismatch
			{
				Exact_Match_Forward[Start+Tag.Level]=Tag;
				return;
			}

		}
	}
}

void Search_Forwards(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSStack_Top=0;
	FSSStack[0]=Tag;
	struct SARange Range,Temp_Range;
	while(FSStack_Top!=-1)//While Stack non-empty....
	{
		Range=FSSStack[FSStack_Top];
		FSStack_Top--;		//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Forwards_OneSA(Range,Count,Start,StringLength,revfmi);
			if(MAXHITS==Hits) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start);//One_Branch(Range,revfmi);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{

						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;

					}


					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if (Temp_Range.Mismatches == Count) //dont print exact matches
							{
								Print_LocationX(Temp_Range);
								if (MAXHITS==Hits) return;
							}
							else continue;
						}
						else 
						{

							FSStack_Top++;//Push range
							Temp_Range.Level++;
							FSSStack[FSStack_Top]=Temp_Range;
						}
					}
					else
					{
						if(5 >Count)//store only for one mismatch... and last node will not branch
						{
							if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							else //2 mismatches with the last at the end...
							{
								if(Two_Mismatches_At_End_Forward_Pointer < END_BOUND)
								{
									Two_Mismatches_At_End_Forward[Two_Mismatches_At_End_Forward_Pointer]=Temp_Range;
									Two_Mismatches_At_End_Forward_Pointer++;
								}
								continue;
							}
							if(Mismatches_Forward_Pointer < ARRAY_BOUND)
							{
								Mismatches_Forward[Mismatches_Forward_Pointer]=Temp_Range;
								Mismatches_Forward_Pointer++;
							}
							if (HIGHSCAN && !Hits && Mismatches_Forward_Pointer+Mismatches_Backward_Pointer+Left_Mishits_Pointer+Right_Mishits_Pointer>EXTCUT) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
						}
						continue;
					}
				} 
			}
		}
	}
	return;
}

void Search_Forwards_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
		if (!Do_Branch[Tag.Level+Start] && Current_Tag[Tag.Level+Start]!=Now) return;  
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Tag.Level+Start]!=Now)
		{

			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
			
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches==Count)//avoid printing exact matches...
				{
					//if (Tag.Skip) Tag.Start= Tag.End; else Tag.End=Tag.Start;
					if (!Tag.Skip) Tag.End=Tag.Start;
					Print_LocationX(Tag);
				}
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else//log 2 mismatches 
		{
			if(5 > Count)//store only for one mismatch... later report these on a seperate stack..
			{
				if (!Tag.Skip) Tag.End=Tag.Start;//possibly two mismatch exists..
				if (Tag.Level != StringLength) Tag.Level++; 
				else //2 mismatches occuring in last position...
				{
					if(Two_Mismatches_At_End_Forward_Pointer < END_BOUND)
					{
						//if(Tag.Skip) Tag.Start=Tag.End;
						Two_Mismatches_At_End_Forward[Two_Mismatches_At_End_Forward_Pointer]=Tag;
						Two_Mismatches_At_End_Forward_Pointer++;
					}
					return;
				}
				if(Mismatches_Forward_Pointer < ARRAY_BOUND)
				{
					Mismatches_Forward[Mismatches_Forward_Pointer]=Tag;
					Mismatches_Forward_Pointer++;
				}
				if (HIGHSCAN && !Hits && Mismatches_Forward_Pointer+Mismatches_Backward_Pointer+Left_Mishits_Pointer+Right_Mishits_Pointer>EXTCUT) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
			}
			return;
		}
	}
}

void Search_Forwards_Indel( SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	if (Tag.Skip) {Tag.End=Tag.Start;Tag.Skip=0;}// maybe optimise later...
	Start=Start-2;//Adjust for offset difference
	SARange TRange;
	int FSStack_Top=0;
	FSSStack[0]=Tag;
	struct SARange Range,Temp_Range;
	while(FSStack_Top!=-1)//While Stack non-empty....
	{
		Range=FSSStack[FSStack_Top];
		FSStack_Top--;		//Pop the range

		if(Range.Mismatch_Pos[MAX_MISMATCHES_BOUND-1]!=INDELMARK)//No indels?
		{
			TRange=Range;
			TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-1]=INDELMARK;
			TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-2]=DELETEMARK;
			TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-3]=(Start+TRange.Level);//indicate position
			FSStack_Top++;
			TRange.Level=TRange.Level+1;
			FSSStack[FSStack_Top]=TRange;
		}

		Branch_Detect(Range,revfmi,Start);//One_Branch(Range,revfmi);

		for(int Branch=0;Branch<4;Branch++)
		{
			if (Branch_Characters[Branch])//This character actually branches
			{
				Temp_Range=Range;
				Temp_Range.Start = Branch_Ranges[Branch].Start;
				Temp_Range.End = Branch_Ranges[Branch].End;
				if(Temp_Range.Mismatch_Pos[MAX_MISMATCHES_BOUND-1]!=INDELMARK)//No indels?
				{
					TRange=Temp_Range;
					TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-1]=INDELMARK;
					TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-2]=INSERTMARK;
					TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-3]=(Start+TRange.Level);//indicate position
					TRange.Mismatch_Char=TRange.Mismatch_Char | (Branch<<2*2);
					FSStack_Top++;
					FSSStack[FSStack_Top]=TRange;
				}
				if (Current_Tag[Temp_Range.Level+Start]!=Branch)
				{

					Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
					Temp_Range.Mismatch_Pos[Temp_Range.Mismatches] =(Start+Temp_Range.Level);
					Temp_Range.Mismatches++;

				}

				if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
				{
					if(Temp_Range.Level== StringLength && Temp_Range.Mismatch_Pos[MAX_MISMATCHES_BOUND-1]==INDELMARK)
					{
						Print_LocationX(Temp_Range);
						if (MAXHITS==Hits) return;
					}
					else 
					{
						

						FSStack_Top++;//Push range
						Temp_Range.Level++;
						FSSStack[FSStack_Top]=Temp_Range;
					}
				}
			} 
		}
	}
	return;
}


void One_Branch(struct SARange Tag,BWT *fmi)
{
	unsigned Last, First;
	char Now;
	Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

	if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
	for (unsigned long Pos=First;Pos<=Last;Pos++)
	{
		Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
		Branch_Characters[Now]++;	
	}

}


void Branch_Detect (const struct SARange Tag,BWT *fmi,int Start)
{

	Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;
	if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
	{
		unsigned Last, First;
		char Now;

		if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 

		for (unsigned long Pos=First;Pos<=Last;Pos++)
		{
			Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
			Branch_Characters[Now]++;	
		}

		for (int Branch=0;Branch<4;Branch++)
		{
			if ( !Do_Branch[Tag.Level+Start] && Branch != Current_Tag[Start+Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else if (Branch_Characters[Branch])
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = Branch_Ranges[Branch].Start + Branch_Characters[Branch]-1;// Calculate SAranges
			}
		}
	}
	else
	{
		for (int Branch=0;Branch<4;Branch++)
		{
			if ( !Do_Branch[Tag.Level+Start] && Branch != Current_Tag[Start+Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.End+1, Branch);
				if(!(Branch_Ranges[Branch].End<Branch_Ranges[Branch].Start)) Branch_Characters[Branch]=1;
			}
		}

	}
}

void Branch_Detect_Backwards (const struct SARange Tag,BWT *fmi,int Start)
{

	Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;
	if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
	{
		unsigned Last, First;
		char Now;

		if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 

		for (unsigned long Pos=First;Pos<=Last;Pos++)
		{
			Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
			Branch_Characters[Now]++;	
		}

		for (int Branch=0;Branch<4;Branch++)
		{
			if ( !Do_Branch[Start-Tag.Level] && Branch != Current_Tag[Start-Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else if (Branch_Characters[Branch])
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = Branch_Ranges[Branch].Start + Branch_Characters[Branch]-1;// Calculate SAranges
			}
		}
	}
	else
	{
		for (int Branch=0;Branch<4;Branch++)
		{
			if ( !Do_Branch[Start-Tag.Level] && Branch != Current_Tag[Start-Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.End+1, Branch);
				if(!(Branch_Ranges[Branch].End<Branch_Ranges[Branch].Start)) Branch_Characters[Branch]=1;
			}
		}

	}
}
//}-----------------------------  FORWARD SEARCH ROUTINE  -------------------------------------------------/

//{-----------------------------  BACKWARD SEARCH ROUTINE  -------------------------------------------------/

void Search_Backwards_Exact(struct SARange & Tag,int Start,int StringLength,BWT *fmi)
{

	int Level;
	unsigned long Index,Now,First,Last;
	if (!Tag.Start) return;

	for(;;)	
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			Level=Tag.Level;
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
			{
				Tag.Skip++;
				Tag.End=Tag.Start;
			}
			for(;;)
			{
				Index=Tag.Start;
				if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
				Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
				if (Current_Tag[Start-Level] == Now)
				{
					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0) 
					{
						Tag.Skip++;Tag.End=Tag.Start;
					}

					//Exact_Match_Backward[Level]=Tag;
					//if (!Tag.Skip) Exact_Match_Backward[Level].End=Tag.Start;
					if(Level== StringLength)//no need to print as we are still halfway..
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						return;	
					}
					else {Level++;continue;}
				} 
				else
				{
					Tag.Start=0;//Tag.End=0;
					return;	
				}
			}
		}
		else 
		{
			Exact_Match_Backward[Tag.Level]=Tag;
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
				for (unsigned long Pos=First;Pos<=Last;Pos++)
				{
					Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;	
				}

				Now=Current_Tag[Start-Tag.Level];
				if (Branch_Characters[Now])//we have a match... 
				{

					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges
				
				}
				else
				{
					Tag.Start=0;//Tag.End=0;
				}
			}
			else Get_SARange_Fast(Current_Tag[Start-Tag.Level],Tag,fmi);

			if (Tag.Start!=0)
			{
				if(Tag.Level== StringLength)
				{
					return;
				}
				else {Tag.Level++;continue;}
			} 
			else
			{
				return;
			}

		}
	}
}

//}-----------------------------  BACKWARD SEARCH ROUTINE  -------------------------------------------------/

void Search_Backwards(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	int BMStack_Top=0;
	BMStack[0]=Tag;
	struct SARange Range,Temp_Range;
	while(BMStack_Top!=-1)//While Stack non-empty....
	{
		Range=BMStack[BMStack_Top];
		BMStack_Top--;	//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Backwards_OneSA(Range,Count,Start,StringLength,fmi);
			if(MAXHITS==Hits) return;
		}
		else
		{
			Branch_Detect_Backwards(Range,fmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=Start-Temp_Range.Level;
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches==Count)
							{
								Print_LocationX(Temp_Range);
								if(MAXHITS==Hits) return;
							}
							else continue;
						}
						else
						{
							BMStack_Top++;//Push range
							Temp_Range.Level++;
							BMStack[BMStack_Top]=Temp_Range;
						}
					}
					else 
					{
						if(5 > Count)// 2 mismatches...
						{
							if(Temp_Range.Level != StringLength) Temp_Range.Level++;
							//if((Start-Temp_Range.Level) != 0) Temp_Range.Level++;
							else // 2 mismatches with the last at the end?
							{
								if(Two_Mismatches_At_End_Pointer < END_BOUND)
								{
									Two_Mismatches_At_End[Two_Mismatches_At_End_Pointer]=Temp_Range;
									Two_Mismatches_At_End_Pointer++;
								}
								continue;
							}
							if (Mismatches_Backward_Pointer < ARRAY_BOUND)
							{
								Mismatches_Backward[Mismatches_Backward_Pointer]=Temp_Range;
								Mismatches_Backward_Pointer++;
							}
							if (HIGHSCAN && !Hits &&  Mismatches_Forward_Pointer+Mismatches_Backward_Pointer+Left_Mishits_Pointer+Right_Mishits_Pointer>EXTCUT) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
						}
						continue;
					}
				} 
			}
		}
	}
	return;
}

void Search_Backwards_Indel(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{

	if (!Tag.Start) return;
	int BMStack_Top=0;
	int Branch;
	BMStack[0]=Tag;
	struct SARange Range,Temp_Range,TRange;
	while(BMStack_Top!=-1)//While Stack non-empty....
	{

		Range=BMStack[BMStack_Top];
		BMStack_Top--;	//Pop the range

		if(Range.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] !=INDELMARK)//No indels?
		{
			TRange=Range;
			TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-1]=INDELMARK;
			TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-2]=DELETEMARK;
			TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-3]=(Start-TRange.Level);//indicate position
			BMStack_Top++;
			TRange.Level=TRange.Level+1;
			BMStack[BMStack_Top]=TRange;
		}

		Branch_Detect_Backwards(Range,fmi,Start);
		{
			for (Branch=0;Branch<4;Branch++)
			{
				if(Branch_Characters[Branch] && (Range.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] !=INDELMARK))//No indels?
				{
					TRange=Range;
					TRange.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Range.Start, Branch) + 1;
					TRange.End = Branch_Ranges[Branch].End;//Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges
					TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-1]=INDELMARK;//indicate indel
					TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-2]=INSERTMARK;//indicate insertl
					TRange.Mismatch_Pos[MAX_MISMATCHES_BOUND-3]=(Start-TRange.Level+1);//indicate position
					TRange.Mismatch_Char=TRange.Mismatch_Char | (Branch<<2*2);
					BMStack_Top++;
					BMStack[BMStack_Top]=TRange;
				}
			}

			Branch=Current_Tag[Start-Range.Level];
			if (Branch_Characters[Branch])//This character actually branches
			{
				
				Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Range.Start, Branch) + 1;
				Range.End = Branch_Ranges[Branch].End;//Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

				

				if(Range.Level== StringLength && Range.Mismatch_Pos[MAX_MISMATCHES_BOUND-1]==INDELMARK)
				{
					Print_LocationX(Range);
					if(MAXHITS==Hits) return;
				}
				else
				{
					BMStack_Top++;//Push range
					Range.Level++;
					BMStack[BMStack_Top]=Range;
				}
			} 
		}
	}
	return;

}

void Search_Backwards_OneSA(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if (!Do_Branch[Start-Tag.Level] && Current_Tag[Start-Tag.Level]!=Now) return;  
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches==Count)
				{
					//if (Tag.Skip) Tag.Start= Tag.End; else Tag.End=Tag.Start;
					if (!Tag.Skip) Tag.End=Tag.Start;
					Print_LocationX(Tag);
				}
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else 
		{
			if(5 >= Tag.Mismatches && 5 > Count)// 2 mismatches
			{
				if(!Tag.Skip) Tag.End=Tag.Start;//possibly two mismatch exists..
				if (Tag.Level != StringLength) Tag.Level++; 
				else//two mismatches with the last at the end ... 
				{
					//if(Tag.Skip) Tag.Start=Tag.End;
					if(Two_Mismatches_At_End_Pointer < END_BOUND)
					{
						Two_Mismatches_At_End[Two_Mismatches_At_End_Pointer]=Tag;
						Two_Mismatches_At_End_Pointer++;
					}
					return;
				}
				if(Mismatches_Backward_Pointer < ARRAY_BOUND)
				{
					Mismatches_Backward[Mismatches_Backward_Pointer]=Tag;
					Mismatches_Backward_Pointer++;
				}
				if (HIGHSCAN && !Hits && Mismatches_Forward_Pointer+Mismatches_Backward_Pointer+Left_Mishits_Pointer+Right_Mishits_Pointer>EXTCUT) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
			}
			return;
		} 
	}
}


//{----------------------------------- FILE HANDLING ---------------------------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  File_Open
 *  Description:  Open a file:
 *  Mode - "w" create/Write text from top    - "wb" Write/Binary  -"w+" open text read/write            -"w+b" same as prev/Binary
 *         "r" Read text from top            - "rb" Read/Binary   -"r+" open text read/write/nocreate   -"r+b" same as prev/binary
 *       - "a" text append/write                                  -"a+" text open file read/append      -"a+b" open binary read/append
 *
 * =====================================================================================
 */
FILE* File_Open(const char* File_Name,const char* Mode)
{
	FILE* Handle;
	Handle=fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File %s Cannot be opened ....",File_Name);
		exit(1);
	}
	else return Handle;
}

void File_OpenZ(const char* File_Name,const char* Mode,gzFile & Handle)
{
	Handle=gzopen(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File %s Cannot be opened ....",File_Name);
		exit(1);
	}
}

size_t fwriteX ( const void * ptr, size_t size, size_t count, void* stream )
{
	if (OUTPUT_ZIPPED)
	{
		return gzwrite ( (gzFile*) stream ,ptr, size*count);
	}
	else
	{
		return fwrite ( ptr, size, count, (FILE*) stream );
	}
}


void fprintfX(void* Handle,char* Format, char* String)
{

	if (OUTPUT_ZIPPED)
	{
		gzprintf((gzFile*) Handle, Format, String);
	}
	else
	{
		fprintf((FILE*) Handle, Format, String);
	}
}

//}----------------------------------- FILE HANDLING ---------------------------------------------------------

//{----------------------------------- FM INDEX ROUTINES ---------------------------------------------------------
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initFMI
 *  Description:  Opens FM index fmiFile
 * =====================================================================================
 */
BWT* initFMI(const char* BWTCodeFileName,const char* BWTOccValueFileName ) 

{
	BWT *fmi;
        int PoolSize = 524288;
	MMMasterInitialize(3, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);

	fmi = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, NULL, NULL, NULL, NULL);//Load FM index
	return fmi;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_SARange
 *  Description:  gets the SA range of strings having prefix [New_Char][Range]
 * =====================================================================================
 */
SARange Get_SARange( long New_Char,struct SARange Range,BWT *fmi)
{

	Range.Start = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.Start, New_Char) + 1;
	Range.End = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.End+1, New_Char);
	if (Range.End<Range.Start) 
	{
		Range.Start=0;
	}
	return Range;

}

void Get_SARange_Fast( long New_Char,struct SARange & Range,BWT *fmi)
{

	Range.Start = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.Start, New_Char) + 1;
	Range.End = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.End+1, New_Char);
	if (Range.End<Range.Start) 
	{
		//Range.End=0;
		Range.Start=0;
	}
}

void Get_SARange_Fast_2( long New_Char,struct SARange & Start_Range, struct SARange & Dest_Range,BWT *fmi)
{

	Dest_Range.Start = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Start_Range.Start, New_Char) + 1;
	Dest_Range.End = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Start_Range.End+1, New_Char);
	if (Dest_Range.End<Dest_Range.Start) 
	{
		//Range.End=0;
		Dest_Range.Start=0;
	}
}
//}----------------------------------- FM INDEX ROUTINES ---------------------------------------------------------

//{-----------------------------------PRINT ROUTINE---------------------------------------------------------
void Convert_To_Reverse(SARange &Tag)
{
	char New_Char,Temp_One,Temp_Two;
	unsigned pos;
	unsigned Gap=Tag.End-Tag.Start-1;

	if(Tag.Mismatches)
	{
		for( int i=0;i<Tag.Mismatches;i++)
		{
			pos=Tag.Mismatch_Pos[i];
			Temp_Char_Array[i]=Current_Tag[pos];
			Current_Tag[pos]=Tag.Mismatch_Char>>(2*i) & 3;
		}

	}
	if(LOOKUPSIZE == 3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	Tag.Start=Forward_Start_LookupX[c];Tag.Level=LOOKUPSIZE+1;
	Tag.Skip=0;
	if(!Tag.Skip)
	{
		while (Tag.Level <= STRINGLENGTH)
		{
			New_Char=Current_Tag[1-2+Tag.Level];
			Tag.Start = revfmi->cumulativeFreq[New_Char] + BWTOccValue(revfmi, Tag.Start, New_Char) + 1;
			Tag.Level++;
		}
		Tag.End=Tag.Start+Gap;
	}
	else
	{
		Search_Forwards_Exact(Tag,-1,STRINGLENGTH,revfmi);//Find exact matches and report... if not found get the range for 0|?
		/*Tag.Skip=0;
		while (Tag.Level <= STRINGLENGTH)
		{
			Get_SARange_Fast(Current_Tag[Tag.Level-1],Tag,revfmi);
			if (Tag.Skip) Tag.Skip++;
			else if(Tag.End==Tag.Start && Tag.Start % SAINTERVAL == 0) 
			{
				Tag.Skip++;Tag.End=Tag.Start;
			}
		}*/
	}

	if(Tag.Mismatches)
	{
		for( int i=0;i<Tag.Mismatches;i++)
		{
			Current_Tag[Tag.Mismatch_Pos[i]]=Temp_Char_Array[i];
		}
	} 
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Print_Location
 *  Description:  Prints the hex location of start and end of Range
 * =====================================================================================
 */
void Print_LocationX (struct SARange & Tag)
{
	unsigned Gap;
	unsigned Temp_BCX1[4];
	SARange Temp_Branch_RangesFX[4];
	char *C=Current_Tag;
#define MISM_IN_LOW 2
	if (GIS)//new output format...
	{
		char New_Record[4];
		char Lcount=0;char Rcount=0;

//{----------------------------------------------------- HIGH MISMATCH ------------------------------------------
		if (Extending_Tag)
		{

			if(Extending_Tag==LEFTPASS1 ||Extending_Tag==LEFTPASS1Y )
			{
				//Extend left half....
				int FM=FMIndex;
				int st=STRINGLENGTH;
				SARange Temp_Branch_RangesFX[4];
				unsigned Temp_BCX1[4];
				int Brute_Large = MAX_MISMATCHES/2;
				Brute_Large += MAX_MISMATCHES % 2;

				SARange TagT;
				TagT=Tag;
				Temp_BCX1[0]=Branch_Characters[0];Temp_BCX1[1]=Branch_Characters[1];Temp_BCX1[2]=Branch_Characters[2];Temp_BCX1[3]=Branch_Characters[3];
				memcpy(Temp_Branch_RangesFX,Branch_Ranges,4*sizeof(SARange));
				if(FMIndex!=REVERSE)
				{
					Backwards(TagT,1,STRINGLENGTH);
				}
				FMIndex=REVERSE;
				char E=Extending_Tag;
				if(Extending_Tag==LEFTPASS1) Extending_Tag=LEFTPASSX; else Extending_Tag=LEFTPASSXY;
				if (SEED) TagT.Level=SEEDSIZE+1; else TagT.Level=LHt1+1;
				STRINGLENGTH=STRINGLENGTHt1;
				if(SEED)Search_ForwardsX(TagT,10,1,STRINGLENGTHt1,revfmi);
				else Search_ForwardsX(TagT,Brute_Large,1,STRINGLENGTHt1,revfmi);
				memcpy(Branch_Ranges,Temp_Branch_RangesFX,4*sizeof(SARange));
				Extending_Tag=E;
				//Extending_Tag=LEFTPASS1;
				Branch_Characters[0]=Temp_BCX1[0];Branch_Characters[1]=Temp_BCX1[1];Branch_Characters[2]=Temp_BCX1[2];Branch_Characters[3]=Temp_BCX1[3];
				STRINGLENGTH=st;
				FMIndex=FM;
				return;

			}
			else if(Extending_Tag==RIGHTPASS1 ||Extending_Tag==RIGHTPASS1Y)
			{
				//Extend right half....
				int FM=FMIndex;
				int st=STRINGLENGTH;
				int rh=RH;
				int LHtt=LHt;
				int Brute_Large = MAX_MISMATCHES/2;
				Brute_Large += MAX_MISMATCHES % 2;
				SARange Temp_Branch_RangesBX[4];
				unsigned Temp_BCX1[4];

				Temp_BCX1[0]=Branch_Characters[0];Temp_BCX1[1]=Branch_Characters[1];Temp_BCX1[2]=Branch_Characters[2];Temp_BCX1[3]=Branch_Characters[3];
				memcpy(Temp_Branch_RangesBX,Branch_Ranges,4*sizeof(SARange));
				STRINGLENGTH=RHt1;RH=RHt1;
				SARange TagT=Tag;

				if(FMIndex==REVERSE)
				{
					Reverse(TagT,1,STRINGLENGTH);
				}

				if (SEED){LHt=STRINGLENGTHt-SEEDSIZE;}
				RH=rh;

				for( int i=0;i<TagT.Mismatches;i++)
				{
					TagT.Mismatch_Pos[i]+= LHt1;//LHr;
				}
				STRINGLENGTH=STRINGLENGTHt1;//LHr+RHr;
				FMIndex=FORWARD;
				char E=Extending_Tag;
				if(Extending_Tag==RIGHTPASS1) Extending_Tag=RIGHTPASSX; else Extending_Tag=RIGHTPASSXY;
				//Extending_Tag=RIGHTPASSX;
				TagT.Level=1;
				Current_Tag -= LHt1;
				if(SEED) Search_BackwardsX(TagT,10,LHt1,LHt1,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
				else Search_BackwardsX(TagT,Brute_Large,LHt1,LHt1,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
				memcpy(Branch_Ranges,Temp_Branch_RangesBX,4*sizeof(SARange));
				Extending_Tag=E;
				//if(Extending_Tag==RIGHTPASSX) Extending_Tag=RIGHTPASS1; else Extending_Tag=RIGHTPASS1Y;
				//Extending_Tag=RIGHTPASS;
				FMIndex=FM;
				Current_Tag += LHt1;
				Branch_Characters[0]=Temp_BCX1[0];Branch_Characters[1]=Temp_BCX1[1];Branch_Characters[2]=Temp_BCX1[2];Branch_Characters[3]=Temp_BCX1[3];
				STRINGLENGTH=st;
				LHt=LHtt;
				return;
			}
			else if(Extending_Tag==LEFTPASS || Extending_Tag==LEFTPASSX)
			{
				//Extend left half....
				int FM=FMIndex;
				int st=STRINGLENGTH;
				SARange Temp_Branch_RangesFX[4];
				unsigned Temp_BCX[4];

				SARange Temp_Branch_RangesBX[4];
				SARange TagT;
				TagT=Tag;
				Temp_BCX[0]=Branch_Characters[0];Temp_BCX[1]=Branch_Characters[1];Temp_BCX[2]=Branch_Characters[2];Temp_BCX[3]=Branch_Characters[3];
				memcpy(Temp_Branch_RangesFX,Branch_Ranges,4*sizeof(SARange));
				if(FMIndex!=REVERSE)
				{
					Backwards(TagT,1,STRINGLENGTH);
				}
				FMIndex=REVERSE;
				char E=Extending_Tag;
				Extending_Tag=FALSE;
				if (SEED) TagT.Level=SEEDSIZE+1; else TagT.Level=LHt+1;
				STRINGLENGTH=STRINGLENGTHt;
				if(SEED)Search_ForwardsX(TagT,10,1,STRINGLENGTHt,revfmi);
				else Search_ForwardsX(TagT,MAX_MISMATCHES,1,STRINGLENGTHt,revfmi);
				memcpy(Branch_Ranges,Temp_Branch_RangesFX,4*sizeof(SARange));
				Extending_Tag=E;
				Branch_Characters[0]=Temp_BCX[0];Branch_Characters[1]=Temp_BCX[1];Branch_Characters[2]=Temp_BCX[2];Branch_Characters[3]=Temp_BCX[3];
				STRINGLENGTH=st;
				FMIndex=FM;
				return;

			}
			else if (Extending_Tag==LEFTPASSY)//left half of right half..
			{
				//Extend right half....
				int FM=FMIndex;
				int st=STRINGLENGTH;
				int rh=RH;
				int LHtt=LHt;
				SARange Temp_Branch_RangesBX[4];
				unsigned Temp_BCX[4];

				Temp_BCX[0]=Branch_Characters[0];Temp_BCX[1]=Branch_Characters[1];Temp_BCX[2]=Branch_Characters[2];Temp_BCX[3]=Branch_Characters[3];
				memcpy(Temp_Branch_RangesBX,Branch_Ranges,4*sizeof(SARange));
				STRINGLENGTH=RHt;RH=RHt;
				SARange TagT=Tag;

				if(FMIndex==REVERSE)
				{
					Reverse(TagT,1,STRINGLENGTH);
				}

				if (SEED){LHt=STRINGLENGTHt-SEEDSIZE;}
				RH=rh;

				for( int i=0;i<TagT.Mismatches;i++)
				{
					TagT.Mismatch_Pos[i]+= LHt;
				}
				STRINGLENGTH=STRINGLENGTHt;
				FMIndex=FORWARD;
				char E=Extending_Tag;
				Extending_Tag=FALSE;
				TagT.Level=1;
				Current_Tag -= LHt;
				if(SEED) Search_BackwardsX(TagT,10,LHt,LHt,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
				else Search_BackwardsX(TagT,MAX_MISMATCHES,LHt,LHt,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
				memcpy(Branch_Ranges,Temp_Branch_RangesBX,4*sizeof(SARange));
				Extending_Tag=E;
				//Extending_Tag=LEFTPASSY;
				FMIndex=FM;
				Current_Tag += LHt;
				Branch_Characters[0]=Temp_BCX[0];Branch_Characters[1]=Temp_BCX[1];Branch_Characters[2]=Temp_BCX[2];Branch_Characters[3]=Temp_BCX[3];
				STRINGLENGTH=st;
				LHt=LHtt;
				return;
			}
			else
			{
				//Extend right half....
				int FM=FMIndex;
				int st=STRINGLENGTH;
				int rh=RH;
				int LHtt=LHt;
				SARange Temp_Branch_RangesBX[4];
				unsigned Temp_BCX[4];

				Temp_BCX[0]=Branch_Characters[0];Temp_BCX[1]=Branch_Characters[1];Temp_BCX[2]=Branch_Characters[2];Temp_BCX[3]=Branch_Characters[3];
				memcpy(Temp_Branch_RangesBX,Branch_Ranges,4*sizeof(SARange));
				STRINGLENGTH=RHt;RH=RHt;
				SARange TagT=Tag;

				if(FMIndex==REVERSE)
				{
					Reverse(TagT,1,STRINGLENGTH);
				}

				if (SEED){LHt=STRINGLENGTHt-SEEDSIZE;}
				RH=rh;

				for( int i=0;i<TagT.Mismatches;i++)
				{
					TagT.Mismatch_Pos[i]+= LHt;
				}
				STRINGLENGTH=STRINGLENGTHt;
				FMIndex=FORWARD;
				char E=Extending_Tag;
				Extending_Tag=FALSE;
				TagT.Level=1;
				Current_Tag -= LHt;
				if(SEED) Search_BackwardsX(TagT,10,LHt,LHt,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
				else Search_BackwardsX(TagT,MAX_MISMATCHES,LHt,LHt,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
				memcpy(Branch_Ranges,Temp_Branch_RangesBX,4*sizeof(SARange));
				Extending_Tag=E;//RIGHTPASS;
				FMIndex=FM;
				Current_Tag += LHt;
				Branch_Characters[0]=Temp_BCX[0];Branch_Characters[1]=Temp_BCX[1];Branch_Characters[2]=Temp_BCX[2];Branch_Characters[3]=Temp_BCX[3];
				STRINGLENGTH=st;
				LHt=LHtt;
				return;
			}
		}
//}----------------------------------------------------- HIGH MISMATCH ------------------------------------------

		if (NPOLICY && NCount && Tag.Mismatches)
		{
			//search for places where N's and mismatches overlap
			int Mis_InN=0,Mis_NotinN=0;
			for (int i=0;i<Tag.Mismatches;i++) 
				if (NLocations[Tag.Mismatch_Pos[i]]) {Mismatches_InN[Mis_InN++]=Current_Tag[Tag.Mismatch_Pos[i]];}
			//If mismatches fall on an N, discard it...
			Mis_NotinN=Tag.Mismatches-Mis_InN;
			if(NISMISMATCH){if (NCount+Mis_NotinN > MAX_MISMATCHES ) return;}
			else if (Mis_NotinN > MAX_MISMATCHES) return;
		}

		if (HIGHSCAN && Last_Mismatch_Written )//&& !FILTERUNIQUEHITS)//choose best hits...
		{
			if (Tag.Mismatches > Last_Mismatch_Written) return;//we already have a better hit..
			if (Tag.Mismatches == Last_Mismatch_Written)//we already have a better hit..
			{
				if(!FILTERUNIQUEHITS) return;
				else
				{
					char Same_Hit=TRUE;
					if(MAXSPECIFIED)
					{
						if ( FORWARD==FMIndex || Tag.Skip )// forward search index...
						{
							int S=STRINGLENGTH;
							int Gap=Tag.Skip? 0:(Tag.End-Tag.Start);
							STRINGLENGTH=STRINGLENGTHt;
							Convert_To_Reverse(Tag);
							STRINGLENGTH=S;
							Tag.End=Gap+Tag.Start;
							FMIndex=REVERSE;
						}
						if(Multi_Hit.find(Tag.Start)!=Multi_Hit.end()) return;
						Same_Hit=FALSE;
					}
					else
					{
						for( int i=0;i<Tag.Mismatches;i++)
						{
							if (Mismatch_Check[Tag.Mismatch_Pos[i]]!=(Tag.Mismatch_Char>>(2*i) & 3)+'0')
							{
								Same_Hit = FALSE;break;
							}
						}
					}
					if (Same_Hit) return;
				}
			}
			else //we have a better hit... 
			{
				if(MAXSPECIFIED && !LEAST_MISMATCH)
				{
					//Need all hits...
				}
				else
				{
					Write_Buffer_Ptr=Last_Write_Buffer;
					Hits=0;Stats[6]=0;
				}
			}
		}

		if(Tag.Skip) 
		{
			Gap=0;Tag.Start=Tag.End;
		} 
		else 
		{
			Gap=Tag.End-Tag.Start;
		}
		Record.Gap=Gap;Gap++;
		Stats[In_Mismatch] +=Gap;

		if(STATMODE)
		{
			if (Stats[0] +Stats[1] > THRESHOLD01 ||Stats[0]+Stats[1]+Stats[2]>THRESHOLD012 )//bad tag....
			{
				Write_Buffer_Ptr=Write_Buffer;
				Total_Tags -= Tags_From_Head;Total_Hits -=Hits1;
				Hits=MAXHITS;
				Rollover_Step=TRUE;
				Tag_Stat_Bad=TRUE;
				return;
			}
		}

		if ( ONEFMINDEX && FORWARD==FMIndex )// forward search index...
		{
			int S=STRINGLENGTH;
			STRINGLENGTH=STRINGLENGTHt;
			Convert_To_Reverse(Tag);
			STRINGLENGTH=S;
			if (Tag.Skip) Tag.Start=Tag.End;
		}
		
		if (!Hits) // First Hit, so write header...
		{
			if(!Tag_Printed) Tag_Printed=TRUE;//A hit for this tag has been found...
			if (SCANBOTH && First_Pass_Hits) {Total_Tags--;Tags_From_Head--;}
			Total_Tags++;Tags_From_Head++;
		}


		if (COUNT_ALLHITS)
		{
			if(Hits + Gap > MAXHITS ) Hits = MAXHITS; else Hits=Hits+Gap;//COUNT_ALLHITS =1 => count all the hits, COUNT_ALLHITS =0 => count saranges... i.e. unique hits...
		}

		if ((!MAXSPECIFIED && FILTERUNIQUEHITS && Hits !=1))
		{
			int Reject=0;
			if (In_Mismatch<6) Reject=TRUE;
			if(HIGHSCAN) Write_Buffer_Ptr=Last_Write_Buffer;
			/*else//high mismatch situation...
			{
				if (Tag.Mismatches > Last_Mismatch_Written)//we already have a better hit..
				{
					Hits = 1;
					return;
				}
				if (Gap==1)
				{
					Reject=TRUE;//if mismatches are the same as before...
					if (Tag.Mismatches < Last_Mismatch_Written)//we have a better hit... 
					{
						Write_Buffer_Ptr=Last_Write_Buffer;
						Hits=1;Reject=FALSE;
					}
					else if(Tag.Mismatches>Last_Mismatch_Written) {Hits--;return;}//worse hit. Ignore...
				}
				else//current hit has many occurances...
				{
					Reject=TRUE;
				}

			}*/
			if (Reject)
			{

				Hits=MAXHITS;
				In_Mismatch=7;
				Write_Buffer_Ptr=Last_Write_Buffer;
				return;
			}
		}

		New_Record[0]='%';
		New_Record[1]=Tag_Number;
		New_Record[2]=Current_Tag[STRINGLENGTH];
		New_Record[3]=NCount;
		memcpy(Write_Buffer_Ptr,New_Record,4);Write_Buffer_Ptr += 4;
		if (NCount)
		{
			memcpy(Write_Buffer_Ptr,N,2*NCount);Write_Buffer_Ptr += 2*NCount;
		}
		memcpy(MismatchesGIS.Mismatch_Pos,Tag.Mismatch_Pos,Tag.Mismatches);//MAX_MISMATCHES_BOUND);
		MismatchesGIS.Mismatch_Char=Tag.Mismatch_Char;
		Record.Start=Tag.Start;
		Record.Tag=Tag.Tag;
		Record.Skip=Tag.Skip;
		Record.Mismatches=Tag.Mismatches;
		if (ONEFMINDEX) Record.Index=REVERSE; else Record.Index=FMIndex;
		int Skip_Length=Tag.Mismatches+sizeof(unsigned);
		memcpy(Write_Buffer_Ptr,&Record,sizeof(Record));Write_Buffer_Ptr += sizeof(Record);
		memcpy(Write_Buffer_Ptr,&MismatchesGIS,Skip_Length);Write_Buffer_Ptr += Skip_Length;//+1;
		if (Write_Buffer_Ptr >= Write_Buffer + Write_Buf_Size ) 
		{
			printf("Write_Buffer small... \n");exit(100);
		}
		Last_Mismatch_Written=Tag.Mismatches;Last_In_Mis=In_MismatchX;
		if(MAXSPECIFIED && HIGHSCAN)
		{
			if ( FORWARD==FMIndex || Tag.Skip)// forward search index...
			{
				int S=STRINGLENGTH;
				STRINGLENGTH=STRINGLENGTHt;
				Convert_To_Reverse(Tag);
				STRINGLENGTH=S;
			}
			Multi_Hit[Tag.Start]=true;
			if(LEAST_MISMATCH){MAX_MISMATCHES=Last_Mismatch_Written;Mismatches_Reduced=TRUE;}
			FILTERUNIQUEHITS=TRUE;
			return;
		}
		if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 6) {Hits=MAXHITS;}
		if (Larger_Than_Ten && Last_Mismatch_Written == 12) {Hits=MAXHITS;}
		if (Last_Mismatch_Written >=6)
		{
			if (FILTERUNIQUEHITS)
			{
				if (Hits >= 2)
				{
					Write_Buffer_Ptr=Last_Write_Buffer;
					if (Last_Mismatch_Written == 6) {Hits=MAXHITS;return;}
					Hits=0;Stats[6]=0;
					MAX_MISMATCHES=Last_Mismatch_Written-1;
					Mismatches_Reduced=TRUE;
				}
				else 
				{
					MAX_MISMATCHES=Last_Mismatch_Written;Mismatches_Reduced=TRUE;
					memset(Mismatch_Check,0,MAXTAG);
					for( int i=0;i<Tag.Mismatches;i++)
					{
						Mismatch_Check[Tag.Mismatch_Pos[i]]=(Tag.Mismatch_Char>>(2*i) & 3)+'0';
					}

				}
			}
			else {MAX_MISMATCHES=Last_Mismatch_Written-1;Mismatches_Reduced=TRUE;}
		}
		return;
	}
	
	switch(HITMODE)
	{
		case(DEFAULT):// lowest mismatch output of "best hit"
			{
				//Total_Hits++;
				Total_Tags++;
				if(Tag.Skip) {Gap=0;Tag.Start=Tag.End;} else Gap=Tag.End-Tag.Start;
				Record.Tag=Tag.Tag;
				if(ONEFMINDEX && FORWARD == FMIndex) 
				{
					Convert_To_Reverse(Tag);
					if (Tag.Skip) Tag.Start=Tag.End;
					Record.Index=REVERSE;
				} 
				else 
				{
					Record.Index=FMIndex;
				}
				Record.Start=Tag.Start;
				//Record.Index=FMIndex;
				if (Current_Tag[STRINGLENGTH]=='-') Tag.Mismatches=Tag.Mismatches+100;//code - strand....
				/*if(BITMASK==(BITMASK & (Tag.Mismatch_Pos>>(BTS_PER_LOC)))) 
				{
					if(BITMASK==(BITMASK & (Tag.Mismatch_Pos>>(BTS_PER_LOC*2))))
					{
						Tag.Mismatches = Tag.Mismatches+75;//insert...
					}
					else
					{
						Tag.Mismatches = Tag.Mismatches+50;//del...
					}
				}*/
				Record.Mismatches=Tag.Mismatches;
				Record.Gap=Gap;Gap++;
				fwriteX(&Record,sizeof(Record),1,Output_File);
				if (COUNT_ALLHITS)
				{
					if(Hits + Gap > MAXHITS ) Hits = MAXHITS; else Hits=Hits+Gap;//COUNT_ALLHITS =1 => count all the hits, COUNT_ALLHITS =0 => count saranges... i.e. unique hits...
				}
				else
				{
					Hits++;
				}
			//	return;
			}
			return;
		case(DEEP)://output first MAXHITS
			{
				if (Extending_Tag)
				{
//Gooj
					if(Extending_Tag==LEFTPASS)
					{
//Extend left half....
						int FM=FMIndex;
						int st=STRINGLENGTH;
						SARange TagT;
						TagT=Tag;
						Temp_BCX[0]=Branch_Characters[0];Temp_BCX[1]=Branch_Characters[1];Temp_BCX[2]=Branch_Characters[2];Temp_BCX[3]=Branch_Characters[3];
						memcpy(Temp_Branch_RangesFX,Branch_Ranges,4*sizeof(SARange));
						if(FMIndex!=REVERSE)
						{
							Backwards(TagT,1,STRINGLENGTH);
						}
						FMIndex=REVERSE;
						Extending_Tag=FALSE;
						TagT.Level=RHt+1;
						STRINGLENGTH=STRINGLENGTHt;
						Search_ForwardsX(TagT,MAX_MISMATCHES,1,STRINGLENGTHt,revfmi);
						memcpy(Branch_Ranges,Temp_Branch_RangesFX,4*sizeof(SARange));
						Extending_Tag=LEFTPASS;
						Branch_Characters[0]=Temp_BCX[0];Branch_Characters[1]=Temp_BCX[1];Branch_Characters[2]=Temp_BCX[2];Branch_Characters[3]=Temp_BCX[3];
						STRINGLENGTH=st;
						FMIndex=FM;
						return;

					}
					else
					{
						//Extend right half....
						int FM=FMIndex;
						int st=STRINGLENGTH;
						int rh=RH;
						Temp_BCX[0]=Branch_Characters[0];Temp_BCX[1]=Branch_Characters[1];Temp_BCX[2]=Branch_Characters[2];Temp_BCX[3]=Branch_Characters[3];
						memcpy(Temp_Branch_RangesBX,Branch_Ranges,4*sizeof(SARange));
						STRINGLENGTH=LHt;RH=LHt;
						SARange TagT=Tag;

						if(FMIndex==REVERSE)
						{
							Reverse(TagT,1,STRINGLENGTH);
						}
						for( int i=0;i<TagT.Mismatches;i++)
						{
							TagT.Mismatch_Pos[i]+= LHt;
						}
						RH=rh;STRINGLENGTH=STRINGLENGTHt;
						FMIndex=FORWARD;
						Extending_Tag=FALSE;
						TagT.Level=1;
						Current_Tag -= LHt;
						Search_BackwardsX(TagT,MAX_MISMATCHES,LHt,LHt,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
						memcpy(Branch_Ranges,Temp_Branch_RangesBX,4*sizeof(SARange));
						Extending_Tag=RIGHTPASS;
						FMIndex=FM;
						Current_Tag += LHt;
						Branch_Characters[0]=Temp_BCX[0];Branch_Characters[1]=Temp_BCX[1];Branch_Characters[2]=Temp_BCX[2];Branch_Characters[3]=Temp_BCX[3];
						STRINGLENGTH=st;
						return;
					}
				}
				if(Tag.Skip) {Gap=0;Tag.Start=Tag.End;} else Gap=Tag.End-Tag.Start;
				Stats[In_Mismatch] +=Gap+1;
				/*if (STATMODE)
				{
					if (Stats[0] +Stats[1] > THRESHOLD01 ||Stats[0]+Stats[1]+Stats[2]>THRESHOLD012 )//bad tag....
					{
						Hits=MAXHITS;
						Rollover_Step=TRUE;
						Tag_Stat_Bad=TRUE;
						return;
					}
				}*/
				if ( ONEFMINDEX && FORWARD==FMIndex )// forward search index...
				{
					Convert_To_Reverse(Tag);
					if (Tag.Skip) Tag.Start=Tag.End;
				}
				if(UNIQUEHITS && !Gap)//a unique hit...
				{
					Record.Start=Tag.Start;
					Record.Tag=Tag.Tag;
					if(ONEFMINDEX) Record.Index=REVERSE; else Record.Index=FMIndex;
					Record.Mismatches=Tag.Mismatches;
					fwriteX(&Record,sizeof(Record),1,Unique_File);
					//if (! ALLHITS) {Hits++;Total_Hits++;Total_Tags++;return;}//write only the unique hits...
					if (! ALLHITS) {Hits++;Total_Tags++;return;}//write only the unique hits...
				}
				if(ALLHITS)//do we need to report all hits?
				{
					if(!Hits || Print_Header) //first hit...
					{
						New_Record='@';
						fwriteX(&New_Record,1,1,Output_File);//write new record marker... later make gap global
						
						if(PRINT_DESC) fprintfX(Output_File,"%s",Description);
						for( int i=0;i<STRINGLENGTH; i++) Translated_String[i]=Code_To_Char[Current_Tag[i]];
						Translated_String[STRINGLENGTH]=Current_Tag[STRINGLENGTH];
						fwriteX(Translated_String,1,STRINGLENGTH+1,Output_File);//write tag...
						if(!COUNT_ALLHITS) Hits++;//count distinct saranges
						Total_Tags++;
						Print_Header=FALSE;
						Tag_Printed=TRUE;
					}
					//else
					{
						New_Record='%';
						fwriteX(&New_Record,1,1,Output_File);//write new record marker... later make gap global
					}
					memcpy(Mismatches.Mismatch_Pos,Tag.Mismatch_Pos,MAX_MISMATCHES_BOUND);
					Mismatches.Mismatch_Char=Tag.Mismatch_Char;

					Mismatches.Gap=Gap;Gap++;

					//fwriteX(&Mismatches,sizeof(Mismatches),1,Output_File);
					fwriteX(&Mismatches,sizeof(Mismatches),1,Output_File);
					if (COUNT_ALLHITS)
					{
						if(Hits + Gap > MAXHITS ) Hits = MAXHITS; else Hits=Hits+Gap;//COUNT_ALLHITS =1 => count all the hits, COUNT_ALLHITS =0 => count saranges... i.e. unique hits...
					}
					Record.Start=Tag.Start;
					Record.Tag=Tag.Tag;
					Record.Skip=Tag.Skip;
					if (ONEFMINDEX) Record.Index=REVERSE; else Record.Index=FMIndex;
					Record.Mismatches=Tag.Mismatches;
					fwriteX(&Record,sizeof(Record),1,Output_File);
					//Total_Hits++;
				}
			}
			return;
		case(PAIREND)://output first MAXHITS
			{

//{----------------------------------------------------- HIGH MISMATCH ------------------------------------------
				if (Extending_Tag)
				{
					if(Extending_Tag==LEFTPASS)
					{
						//Extend left half....
						int FM=FMIndex;
						int st=STRINGLENGTH;
						SARange TagT;
						TagT=Tag;
						Temp_BCX[0]=Branch_Characters[0];Temp_BCX[1]=Branch_Characters[1];Temp_BCX[2]=Branch_Characters[2];Temp_BCX[3]=Branch_Characters[3];
						memcpy(Temp_Branch_RangesFX,Branch_Ranges,4*sizeof(SARange));
						if(FMIndex!=REVERSE)
						{
							Backwards(TagT,1,STRINGLENGTH);
						}
						FMIndex=REVERSE;
						Extending_Tag=FALSE;
						TagT.Level=RHt+1;
						STRINGLENGTH=STRINGLENGTHt;
						Search_ForwardsX(TagT,MAX_MISMATCHES,1,STRINGLENGTHt,revfmi);
						memcpy(Branch_Ranges,Temp_Branch_RangesFX,4*sizeof(SARange));
						Extending_Tag=LEFTPASS;
						Branch_Characters[0]=Temp_BCX[0];Branch_Characters[1]=Temp_BCX[1];Branch_Characters[2]=Temp_BCX[2];Branch_Characters[3]=Temp_BCX[3];
						STRINGLENGTH=st;
						FMIndex=FM;
						return;

					}
					else
					{
						//Extend right half....
						int FM=FMIndex;
						int st=STRINGLENGTH;
						int rh=RH;
						Temp_BCX[0]=Branch_Characters[0];Temp_BCX[1]=Branch_Characters[1];Temp_BCX[2]=Branch_Characters[2];Temp_BCX[3]=Branch_Characters[3];
						memcpy(Temp_Branch_RangesBX,Branch_Ranges,4*sizeof(SARange));
						STRINGLENGTH=LHt;RH=LHt;
						SARange TagT=Tag;

						if(FMIndex==REVERSE)
						{
							Reverse(TagT,1,STRINGLENGTH);
						}
						for( int i=0;i<TagT.Mismatches;i++)
						{
							TagT.Mismatch_Pos[i]+= LHt;
						}
						RH=rh;STRINGLENGTH=STRINGLENGTHt;
						FMIndex=FORWARD;
						Extending_Tag=FALSE;
						TagT.Level=1;
						Current_Tag -= LHt;
						Search_BackwardsX(TagT,MAX_MISMATCHES,LHt,LHt,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
						memcpy(Branch_Ranges,Temp_Branch_RangesBX,4*sizeof(SARange));
						Extending_Tag=RIGHTPASS;
						FMIndex=FM;
						Current_Tag += LHt;
						Branch_Characters[0]=Temp_BCX[0];Branch_Characters[1]=Temp_BCX[1];Branch_Characters[2]=Temp_BCX[2];Branch_Characters[3]=Temp_BCX[3];
						STRINGLENGTH=st;
						return;
					}
				}
//}----------------------------------------------------- HIGH MISMATCH ------------------------------------------
				if(Tag.Skip) {Gap=0;Tag.Start=Tag.End;} else Gap=Tag.End-Tag.Start;
				if ( ONEFMINDEX && FORWARD==FMIndex )// forward search index...
				{
					Convert_To_Reverse(Tag);
					//if (Tag.Skip) Tag.Start=Tag.End;
				}

				if(!Hits) //first hit...
				{
					if(!First_Pass_Printed)//first pass gets switched ....
					{
						New_Record='@';
						fwriteX(&New_Record,1,1,Output_File);//write new record marker... later make gap global

						if(PRINT_DESC) fprintfX(Output_File,"%s",Description);
						fprintfX(Output_File,"%s",Tag_Copy);
						First_Pass_Printed=TRUE;
					}
					if(!COUNT_ALLHITS) Hits++;//count distinct saranges
					Total_Tags++;
					Tag_Printed=TRUE;
				}
				//else
				{
					New_Record='%';
					fwriteX(&New_Record,1,1,Output_File);//write new record marker... later make gap global
					if(!First_Pass) Translated_String[0]='H'; else Translated_String[0]='T';
					Translated_String[1]=Current_Tag[STRINGLENGTH];//Translated_String[2]=0;
					fwriteX(&Translated_String,2,1,Output_File);//write head/tail and orientation...
				}
				memcpy(Mismatches.Mismatch_Pos,Tag.Mismatch_Pos,MAX_MISMATCHES_BOUND);
				Mismatches.Mismatch_Char=Tag.Mismatch_Char;
				Mismatches.Gap=Gap;
				Gap++;

				fwriteX(&Mismatches,sizeof(Mismatches),1,Output_File);
				if (COUNT_ALLHITS)
				{
					if(Hits + Gap > MAXHITS ) Hits = MAXHITS; else Hits=Hits+Gap;//COUNT_ALLHITS =1 => count all the hits, COUNT_ALLHITS =0 => count saranges... i.e. unique hits...
				}
				Record.Start=Tag.Start;
				Record.Tag=Tag.Tag;
				Record.Skip=Tag.Skip;
				if (ONEFMINDEX) Record.Index=REVERSE; else Record.Index=FMIndex;
				Record.Mismatches=Tag.Mismatches;
				fwriteX(&Record,sizeof(Record),1,Output_File);
				//Total_Hits++;
			}
			return;
	}
}
//}-----------------------------------PRINT ROUTINE---------------------------------------------------------

//{-----------------------------  Parse Command Line  -------------------------------------------------
void Parse_Command_line(int argc, char* argv[])
{
	int Current_Option=0;
	char* Short_Options ="a:A:dY:yl:Lw:W:t:F:hq:o:b:m:f:Gg:sXx::n:pcC:iI:rRBMzZSe:E:Uuk";//allowed options....
	char* This_Program = argv[0];//Current program name....
	char* Help_String=
"Parameters:\n"
" --help | -h\t\t\t\t Print help\n"
" --query | -q <filename>\t\t Query file(File of Tags)\n"
" --output | -o <filename>\t\t Name of output file\n"
" --usequality | -p \t\t use quality score..\n"
" --genome | -g <filename>\t\t Name of the reference genome\n"
" --buffersize | -b <integer> \t\t Size of disk buffers\n"
" --maxhits | -m <integer> \t\t Maximum number of hits to output (-mall for all hits)...\n"
" --singleindex | -s <integer> \t\t Convert output to be read by reverse index only...\n"
" --statmode | -S \t\t\t Filter reads if too many low mismatches...\n"
" --threshold1 | -w <number> \t\t parameters for statmode..\n"
" --threshold2 | -W <number> \t\t parameters for statmode..\n"
" --mishits | -x <file name> \t\t Print unprocessed hits to a file..\n"
" --maxmismatches | -n <number> \t\t maximum mismatches allowed..\n"
" --scanboth | M \t\t scan both directions of a tag, even if default yields hit..\n"
" --nocomplement | -c <filename> \t\t Do not scan for the reverse complement of tag..\n"
" --indelscan | i \t\t  scan for indels ...\n"
" --ignorehead | I <number> \t\t Ignore specified number of bases from start...\n"
" --forcelength | F <number> \t\t Force the string lengths to be that specified by number...\n"
" --rollover | -r \t\t when MAXHITS hits are not found, scan the reverse complement..\n"
" --printblanks | -B \t\t write an empty record for unprocessed tags..\n"
" --maxtags | -t <number>\t\t Scan only <number> of tags..\n"
" --zip  | -z \t\t\t Compress the output file..\n"
" --zigzag | -Z \t\t\t Scan in zigzag mode...\n"
" --guess | -G \t\t\t Guess strand direction...\n"
" --filteruniquehits | U \t\t\t Write a hit if it is the uniquehit...\n" 
" --nbound | C \t\t\t Maximum number of N's to allow in a read..\n" 
" --randomletter | R \t\t\t Substitute a random letter in place of N in read..\n"
" --nismismatch | X \t\t\t Treat each N as a  mismatch..\n"
" --ignoren | y \t\t\t maxmismatches + number of N mismatches are allowed..\n"
" --email | -e <number> \t\t Notify by email when finished (0=User,1=mail,2=sendmail)..\n"
" --forcesolid| -L \t\t force file to be considered as SOLiD. Use when using dibase encoded files..\n"
" --dust | -d \t\t Filter low complexity reads ...\n"
" --seed | Y \t\t Seed Mode..\n"
" --slider | a <num>\t\t Slider Mode..\n"
" --blast | A <num>\t\t Slider blast Mode..\n"
" --maxext | -E \t\t Limit of extension to use in heuristics...\n"
" --logfile | -l <filename>\t\t Log batman output to this file...\n"
" --leastmismatch | -k \t\t Restrict maxhit output to least mismatches...\n"
" --chopwhitespace |-B \t\t remove white spaces from description...\n"
;
	if(argc == 1) {printf("%s \n",Help_String);exit(0);}
	Source=(char*)malloc(sizeof(char)*5000);//create space for file names...
	char *options, *value; 
	char* Name;int Last_Dash;char* Genome_Name;

	PATTERNFILE=PATTERNFILE_DEFAULT;HITSFILE=HITSFILE_DEFAULT;GENOMEFILE=GENOMEFILE_DEFAULT;
	MISHITFILE=MISHITFILE_DEFAULT;AMBIGUOUSFILE=AMBIGUOUSFILE_DEFAULT;
	BWTFILE = BWTFILE_DEFAULT; 
	OCCFILE = OCCFILE_DEFAULT;
	REVBWTINDEX=REVBWTINDEX_DEFAULT;
	REVOCCFILE=REVOCCFILE_DEFAULT;

	Translated_String[STRINGLENGTH]='\0';//temporary buffer for translation of string...

	for(;;)	
	{
		Current_Option=getopt_long(argc, argv, Short_Options, Long_Options, NULL);
		if (Current_Option == -1 ) break;
		switch(Current_Option)
		{
			case 'h':
				printf("%s \n",Help_String);exit(0);
			case 'u':
				CHOPWHITE=TRUE;
				break;
			case 'd':
				DUST=TRUE;
				break;
			case 'R':
				RANDOMLETTER=TRUE;
				break;
			case 'k':
				LEAST_MISMATCH=TRUE;
				break;
			case 'Y':
				SEED=TRUE;
				SEEDSIZE=atoi(optarg);
				break;
			case 'L':
				FORCESOLID=TRUE;
				break;
			case 'X':
				NPOLICY=TRUE;NISMISMATCH=TRUE;
				break;
			case 'y':
				NPOLICY=TRUE;
				break;
			case 't':
				MAX_TAGS_TO_PROCESS=atoi(optarg);
				break;
			case 'a':
				SLIDER=atoi(optarg);
				break;
			case 'A':
				SLIDER=atoi(optarg);
				SLIDERBLAST=TRUE;
				break;
			case 'U':
				FILTERUNIQUEHITS=TRUE;
				ZIGZAGMODE=TRUE;TAG_GUESS=FALSE;
				break;
			case 'S':
				STATMODE=TRUE;
				break;
			case 'z':
				OUTPUT_ZIPPED=TRUE;
				break;
			case 'Z':
				ZIGZAGMODE=TRUE;TAG_GUESS=FALSE;
				break;
			case 'G':
				ZIGZAGMODE=FALSE;TAG_GUESS=TRUE;
				break;
			case 'i':
				INDELSCAN=TRUE;
				break;
			case 'r':
				ROLLOVER=TRUE;
				break;
			case 'B':
				PRINTBLANKS=TRUE;
				break;
			case 'M':
				SCANBOTH=TRUE;
				break;
			case 'C':
				NBOUND=atoi(optarg);
				break;
			case 'w':
				THRESHOLD01=atoi(optarg);
				break;
			case 'W':
				THRESHOLD012=atoi(optarg);
				break;

			case 'e':
				MAILCLIENT=atoi(optarg);
				EMAIL=TRUE;
				break;
			case 'I':
				IGNOREHEAD=atoi(optarg);
				break;
			case 'F':
				FORCELENGTH=atoi(optarg);
				break;
			case 'c':
				DIRECTIONS=1;
				TAG_GUESS=FALSE;
				SCANCOMPONLY=TRUE;
				break;

			case 'p':
				USEQUALITY=TRUE;
				break;
			case 'n':
				MAX_MISMATCHES = atoi(optarg);
				if (MAX_MISMATCHES <0 or MAX_MISMATCHES >MAX_MISMATCHES_BOUND) MAX_MISMATCHES=5;
				break;
			case 'x':
				PRINT_MISHITS=TRUE;
				if(optarg) MISHITFILE=optarg;
				//Mishit_File=File_Open("mishits.fq","w");
				break;
			case 'l':
				LOG=TRUE;
				if(optarg) LOGFILE=optarg;
				break;
			case 's':
				ONEFMINDEX =TRUE;
				break;
			case 'q':
				if(!Patternfile_Count){PATTERNFILE=optarg;}
				else PATTERNFILE1=optarg;
				Patternfile_Count++;
				break;
			case 'o':
				HITSFILE=optarg;
				break;
			case 'b':
				DISKBUFFERSIZE=atol(optarg);
				break;
			case 'm':
				if(!strcmp(optarg,"all")) MAXHITS=UINT_MAX;
				else if(!(MAXHITS=atoi(optarg))) {printf("Maximum hits defaulted to 1\n");MAXHITS=1;};
				FILTERUNIQUEHITS=FALSE;
				if(MAXHITS!=1) MAXSPECIFIED=TRUE;
				break;
			/*case 'l':
				STRINGLENGTH=atoi(optarg);
				break;*/
			case 'E':
				EXTCUT=atoi(optarg);
				break;
			case 'f':
				printf("WARNING: THESE MODES ARE DEPRECATED...\n");
				HITMODE=atoi(optarg);
				GIS=FALSE;
				if(HITMODE ==0) { UNIQUEHITS =FALSE;ALLHITS = FALSE; }
				else if(HITMODE ==1) { HITMODE=DEEP;UNIQUEHITS =TRUE;ALLHITS = FALSE; }
				else if(HITMODE ==2) { HITMODE=DEEP;UNIQUEHITS =TRUE;ALLHITS = TRUE; }
				else if(HITMODE ==3) { HITMODE=DEEP;UNIQUEHITS =FALSE;ALLHITS = TRUE;}
				else if(HITMODE ==4) { HITMODE=DEEP;UNIQUEHITS =FALSE;ALLHITS = TRUE;PRINT_DESC=FALSE;}
				break;
			case 'g':
				Name=optarg;Last_Dash=0;Genome_Name=optarg;
				for(;Name[0]!=0;Name++)
				{
					if (Name[0]=='/') 
					{
						Last_Dash++;Genome_Name=Name;
					}
				}

				REVBWTINDEX = (char*)Source;
				if(Last_Dash) Last_Dash=Genome_Name-optarg+1; else Genome_Name--;
				strncpy(REVBWTINDEX,optarg,Last_Dash);
				REVBWTINDEX[Last_Dash+0]='r';REVBWTINDEX[Last_Dash+1]='e';REVBWTINDEX[Last_Dash+2]='v';
				strcpy(REVBWTINDEX+Last_Dash+3,Genome_Name+1);
				strcat(REVBWTINDEX+Last_Dash+3,".bwt"); 

				BWTFILE=REVBWTINDEX+600;
				strncpy(BWTFILE,optarg,Last_Dash);
				strcpy(BWTFILE+Last_Dash,Genome_Name+1);
				strcat(BWTFILE+Last_Dash,".bwt"); 


				REVOCCFILE = BWTFILE+600;
				strncpy(REVOCCFILE,optarg,Last_Dash);
				REVOCCFILE[Last_Dash+0]='r';REVOCCFILE[Last_Dash+1]='e';REVOCCFILE[Last_Dash+2]='v';
				strcpy(REVOCCFILE+Last_Dash+3,Genome_Name+1);
				strcat(REVOCCFILE+Last_Dash+3,".fmv"); 


				OCCFILE=REVOCCFILE+600;			
				strncpy(OCCFILE,optarg,Last_Dash);
				strcpy(OCCFILE+Last_Dash,Genome_Name+1);
				strcat(OCCFILE+Last_Dash,".fmv"); 


				break;
			default:
				printf("%s \n",Help_String);
				exit(0);
		}
	}	
	if (SCANBOTH) ROLLOVER=TRUE;
	if (FILTERUNIQUEHITS) MAXHITS=2;
	Patternfile_Count--;

}

//}-----------------------------  Parse Command Line  -------------------------------------------------
void Build_Tables()
{
	if (NORMAL_TAGS)
	{
		if (STRINGLENGTH < 28) LOOKUPSIZE =3;
		if (STRINGLENGTH <=51 && MAX_MISMATCHES >5) LOOKUPSIZE = 3;
	}
	else
	{
		if (PAIR_LENGTH_LEFT< 28) LOOKUPSIZE =3;if (PAIR_LENGTH_RIGHT< 28) LOOKUPSIZE =3;
		if (PAIR_LENGTH_LEFT <=51 && MAX_MISMATCHES >5) LOOKUPSIZE = 3;if (PAIR_LENGTH_RIGHT <=51 && MAX_MISMATCHES >5) LOOKUPSIZE = 3;
	}
	Range LRange;
	for( int i=0;i<4;i++)//fill lookup tables for first character...
	{
		Forward_Start_Lookup[i]=revfmi->cumulativeFreq[i] + 1;
		Forward_End_Lookup[i]=revfmi->cumulativeFreq[i + 1];
		Backward_Start_Lookup[i]=fwfmi->cumulativeFreq[i] + 1;
		Backward_End_Lookup[i]=fwfmi->cumulativeFreq[i + 1];//
		LRange.Start=revfmi->cumulativeFreq[i] + 1;
		LRange.End=revfmi->cumulativeFreq[i + 1];
		if (LRange.Start>LRange.End) LRange.Start=0;
		LRange.Label=i;
		for(int j=0;j<4;j++) Build_Preindex_Forward(LRange, 2, j);
		LRange.Start=fwfmi->cumulativeFreq[i] + 1;
		LRange.End=fwfmi->cumulativeFreq[i + 1];
		if (LRange.Start>LRange.End) LRange.Start=0;
		LRange.Label=i;
		for(int j=0;j<4;j++) Build_Preindex_Backward(LRange, 2, j);

	}
}

void Allocate_Memory()
{
	int STRINGLENGTH =36;
	int MAXSTRINGLEN=255;
	int Max_Allocate=1;
	int Max_Limit=5;
	for (int i=0;i<Max_Limit-1;i++) Max_Allocate=Max_Allocate*STRINGLENGTH;
	Forward_Start_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*LOOKUPSIZE));
	Forward_End_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*LOOKUPSIZE));
	Backward_Start_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*LOOKUPSIZE));
	Backward_End_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*LOOKUPSIZE));
	BMHStack=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	FSHStack=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	FSHStackX0X=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	FSSStack=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	FSSStackX=(SARange*)malloc(sizeof(SARange)*42*MAXSTRINGLEN);	
	BMStack=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	BMStackX=(SARange*)malloc(sizeof(SARange)*42*MAXSTRINGLEN);	
	BMStack_X11=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	BMStack_X11H=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	PSBStack=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	Exact_Match_ForwardF=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	Exact_Match_ForwardC=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	//Exact_Match_Forward=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);
	Exact_Match_Backward=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);
	//Exact_Match_BackwardF=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	//Exact_Match_BackwardC=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	//Exact_Match_Backward=Exact_Match_BackwardF;
	ARRAY_BOUND =sizeof(SARange)*4*Max_Allocate;
	Left_MishitsF=(SARange*)malloc(ARRAY_BOUND);	
	Left_Mishits=Left_MishitsF;
	Right_MishitsF=(SARange*)malloc(ARRAY_BOUND);	
	Right_Mishits=Right_MishitsF;
	Mismatches_BackwardF=(SARange*)malloc(ARRAY_BOUND);//STRINGLENGTH*STRINGLENGTH);//*STRINGLENGTH);	
	Mismatches_Backward=Mismatches_BackwardF;
	Mismatches_ForwardF=(SARange*)malloc(ARRAY_BOUND);//STRINGLENGTH*STRINGLENGTH*STRINGLENGTH);	
	Mismatches_Forward=Mismatches_ForwardF;
        END_BOUND=sizeof(SARange)*2*STRINGLENGTH*STRINGLENGTH*STRINGLENGTH;
	Two_Mismatches_At_End_ForwardF=(SARange*)malloc(END_BOUND);	
	Two_Mismatches_At_End_Forward=Two_Mismatches_At_End_ForwardF;
	Two_Mismatches_At_EndF=(SARange*)malloc(END_BOUND);	
	Two_Mismatches_At_End=Two_Mismatches_At_EndF;
	Possible_20F=(SARange*)malloc(END_BOUND);///no need extra stringlength.	
	Possible_20=Possible_20F;
	Possible_02F=(SARange*)malloc(sizeof(SARange)*2*STRINGLENGTH*STRINGLENGTH*STRINGLENGTH);///no need..	
	Possible_02=Possible_02F;
	END_BOUND=END_BOUND/sizeof(SARange);
	ARRAY_BOUND=ARRAY_BOUND/sizeof(SARange);

	Write_Buf_Size=((sizeof(Mismatches)+sizeof(Record)+3+100)*1500);//(MAXHITS+1));
	if (MAXHITS==UINT_MAX) Write_Buf_Size=((sizeof(Mismatches)+sizeof(Record)+3+100)*150000);//(MAXHITS+1));
	if(!(Write_Buffer=(char*)malloc(Write_Buf_Size))) 
	{
		printf("out of memory\n");exit(1);
	}
	Write_Buffer_Ptr=Write_Buffer;
	Last_Write_Buffer=Write_Buffer_Ptr;

	if(ZIGZAGMODE)
	{
		Left_MishitsC=(SARange*)malloc(sizeof(SARange)*4*Max_Allocate);	
		Right_MishitsC=(SARange*)malloc(sizeof(SARange)*4*Max_Allocate);	
		Mismatches_BackwardC=(SARange*)malloc(sizeof(SARange)*4*Max_Allocate);//STRINGLENGTH*STRINGLENGTH);//*STRINGLENGTH);	
		Mismatches_ForwardC=(SARange*)malloc(sizeof(SARange)*4*Max_Allocate);//STRINGLENGTH*STRINGLENGTH*STRINGLENGTH);	
		Two_Mismatches_At_End_ForwardC=(SARange*)malloc(sizeof(SARange)*2*STRINGLENGTH*STRINGLENGTH*STRINGLENGTH);	
		Possible_20C=(SARange*)malloc(sizeof(SARange)*2*STRINGLENGTH*STRINGLENGTH*STRINGLENGTH);	
		Two_Mismatches_At_EndC=(SARange*)malloc(sizeof(SARange)*2*STRINGLENGTH*STRINGLENGTH*STRINGLENGTH);	
		Possible_02C=(SARange*)malloc(sizeof(SARange)*2*STRINGLENGTH*STRINGLENGTH*STRINGLENGTH);	
		ORG_MISMATCHES=MAX_MISMATCHES;

		if (NULL==Possible_20C||NULL==Left_MishitsC||NULL== Possible_02C||NULL==Mismatches_BackwardC|NULL==Two_Mismatches_At_End_ForwardC||NULL==Two_Mismatches_At_EndC||NULL==BMHStack||NULL==FSHStackX0X||NULL==FSHStack||NULL==FSSStack||NULL==BMStack_X11H||NULL==BMStack_X11||NULL==BMStack||NULL==Exact_Match_Backward||NULL==Exact_Match_ForwardF||NULL==Exact_Match_ForwardC||NULL==Mismatches_ForwardC||NULL==PSBStack||NULL==Forward_End_LookupX||NULL==Forward_Start_LookupX||NULL==Source)
		{
			printf("out of memory");
			exit(1);
		}
	}

	if (NULL== Possible_02||NULL==Mismatches_Backward||NULL==Two_Mismatches_At_End_Forward||NULL==Two_Mismatches_At_End||NULL==BMHStack||NULL==FSHStackX0X||NULL==FSHStack||NULL==FSSStack||NULL==FSSStackX||NULL==BMStack_X11H||NULL==BMStack_X11||NULL==BMStack||NULL==BMStackX||NULL==Exact_Match_Backward||NULL==Exact_Match_ForwardF||NULL==Exact_Match_ForwardC||NULL==Mismatches_Forward||NULL==PSBStack||NULL==Forward_End_LookupX||NULL==Forward_Start_LookupX||NULL==Source)
	{
		printf("out of memory");
		exit(1);
	}
	if (NULL==Possible_20F||NULL==Left_MishitsF||NULL== Possible_02F||NULL==Mismatches_BackwardF||NULL==Two_Mismatches_At_End_ForwardF||NULL==Two_Mismatches_At_EndF||NULL==BMHStack||NULL==FSHStackX0X||NULL==FSHStack||NULL==FSSStack||NULL==BMStack_X11H||NULL==BMStack_X11||NULL==BMStack||NULL==Exact_Match_Backward||NULL==Exact_Match_ForwardF||NULL==Exact_Match_ForwardC||NULL==Mismatches_ForwardF||NULL==PSBStack||NULL==Forward_End_LookupX||NULL==Forward_Start_LookupX||NULL==Source)
	{
		printf("out of memory");
		exit(1);
	}


	Exact_Match_Forward=Exact_Match_ForwardF;
}

void Open_Files()
{

	//File_OpenZ(PATTERNFILE,"r",Input_File);//Load tags
	Input_File=File_Open(PATTERNFILE,"r");//Load tags
	//if(Patternfile_Count) File_OpenZ(PATTERNFILE1,"r",Mate_File);//Load tags
	if(Patternfile_Count) Mate_File=File_Open(PATTERNFILE1,"r");//Load tags
	
	if (FILTER_AMBIGUOUS) Ambiguous_File=File_Open(AMBIGUOUSFILE,"w");
	if (SLIDERBLAST) BLAST_File=File_Open(BLASTFILE,"w");
	if (PRINT_MISHITS) Mishit_File=File_Open(MISHITFILE,"w"); 
	Log_File=File_Open(LOGFILE,"w"); 
	//fprintf(Log_File,"%d\n",argc);
	for (int x=0; x<Argc; x++) fprintf(Log_File,"%s ",Argv[x]); fprintf(Log_File,"\n "); 
	if(OUTPUT_ZIPPED)
	{
		File_OpenZ(HITSFILE,"w",Output_FileG);//Open output file...
		Output_File=Output_FileG;
	}
	else
	{
		Output_FileF=File_Open(HITSFILE,"w");//Open output file...
		Output_File=Output_FileF;
		if(setvbuf((FILE*)Output_File,NULL,_IOFBF,DISKBUFFERSIZE*sizeof(long)))
		{
			printf("Allocating disk buffers failed...\n");
			exit(1);
		}
	}
	
	//gz_stream *s=(gz_stream*)Input_File;
	//if (s->transparent) INPUT_ZIPPED=FALSE; else INPUT_ZIPPED=TRUE;
	//Input_FileO=s->file;


}

void Init_Variables()
{
	for(;;)//ignore comments...
	{
		//gzgets(Input_File,Description,MAXDES);
		fgets(Description,MAXDES,Input_File);
		if (Description[0] != '#') break;
	}
	//if (gzgets(Input_File,Description,MAXDES)!=0)//Measure tag length
	{
		//gzgets(Input_File,Current_Tag,MAXTAG);
		fgets(Current_Tag,MAXTAG,Input_File);
		if(Current_Tag[2] == '.' || (Current_Tag[2]>='0' && Current_Tag[2] <='3')) {SOLID=TRUE;IGNOREHEAD +=2;Random_Array=Random_ArrayC;}
		if (FORCESOLID) {SOLID=TRUE;IGNOREHEAD +=2;}//dibase...
		for(TAG_COPY_LEN=0;Current_Tag[TAG_COPY_LEN]!='\n' && Current_Tag[TAG_COPY_LEN]!='\r' && Current_Tag[TAG_COPY_LEN]!=0;TAG_COPY_LEN++);//TAG_COPY_LEN++;
		if(Patternfile_Count)
		{
			//gzgets(Mate_File,Description,MAXDES);
			fgets(Description,MAXDES,Mate_File);
			//Current_Tag[TAG_COPY_LEN++]='\t';gzgets(Mate_File,Current_Tag+TAG_COPY_LEN,MAXTAG);
			Current_Tag[TAG_COPY_LEN++]='\t';fgets(Current_Tag+TAG_COPY_LEN,MAXTAG,Mate_File);
			for(TAG_COPY_LEN=0;Current_Tag[TAG_COPY_LEN]!='\n' && Current_Tag[TAG_COPY_LEN]!='\r' && Current_Tag[TAG_COPY_LEN]!=0;TAG_COPY_LEN++);//TAG_COPY_LEN++;
		}
		for(STRINGLENGTH=0;Current_Tag[STRINGLENGTH]!='\n' && Current_Tag[STRINGLENGTH]!='\r' && Current_Tag[STRINGLENGTH]!=0 && Current_Tag[STRINGLENGTH]!=PAIR_END_SEPERATOR;STRINGLENGTH++);
		if(Current_Tag[STRINGLENGTH]==PAIR_END_SEPERATOR) 
		{
			NORMAL_TAGS=FALSE;//we have pair ended tags..
			if(Patternfile_Count) {PAIRING_TYPE=TWOFILE;}else {PAIRING_TYPE=TAB;}
			if (FORCELENGTH) STRINGLENGTH=FORCELENGTH;
			if (IGNOREHEAD) STRINGLENGTH=STRINGLENGTH-IGNOREHEAD;	
			PAIR_LENGTH_LEFT=STRINGLENGTH;
			if (PAIR_LENGTH_LEFT < 28) PLOOKUPSIZE=3;
			PLH=PAIR_LENGTH_LEFT/2;//calculate tag portions...
			PLHQL=PLH/2;
			if ((PAIR_LENGTH_LEFT % 2)) {PLH++;PLHQL++;}	
			PLHQR=PLH-PLHQL;
			PRH=PAIR_LENGTH_LEFT-PLH;
			PRHQL=PRH/2;PRHQR=PRH-PRHQL;

			for(PAIR_LENGTH_RIGHT=0;Current_Tag[STRINGLENGTH+1+PAIR_LENGTH_RIGHT+IGNOREHEAD]!='\n' && Current_Tag[STRINGLENGTH+1+PAIR_LENGTH_RIGHT+IGNOREHEAD]!='\r' && Current_Tag[STRINGLENGTH+1+PAIR_LENGTH_RIGHT+IGNOREHEAD]!=0;PAIR_LENGTH_RIGHT++);
			if (FORCELENGTH) PAIR_LENGTH_RIGHT=FORCELENGTH;
			if (IGNOREHEAD) PAIR_LENGTH_RIGHT -= IGNOREHEAD;	
			if (PAIR_LENGTH_RIGHT < 28) PRLOOKUPSIZE=3;
			PRLH=PAIR_LENGTH_RIGHT/2;//calculate tag portions...
			PRLHQL=PRLH/2;
			if ((PAIR_LENGTH_RIGHT % 2)) {PRLH++;PRLHQL++;}	
			PRLHQR=PRLH-PRLHQL;
			PRRH=PAIR_LENGTH_RIGHT-PRLH;
			PRRHQL=PRRH/2;PRRHQR=PRRH-PRRHQL;
			
			if (PAIR_LENGTH_RIGHT> PAIR_LENGTH_LEFT) STRINGLENGTH=PAIR_LENGTH_RIGHT;
			//gzgets(Input_File,Quality,MAXTAG);//plus
			fgets(Quality,MAXTAG,Input_File);//plus
			if (Quality[0]=='>') FILETYPE=FA;else FILETYPE=FQ;
			if (FILETYPE == FQ && Quality[0] != '+' && Description[0] != '@') {printf("Init_Variables: Cannot determine file type ...\n");exit(1);}
		}
		else
		{
			//gzgets(Input_File,Quality,MAXTAG);//plus
			if(fgets(Quality,MAXTAG,Input_File)==NULL) Quality[0]='>';//plus
			if (Quality[0]=='>') FILETYPE=FA;else FILETYPE=FQ;
			if (FILETYPE == FQ && Quality[0] != '+' && Description[0] != '@') {printf("Init_Variables: Cannot determine file type ...\n");exit(1);}
			//gzgets(Input_File,Quality,MAXTAG);//phred
			fgets(Quality,MAXTAG,Input_File);//phred

			if (FORCELENGTH) STRINGLENGTH=FORCELENGTH;
			if (IGNOREHEAD) STRINGLENGTH=STRINGLENGTH-IGNOREHEAD;	

		}

		if ((PAIR_LENGTH_RIGHT) && (PAIR_LENGTH_RIGHT < 28 || PAIR_LENGTH_LEFT < 28))LOOKUPSIZE=3;
		if (MAX_MISMATCHES > 10) LOOKUPSIZE=3;

		LH=STRINGLENGTH/2;//calculate tag portions...
		LHQL=LH/2;
		if ((STRINGLENGTH % 2)) {LH++;LHQL++;}	
		LHQR=LH-LHQL;
		RH=STRINGLENGTH-LH;
		RHQL=RH/2;RHQR=RH-RHQL;
//For 10 mismatch extension...
		STRINGLENGTHl=RH;//calculate tag portions...
		LHl=STRINGLENGTHl/2;//calculate tag portions...
		LHQLl=LHl/2;
		if ((STRINGLENGTHl % 2)) {LHl++;LHQLl++;}	
		LHQRl=LHl-LHQLl;
		RHl=STRINGLENGTHl-LHl;
		RHQLl=RHl/2;RHQRl=RHl-RHQLl;

		STRINGLENGTHr=LH;//calculate tag portions...
		LHr=STRINGLENGTHr/2;//calculate tag portions...
		LHQLr=LHr/2;
		if ((STRINGLENGTHr % 2)) {LHr++;LHQLr++;}	
		LHQRr=LHr-LHQLr;
		RHr=STRINGLENGTHr-LHr;
		RHQLr=RHr/2;RHQRr=RHr-RHQLr;

	}
	//else {printf("Tag file error : cannot get string length\n");exit(1);}

	for(int i=0;i<4;i++)//sample 4 reads for length...
	{
		//if (gzgets(Input_File,Description,MAXDES)!=0)//des .... read a tag
		if (fgets(Description,MAXDES,Input_File)!=0)//des .... read a tag
		{
			//gzgets(Input_File,Current_Tag,MAXDES);//tag
			fgets(Current_Tag,MAXDES,Input_File);//tag
			if (NORMAL_TAGS)
			{
				//gzgets(Input_File,Quality,MAXTAG);//plus
				fgets(Quality,MAXTAG,Input_File);//plus
				//gzgets(Input_File,Quality,MAXTAG);//phred
				fgets(Quality,MAXTAG,Input_File);//phred
			}
		}
	}

	//Average_Length=ftello64(Input_FileO)/5;
	Average_Length=ftello64(Input_File)/5;
	//fseek(Input_FileO, 0L, SEEK_END);
	fseek(Input_File, 0L, SEEK_END);
	//File_Size = ftello64(Input_FileO);
	File_Size = ftello64(Input_File);
	Number_of_Tags=(File_Size/Average_Length)/20;if(Number_of_Tags <10) Number_of_Tags=10;
	if (MAX_TAGS_TO_PROCESS) Number_of_Tags=1000;

	//gzseek(Input_File,0,SEEK_SET);//go top
	fseek(Input_File,0,SEEK_SET);//go top
	//gzseek(Mate_File,0,SEEK_SET);//go top
	if(Patternfile_Count) fseek(Mate_File,0,SEEK_SET);//go top


	Low_QualityC=(char*)malloc(STRINGLENGTH+1);
	Low_QualityF=(char*)malloc(STRINGLENGTH+1);
	Low_Quality=Low_QualityF;

	NLocations=(char*)malloc(STRINGLENGTH+1);

	Do_All=(char*)malloc(STRINGLENGTH+1);
	

	Char_To_Code['N']=0;Char_To_Code['n']=0;Char_To_Code['A']=0;Char_To_Code['C']=1;Char_To_Code['G']=2;Char_To_Code['T']=3;Char_To_Code['a']=0;Char_To_Code['c']=1;Char_To_Code['g']=2;Char_To_Code['t']=3;Char_To_Code['+']='+';Char_To_Code['-']='-';//we are using character count to store the fmicode for acgt
	Char_To_Code['0']=0;Char_To_Code['1']=1;Char_To_Code['2']=2;Char_To_Code['3']=3;Char_To_Code['.']=0;//for SOLiD
	Char_To_Code[0]=0;Char_To_Code[1]=1;Char_To_Code[2]=2;Char_To_Code[3]=3;//identitiy..
	Char_To_CodeC['N']=3;Char_To_CodeC['n']=3;Char_To_CodeC['.']=3;Char_To_CodeC[0]=3;Char_To_CodeC[1]=2;Char_To_CodeC[2]=1;Char_To_CodeC[3]=0;Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;Char_To_CodeC['-']='-';Char_To_CodeC['+']='+';//we are using character count to store the fmicode for acgt
	for (int i=0;i<MAXTAG+1;i++) Quality_Copy[i]=0;

	if(!MAX_MISMATCHES_H) MAX_MISMATCHES_H=MAX_MISMATCHES;
	if (!MAX_HIGH_QUALITY_LENGTH) MAX_HIGH_QUALITY_LENGTH=STRINGLENGTH-MAX_MISMATCHES_H;
	if (MAX_MISMATCHES >5) Stat_Size=7; else Stat_Size=MAX_MISMATCHES+1;
	if (MAX_MISMATCHES <=5 || FILETYPE != FQ) USEQUALITYH = FALSE;//use quality in heuristics.. 
	if (FILETYPE == FA) {Quality[0]='*';Quality[1]=0;}
	fa.seq=Translated_String;
	if (!NBOUND) 
	{
		if (MAX_MISMATCHES <=2) NBOUND=2; else NBOUND=200;
	}
	if (!NORMAL_TAGS)
	{
		ZIGZAGMODE=TRUE;TAG_GUESS=FALSE;
	}

}

void Load_Indexes()
{
	fwfmi=initFMI(BWTFILE,OCCFILE);//Load FM indexes
	revfmi=initFMI(REVBWTINDEX,REVOCCFILE);
	SOURCELENGTH = fwfmi->textLength;
	if (SOURCELENGTH!=revfmi->textLength)
	{ 
		printf("FM index load error \n"); 
		exit(1);
	}
	FWDInverseSA0=fwfmi->inverseSa0;
}

void Verbose(FILE* out)
{
	fprintf(out,"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n");
	fprintf(out,"Using the genome files\n %s\t %s\n %s\t %s\n", BWTFILE,OCCFILE,REVBWTINDEX,REVOCCFILE); 
	fprintf(out,"Query File : %s \t\t Output file: %s\n",PATTERNFILE,HITSFILE);
	if(Patternfile_Count) fprintf(out,"Mate File : %s \n ",PATTERNFILE1);
	fprintf(out,"Length of Tags: %d\t", STRINGLENGTH);
	fprintf(out,"Mismatches allowed : %d\n",MAX_MISMATCHES);
	if (MAX_MISMATCHES >5){if (EXTCUT == INT_MAX) fprintf(out,"Extentions: Unlimited\n"); else fprintf(out,"Extensions: %u\n",EXTCUT);}
	if (SOLID) 
	{
		if (FORCESOLID) fprintf (out,"DIBASE-SOLiD reads...\n");
		else fprintf (out,"SOLiD reads...\n");
	}
	else {if (FILETYPE == FQ) fprintf(out,"FASTQ file..\n"); else fprintf(out,"FASTA file..\n");}
	if (!NORMAL_TAGS) fprintf(out,"Pair end data ... Left %d  Right %d\n",PAIR_LENGTH_LEFT,PAIR_LENGTH_RIGHT); 
	fprintf(out,"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n");
	if(USEQUALITY) fprintf(out,"Using Quality info ...\n"); else fprintf(out,"Not using Quality info...\n");
	if (DIRECTIONS ==2 ) fprintf(out,"Scanning complement...\n"); else fprintf(out,"Not scanning complement...\n");
	if(PRINT_MISHITS) fprintf(out,"Mishits file: %s\n",MISHITFILE);
	if(FILTER_AMBIGUOUS) fprintf(out,"Ambiguous file: %s\n",AMBIGUOUSFILE);
	if(INDELSCAN && MAX_MISMATCHES >=3 ) fprintf(out,"Scanning indels...\n"); else fprintf(out,"Not scanning indels...\n");
	if (SOLID) {if( IGNOREHEAD>2) fprintf (out,"Ignoring %d bases from start ..\n", IGNOREHEAD-2);}
	else if (IGNOREHEAD) fprintf (out,"Ignoring %d bases from start ..\n", IGNOREHEAD);
	if(INPUT_ZIPPED) fprintf (out,"==< Reading compressed Tags >==\n");
	if(OUTPUT_ZIPPED) fprintf (out,"==|Compressing hits|==\n");
	if(ZIGZAGMODE) fprintf (out,"^^Scanning in zigzag mode^^\n");
	if(STATMODE) fprintf (out,"Statmode...\n");
	if (MAXHITS == UINT_MAX) fprintf(out,"Output all hits..\n");
	if(FILTERUNIQUEHITS) fprintf(out,"Writing only unique hits...\n");
	if(EMAIL) fprintf (out,"Confirm by email..\n");
	if(DUST) fprintf(out,"Filtering low complexity reads...\n");
	if(HEURISTIC)
	{
		fprintf (out,"Heuristics on with length = %d, mismatches %d\n",MAX_HIGH_QUALITY_LENGTH,MAX_MISMATCHES_H);
	}
	if(SCANBOTH) fprintf (out,"Scanning both directions of the tag...\n");
}


void Show_Progress(unsigned Percentage)
{
	if (Percentage >=97) return;
	printf("+%d%\b\b\b",Percentage);
	fprintf(Log_File,"%d%%->",Percentage);
	fflush(stdout);
}

void Read_INI()
{

	dictionary* Dictionary;
	Dictionary=iniparser_load("batman.ini",FALSE);
	if (!Dictionary) Dictionary=iniparser_load("~/batman.ini",FALSE);
	if (Dictionary)
	{
		MAX_MISMATCHES = iniparser_getint(Dictionary,"settings:mismatches",4);
		INDELSCAN = iniparser_getint(Dictionary,"settings:indels",1);
		USEQUALITY = iniparser_getint(Dictionary,"settings:quality",0);
		ONEFMINDEX=iniparser_getint(Dictionary,"settings:singleindex",0);
		SUPERACCURATE=iniparser_getint(Dictionary,"settings:singleindex",0);
		HEURISTIC=iniparser_getint(Dictionary,"settings:heuristics",0);
		IGNOREHEAD=iniparser_getint(Dictionary,"settings:ignorehead",0);

		PRINTBLANKS=iniparser_getint(Dictionary,"scan:printblanks",0);
		ROLLOVER=iniparser_getint(Dictionary,"scan:rollover",0);
		SCANBOTH=iniparser_getint(Dictionary,"scan:scanboth",0);
		if(SCANBOTH) ROLLOVER=TRUE;
		HITMODE=iniparser_getint(Dictionary,"output:hitmode",0);
		OUTPUT_ZIPPED=iniparser_getint(Dictionary,"output:zipped",0);
		if(HITMODE ==0) { UNIQUEHITS =FALSE;ALLHITS = FALSE; }
		else if(HITMODE ==1) { HITMODE=DEEP;UNIQUEHITS =TRUE;ALLHITS = FALSE; }
		else if(HITMODE ==2) { HITMODE=DEEP;UNIQUEHITS =TRUE;ALLHITS = TRUE; }
		else if(HITMODE ==3) { HITMODE=DEEP;UNIQUEHITS =FALSE;ALLHITS = TRUE;}
		else if(HITMODE ==4) { HITMODE=DEEP;UNIQUEHITS =FALSE;ALLHITS = TRUE;PRINT_DESC=FALSE;}


	}
	iniparser_freedict(Dictionary);
}

//{-----------------------------  Extend_Left  -------------------------------------------------
void Extend_Left_Half()
{

	struct SARange Range,TRange;
	int Start= -1;
	char MAX_MISMATCHES = MAX_MISMATCHES_IN_EXTENSION;
	Left_Mishits_Pointer=0;
	Right_Mishits_Pointer=0;
	Mismatches_Forward_Pointer=0;//first node where SA range was not found, all other nodes will not have matches..
	Mismatches_Backward_Pointer=0;
	Two_Mismatches_At_End_Pointer=0;
	Two_Mismatches_At_End_Forward_Pointer=0;
	Possible_20_Pointer=0;
	Possible_02_Pointer=0;
	int Possible_03_Pointer;
	int Possible_30_Pointer;
	int Possible_04_Pointer;
	int Possible_40_Pointer;
	int Possible_05_Pointer;int Possible_05_PointerF=0;int Possible_05_PointerC=0;
	int Possible_50_Pointer;int Possible_50_PointerF=0;int Possible_50_PointerC=0;
	int Mismatches_Forward_Pointer_Last4;int Mismatches_Forward_Pointer_Last4F=0;int Mismatches_Forward_Pointer_Last4C=0;
	int Mismatches_Forward_Pointer_Last5;int Mismatches_Forward_Pointer_Last5F=0;int Mismatches_Forward_Pointer_Last5C=0;
	int Left_Mishits_Pointer_1;int Left_Mishits_Pointer_1F=0;int Left_Mishits_Pointer_1C=0;

	Start=1-2;//Adjust for offsets...
	SATot=0,S=0;

	FMIndex=REVERSE;
	In_MismatchX=0;
	if(LOOKUPSIZE ==3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Tag=Actual_Tag;
	Range.Level=LOOKUPSIZE+1;Range.Mismatches=0;Range.Skip=0;
	memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
	Search_Forwards_Exact(Range,Start,STRINGLENGTH,revfmi);//Find exact matches and report... if not found get the range for 0|?
	if(MAXHITS==Hits) return;


//{------------------------------------------- ONE MISMATCH ---------------------------------------------------------------------------------------------
	In_MismatchX=1;
//One mismatches...
	FMIndex=REVERSE;
	Range=Exact_Match_Forward[Start+LH];
	if(Range.Start && Range.Tag == Actual_Tag)//if there are hits of the form 0|?
	{
		Range.Level=1;
		if (Range.Skip) SATot++; else SATot += Range.End-Range.Start+1;
		S++;
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			Search_Forwards(Range,1,LH+1,RH,revfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		Search_Forwards(Range,1,LH+1,RH,revfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
		if(MAXHITS==Hits)  return;


	}		
	FMIndex=FORWARD;
	if(LOOKUPSIZE ==3)
	{
		c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
	}
	Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];Range.Level=LOOKUPSIZE+1;
	Range.Mismatches=0;Range.Tag=Actual_Tag;Range.Skip=0;
	memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
	Search_Backwards_Exact( Range,STRINGLENGTH,RH,fwfmi);//Backward scan for ?|0
	if(Range.Start)//if there are possible hits of the form ?|0
	{
		if (Range.Skip) SATot++; else SATot += Range.End-Range.Start+1;
		S++;
		Range.Level=1;
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			Search_Backwards(Range,1,LH,LH,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		Search_Backwards(Range,1,LH,LH,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
		if(MAXHITS==Hits) return;
	}
	if (MAX_MISMATCHES == 1 ) return;
	//}------------------------------------------- ONE MISMATCH ---------------------------------------------------------------------------------------------

	//{------------------------------------------- TWO MISMATCH ---------------------------------------------------------------------------------------------
	FMIndex=REVERSE;
	In_MismatchX=2;
	//if (Last_Mismatch_Written && Last_In_Mis< In_MismatchX) return;
	if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			Two_Mismatches_At_End_Forward[i].Mismatch_Pos[1]= (STRINGLENGTH-1);//mismatches of the form 0|2, with last mismatch at the end...
			Print_LocationX(Two_Mismatches_At_End_Forward[i]);
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}
	Two_Mismatches_At_End_Forward_Pointer=0;

	FMIndex=FORWARD;
	if(Two_Mismatches_At_End_Pointer)
	{
		Do_Branch=Do_All;
		for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
		{
			Print_LocationX(Two_Mismatches_At_End[i]);//Mismatches of the form 2|0, with one mismatch at the first position
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}

	Two_Mismatches_At_End_Pointer=0;
	FMIndex=REVERSE;
	Possible_03_Pointer=Mismatches_Forward_Pointer;
	if(Mismatches_Forward_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Possible_03_Pointer-1;i>=0;i--)
			{
				Search_Forwards(Mismatches_Forward[i],2,LH+1,RH,revfmi);//scan for possible two mismatches of the form 0|2, and store candidates for 0|3
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=Possible_03_Pointer-1;i>=0;i--)
		{
			Search_Forwards(Mismatches_Forward[i],2,LH+1,RH,revfmi);//scan for possible two mismatches of the form 0|2, and store candidates for 0|3
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}

	FMIndex=FORWARD;
	Possible_30_Pointer=Mismatches_Backward_Pointer;
	if(Mismatches_Backward_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Possible_30_Pointer-1;i>=0;i--)
			{
				Search_Backwards(Mismatches_Backward[i],2,LH,LH,fwfmi);//scan for possible two mismatches of the form 2|0, and stores the candidates for 3|0
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=Possible_30_Pointer-1;i>=0;i--)
		{
			Search_Backwards(Mismatches_Backward[i],2,LH,LH,fwfmi);//scan for possible two mismatches of the form 2|0, and stores the candidates for 3|0
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}

	//----------------------------------------------------------------------------------------------------------------------------------------
	if(LOOKUPSIZE==3)
	{
		c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
	}
	Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];
	Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;//Range.Tag=Actual_Tag;
	memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;

	Search_Backwards_Exact_X0( Range,STRINGLENGTH,RHQR,fwfmi);// ?|?|0
	Range.Level=1;

	if(USEQUALITY)
	{
		Do_Branch=Low_Quality;
		Search_Backwards_X10(Range,1,LH + RHQL, RHQL,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if(MAXHITS==Hits) return;
	}

	Do_Branch=Do_All;
	Search_Backwards_X10(Range,1,LH + RHQL, RHQL,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
	if(MAXHITS==Hits) return;
	//----------------------------------------------------------------------------------------------------------------------------------------
	if(LOOKUPSIZE==3)
	{
		c=Current_Tag[LH+0] | (Current_Tag[LH+1]<<2) | (Current_Tag[LH+2]<<4);// | (Current_Tag[LH+3]<<6) | Current_Tag[LH+4]<<8 | (Current_Tag[LH+5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[LH+0] | (Current_Tag[LH+1]<<2) | (Current_Tag[LH+2]<<4) | (Current_Tag[LH+3]<<6) | Current_Tag[LH+4]<<8 | (Current_Tag[LH+5]<<10);//Use lookup table...
	}
	Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];
	Range.Mismatches=0;Range.Level=LOOKUPSIZE+1; Range.Skip=0;
	memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
	Search_Forwards_0X(Range,LH+1,RHQL,revfmi);
	Range.Level=1;TRange=Range;
	if(USEQUALITY)
	{
		Do_Branch=Low_Quality;

		Search_X01(Range,1,LH + RHQL +1,RHQR,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if(MAXHITS==Hits) return;
	}

	Do_Branch=Do_All;
	Search_X01(TRange,1,LH + RHQL +1,RHQR,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
	if(MAXHITS==Hits) return;
	if( MAX_MISMATCHES ==2) return;

	//}------------------------------------------- TWO MISMATCH ---------------------------------------------------------------------------------------------

	//{------------------------------------------- THREE MISMATCH ---------------------------------------------------------------------------------------------
	//Find three mismatches....
	In_MismatchX=3;
	//if (Last_Mismatch_Written && Last_In_Mis< In_MismatchX) return;
	FMIndex=REVERSE;
	if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			Print_LocationX(Two_Mismatches_At_End_Forward[i]);//mismatches of the form 0|3, with last mismatch at the end...
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}
	Two_Mismatches_At_End_Forward_Pointer=0;

	FMIndex=FORWARD;
	if(Two_Mismatches_At_End_Pointer)
	{
		for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
		{
			Print_LocationX(Two_Mismatches_At_End[i]);//Mismatches of the form 3|0, with one mismatch at the first position
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}
	Two_Mismatches_At_End_Pointer=0;

	FMIndex=REVERSE;
	Possible_04_Pointer=Mismatches_Forward_Pointer;
	if(Mismatches_Forward_Pointer!=Possible_03_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Possible_04_Pointer-1;i>=Possible_03_Pointer;i--)
			{
				Search_Forwards(Mismatches_Forward[i],3,LH+1,RH,revfmi);//scan for possible three mismatches of the form 0|3, and finds mismatches of the form 1|2, stores possibles in the form 1|3
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=Possible_04_Pointer-1;i>=Possible_03_Pointer;i--)
		{
			Search_Forwards(Mismatches_Forward[i],3,LH+1,RH,revfmi);//scan for possible three mismatches of the form 0|3, and finds mismatches of the form 1|2, stores possibles in the form 1|3
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}

	FMIndex=FORWARD;
	Possible_40_Pointer=Mismatches_Backward_Pointer;
	if(Mismatches_Backward_Pointer!=Possible_30_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Possible_40_Pointer-1;i>=Possible_30_Pointer;i--)
			{
				Search_Backwards(Mismatches_Backward[i],3,LH,LH,fwfmi);//scan for possible mismatches of the form 3|0, 2|1 and sotres the candidates for 4|0, 3|1
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=Possible_40_Pointer-1;i>=Possible_30_Pointer;i--)
		{
			Search_Backwards(Mismatches_Backward[i],3,LH,LH,fwfmi);//scan for possible mismatches of the form 3|0, 2|1 and sotres the candidates for 4|0, 3|1
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;

	}

	//----------------------------------------------------------------------------------------------------------------------------------------
	FMIndex=REVERSE;
	if(LOOKUPSIZE==3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;
	memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0; Range.Skip=0;
	Search_Forwards_0X(Range,1,LHQL,revfmi);
	Range.Level=1;
	if(USEQUALITY)
	{
		Do_Branch=Low_Quality;
		Search_01X(Range,1,LHQL +1,LHQR,revfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
		if(MAXHITS==Hits) return;
	}

	Do_Branch=Do_All;
	Search_01X(Range,1,LHQL +1,LHQR,revfmi);
	if(MAXHITS==Hits) return;
	//----------------------------------------------------------------------------------------------------------------------------------------
	if(LOOKUPSIZE==3)
	{
		c=Current_Tag[LH-1-0] | (Current_Tag[LH-1-1]<<2) | (Current_Tag[LH-1-2]<<4);// | (Current_Tag[LH-1-3]<<6) | Current_Tag[LH-1-4]<<8 | (Current_Tag[LH-1-5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[LH-1-0] | (Current_Tag[LH-1-1]<<2) | (Current_Tag[LH-1-2]<<4) | (Current_Tag[LH-1-3]<<6) | Current_Tag[LH-1-4]<<8 | (Current_Tag[LH-1-5]<<10);//Use lookup table...
	}

	Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;//Range.Tag=Actual_Tag;
	memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
	Search_Backwards_Exact_X0( Range,LH,LHQR,fwfmi);// ?|0|?
	Range.Level=1;
	if(USEQUALITY)
	{
		TRange=Range;
		Do_Branch=Low_Quality;
		Search_10X(TRange,1,LHQL, LHQL,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3 
		if(MAXHITS==Hits) return;
	}
	Do_Branch=Do_All;
	Search_10X(Range,1,LHQL, LHQL,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
	if(MAXHITS==Hits) return;
	if( MAX_MISMATCHES ==3 && !INDELSCAN) return;
	//}------------------------------------------- THREE MISMATCH ---------------------------------------------------------------------------------------------
	//{-------------------------------------------  INDEL  ---------------------------------------------------------------------------------------------
	if (INDELSCAN)
	{
		FMIndex=REVERSE;
		Range=Exact_Match_Forward[Start+LH];
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
		if(Range.Start && Range.Tag == Actual_Tag)//if there are hits of the form 0|?
		{
			Range.Level=1;
			Do_Branch=Do_All;
			Search_Forwards_Indel(Range,1,LH+1,RH,revfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
			if(MAXHITS==Hits) return;
		}

		FMIndex=FORWARD;
		if(LOOKUPSIZE ==3)
		{
			c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;Range.Tag=Actual_Tag;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		Search_Backwards_Exact( Range,STRINGLENGTH,RH,fwfmi);//Backward scan for ?|0
		if(Range.Start)//if there are possible hits of the form ?|0
		{
			Range.Level=1;

			Do_Branch=Do_All;
			Search_Backwards_Indel(Range,0,LH,LH,fwfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
			if(MAXHITS==Hits) return;
		}
		if( MAX_MISMATCHES ==3 ) return;
	}


	//}-------------------------------------------  INDEL  ---------------------------------------------------------------------------------------------

	//{------------------------------------------- FOUR MISMATCH ---------------------------------------------------------------------------------------------
	In_MismatchX=4;
if (Larger_Than_Ten) return;
if (!FILTERUNIQUEHITS && Last_Mismatch_Written && Last_Mismatch_Written <=8) return;
	//if (Last_Mismatch_Written && Last_In_Mis< In_MismatchX) return;
	FMIndex=REVERSE;
	if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			Print_LocationX(Two_Mismatches_At_End_Forward[i]);//mismatches of the form 0|4, with last mismatch at the end...
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}

	Two_Mismatches_At_End_Forward_Pointer=0;
	FMIndex=FORWARD;
	if(Two_Mismatches_At_End_Pointer)
	{
		for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
		{
			Print_LocationX(Two_Mismatches_At_End[i]);//mismatches of the form 0|4, with one mismatch at the start...
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}
	Two_Mismatches_At_End_Pointer=0;

	FMIndex=REVERSE;
	Possible_05_Pointer=Mismatches_Forward_Pointer;
	int Wrap=FALSE;
	if(Mismatches_Forward_Pointer)
	{
		if (Possible_04_Pointer > 46000) {Mismatches_Forward_Pointer=0;Wrap=TRUE;}
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Possible_04_Pointer;i<Possible_05_Pointer;i++)//Mismatches_Forward_Pointer;i++)
			{
				Search_Forwards(Mismatches_Forward[i],4,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|4, and finds mismatches of the form 1|3, stores possibles in the form 1|4
				if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=Possible_04_Pointer;i<Possible_05_Pointer;i++)//Mismatches_Forward_Pointer;i++)
		{
			Search_Forwards(Mismatches_Forward[i],4,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|4, and finds mismatches of the form 1|3, stores possibles in the form 1|4
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}
	if(Wrap) Possible_05_Pointer=0;
	Mismatches_Forward_Pointer_Last4=Mismatches_Forward_Pointer;

	FMIndex=FORWARD;
	Possible_50_Pointer=Mismatches_Backward_Pointer;
	if(Mismatches_Backward_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Possible_50_Pointer-1;i>=Possible_40_Pointer;i--)//Mismatches_Backward_Pointer-1;i>=0;i--)
			{
				Search_Backwards(Mismatches_Backward[i],4,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
				if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=Possible_50_Pointer-1;i>=Possible_40_Pointer;i--)//Mismatches_Backward_Pointer-1;i>=0;i--)
		{
			Search_Backwards(Mismatches_Backward[i],4,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;

	}

	FMIndex=REVERSE;
	Left_Mishits_Pointer_1=Left_Mishits_Pointer;

	if(Left_Mishits_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=0;i<Left_Mishits_Pointer;i++)
			{
				TRange=Left_Mishits[i];
				if (Left_Mishits[i].Level != LH+1)
				{
					Search_Exact(Left_Mishits[i],-1,LH,revfmi);
					Left_Mishits[i].Level=LH+1;//search exact does not update level..
				}
				if (Left_Mishits[i].Start)
				{
					SATot += Left_Mishits[i].End-Left_Mishits[i].Start+1;
					S++;
					if (S > SACUT && HIGHSCAN) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
					Search_Forwards(Left_Mishits[i],4,1,STRINGLENGTH,revfmi);//find mismatches of the form 022 form, stores possibles of the form 023
				}
				Left_Mishits[i]=TRange;
				if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=0;i<Left_Mishits_Pointer;i++)
		{
			if (Left_Mishits[i].Level != LH+1)
			{
				Search_Exact(Left_Mishits[i],-1,LH,revfmi);
				Left_Mishits[i].Level=LH+1;//search exact does not update level..
			}
			if (Left_Mishits[i].Start)
			{
				if (Left_Mishits[i].Skip) SATot++; else SATot += Left_Mishits[i].End-Left_Mishits[i].Start+1;
				S++;
				if (S > SACUT && HIGHSCAN) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
				Search_Forwards(Left_Mishits[i],4,1,STRINGLENGTH,revfmi);//find mismatches of the form 022 form, stores possibles of the form 023
				if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
			}
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;

	}


	Mismatches_Forward_Pointer_Last5=Mismatches_Forward_Pointer;
	if( Right_Mishits_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=0;i<Right_Mishits_Pointer;i++)
			{
				TRange=Right_Mishits[i];
				if(Right_Mishits[i].Level!=LHQL) 
				{
					Search_Backwards_Exact( Right_Mishits[i],LHQL,LHQL,fwfmi);//finds mismatches of the form 202, stores possibles of the form 203
				}
				if(Right_Mishits[i].Start)
				{	
					Backwards(Right_Mishits[i],1,LH);
					if(Right_Mishits[i].Start)
					{
						Right_Mishits[i].Level=1;
						Search_Forwards(Right_Mishits[i],4,LH+1,RH,revfmi);
						if(MAXHITS==Hits) break;
					}
				}
				if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
				Right_Mishits[i]=TRange;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=0;i<Right_Mishits_Pointer;i++)
		{
			if(Right_Mishits[i].Level!=LHQL) 
			{
				Search_Backwards_Exact( Right_Mishits[i],LHQL,LHQL,fwfmi);//finds mismatches of the form 202, stores possibles of the form 203
			}
			if(Right_Mishits[i].Start)
			{	
				if (Right_Mishits[i].Skip) SATot++; else SATot+=Right_Mishits[i].End-Right_Mishits[i].Start+1;
				S++;
				if (S > SACUT && HIGHSCAN) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
				Backwards(Right_Mishits[i],1,LH);
				if(Right_Mishits[i].Start)
				{
					Right_Mishits[i].Level=1;
					Search_Forwards(Right_Mishits[i],4,LH+1,RH,revfmi);
					if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
					if(MAXHITS==Hits) break;
				}
			}
		}
		if(MAXHITS==Hits) return;
	}


	int LHQLrx=LHQL/2;
	if (LHQL % 2) LHQLrx++; int LHQRrx=LHQL-LHQLrx;

	FMIndex=REVERSE;
	if (LOOKUPSIZE >= LHQLrx)
	{
		Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
	}
	else
	{
		if(LOOKUPSIZE==3)
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0; Range.Skip=0;
	}
	Search_Forwards_0X(Range,1,LHQLrx,revfmi);
	Range.Level=1;
	if(USEQUALITY)
	{
		TRange=Range;
		Do_Branch=Low_Quality;
		Search_01LX(Range,1,LHQLrx +1,LHQRrx,revfmi);
		Range=TRange;
		if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
		if(MAXHITS==Hits) return;
	}

	Do_Branch=Do_All;
	Search_01LX(Range,1,LHQLrx +1,LHQRrx,revfmi);
	if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
	if(MAXHITS==Hits) return;
	//--------------------------------
	if (LOOKUPSIZE >= LHQRrx)
	{
		Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
	}
	else
	{
		if(LOOKUPSIZE==3)
		{
			c=Current_Tag[LHQL-1-0] | (Current_Tag[LHQL-1-1]<<2) | (Current_Tag[LHQL-1-2]<<4);// | (Current_Tag[LH-1-3]<<6) | Current_Tag[LH-1-4]<<8 | (Current_Tag[LH-1-5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[LHQL-1-0] | (Current_Tag[LHQL-1-1]<<2) | (Current_Tag[LHQL-1-2]<<4) | (Current_Tag[LHQL-1-3]<<6) | Current_Tag[LHQL-1-4]<<8 | (Current_Tag[LHQL-1-5]<<10);//Use lookup table...
		}

		Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;//Range.Tag=Actual_Tag;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
	}
	Search_Backwards_Exact_X0( Range,LHQL,LHQRrx,fwfmi);// ?|0|?
	Range.Level=1;
	if(USEQUALITY)
	{
		TRange=Range;
		Do_Branch=Low_Quality;
		Search_10LX(TRange,1,LHQL, LHQL,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3 
		if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==9 ) return;
		if(MAXHITS==Hits) return;
	}
	Do_Branch=Do_All;
	Search_10LX(Range,1,LHQLrx, LHQLrx,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
	if (!FILTERUNIQUEHITS && Last_Mismatch_Written ==8 ) return;
	if(MAXHITS==Hits) return;

	if( MAX_MISMATCHES ==4) return;
	//}------------------------------------------- FOUR MISMATCH ---------------------------------------------------------------------------------------------
	//{------------------------------------------- FIVE MISMATCH ---------------------------------------------------------------------------------------------
	In_MismatchX=5;
	if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
	if (MAX_ORG ==10 && Ext_Scan ==2) return;

	FMIndex=REVERSE;
	if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			Print_LocationX(Two_Mismatches_At_End_Forward[i]);//mismatches of the form 0|4, with last mismatch at the end...
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}
	Two_Mismatches_At_End_Forward_Pointer=0;


	FMIndex=FORWARD;
	if(Two_Mismatches_At_End_Pointer)
	{
		for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
		{
			Print_LocationX(Two_Mismatches_At_End[i]);//mismatches of the form 0|5, with one mismatch at the start...
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}
	Two_Mismatches_At_End_Pointer=0;

	FMIndex=REVERSE;
	if(Mismatches_Forward_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Possible_05_Pointer;i<Mismatches_Forward_Pointer_Last4;i++)//Mismatches_Forward_Pointer;i++)
			{
				Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
				if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=Possible_05_Pointer;i<Mismatches_Forward_Pointer_Last4;i++)//Mismatches_Forward_Pointer;i++)
		{
			Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;

	}


	FMIndex=REVERSE;
	if(Mismatches_Forward_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Mismatches_Forward_Pointer_Last4;i<Mismatches_Forward_Pointer_Last5;i++)//Mismatches_Forward_Pointer;i++)
			{
				Search_Forwards(Mismatches_Forward[i],5,1,STRINGLENGTH,revfmi);//scan for possible five mismatches of the form 0|5, and finds mismatches of the form 1|4,2|3 
				if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=Mismatches_Forward_Pointer_Last4;i<Mismatches_Forward_Pointer_Last5;i++)//Mismatches_Forward_Pointer;i++)
		{
			Search_Forwards(Mismatches_Forward[i],5,1,STRINGLENGTH,revfmi);//scan for possible five mismatches of the form 0|5, and finds mismatches of the form 1|4,2|3 
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}

	FMIndex=REVERSE;
	if(Mismatches_Forward_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Mismatches_Forward_Pointer_Last5;i<Mismatches_Forward_Pointer;i++)//Mismatches_Forward_Pointer;i++)
			{
				Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
				if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=Mismatches_Forward_Pointer_Last5;i<Mismatches_Forward_Pointer;i++)//Mismatches_Forward_Pointer;i++)
		{
			Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;
	}

	FMIndex=FORWARD;
	if(Mismatches_Backward_Pointer!=Possible_50_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Mismatches_Backward_Pointer-1;i>=Possible_50_Pointer;i--)
			{
				Search_Backwards(Mismatches_Backward[i],5,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
				if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=Mismatches_Backward_Pointer-1;i>=Possible_50_Pointer;i--)
		{
			Search_Backwards(Mismatches_Backward[i],5,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
			if(MAXHITS==Hits) break;
		}
		if(MAXHITS==Hits) return;

	}
	FMIndex=FORWARD;
	if (Possible_20_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Possible_20_Pointer-1;i>=0;i--)
			{
				if(Possible_20[i].Level!=RHQL) 
				{
					Possible_20[i].Level++;
					Search_Backwards_Exact( Possible_20[i],LH+RHQL,RHQL,fwfmi);
				}

				if(Possible_20[i].Start)
				{	
					S++;
					SATot += Possible_20[i].End-Possible_20[i].Start+1;
					Possible_20[i].Level=1;
					Search_Backwards(Possible_20[i],5,LH,LH,fwfmi);
					if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
					if(MAXHITS==Hits) break;
				}
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=Possible_20_Pointer-1;i>=0;i--)
		{
			if(Possible_20[i].Level!=RHQL) 
			{
				Possible_20[i].Level++; 
				Search_Backwards_Exact( Possible_20[i],LH+RHQL,RHQL,fwfmi);
			}

			if(Possible_20[i].Start)
			{	
				S++;
				if (Possible_20[i].Skip) SATot++; else SATot += Possible_20[i].End-Possible_20[i].Start+1;
				if (S > SACUT && HIGHSCAN) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
				Possible_20[i].Level=1;
				Search_Backwards(Possible_20[i],5,LH,LH,fwfmi);
				if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
				if(MAXHITS==Hits) break;
			}
		}
		if(MAXHITS==Hits) return;

	}

	if(Possible_02_Pointer)
	{
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=0;i<Possible_02_Pointer;i++)
			{
				TRange=Possible_02[i];
				if(Possible_02[i].Level!=RHQR) 
				{
					Search_Exact( Possible_02[i],LH + RHQL-1 ,RHQR,revfmi);//finds mismatches of the form 202, stores possibles of the form 203
				}
				if(Possible_02[i].Start) 
				{
					Reverse(Possible_02[i],STRINGLENGTH,RH);
					if(Possible_02[i].Start)
					{
						Possible_02[i].Level=1;
						Search_Backwards(Possible_02[i],5,LH,LH,fwfmi);
					}
				}
				if(MAXHITS==Hits) break;
				Possible_02[i]=TRange;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		for(int i=0;i<Possible_02_Pointer;i++)
		{

			if(Possible_02[i].Level!=RHQR) 
			{
				Search_Exact( Possible_02[i],LH + RHQL-1 ,RHQR,revfmi);//finds mismatches of the form 202, stores possibles of the form 203
			}
			if(Possible_02[i].Start) 
			{
				S++;
				if (Possible_02[i].Skip) SATot++; else SATot+=Possible_02[i].End-Possible_02[i].Start+1;
				if (S > SACUT && HIGHSCAN) {Hits=MAXHITS;Tag_Stat_Bad=TRUE;return;}
				Reverse(Possible_02[i],STRINGLENGTH,RH);
				if(Possible_02[i].Start)
				{
					Possible_02[i].Level=1;
					Search_Backwards(Possible_02[i],5,LH,LH,fwfmi);
					if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
				}
			}
			if(MAXHITS==Hits) break;
		}

		if(MAXHITS==Hits) return;
	}

	int RHQLlx=RHQR/2;
	if (RHQR % 2) RHQLlx++; int RHQRlx=RHQR-RHQLlx;
	if (LOOKUPSIZE >= RHQRlx)
	{
		Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
	}
	else
	{
		if(LOOKUPSIZE==3)
		{
			c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];
		Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;//Range.Tag=Actual_Tag;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
	}

	Search_Backwards_Exact_X0( Range,STRINGLENGTH,RHQRlx,fwfmi);// ?|?|0
	Range.Level=1;

	if(USEQUALITY)
	{
		Do_Branch=Low_Quality;
		TRange=Range;
		Search_Backwards_XL10(Range,1,LH + RHQL+RHQLlx, RHQLlx,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
		Range=TRange;
		if(MAXHITS==Hits) return;
	}

	Do_Branch=Do_All;
	Search_Backwards_XL10(Range,1,LH + RHQL+RHQLlx, RHQLlx,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
	if(MAXHITS==Hits) return;



	if (LOOKUPSIZE >= RHQLlx)
	{
		Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
	}
	else
	{
		if(LOOKUPSIZE==3)
		{
			c=Current_Tag[LH+RHQL+0] | (Current_Tag[LH+RHQL+1]<<2) | (Current_Tag[LH+RHQL+2]<<4);// | (Current_Tag[LH+3]<<6) | Current_Tag[LH+4]<<8 | (Current_Tag[LH+5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[LH+RHQL+0] | (Current_Tag[LH+RHQL+1]<<2) | (Current_Tag[LH+RHQL+2]<<4) | (Current_Tag[LH+RHQL+3]<<6) | Current_Tag[LH+RHQL+4]<<8 | (Current_Tag[LH+RHQL+5]<<10);//Use lookup table...
		}
		Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];
		Range.Mismatches=0;Range.Level=LOOKUPSIZE+1; Range.Skip=0;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
	}

	Search_Forwards_0X(Range,LH+RHQL+1,RHQLlx,revfmi);
	Range.Level=1;
	TRange=Range;
	if(USEQUALITY)
	{
		Do_Branch=Low_Quality;
		Search_XL01(Range,1,LH + RHQL+RHQLlx +1,RHQRlx,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if(MAXHITS==Hits) return;
	}

	Do_Branch=Do_All;
	Search_XL01(TRange,1,LH + RHQL+RHQLlx +1,RHQRlx,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
	if (!FILTERUNIQUEHITS && Last_Mismatch_Written == 10 ) return;
	if(MAXHITS==Hits) return;

	//}------------------------------------------- FIVE MISMATCH ---------------------------------------------------------------------------------------------


}

void Search_Backwards_OneSAX(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if (!Do_Branch[Start-Tag.Level] && Current_Tag[Start-Tag.Level]!=Now) return;  
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written && Tag.Mismatches >= Last_Mismatch_Written) return;

		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if (Extending_Tag != RIGHTPASSX && Extending_Tag != RIGHTPASSXY) 
				{
					if (Tag.Skip) Tag.Start= Tag.End; else Tag.End=Tag.Start;
					Print_LocationX(Tag);
				}
				else
				{
					if(!Tag.Skip) Tag.End=Tag.Start;
char E=Extending_Tag;
					if(Extending_Tag==RIGHTPASSX) Extending_Tag=LEFTPASS; else Extending_Tag=LEFTPASSY;
					Print_LocationX(Tag);
Extending_Tag=E;
					//if(Extending_Tag==LEFTPASS) Extending_Tag=RIGHTPASSX; else Extending_Tag=RIGHTPASSXY;
				}
/*
				if (Extending_Tag != RIGHTPASSX) 
				{
					if (Tag.Skip) Tag.Start= Tag.End; else Tag.End=Tag.Start;
					Print_LocationX(Tag);
				}
				else
				{
					Extending_Tag = LEFTPASS; 
					Print_LocationX(Tag);
					Extending_Tag = RIGHTPASSX; 
				}*/
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else 
		{
			return;
		} 
	}
}

class CompSA
{
	public:
		bool operator()(SARange& Tag1,SARange& Tag2)
		{
			if (Tag1.Mismatches >= Tag2.Mismatches) return true; else return false;
		}
};

void Search_BackwardsX(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	if (!Tag.Start) return;
	int BMStack_Top=0;
	//BMStackX[0]=Tag;
	priority_queue <SARange,vector <SARange>,CompSA> q;
	q.push(Tag);
	struct SARange Range,Temp_Range;
	while(!q.empty())//BMStack_Top!=-1)//While Stack non-empty....
	{
		Range=q.top();q.pop();//BMStackX[BMStack_Top];
		//BMStack_Top--;	//Pop the range
		if (Range.End==Range.Start || Tag.Skip)//does this SArange have only one branch?
		{
			Search_Backwards_OneSAX(Range,Count,Start,StringLength,fmi);
			if (Mismatches_Reduced) {Count=MAX_MISMATCHES;Mismatches_Reduced=FALSE;}
			if(MAXHITS==Hits) return;
		}
		else
		{
			Branch_Detect_Backwards(Range,fmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=Start-Temp_Range.Level;
						Temp_Range.Mismatches++;
						if (!FILTERUNIQUEHITS && Last_Mismatch_Written && Temp_Range.Mismatches >= Last_Mismatch_Written) continue;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if (Extending_Tag != RIGHTPASSX && Extending_Tag != RIGHTPASSXY) 
							{
								Print_LocationX(Temp_Range);
							}
							else
							{
								char E=Extending_Tag;
								if(Extending_Tag==RIGHTPASSX) Extending_Tag=LEFTPASS; else Extending_Tag=LEFTPASSY;
								Print_LocationX(Temp_Range);
								Extending_Tag=E;
							}
							if (Mismatches_Reduced) {Count=MAX_MISMATCHES;Mismatches_Reduced=FALSE;}
							if(MAXHITS==Hits) return;
							else continue;
						}
						else
						{
							//BMStack_Top++;//Push range
							Temp_Range.Level++;
							q.push(Temp_Range);
							//BMStackX[BMStack_Top]=Temp_Range;
						}
					}
					else 
					{
						continue;
					}
				} 
			}
		}
	}
	return;
}


void Search_ForwardsX(const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{

	if (!Tag.Start) return;
	priority_queue <SARange,vector <SARange>,CompSA> q;
	Start=Start-2;//Adjust for offset difference
	int FSStack_Top=0;
	q.push(Tag);//FSSStackX[0]=Tag;
	struct SARange Range,Temp_Range;
	while(!q.empty())//FSStack_Top!=-1)//While Stack non-empty....
	{
		Range=q.top();q.pop();//FSSStackX[FSStack_Top];
		FSStack_Top--;		//Pop the range
		if (Range.End==Range.Start || Tag.Skip)//does this SArange have only one branch?
		{
			Search_Forwards_OneSAX(Range,Count,Start,StringLength,revfmi);
			if (Mismatches_Reduced) {Count=MAX_MISMATCHES;Mismatches_Reduced=FALSE;}
			if(MAXHITS==Hits) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start);//One_Branch(Range,revfmi);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{

						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
						if (!FILTERUNIQUEHITS && Last_Mismatch_Written && Temp_Range.Mismatches >= Last_Mismatch_Written) continue;

					}


					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if (!SEED && Temp_Range.Mismatches < 5) continue;//dont print low  matches
							if (Extending_Tag != LEFTPASSX && Extending_Tag != LEFTPASSXY) 
							{
								Print_LocationX(Temp_Range);
							}
							else
							{
								char E=Extending_Tag;
								if(Extending_Tag==LEFTPASSX) Extending_Tag=LEFTPASS; else Extending_Tag=LEFTPASSY;
								Print_LocationX(Temp_Range);
								Extending_Tag=E;
							}
							if (Mismatches_Reduced) {Count=MAX_MISMATCHES;Mismatches_Reduced=FALSE;}
							if (MAXHITS==Hits) return;
							continue;
						}
						else 
						{

							//FSStack_Top++;//Push range
							Temp_Range.Level++;
							q.push(Temp_Range);
							//FSSStackX[FSStack_Top]=Temp_Range;
						}
					}
					else
					{
						continue;
					}
				} 
			}
		}
	}
	return;
}

void Search_Forwards_OneSAX(struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
		if (!Do_Branch[Tag.Level+Start] && Current_Tag[Tag.Level+Start]!=Now) return;  
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Tag.Level+Start]!=Now)
		{

			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
			if (!FILTERUNIQUEHITS && Last_Mismatch_Written && Tag.Mismatches >= Last_Mismatch_Written) return;

		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if (!SEED && Tag.Mismatches < 5) return;//dont print low  matches
				if (Extending_Tag != LEFTPASSX && Extending_Tag != LEFTPASSXY) 
				{
					if (Tag.Skip) Tag.Start= Tag.End; else Tag.End=Tag.Start;
					Print_LocationX(Tag);
				}
				else
				{
					if(!Tag.Skip) Tag.End=Tag.Start;
char E=Extending_Tag;
					if(Extending_Tag==LEFTPASSX) Extending_Tag=LEFTPASS; else Extending_Tag=LEFTPASSY;
					Print_LocationX(Tag);
Extending_Tag=E;
					//if(Extending_Tag==LEFTPASS) Extending_Tag=LEFTPASSX; else Extending_Tag=LEFTPASSXY;
				}
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else//log 2 mismatches 
		{
			return;
		}
	}
}

//}-----------------------------  Extend_Left  -------------------------------------------------

//{-----------------------------  Zig_Zag  -------------------------------------------------
void Zig_Zag()
{
	struct SARange Range,TRange;
	int Start= -1;
	Left_Mishits_Pointer=0;
	int Left_Mishits_PointerF=0;
	int Left_Mishits_PointerC=0;
	Right_Mishits_Pointer=0;
	int Right_Mishits_PointerF=0;
	int Right_Mishits_PointerC=0;
	//int Possible_20_Pointer=0;
	int Possible_20_PointerF=0;
	int Possible_20_PointerC=0;
	//int Possible_02_Pointer=0;
	int Possible_02_PointerF=0;
	int Possible_02_PointerC=0;
        Mismatches_Forward_Pointer=0;//first node where SA range was not found, all other nodes will not have matches..
	int Mismatches_Forward_PointerF=0;//first node where SA range was not found, all other nodes will not have matches..
	int Mismatches_Forward_PointerC=0;//first node where SA range was not found, all other nodes will not have matches..
	Mismatches_Backward_Pointer=0;
	int Mismatches_Backward_PointerF=0;
	int Mismatches_Backward_PointerC=0;
	Two_Mismatches_At_End_Pointer=0;
	int Two_Mismatches_At_End_PointerF=0;
	int Two_Mismatches_At_End_PointerC=0;
	Two_Mismatches_At_End_Forward_Pointer=0;
	int Two_Mismatches_At_End_Forward_PointerF=0;
	int Two_Mismatches_At_End_Forward_PointerC=0;
	int Possible_03_Pointer;
	int Possible_03_PointerF=0;int Possible_03_PointerC=0;
	int Possible_30_Pointer;
	int Possible_30_PointerF=0;int Possible_30_PointerC=0;
	int Possible_04_Pointer;
	int Possible_04_PointerF=0;
	int Possible_04_PointerC=0;
	int Possible_40_Pointer;
	int Possible_40_PointerF=0;
	int Possible_40_PointerC=0;
	int Possible_05_Pointer;int Possible_05_PointerF=0;int Possible_05_PointerC=0;
	int Possible_50_Pointer;int Possible_50_PointerF=0;int Possible_50_PointerC=0;
	int Mismatches_Forward_Pointer_Last4;int Mismatches_Forward_Pointer_Last4F=0;int Mismatches_Forward_Pointer_Last4C=0;
	int Mismatches_Forward_Pointer_Last5;int Mismatches_Forward_Pointer_Last5F=0;int Mismatches_Forward_Pointer_Last5C=0;
	int Left_Mishits_Pointer_1;int Left_Mishits_Pointer_1F=0;int Left_Mishits_Pointer_1C=0;

	Start=1-2;//Adjust for offsets...

	ZigZag=FORWARD;
	if (!NORMAL_TAGS)
	{
		MAX_MISMATCHES=ORG_MISMATCHES;
		if (MAXHITS == 1) {MAXHITS= 100000;}
	}
	char Best_Hit_Found=FALSE;

	In_Mismatch=0;
	Exact_Match_Forward=Exact_Match_ForwardF;
	Current_Tag=Original+IGNOREHEAD;
	int Org_Length,Complement_Length;
	for(;;)
	{
		FMIndex=REVERSE;
		if(LOOKUPSIZE ==3)
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Tag=Actual_Tag;
		Range.Level=LOOKUPSIZE+1;Range.Mismatches=0;Range.Skip=0;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		Search_Forwards_Exact(Range,Start,STRINGLENGTH,revfmi);//Find exact matches and report... if not found get the range for 0|?
		Complement_Length=Range.Level;
		if(MAXHITS==Hits) return;
		if(ZigZag) break;
		ZigZag++;
		Org_Length=Range.Level;
		Current_Tag=Complement;
		Exact_Match_Forward=Exact_Match_ForwardC;
	}
	if (FILTERUNIQUEHITS && Hits) return;
	if (!NORMAL_TAGS && !Best_Hit_Found && Hits && ORG_MISMATCHES !=0) {MAX_MISMATCHES=1;Best_Hit_Found=TRUE;} //Go opthit + 1
	if (MAX_MISMATCHES == 0) return;
	if (LEAST_MISMATCH && Hits) return;

	if (Complement_Length > Org_Length)//first try the complement
	{
		Guessed=Complement;Guessed_Quality=Low_QualityC+IGNOREHEAD;
		Guess_Complement=Original+IGNOREHEAD;Guessed_Complement_Quality=Low_QualityF+IGNOREHEAD;
		Exact_Match_ForwardC=Exact_Match_ForwardF;Exact_Match_ForwardF=Exact_Match_Forward;
	}
	else
	{

		Guessed=Original+IGNOREHEAD;Guessed_Quality=Low_QualityF+IGNOREHEAD;
		Guess_Complement=Complement;Guessed_Complement_Quality=Low_QualityC+IGNOREHEAD;
	}

//{------------------------------------------- ONE MISMATCH ---------------------------------------------------------------------------------------------
	//One mismatches...
	In_Mismatch=1;
	ZigZag=FORWARD;
	Current_Tag=Guessed;//Original+IGNOREHEAD;
	Exact_Match_Forward=Exact_Match_ForwardF;
	Mismatches_Forward=Mismatches_ForwardF;
	Mismatches_Backward=Mismatches_BackwardF;
	Two_Mismatches_At_End_Forward=Two_Mismatches_At_End_ForwardF;
	Two_Mismatches_At_End=Two_Mismatches_At_EndF;

	Mismatches_Forward_Pointer=Mismatches_Forward_PointerF;
	Mismatches_Backward_Pointer=Mismatches_Backward_PointerF;
	Two_Mismatches_At_End_Pointer=Two_Mismatches_At_End_PointerF;
	Two_Mismatches_At_End_Forward_Pointer=Two_Mismatches_At_End_Forward_PointerF;

	for(;;)
	{
		FMIndex=REVERSE;
		Range=Exact_Match_Forward[Start+LH];
		if(Range.Start && Range.Tag == Actual_Tag)//if there are hits of the form 0|?
		{
			Range.Level=1;
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				Search_Forwards(Range,1,LH+1,RH,revfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			Search_Forwards(Range,1,LH+1,RH,revfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
			if(MAXHITS==Hits) return;


		}		
		FMIndex=FORWARD;
		if(LOOKUPSIZE ==3)
		{
			c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];Range.Level=LOOKUPSIZE+1;
		Range.Mismatches=0;Range.Tag=Actual_Tag;Range.Skip=0;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		Search_Backwards_Exact( Range,STRINGLENGTH,RH,fwfmi);//Backward scan for ?|0
		if(Range.Start)//if there are possible hits of the form ?|0
		{
			Range.Level=1;
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				Search_Backwards(Range,1,LH,LH,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			Search_Backwards(Range,1,LH,LH,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
			if(MAXHITS==Hits) return;
		}
		if(ZigZag) break;
		ZigZag++;
		Current_Tag=Guess_Complement;
		Exact_Match_Forward=Exact_Match_ForwardC;
		Mismatches_Forward=Mismatches_ForwardC;
		Mismatches_Backward=Mismatches_BackwardC;
		Two_Mismatches_At_End_Forward=Two_Mismatches_At_End_ForwardC;
		Two_Mismatches_At_End=Two_Mismatches_At_EndC;

		Mismatches_Forward_PointerF=Mismatches_Forward_Pointer;
		Mismatches_Backward_PointerF=Mismatches_Backward_Pointer;
		Two_Mismatches_At_End_PointerF=Two_Mismatches_At_End_Pointer;
		Two_Mismatches_At_End_Forward_PointerF=Two_Mismatches_At_End_Forward_Pointer;

		Mismatches_Forward_Pointer=Mismatches_Forward_PointerC;
		Mismatches_Backward_Pointer=Mismatches_Backward_PointerC;
		Two_Mismatches_At_End_Forward_Pointer=Two_Mismatches_At_End_Forward_PointerC;
		Two_Mismatches_At_End_Pointer=Two_Mismatches_At_End_PointerC;
	}
	Mismatches_Forward_PointerC=Mismatches_Forward_Pointer;
	Mismatches_Backward_PointerC=Mismatches_Backward_Pointer;
	Two_Mismatches_At_End_Forward_PointerC=Two_Mismatches_At_End_Forward_Pointer;
	Two_Mismatches_At_End_PointerC=Two_Mismatches_At_End_Pointer;

	if (FILTERUNIQUEHITS && Hits) return;
	if (!NORMAL_TAGS && !Best_Hit_Found && Hits && ORG_MISMATCHES !=1) {MAX_MISMATCHES=2;Best_Hit_Found=TRUE;} //Go opthit + 1
	if (NPOLICY && NCount){if (!NISMISMATCH && MAX_MISMATCHES + NCount == 1) return;else if (MAX_MISMATCHES==NCount) return;}
	else if (MAX_MISMATCHES == 1 ) return;
	if (LEAST_MISMATCH && Hits) return;
//}------------------------------------------- ONE MISMATCH ---------------------------------------------------------------------------------------------

//{------------------------------------------- TWO MISMATCH ---------------------------------------------------------------------------------------------
	In_Mismatch=2;
	ZigZag=FORWARD;
	Current_Tag=Guessed;//Original+IGNOREHEAD;
	Exact_Match_Forward=Exact_Match_ForwardF;
	Mismatches_Forward=Mismatches_ForwardF;
	Mismatches_Backward=Mismatches_BackwardF;
	Two_Mismatches_At_End_Forward=Two_Mismatches_At_End_ForwardF;
	Two_Mismatches_At_End=Two_Mismatches_At_EndF;
	Possible_20=Possible_20F;
	Possible_02=Possible_02F;

	Mismatches_Forward_Pointer=Mismatches_Forward_PointerF;
	Mismatches_Backward_Pointer=Mismatches_Backward_PointerF;
	Two_Mismatches_At_End_Pointer=Two_Mismatches_At_End_PointerF;
	Two_Mismatches_At_End_Forward_Pointer=Two_Mismatches_At_End_Forward_PointerF;

	Possible_20_Pointer=Possible_20_PointerF;
	Possible_02_Pointer=Possible_02_PointerF;
	Possible_03_Pointer=Possible_03_PointerF;
	Possible_30_Pointer=Possible_30_PointerF;

	for(;;)
	{
		FMIndex=REVERSE;
		if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
		{
			for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
			{
				Two_Mismatches_At_End_Forward[i].Mismatch_Pos[1]= (STRINGLENGTH-1);//mismatches of the form 0|2, with last mismatch at the end...
				Print_LocationX(Two_Mismatches_At_End_Forward[i]);
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}
		Two_Mismatches_At_End_Forward_Pointer=0;

		FMIndex=FORWARD;
		if(Two_Mismatches_At_End_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
				{
					Print_LocationX(Two_Mismatches_At_End[i]);//Mismatches of the form 2|0, with one mismatch at the first position
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
			{
				Print_LocationX(Two_Mismatches_At_End[i]);//Mismatches of the form 2|0, with one mismatch at the first position
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;

		}

		Two_Mismatches_At_End_Pointer=0;
		FMIndex=REVERSE;
		Possible_03_Pointer=Mismatches_Forward_Pointer;
		if(Mismatches_Forward_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=Possible_03_Pointer-1;i>=0;i--)
				{
					Search_Forwards(Mismatches_Forward[i],2,LH+1,RH,revfmi);//scan for possible two mismatches of the form 0|2, and store candidates for 0|3
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=Possible_03_Pointer-1;i>=0;i--)
			{
				Search_Forwards(Mismatches_Forward[i],2,LH+1,RH,revfmi);//scan for possible two mismatches of the form 0|2, and store candidates for 0|3
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		FMIndex=FORWARD;
		Possible_30_Pointer=Mismatches_Backward_Pointer;
		if(Mismatches_Backward_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=Possible_30_Pointer-1;i>=0;i--)
				{
					Search_Backwards(Mismatches_Backward[i],2,LH,LH,fwfmi);//scan for possible two mismatches of the form 2|0, and stores the candidates for 3|0
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=Possible_30_Pointer-1;i>=0;i--)
			{
				Search_Backwards(Mismatches_Backward[i],2,LH,LH,fwfmi);//scan for possible two mismatches of the form 2|0, and stores the candidates for 3|0
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		//----------------------------------------------------------------------------------------------------------------------------------------
		if(LOOKUPSIZE==3)
		{
			c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];
		Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;Range.Skip=0;//Range.Tag=Actual_Tag;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;

		Search_Backwards_Exact_X0( Range,STRINGLENGTH,RHQR,fwfmi);// ?|?|0
		Range.Level=1;

		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			Search_Backwards_X10(Range,1,LH + RHQL, RHQL,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		Search_Backwards_X10(Range,1,LH + RHQL, RHQL,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if(MAXHITS==Hits) return;

		//----------------------------------------------------------------------------------------------------------------------------------------
		if(LOOKUPSIZE==3)
		{
			c=Current_Tag[LH+0] | (Current_Tag[LH+1]<<2) | (Current_Tag[LH+2]<<4);// | (Current_Tag[LH+3]<<6) | Current_Tag[LH+4]<<8 | (Current_Tag[LH+5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[LH+0] | (Current_Tag[LH+1]<<2) | (Current_Tag[LH+2]<<4) | (Current_Tag[LH+3]<<6) | Current_Tag[LH+4]<<8 | (Current_Tag[LH+5]<<10);//Use lookup table...
		}
		Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];
		Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;Range.Skip=0;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		Search_Forwards_0X(Range,LH+1,RHQL,revfmi);
		Range.Level=1;TRange=Range;
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;

			Search_X01(Range,1,LH + RHQL +1,RHQR,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		Search_X01(TRange,1,LH + RHQL +1,RHQR,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if(MAXHITS==Hits) return;

		if(ZigZag) break;
		ZigZag++;
		Current_Tag=Guess_Complement;
		Exact_Match_Forward=Exact_Match_ForwardC;
		Mismatches_Forward=Mismatches_ForwardC;
		Mismatches_Backward=Mismatches_BackwardC;
		Two_Mismatches_At_End_Forward=Two_Mismatches_At_End_ForwardC;
		Two_Mismatches_At_End=Two_Mismatches_At_EndC;
		Possible_20=Possible_20C;
		Possible_02=Possible_02C;

		Mismatches_Forward_PointerF=Mismatches_Forward_Pointer;
		Mismatches_Backward_PointerF=Mismatches_Backward_Pointer;
		Two_Mismatches_At_End_Forward_PointerF=Two_Mismatches_At_End_Forward_Pointer;
		Two_Mismatches_At_End_PointerF=Two_Mismatches_At_End_Pointer;

		Possible_20_PointerF=Possible_20_Pointer;
		Possible_02_PointerF=Possible_02_Pointer;
		Possible_03_PointerF=Possible_03_Pointer;
		Possible_30_PointerF=Possible_30_Pointer;

		Mismatches_Forward_Pointer=Mismatches_Forward_PointerC;
		Mismatches_Backward_Pointer=Mismatches_Backward_PointerC;
		Two_Mismatches_At_End_Forward_Pointer=Two_Mismatches_At_End_Forward_PointerC;
		Two_Mismatches_At_End_Pointer=Two_Mismatches_At_End_PointerC;

		Possible_20_Pointer=Possible_20_PointerC;
		Possible_02_Pointer=Possible_02_PointerC;
		Possible_03_Pointer=Possible_03_PointerC;
		Possible_30_Pointer=Possible_30_PointerC;

	}
	Mismatches_Forward_PointerC=Mismatches_Forward_Pointer;
	Mismatches_Backward_PointerC=Mismatches_Backward_Pointer;
	Two_Mismatches_At_End_Forward_PointerC=Two_Mismatches_At_End_Forward_Pointer;
	Two_Mismatches_At_End_PointerC=Two_Mismatches_At_End_Pointer;

	Possible_20_PointerC=Possible_20_Pointer;
	Possible_02_PointerC=Possible_02_Pointer;
	Possible_03_PointerC=Possible_03_Pointer;
	Possible_30_PointerC=Possible_30_Pointer;


	if (FILTERUNIQUEHITS && Hits) return;
	if (!NORMAL_TAGS && !Best_Hit_Found && Hits && ORG_MISMATCHES !=2) {MAX_MISMATCHES=3;Best_Hit_Found=TRUE;} //Go opthit + 1
	if (NPOLICY && NCount){if (!NISMISMATCH && MAX_MISMATCHES + NCount == 2) return;else if (MAX_MISMATCHES==NCount) return;}
	else if( MAX_MISMATCHES ==2) return;
	if (LEAST_MISMATCH && Hits) return;

//}------------------------------------------- TWO MISMATCH ---------------------------------------------------------------------------------------------

//{------------------------------------------- THREE MISMATCH ---------------------------------------------------------------------------------------------
	//Find three mismatches....
	In_Mismatch=3;
	ZigZag=FORWARD;
	Current_Tag=Guessed;//Original+IGNOREHEAD;
	Exact_Match_Forward=Exact_Match_ForwardF;
	Mismatches_Forward=Mismatches_ForwardF;
	Mismatches_Backward=Mismatches_BackwardF;
	Two_Mismatches_At_End_Forward=Two_Mismatches_At_End_ForwardF;
	Two_Mismatches_At_End=Two_Mismatches_At_EndF;
	Possible_20=Possible_20F;
	Possible_02=Possible_02F;
	Right_Mishits=Right_MishitsF;
	Left_Mishits=Left_MishitsF;

	Mismatches_Forward_Pointer=Mismatches_Forward_PointerF;
	Mismatches_Backward_Pointer=Mismatches_Backward_PointerF;
	Two_Mismatches_At_End_Forward_Pointer=Two_Mismatches_At_End_Forward_PointerF;
	Two_Mismatches_At_End_Pointer=Two_Mismatches_At_End_PointerF;

	Possible_20_Pointer=Possible_20_PointerF;
	Possible_02_Pointer=Possible_02_PointerF;
	Possible_03_Pointer=Possible_03_PointerF;
	Possible_30_Pointer=Possible_30_PointerF;
	Possible_04_Pointer=Possible_04_PointerF;
	Possible_40_Pointer=Possible_40_PointerF;
	Right_Mishits_Pointer=Right_Mishits_PointerF;
	Left_Mishits_Pointer=Left_Mishits_PointerF;

	for(;;)
	{
		FMIndex=REVERSE;
		if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
		{
			for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
			{
				Print_LocationX(Two_Mismatches_At_End_Forward[i]);//mismatches of the form 0|3, with last mismatch at the end...
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}
		Two_Mismatches_At_End_Forward_Pointer=0;

		FMIndex=FORWARD;
		if(Two_Mismatches_At_End_Pointer)
		{
			for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
			{
				Print_LocationX(Two_Mismatches_At_End[i]);//Mismatches of the form 3|0, with one mismatch at the first position
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}
		Two_Mismatches_At_End_Pointer=0;

		FMIndex=REVERSE;
		Possible_04_Pointer=Mismatches_Forward_Pointer;
		if(Mismatches_Forward_Pointer!=Possible_03_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=Possible_04_Pointer-1;i>=Possible_03_Pointer;i--)
				{
					Search_Forwards(Mismatches_Forward[i],3,LH+1,RH,revfmi);//scan for possible three mismatches of the form 0|3, and finds mismatches of the form 1|2, stores possibles in the form 1|3
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=Possible_04_Pointer-1;i>=Possible_03_Pointer;i--)
			{
				Search_Forwards(Mismatches_Forward[i],3,LH+1,RH,revfmi);//scan for possible three mismatches of the form 0|3, and finds mismatches of the form 1|2, stores possibles in the form 1|3
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		FMIndex=FORWARD;
		Possible_40_Pointer=Mismatches_Backward_Pointer;
		if(Mismatches_Backward_Pointer!=Possible_30_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=Possible_40_Pointer-1;i>=Possible_30_Pointer;i--)
				{
					Search_Backwards(Mismatches_Backward[i],3,LH,LH,fwfmi);//scan for possible mismatches of the form 3|0, 2|1 and sotres the candidates for 4|0, 3|1
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=Possible_40_Pointer-1;i>=Possible_30_Pointer;i--)
			{
				Search_Backwards(Mismatches_Backward[i],3,LH,LH,fwfmi);//scan for possible mismatches of the form 3|0, 2|1 and sotres the candidates for 4|0, 3|1
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;

		}

		//----------------------------------------------------------------------------------------------------------------------------------------
		FMIndex=REVERSE;
		if(LOOKUPSIZE==3)
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		Search_Forwards_0X(Range,1,LHQL,revfmi);
		Range.Level=1;
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			Search_01X(Range,1,LHQL +1,LHQR,revfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		Search_01X(Range,1,LHQL +1,LHQR,revfmi);
		if(MAXHITS==Hits) return;
		//----------------------------------------------------------------------------------------------------------------------------------------
		if(LOOKUPSIZE==3)
		{
			c=Current_Tag[LH-1-0] | (Current_Tag[LH-1-1]<<2) | (Current_Tag[LH-1-2]<<4);// | (Current_Tag[LH-1-3]<<6) | Current_Tag[LH-1-4]<<8 | (Current_Tag[LH-1-5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[LH-1-0] | (Current_Tag[LH-1-1]<<2) | (Current_Tag[LH-1-2]<<4) | (Current_Tag[LH-1-3]<<6) | Current_Tag[LH-1-4]<<8 | (Current_Tag[LH-1-5]<<10);//Use lookup table...
		}

		Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;//Range.Tag=Actual_Tag;
		memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		Search_Backwards_Exact_X0( Range,LH,LHQR,fwfmi);// ?|0|?
		Range.Level=1;
		if(USEQUALITY)
		{
			TRange=Range;
			Do_Branch=Low_Quality;
			Search_10X(TRange,1,LHQL, LHQL,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
			if(MAXHITS==Hits) return;
		}
		Do_Branch=Do_All;
		Search_10X(Range,1,LHQL, LHQL,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
		if(MAXHITS==Hits) return;
		
		if(ZigZag) break;
		ZigZag++;
		Current_Tag=Guess_Complement;
		Exact_Match_Forward=Exact_Match_ForwardC;
		Mismatches_Forward=Mismatches_ForwardC;
		Two_Mismatches_At_End_Forward=Two_Mismatches_At_End_ForwardC;
		Two_Mismatches_At_End=Two_Mismatches_At_EndC;
		Mismatches_Backward=Mismatches_BackwardC;
		Possible_20=Possible_20C;
		Possible_02=Possible_02C;
		Right_Mishits=Right_MishitsC;
		Left_Mishits=Left_MishitsC;

		Mismatches_Forward_PointerF=Mismatches_Forward_Pointer;
		Mismatches_Backward_PointerF=Mismatches_Backward_Pointer;
		Two_Mismatches_At_End_Forward_PointerF=Two_Mismatches_At_End_Forward_Pointer;
		Two_Mismatches_At_End_PointerF=Two_Mismatches_At_End_Pointer;

		Possible_20_PointerF=Possible_20_Pointer;
		Possible_02_PointerF=Possible_02_Pointer;
		Possible_03_PointerF=Possible_03_Pointer;
		Possible_30_PointerF=Possible_30_Pointer;
		Possible_04_PointerF=Possible_04_Pointer;
		Possible_40_PointerF=Possible_40_Pointer;
		Right_Mishits_PointerF=Right_Mishits_Pointer;
		Left_Mishits_PointerF=Left_Mishits_Pointer;

		Mismatches_Forward_Pointer=Mismatches_Forward_PointerC;
		Mismatches_Backward_Pointer=Mismatches_Backward_PointerC;
		Two_Mismatches_At_End_Forward_Pointer=Two_Mismatches_At_End_Forward_PointerC;
		Two_Mismatches_At_End_Pointer=Two_Mismatches_At_End_PointerC;

		Possible_20_Pointer=Possible_20_PointerC;
		Possible_02_Pointer=Possible_02_PointerC;
		Possible_03_Pointer=Possible_03_PointerC;
		Possible_30_Pointer=Possible_30_PointerC;
		Possible_04_Pointer=Possible_04_PointerC;
		Possible_40_Pointer=Possible_40_PointerC;
		Right_Mishits_Pointer=Right_Mishits_PointerC;
		Left_Mishits_Pointer=Left_Mishits_PointerC;

	}
	Mismatches_Forward_PointerC=Mismatches_Forward_Pointer;
	Mismatches_Backward_PointerC=Mismatches_Backward_Pointer;
	Two_Mismatches_At_End_Forward_PointerC=Two_Mismatches_At_End_Forward_Pointer;
	Two_Mismatches_At_End_PointerC=Two_Mismatches_At_End_Pointer;

	Possible_20_PointerC=Possible_20_Pointer;
	Possible_02_PointerC=Possible_02_Pointer;
	Possible_03_PointerC=Possible_03_Pointer;
	Possible_30_PointerC=Possible_30_Pointer;
	Possible_04_PointerC=Possible_04_Pointer;
	Possible_40_PointerC=Possible_40_Pointer;
	Right_Mishits_PointerC=Right_Mishits_Pointer;
	Left_Mishits_PointerC=Left_Mishits_Pointer;

	if (FILTERUNIQUEHITS && Hits) return;
	if (!NORMAL_TAGS && !Best_Hit_Found && Hits && ORG_MISMATCHES !=3) {MAX_MISMATCHES=4;Best_Hit_Found=TRUE;} //Go opthit + 1
	if (NPOLICY && NCount){if (!NISMISMATCH && MAX_MISMATCHES + NCount == 3) return;else if (MAX_MISMATCHES==NCount) return;}
	else if( MAX_MISMATCHES ==3 && !INDELSCAN) return;
	if (LEAST_MISMATCH && Hits) return;
//}------------------------------------------- THREE MISMATCH ---------------------------------------------------------------------------------------------
//{------------------------------------------- FOUR MISMATCH ---------------------------------------------------------------------------------------------

	In_Mismatch=4;
	ZigZag=FORWARD;
	Current_Tag=Guessed;//Original+IGNOREHEAD;
	Exact_Match_Forward=Exact_Match_ForwardF;
	Mismatches_Forward=Mismatches_ForwardF;
	Mismatches_Backward=Mismatches_BackwardF;
	Two_Mismatches_At_End_Forward=Two_Mismatches_At_End_ForwardF;
	Two_Mismatches_At_End=Two_Mismatches_At_EndF;
	Possible_20=Possible_20F;
	Possible_02=Possible_02F;
	Right_Mishits=Right_MishitsF;
	Left_Mishits=Left_MishitsF;

	Mismatches_Forward_Pointer=Mismatches_Forward_PointerF;
	Mismatches_Backward_Pointer=Mismatches_Backward_PointerF;
	Two_Mismatches_At_End_Pointer=Two_Mismatches_At_End_PointerF;
	Two_Mismatches_At_End_Forward_Pointer=Two_Mismatches_At_End_Forward_PointerF;

	Possible_20_Pointer=Possible_20_PointerF;
	Possible_02_Pointer=Possible_02_PointerF;
	Possible_03_Pointer=Possible_03_PointerF;//redundant?
	Possible_30_Pointer=Possible_30_PointerF;
	Possible_04_Pointer=Possible_04_PointerF;
	Possible_40_Pointer=Possible_40_PointerF;
	Possible_05_Pointer=Possible_05_PointerF;
	Possible_50_Pointer=Possible_50_PointerF;
	Right_Mishits_Pointer=Right_Mishits_PointerF;
	Mismatches_Forward_Pointer_Last4=Mismatches_Forward_Pointer_Last4F;
	Mismatches_Forward_Pointer_Last5=Mismatches_Forward_Pointer_Last5F;
	Left_Mishits_Pointer_1=Left_Mishits_Pointer_1F;
	Left_Mishits_Pointer=Left_Mishits_PointerF;

	for(;;)
	{

		In_Mismatch=4;
		FMIndex=REVERSE;
		if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
		{
			for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
			{
				Print_LocationX(Two_Mismatches_At_End_Forward[i]);//mismatches of the form 0|4, with last mismatch at the end...
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Two_Mismatches_At_End_Forward_Pointer=0;
		FMIndex=FORWARD;
		if(Two_Mismatches_At_End_Pointer)
		{
			for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
			{
				Print_LocationX(Two_Mismatches_At_End[i]);//mismatches of the form 0|4, with one mismatch at the start...
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}
		Two_Mismatches_At_End_Pointer=0;

		FMIndex=REVERSE;
		Possible_05_Pointer=Mismatches_Forward_Pointer;
		int Wrap=FALSE;
		if(Mismatches_Forward_Pointer)
		{
			if (Possible_04_Pointer > 46000) {Mismatches_Forward_Pointer=0;Wrap=TRUE;}
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=Possible_04_Pointer;i<Possible_05_Pointer;i++)//Mismatches_Forward_Pointer;i++)
				{
					Search_Forwards(Mismatches_Forward[i],4,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|4, and finds mismatches of the form 1|3, stores possibles in the form 1|4
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=Possible_04_Pointer;i<Possible_05_Pointer;i++)//Mismatches_Forward_Pointer;i++)
			{
				Search_Forwards(Mismatches_Forward[i],4,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|4, and finds mismatches of the form 1|3, stores possibles in the form 1|4
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}
		if(Wrap) Possible_05_Pointer=0;
		Mismatches_Forward_Pointer_Last4=Mismatches_Forward_Pointer;

		FMIndex=FORWARD;
		Possible_50_Pointer=Mismatches_Backward_Pointer;
		if(Mismatches_Backward_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=Possible_50_Pointer-1;i>=Possible_40_Pointer;i--)//Mismatches_Backward_Pointer-1;i>=0;i--)
				{
					Search_Backwards(Mismatches_Backward[i],4,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=Possible_50_Pointer-1;i>=Possible_40_Pointer;i--)//Mismatches_Backward_Pointer-1;i>=0;i--)
			{
				Search_Backwards(Mismatches_Backward[i],4,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;

		}

		FMIndex=REVERSE;
		Left_Mishits_Pointer_1=Left_Mishits_Pointer;//Polish some more....
		if(Left_Mishits_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=0;i<Left_Mishits_Pointer;i++)
				{
					TRange=Left_Mishits[i];
					if (Left_Mishits[i].Level != LH+1)
					{
						Search_Exact(Left_Mishits[i],-1,LH,revfmi);
						Left_Mishits[i].Level=LH+1;//search exact does not update level..
					}
					if (Left_Mishits[i].Start)
					{
						Search_Forwards(Left_Mishits[i],4,1,STRINGLENGTH,revfmi);//find mismatches of the form 022 form, stores possibles of the form 023
					}
					Left_Mishits[i]=TRange;
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=0;i<Left_Mishits_Pointer;i++)
			{
				if (Left_Mishits[i].Level != LH+1)
				{
					Search_Exact(Left_Mishits[i],-1,LH,revfmi);
					Left_Mishits[i].Level=LH+1;//search exact does not update level..
				}
				if (Left_Mishits[i].Start)
				{
					Search_Forwards(Left_Mishits[i],4,1,STRINGLENGTH,revfmi);//find mismatches of the form 022 form, stores possibles of the form 023
				}
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;

		}


		Mismatches_Forward_Pointer_Last5=Mismatches_Forward_Pointer;
		if( Right_Mishits_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=0;i<Right_Mishits_Pointer;i++)
				{
					TRange=Right_Mishits[i];
					if(Right_Mishits[i].Level!=LHQL) 
					{
						Search_Backwards_Exact( Right_Mishits[i],LHQL,LHQL,fwfmi);//finds mismatches of the form 202, stores possibles of the form 203
					}
					if(Right_Mishits[i].Start)
					{	
						Backwards(Right_Mishits[i],1,LH);
						if(Right_Mishits[i].Start)
						{
							Right_Mishits[i].Level=1;
							Search_Forwards(Right_Mishits[i],4,LH+1,RH,revfmi);
							if(MAXHITS==Hits) break;
						}
					}
					Right_Mishits[i]=TRange;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=0;i<Right_Mishits_Pointer;i++)
			{

				if(Right_Mishits[i].Level!=LHQL) 
				{
					Search_Backwards_Exact( Right_Mishits[i],LHQL,LHQL,fwfmi);//finds mismatches of the form 202, stores possibles of the form 203
				}
				if(Right_Mishits[i].Start)
				{	
					Backwards(Right_Mishits[i],1,LH);
					if(Right_Mishits[i].Start)
					{
						Right_Mishits[i].Level=1;
						Search_Forwards(Right_Mishits[i],4,LH+1,RH,revfmi);
						if(MAXHITS==Hits) break;
					}
				}
			}
			if(MAXHITS==Hits) return;
		}

		FMIndex=REVERSE;
		int LHQLrx=LHQL/2;
		if (LHQL % 2) LHQLrx++; int LHQRrx=LHQL-LHQLrx;

		if (LOOKUPSIZE >= LHQLrx)
		{
			Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		}
		else
		{
			if(LOOKUPSIZE==3)
			{
				c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
			}
			else
			{
				c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
			}
			Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0; Range.Skip=0;
		}
		Search_Forwards_0X(Range,1,LHQLrx,revfmi);
		Range.Level=1;
		if(USEQUALITY)
		{
			TRange=Range;
			Do_Branch=Low_Quality;
			Search_01LX(Range,1,LHQLrx +1,LHQRrx,revfmi);
			if(MAXHITS==Hits) return;
			Range=TRange;
		}

		Do_Branch=Do_All;
		Search_01LX(Range,1,LHQLrx +1,LHQRrx,revfmi);
		if(MAXHITS==Hits) return;
		//--------------------------------
		if (LOOKUPSIZE >= LHQRrx)
		{
			Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		}
		else
		{
			if(LOOKUPSIZE==3)
			{
				c=Current_Tag[LHQL-1-0] | (Current_Tag[LHQL-1-1]<<2) | (Current_Tag[LHQL-1-2]<<4);// | (Current_Tag[LH-1-3]<<6) | Current_Tag[LH-1-4]<<8 | (Current_Tag[LH-1-5]<<10);//Use lookup table...
			}
			else
			{
				c=Current_Tag[LHQL-1-0] | (Current_Tag[LHQL-1-1]<<2) | (Current_Tag[LHQL-1-2]<<4) | (Current_Tag[LHQL-1-3]<<6) | Current_Tag[LHQL-1-4]<<8 | (Current_Tag[LHQL-1-5]<<10);//Use lookup table...
			}

			Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;//Range.Tag=Actual_Tag;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		}
		Search_Backwards_Exact_X0( Range,LHQL,LHQRrx,fwfmi);// ?|0|?
		Range.Level=1;
		if(USEQUALITY)
		{
			TRange=Range;
			Do_Branch=Low_Quality;
			Search_10LX(TRange,1,LHQL, LHQL,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3 
			if(MAXHITS==Hits) return;
		}
		Do_Branch=Do_All;
		Search_10LX(Range,1,LHQLrx, LHQLrx,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
		if(MAXHITS==Hits) return;
		if(ZigZag) break;
		ZigZag++;
		Current_Tag=Guess_Complement;
		Exact_Match_Forward=Exact_Match_ForwardC;
		Mismatches_Forward=Mismatches_ForwardC;
		Two_Mismatches_At_End_Forward=Two_Mismatches_At_End_ForwardC;
		Two_Mismatches_At_End=Two_Mismatches_At_EndC;
		Mismatches_Backward=Mismatches_BackwardC;
		Possible_20=Possible_20C;
		Possible_02=Possible_02C;
		Right_Mishits=Right_MishitsC;
		Left_Mishits=Left_MishitsC;

		Mismatches_Forward_PointerF=Mismatches_Forward_Pointer;
		Mismatches_Backward_PointerF=Mismatches_Backward_Pointer;
		Two_Mismatches_At_End_Forward_PointerF=Two_Mismatches_At_End_Forward_Pointer;
		Two_Mismatches_At_End_PointerF=Two_Mismatches_At_End_Pointer;

		Possible_20_PointerF=Possible_20_Pointer;
		Possible_02_PointerF=Possible_02_Pointer;
		Possible_03_PointerF=Possible_03_Pointer;
		Possible_30_PointerF=Possible_30_Pointer;
		Possible_04_PointerF=Possible_04_Pointer;
		Possible_40_PointerF=Possible_40_Pointer;
		Possible_05_PointerF=Possible_05_Pointer;
		Possible_50_PointerF=Possible_50_Pointer;
		Right_Mishits_PointerF=Right_Mishits_Pointer;
		Mismatches_Forward_Pointer_Last4F=Mismatches_Forward_Pointer_Last4;
		Mismatches_Forward_Pointer_Last5F=Mismatches_Forward_Pointer_Last5;
		Left_Mishits_Pointer_1F=Left_Mishits_Pointer_1;
		Left_Mishits_PointerF=Left_Mishits_Pointer;

		Mismatches_Forward_Pointer=Mismatches_Forward_PointerC;
		Mismatches_Backward_Pointer=Mismatches_Backward_PointerC;
		Two_Mismatches_At_End_Forward_Pointer=Two_Mismatches_At_End_Forward_PointerC;
		Two_Mismatches_At_End_Pointer=Two_Mismatches_At_End_PointerC;

		Possible_20_Pointer=Possible_20_PointerC;
		Possible_02_Pointer=Possible_02_PointerC;
		Possible_03_Pointer=Possible_03_PointerC;
		Possible_30_Pointer=Possible_30_PointerC;
		Possible_04_Pointer=Possible_04_PointerC;
		Possible_40_Pointer=Possible_40_PointerC;
		Right_Mishits_Pointer=Right_Mishits_PointerC;
		Possible_05_Pointer=Possible_05_PointerC;
		Possible_50_Pointer=Possible_50_PointerC;
		Mismatches_Forward_Pointer_Last4=Mismatches_Forward_Pointer_Last4C;
		Mismatches_Forward_Pointer_Last5=Mismatches_Forward_Pointer_Last5C;
		Left_Mishits_Pointer_1=Left_Mishits_Pointer_1C;
		Left_Mishits_Pointer=Left_Mishits_PointerC;

	}
	Mismatches_Forward_PointerC=Mismatches_Forward_Pointer;
	Mismatches_Backward_PointerC=Mismatches_Backward_Pointer;
	Two_Mismatches_At_End_Forward_PointerC=Two_Mismatches_At_End_Forward_Pointer;
	Two_Mismatches_At_End_PointerC=Two_Mismatches_At_End_Pointer;

	Possible_20_PointerC=Possible_20_Pointer;
	Possible_02_PointerC=Possible_02_Pointer;
	Possible_03_PointerC=Possible_03_Pointer;
	Possible_30_PointerC=Possible_30_Pointer;
	Possible_04_PointerC=Possible_04_Pointer;
	Possible_40_PointerC=Possible_40_Pointer;
	Right_Mishits_PointerC=Right_Mishits_Pointer;
	Possible_05_PointerC=Possible_05_Pointer;
	Possible_50_PointerC=Possible_50_Pointer;
	Mismatches_Forward_Pointer_Last4C=Mismatches_Forward_Pointer_Last4;
	Mismatches_Forward_Pointer_Last5C=Mismatches_Forward_Pointer_Last5;
	Left_Mishits_Pointer_1C=Left_Mishits_Pointer_1;
	Left_Mishits_PointerC=Left_Mishits_Pointer;

	if (FILTERUNIQUEHITS && Hits) return;
	if (NPOLICY && NCount){if (!NISMISMATCH && MAX_MISMATCHES + NCount == 4) return;else if (MAX_MISMATCHES==NCount) return;}
	else if( MAX_MISMATCHES ==4) return;
	if (LEAST_MISMATCH && Hits) return;
//}------------------------------------------- FOUR MISMATCH ---------------------------------------------------------------------------------------------
//{------------------------------------------- FIVE MISMATCH ---------------------------------------------------------------------------------------------

	In_Mismatch=5;
	ZigZag=FORWARD;
	Current_Tag=Guessed;//Original+IGNOREHEAD;
	Exact_Match_Forward=Exact_Match_ForwardF;
	Mismatches_Forward=Mismatches_ForwardF;
	Mismatches_Backward=Mismatches_BackwardF;
	Two_Mismatches_At_End_Forward=Two_Mismatches_At_End_ForwardF;
	Two_Mismatches_At_End=Two_Mismatches_At_EndF;
	Possible_20=Possible_20F;
	Possible_02=Possible_02F;
	Right_Mishits=Right_MishitsF;
	Left_Mishits=Left_MishitsF;

	Mismatches_Forward_Pointer=Mismatches_Forward_PointerF;
	Mismatches_Backward_Pointer=Mismatches_Backward_PointerF;
	Two_Mismatches_At_End_Pointer=Two_Mismatches_At_End_PointerF;
	Two_Mismatches_At_End_Forward_Pointer=Two_Mismatches_At_End_Forward_PointerF;

	Possible_20_Pointer=Possible_20_PointerF;
	Possible_02_Pointer=Possible_02_PointerF;
	Possible_04_Pointer=Possible_04_PointerF;
	Possible_40_Pointer=Possible_40_PointerF;
	Right_Mishits_Pointer=Right_Mishits_PointerF;
	Left_Mishits_Pointer=Left_Mishits_PointerF;
	Possible_05_Pointer=Possible_05_PointerF;
	Possible_50_Pointer=Possible_50_PointerF;
	Mismatches_Forward_Pointer_Last4=Mismatches_Forward_Pointer_Last4F;
	Mismatches_Forward_Pointer_Last5=Mismatches_Forward_Pointer_Last5F;
	Left_Mishits_Pointer_1=Left_Mishits_Pointer_1F;

	for(;;)
	{
		FMIndex=REVERSE;
		if(Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
		{
			for(int i=0;i<Two_Mismatches_At_End_Forward_Pointer;i++)
			{
				Print_LocationX(Two_Mismatches_At_End_Forward[i]);//mismatches of the form 0|4, with last mismatch at the end...
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}
		Two_Mismatches_At_End_Forward_Pointer=0;


		FMIndex=FORWARD;
		if(Two_Mismatches_At_End_Pointer)
		{
			for(int i=0;i<Two_Mismatches_At_End_Pointer;i++)
			{
				Print_LocationX(Two_Mismatches_At_End[i]);//mismatches of the form 0|5, with one mismatch at the start...
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}
		Two_Mismatches_At_End_Pointer=0;

		FMIndex=REVERSE;
		if(Mismatches_Forward_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=Possible_05_Pointer;i<Mismatches_Forward_Pointer_Last4;i++)//Mismatches_Forward_Pointer;i++)
				{
					Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=Possible_05_Pointer;i<Mismatches_Forward_Pointer_Last4;i++)//Mismatches_Forward_Pointer;i++)
			{
				Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;

		}


		FMIndex=REVERSE;
		if(Mismatches_Forward_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=Mismatches_Forward_Pointer_Last4;i<Mismatches_Forward_Pointer_Last5;i++)//Mismatches_Forward_Pointer;i++)
				{
					Search_Forwards(Mismatches_Forward[i],5,1,STRINGLENGTH,revfmi);//scan for possible five mismatches of the form 0|5, and finds mismatches of the form 1|4,2|3 
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=Mismatches_Forward_Pointer_Last4;i<Mismatches_Forward_Pointer_Last5;i++)//Mismatches_Forward_Pointer;i++)
			{
				Search_Forwards(Mismatches_Forward[i],5,1,STRINGLENGTH,revfmi);//scan for possible five mismatches of the form 0|5, and finds mismatches of the form 1|4,2|3 
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		FMIndex=REVERSE;
		if(Mismatches_Forward_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=Mismatches_Forward_Pointer_Last5;i<Mismatches_Forward_Pointer;i++)//Mismatches_Forward_Pointer;i++)
				{
					Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=Mismatches_Forward_Pointer_Last5;i<Mismatches_Forward_Pointer;i++)//Mismatches_Forward_Pointer;i++)
			{
				Search_Forwards(Mismatches_Forward[i],5,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|5
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		FMIndex=FORWARD;
		if(Mismatches_Backward_Pointer!=Possible_50_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=Mismatches_Backward_Pointer-1;i>=Possible_50_Pointer;i--)
				{
					Search_Backwards(Mismatches_Backward[i],5,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
					if(MAXHITS==Hits) break;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=Mismatches_Backward_Pointer-1;i>=Possible_50_Pointer;i--)
			{
				Search_Backwards(Mismatches_Backward[i],5,LH,LH,fwfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;

		}

		FMIndex=FORWARD;
		if (Possible_20_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=Possible_20_Pointer-1;i>=0;i--)
				{
					if(Possible_20[i].Level!=RHQL) 
					{
						Possible_20[i].Level++;
						Search_Backwards_Exact( Possible_20[i],LH+RHQL,RHQL,fwfmi);
					}

					if(Possible_20[i].Start)
					{	
						Possible_20[i].Level=1;
						Search_Backwards(Possible_20[i],5,LH,LH,fwfmi);
						if(MAXHITS==Hits) break;
					}
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=Possible_20_Pointer-1;i>=0;i--)
			{
				if(Possible_20[i].Level!=RHQL) 
				{
					Possible_20[i].Level++; 
					Search_Backwards_Exact( Possible_20[i],LH+RHQL,RHQL,fwfmi);
				}

				if(Possible_20[i].Start)
				{	
					Possible_20[i].Level=1;
					Search_Backwards(Possible_20[i],5,LH,LH,fwfmi);
					if(MAXHITS==Hits) break;
				}
			}
			if(MAXHITS==Hits) return;

		}

		if(Possible_02_Pointer)
		{
			if(USEQUALITY)
			{
				Do_Branch=Low_Quality;
				for(int i=0;i<Possible_02_Pointer;i++)
				{
					TRange=Possible_02[i];
					if(Possible_02[i].Level!=RHQR) 
					{
						Search_Exact( Possible_02[i],LH + RHQL-1 ,RHQR,revfmi);//finds mismatches of the form 202, stores possibles of the form 203
					}
					if(Possible_02[i].Start) 
					{
						Reverse(Possible_02[i],STRINGLENGTH,RH);
						if(Possible_02[i].Start)
						{
							Possible_02[i].Level=1;
							Search_Backwards(Possible_02[i],5,LH,LH,fwfmi);
						}
					}
					if(MAXHITS==Hits) break;
					Possible_02[i]=TRange;
				}
				if(MAXHITS==Hits) return;
			}

			Do_Branch=Do_All;
			for(int i=0;i<Possible_02_Pointer;i++)
			{

				if(Possible_02[i].Level!=RHQR) 
				{
					Search_Exact( Possible_02[i],LH + RHQL-1 ,RHQR,revfmi);//finds mismatches of the form 202, stores possibles of the form 203
				}
				if(Possible_02[i].Start) 
				{
					Reverse(Possible_02[i],STRINGLENGTH,RH);
					if(Possible_02[i].Start)
					{
						Possible_02[i].Level=1;
						Search_Backwards(Possible_02[i],5,LH,LH,fwfmi);
					}
				}
				if(MAXHITS==Hits) break;
			}

			if(MAXHITS==Hits) return;
		}

		int RHQLlx=RHQR/2;
		if (RHQR % 2) RHQLlx++; int RHQRlx=RHQR-RHQLlx;
		if (LOOKUPSIZE >= RHQRlx)
		{
			Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		}
		else
		{
			if(LOOKUPSIZE==3)
			{
				c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
			}
			else
			{
				c=Current_Tag[STRINGLENGTH-1-0] | (Current_Tag[STRINGLENGTH-1-1]<<2) | (Current_Tag[STRINGLENGTH-1-2]<<4) | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
			}
			Range.Start=Backward_Start_LookupX[c];Range.End=Backward_End_LookupX[c];
			Range.Mismatches=0;Range.Level=LOOKUPSIZE+1;//Range.Tag=Actual_Tag;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		}

		Search_Backwards_Exact_X0( Range,STRINGLENGTH,RHQRlx,fwfmi);// ?|?|0
		Range.Level=1;

		if(USEQUALITY)
		{
			TRange=Range;
			Do_Branch=Low_Quality;
			Search_Backwards_XL10(Range,1,LH + RHQL+RHQLlx, RHQLlx,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
			if(MAXHITS==Hits) return;
			Range=TRange;
		}

		Do_Branch=Do_All;
		Search_Backwards_XL10(Range,1,LH + RHQL+RHQLlx, RHQLlx,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if(MAXHITS==Hits) return;



		if (LOOKUPSIZE >= RHQLlx)
		{
			Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;Range.Skip=0;
		}
		else
		{
			if(LOOKUPSIZE==3)
			{
				c=Current_Tag[LH+RHQL+0] | (Current_Tag[LH+RHQL+1]<<2) | (Current_Tag[LH+RHQL+2]<<4);// | (Current_Tag[LH+3]<<6) | Current_Tag[LH+4]<<8 | (Current_Tag[LH+5]<<10);//Use lookup table...
			}
			else
			{
				c=Current_Tag[LH+RHQL+0] | (Current_Tag[LH+RHQL+1]<<2) | (Current_Tag[LH+RHQL+2]<<4) | (Current_Tag[LH+RHQL+3]<<6) | Current_Tag[LH+RHQL+4]<<8 | (Current_Tag[LH+RHQL+5]<<10);//Use lookup table...
			}
			Range.Start=Forward_Start_LookupX[c];Range.End=Forward_End_LookupX[c];
			Range.Mismatches=0;Range.Level=LOOKUPSIZE+1; Range.Skip=0;
			memcpy(Range.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Range.Mismatch_Char=0;
		}

		Search_Forwards_0X(Range,LH+RHQL+1,RHQLlx,revfmi);
		Range.Level=1;
		TRange=Range;
		if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			Search_XL01(Range,1,LH + RHQL+RHQLlx +1,RHQRlx,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;
		Search_XL01(TRange,1,LH + RHQL+RHQLlx +1,RHQRlx,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if(MAXHITS==Hits) return;

		if(ZigZag) break;
		ZigZag++;
		Current_Tag=Guess_Complement;
		Exact_Match_Forward=Exact_Match_ForwardC;
		Mismatches_Forward=Mismatches_ForwardC;
		Two_Mismatches_At_End_Forward=Two_Mismatches_At_End_ForwardC;
		Two_Mismatches_At_End=Two_Mismatches_At_EndC;
		Mismatches_Backward=Mismatches_BackwardC;
		Possible_20=Possible_20C;
		Possible_02=Possible_02C;
		Right_Mishits=Right_MishitsC;
		Left_Mishits=Left_MishitsC;

		Mismatches_Forward_PointerF=Mismatches_Forward_Pointer;
		Mismatches_Backward_PointerF=Mismatches_Backward_Pointer;
		Two_Mismatches_At_End_Forward_PointerF=Two_Mismatches_At_End_Forward_Pointer;
		Two_Mismatches_At_End_PointerF=Two_Mismatches_At_End_Pointer;

		Possible_20_PointerF=Possible_20_Pointer;
		Possible_02_PointerF=Possible_02_Pointer;
		Possible_03_PointerF=Possible_03_Pointer;
		Possible_30_PointerF=Possible_30_Pointer;
		Possible_04_PointerF=Possible_04_Pointer;
		Possible_40_PointerF=Possible_40_Pointer;
		Possible_05_PointerF=Possible_05_Pointer;
		Possible_50_PointerF=Possible_50_Pointer;
		Right_Mishits_PointerF=Right_Mishits_Pointer;
		Left_Mishits_PointerF=Left_Mishits_Pointer;
		Mismatches_Forward_Pointer_Last4F=Mismatches_Forward_Pointer_Last4;
		Mismatches_Forward_Pointer_Last5F=Mismatches_Forward_Pointer_Last5;
		Left_Mishits_Pointer_1F=Left_Mishits_Pointer_1;

		Mismatches_Forward_Pointer=Mismatches_Forward_PointerC;
		Mismatches_Backward_Pointer=Mismatches_Backward_PointerC;
		Two_Mismatches_At_End_Forward_Pointer=Two_Mismatches_At_End_Forward_PointerC;
		Two_Mismatches_At_End_Pointer=Two_Mismatches_At_End_PointerC;

		Possible_20_Pointer=Possible_20_PointerC;
		Possible_02_Pointer=Possible_02_PointerC;
		Possible_03_Pointer=Possible_03_PointerC;
		Possible_30_Pointer=Possible_30_PointerC;
		Possible_04_Pointer=Possible_04_PointerC;
		Possible_40_Pointer=Possible_40_PointerC;
		Possible_05_Pointer=Possible_05_PointerC;
		Possible_50_Pointer=Possible_50_PointerC;
		Right_Mishits_Pointer=Right_Mishits_PointerC;
		Left_Mishits_Pointer=Left_Mishits_PointerC;
		Mismatches_Forward_Pointer_Last4=Mismatches_Forward_Pointer_Last4C;
		Mismatches_Forward_Pointer_Last5=Mismatches_Forward_Pointer_Last5C;
		Left_Mishits_Pointer_1=Left_Mishits_Pointer_1C;

	}
	Mismatches_Forward_PointerC=Mismatches_Forward_Pointer;
	Mismatches_Backward_PointerC=Mismatches_Backward_Pointer;
	Two_Mismatches_At_End_Forward_PointerC=Two_Mismatches_At_End_Forward_Pointer;
	Two_Mismatches_At_End_PointerC=Two_Mismatches_At_End_Pointer;

	Possible_20_PointerC=Possible_20_Pointer;
	Possible_02_PointerC=Possible_02_Pointer;
	Possible_03_PointerC=Possible_03_Pointer;
	Possible_30_PointerC=Possible_30_Pointer;
	Possible_04_PointerC=Possible_04_Pointer;
	Possible_40_PointerC=Possible_40_Pointer;
	Possible_05_PointerC=Possible_05_Pointer;
	Possible_50_PointerC=Possible_50_Pointer;
	Right_Mishits_PointerC=Right_Mishits_Pointer;
	Left_Mishits_PointerC=Left_Mishits_Pointer;
	Mismatches_Forward_Pointer_Last4C=Mismatches_Forward_Pointer_Last4;
	Mismatches_Forward_Pointer_Last5C=Mismatches_Forward_Pointer_Last5;
	Left_Mishits_Pointer_1C=Left_Mishits_Pointer_1;
	if (FILTERUNIQUEHITS && Hits) return;
	if( MAX_MISMATCHES ==5) return;
	if (LEAST_MISMATCH && Hits) return;
//}------------------------------------------- FIVE MISMATCH ---------------------------------------------------------------------------------------------

//{------------------------------------------- MULTI ZZ MISMATCH ---------------------------------------------------------------------------------------------
	if (NPOLICY && NCount){if (!NISMISMATCH && MAX_MISMATCHES + NCount == 5) return;else if (MAX_MISMATCHES==NCount) return;}
	if (MAX_MISMATCHES >5 && !Hits)
	{
		MAX_ORG=MAX_MISMATCHES;
		char FILTERUNIQUEHITS_T=FILTERUNIQUEHITS;
		int AB=ARRAY_BOUND;int EB=END_BOUND;
		if (ARRAY_BOUND_BD) {ARRAY_BOUND=ARRAY_BOUND_BD;END_BOUND=ARRAY_BOUND_BD;}
		unsigned TTotal_Hits=Total_Hits,TTotal_Tags=Total_Tags;
		S=SATot=0;
		In_Mismatch=6;
		Ext_Scan=0;
		Current_Tag=Guessed;
		Low_Quality=Guessed_Quality;
		Extentions=0;
		if (MAXSPECIFIED)
		{
			Multi_Hit.clear();
		}
		HIGHSCAN=TRUE;
		Last_Mismatch_Written=0;
		int H=MAXHITS;if(!SUPERACCURATE) MAXHITS=INT_MAX;//10000;
		if (LEAST_MISMATCH) MAXHITS=H;
		Current_Tag=Guessed;
		Low_Quality=Guessed_Quality;
		char L;
		if (Low_Quality[STRINGLENGTH]) L=TRUE; else L=FALSE;
		while(Hits < MAXHITS && Ext_Scan<2) 
		{

			Ext_Scan++;
			if (Ext_Scan==1)
			{
				if (L) 
				{
					Current_Tag=Guessed;
					Low_Quality=Guessed_Quality;
					Left_To_Right();
					if(Hits == MAXHITS && !FILTERUNIQUEHITS) break;
					Current_Tag=Guess_Complement;
					Low_Quality=Guessed_Complement_Quality;
					Right_To_Left();
				}
				else
				{
					Current_Tag=Guessed;
					Low_Quality=Guessed_Quality;
					Right_To_Left();
					if(Hits==MAXHITS && !FILTERUNIQUEHITS) break;
					Current_Tag=Guess_Complement;
					Low_Quality=Guessed_Complement_Quality;
					Left_To_Right();
				}
			}
			else	
			{
				if (L) 
				{
					Current_Tag=Guessed;
					Low_Quality=Guessed_Quality;
					Right_To_Left();
					if(Hits == MAXHITS && !FILTERUNIQUEHITS) break;
					Current_Tag=Guess_Complement;
					Low_Quality=Guessed_Complement_Quality;
					Left_To_Right();
				}
				else
				{
					Current_Tag=Guessed;
					Low_Quality=Guessed_Quality;
					Left_To_Right();
					if(Hits==MAXHITS && !FILTERUNIQUEHITS) break;
					Current_Tag=Guess_Complement;
					Low_Quality=Guessed_Complement_Quality;
					Right_To_Left();
				}
			}


			//Low_Quality=Guessed_Complement_Quality;
		}

		if (Tag_Stat_Bad)
		{
			Total_Hits=TTotal_Hits;Total_Tags=TTotal_Tags;
		}

		Ext_Scan=0;
		//int Ext_Scan=0;
		Current_Tag=Guessed;
		Low_Quality=Guessed_Quality;
		Larger_Than_Ten = TRUE;
		if (!Hits || (Last_Mismatch_Written >12))
		{
			if (!Tag_Stat_Bad && MAX_MISMATCHES >10)
			{
				while(Hits < MAXHITS && Ext_Scan<2) 
				{
					Ext_Scan++;

					Right_To_LeftX();
					if(Hits==MAXHITS) break;
					Left_To_RightX();
					if(Hits==MAXHITS) break;
					Left_To_RightY();
					if(Hits==MAXHITS) break;
					Right_To_LeftY();

					Current_Tag=Guess_Complement;
					Low_Quality=Guessed_Complement_Quality;
				}
			}
		}
		Larger_Than_Ten = FALSE;
		if(ARRAY_BOUND_BD) {ARRAY_BOUND=AB;END_BOUND=EB;}
		MAXHITS=H;MAX_MISMATCHES=MAX_ORG;
		FILTERUNIQUEHITS=FILTERUNIQUEHITS_T;
		if(PRINT_MISHITS)
		{
			if (Tag_Stat_Bad) fprintf(Mishit_File,"%s%sBAD%s%s", Description,Tag_Copy,Plus,Quality);
			else if(!Hits) fprintf(Mishit_File,"%s%s%s%s", Description,Tag_Copy,Plus,Quality);
		}
	}

//}------------------------------------------- MULTI ZZ MISMATCH ---------------------------------------------------------------------------------------------

}
//}-----------------------------  Zig_Zag  -------------------------------------------------


void Search_Brute(int Count1,int Count2,int Start,int StringLength,BWT *fmi)
{
	FMIndex=FORWARD;
	SARange Tag;
	Tag.Start=1;Tag.End=SOURCELENGTH;Tag.Mismatches=0;Tag.Level=1;Tag.Skip=0;
	memcpy(Tag.Mismatch_Pos,Mismatch_Init,MAX_MISMATCHES_BOUND);Tag.Mismatch_Char=0;

	if (!Tag.Start) return;
	int BMStack_Top=0;
	BMStack_X11[0]=Tag;
	struct SARange Range,Temp_Range;
	while(BMStack_Top!=-1)//While Stack non-empty....
	{
		Range=BMStack_X11[BMStack_Top];
		BMStack_Top--;	//Pop the range
		Branch_Detect_Backwards(Range,fmi,Start);

		for(int Branch=0;Branch<4;Branch++)
		{
			if (Branch_Characters[Branch])//This character actually branches
			{
				Temp_Range=Range;//adjust
				Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
				Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

				if (Current_Tag[Start-Temp_Range.Level] != Branch)
				{
					Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
					Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
					Temp_Range.Mismatches++;
				}

				if (Temp_Range.Mismatches<=Count2)//we are guaranteed a valid SA range, check only for mismatches
				{
					if(Temp_Range.Level== StringLength)
					{
						if(Temp_Range.Mismatches>=Count1)
						{
							Temp_Range.Level=1;
							Temp_BC2[0]=Branch_Characters[0];Temp_BC2[1]=Branch_Characters[1];Temp_BC2[2]=Branch_Characters[2];Temp_BC2[3]=Branch_Characters[3];
							memcpy(Temp_Branch_Ranges2,Branch_Ranges,4*sizeof(SARange));//X|8
							Search_BackwardsX(Temp_Range,MAX_MISMATCHES,LH,LH,fwfmi);//LH,fwfmi);
							if(MAXHITS==Hits) return;
							memcpy(Branch_Ranges,Temp_Branch_Ranges2,4*sizeof(SARange));
							Branch_Characters[0]=Temp_BC2[0];Branch_Characters[1]=Temp_BC2[1];Branch_Characters[2]=Temp_BC2[2];Branch_Characters[3]=Temp_BC2[3];

						}
						else continue;
					}
					else
					{
						BMStack_Top++;//Push range
						Temp_Range.Level++;
						BMStack_X11[BMStack_Top]=Temp_Range;
					}
				}
			} 
		}
	}
	return;
}

void Sig_Handle(int Sig)
{
	if (Sig == 11 )
		printf("\nOh No... a Segmantation Fault occured : Guess There is a Bug in the code :-( \n");
	printf("Faulted in\nTag : %d\nDescription: %sRead: %s\nReport this info with the command line + genome used and maybe I can fix it...\n",Actual_Tag,Description,Tag_Copy);

	signal(Sig, SIG_DFL);
} 

//{--------------------------------------- DUST --------------------------------------------------------------------
int dust(FastaSeq* fa) 
{
	int level=28;
      int i, j, l, from, to, a, b, v;
      from = 0;
      to = -1;
      a=0;b=0;

      for (i=0; i < fa->len; i++) Dustiness[i]=0;
      for (i=0; i < fa->len; i += window2) 
      {
            from -= window2;
            to -= window2;
            l = (fa->len > i+window) ? window : fa->len-i;
            v = wo(l, fa->seq+i, &a, &b);
               /* return coordinates from-to*/
               /*fprintf(stderr, "%s : %d - %d \n", fa->id, i+from, i+to); */
	    addDustRgn(fa, i+from, i+to);

	    if (v > level) 
	    {
		    /* return coordinates from a to min(b,window2)*/
		    /* fprintf(stderr, "%s : %d - %d \n", fa->id, 
		       i+a, i+((window2>b)? b : window2)); */
		    addDustRgn(fa, i+a, i + ((window2>b)? b : window2));
		    from = (window2>b)? b : window2;
		    to = b;
	    }
	    else 
	    {
		    from = 0;
		    to = -1;
	    }
      }

      int dust=0;
      for (i=0; i < fa->len; i++) {if(Dustiness[i]) dust++;}
      return dust;
}

inline void addDustRgn(FastaSeq* fa, int f, int t) 
{  
	/* fprintf(stderr, "addDustRgn(, %d, %d)\n", f, t); */
	if (t<0 || f>=t) return; /* not a maskable range */
	for ( int i=f;i<t;i++) Dustiness[i]=1;
	return; 
}

void wo1(int len, char* s, int ivv) 
{
	int i, ii, j, v, t, n, n1, sum;
	static int counts[32*32*32];
	static int iis[32*32*32];
	int js, nis;

	n = 32 * 32 * 32;
	n1 = n - 1;
	nis = 0;
	i = 0;
	ii = 0;
	sum = 0;
	v = 0;
	for (j=0; j < len; j++, s++) 
	{
		ii <<= 5;
		if (isalpha(*s)) 
		{
			if (islower(*s)) {
				ii |= *s - 'a';
			} else {
				ii |= *s - 'A';
			}
		} else {
			i = 0;
			continue;
		}
		ii &= n1;
		i++;
		if (i >= word) 
		{
			for (js=0; js < nis && iis[js] != ii; js++) ;
			if (js == nis) {
				iis[nis] = ii;
				counts[ii] = 0;
				nis++;
			}
			if ((t = counts[ii]) > 0) 
			{
				sum += t;
				v = 10 * sum / j;
				if (mv < v) 
				{
					mv = v;
					iv = ivv;
					jv = j;
				}
			}
			counts[ii]++;
		}
	}
}

int wo(int len, char* s, int* beg, int* end) 
{
	int i, l1;

	l1 = len - word + 1;
	if (l1 < 0) 
	{
		*beg = 0;
		*end = len - 1;
		return 0;
	}
	mv = 0;
	iv = 0;
	jv = 0;
	for (i=0; i < l1; i++) 
	{
		wo1(len-i, s+i, i);
	}
	*beg = iv;
	*end = iv + jv;
	return mv;
}
//}--------------------------------------- DUST --------------------------------------------------------------------

void Quality_Calc()
{
	Minimum_Quality=Quality[0];
	Maximum_Quality=Quality[0];
	memcpy(Quality_Count,All_Zero,256);
	int LSum=0;
	for (int i=0;i<LH;i++)//filter high/low quality scores
	{
		LSum+=Quality[i];
		Quality_Count[Quality[i]]++;
		if (Minimum_Quality>Quality[i]) { Minimum_Quality=Quality[i];}
		if (Maximum_Quality<Quality[i]) { Maximum_Quality=Quality[i];}
	}

	char Sum=0;
	for(int i=Minimum_Quality;i<Maximum_Quality;i++)
	{
		Sum=Sum+Quality_Count[i];
		if (Sum >= MISMATCHES_TO_TRY_FIRST_LEFT){Quality_Bound=i;break;}
	}

	for (int i=0;i<LH;i++)//filter high/low quality scores
	{
		if(Quality[i] > Quality_Bound) Low_QualityF[i]=FALSE;
		else Low_QualityF[i]=TRUE;
	}
	Maximum_Quality=Quality[LH];

	int RSum=0;
	for (int i=LH;i<STRINGLENGTH;i++)//filter high/low quality scores
	{
		RSum+=Quality[i];
		Quality_Count[Quality[i]]++;
		if (Minimum_Quality>Quality[i]) { Minimum_Quality=Quality[i];}
		if (Maximum_Quality<Quality[i]) { Maximum_Quality=Quality[i];}
	}

	Sum=0;
	for(int i=Minimum_Quality;i<Maximum_Quality;i++)
	{
		Sum=Sum+Quality_Count[i];
		if (Sum >= MISMATCHES_TO_TRY_FIRST_RIGHT) {Quality_Bound=i;break;}
	}

	for (int i=LH;i<STRINGLENGTH;i++)//filter high/low quality scores
	{
		if(Quality[i] > Quality_Bound) Low_QualityF[i]=FALSE;
		else Low_QualityF[i]=TRUE;
	}
	for (int i=0;i<STRINGLENGTH;i++) {Low_QualityC[STRINGLENGTH-i-1]=Low_QualityF[i];}
	Low_Quality=Low_QualityF;
//if Lo_quality[STRINGLENGTH]=1 left has highest quality...
	if(LSum >RSum) {Low_QualityF[STRINGLENGTH]=1;Low_QualityC[STRINGLENGTH]=0;} else {Low_QualityF[STRINGLENGTH]=0;Low_QualityC[STRINGLENGTH]=0;}
}
