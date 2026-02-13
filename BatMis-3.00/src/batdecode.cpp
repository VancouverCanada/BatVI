//{-----------------------------  DEFINES  -------------------------------------------------/

using namespace std;
#define TAB	1
#define FQ	2
#define FA	3
#define TWOFILE	4
#define GISMODE 100
#define BTS_PER_LOC 8
#define BITMASK 255 //2^BTS_PER_LOC -1
#define MAXDES 400
#define DEFAULT 0
#define DEEP 1
#define PAIREND 2
#define MAX_MISMATCHES_BOUND 16 //upper boud for mismatch number....

#define DEBUG 0
#define BUFFERSIZE 1000000

/* 	#define BWTFILE "chr11.bwt"// "genome.bwt"// 
 * 	#define OCCFILE "chr11.fmv"// "genome.fmv"//
 * 	#define SAFILE "chr11.sa"
 * 	#define REVBWTFILE "revchr11.bwt"//"revgenome.bwt"//
 * 	#define REVOCCFILE "revchr11.fmv"//"revgenome.fmv"//
 * 	#define REVSAFILE "revchr11.sa"
 */

#define MAXSTRINGLENGTH 256 //36//36//6+EXTRA//6+EXTRA//36

#define MAXTAG 256
#define INDELMARK (MAXTAG-1)
#define INSERTMARK 63
#define DELETEMARK 70

//}-----------------------------  DEFINES  -------------------------------------------------/




//{-----------------------------  INCLUDE FILES  -------------------------------------------------/
#include <stdarg.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>
//#include "fmsearch.h"
#include "zlib.h"
#include <map>
extern "C"
{
#include <time.h>
#include "MemManager.h"
#include "MiscUtilities.h"
#include "TextConverter.h"
#include "iniparser.h"
#include "BWT.h"
}
//}-----------------------------  INCLUDE FILES  -------------------------------------------------/
//#pragma pack(push)
//#pragma pack(1)
struct Offset_Record
{
	char Genome[40];
	unsigned Offset;
};

struct Ann_Info
{
	unsigned Size;
	unsigned Cumulative_Size;
	char* Name;
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


struct Output_Record
{
	//char Name[36];
	unsigned Tag;
	unsigned  Start;
	char Index;
	unsigned char Skip;
	char Mismatches;
	int Gap;
}__attribute__((__packed__));


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


/*struct Mismatches_Record
{

	int 	Gap;
	unsigned Mismatch_Pos;//BTS_PER_LOC|BTS_PER_LOC|...
	unsigned Mismatch_Char;//2|2|...

}__attribute__((__packed__));*/
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
	int Tag;//Tag number
	unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];//BTS_PER_LOC|BTS_PER_LOC|...
	unsigned Mismatch_Char;//2|2|...
	
};
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

//#pragma pack(pop)
//{-----------------------------  GLOBALS  -------------------------------------------------/
void Print_Location (unsigned Range,BWT *fmi);
void Location_To_Genome( unsigned & Location);
void Location_To_Genome(unsigned & Location,Ann_Info & A);
void Process_Default();
void Process_Deep();
void Process_GIS();
void Verify_GIS();
void Extend_GIS();
void Process_SAM();
void Process_Pairend();
void Process_Pairend_Formatted();
void Process_Deep_Formatted();
void Read_INI();
unsigned Get_File_Size(FILE* File);
void Get_Bases (unsigned Location,int StringLength,char Strand);
void cs2nt_DP(int size, const char *nt_ref, char *cs_read,char* cs_qual, char *nt_read, char *btarray);
void cs2nt_JQ(int size, const char *nt_ref, char *cs_read,char* cs_qual, char *nt_read, char *btarray);
char *cs2nt_nt_qual(int size, const char *nt_read, const char *cs_read, char* cs_qual, char *tarray);
void fprintfX(void* Handle, const char* Format_String, ...);
void printfL (char* Format, ...);
FILE* File_Open(const char* File_Name,const char* Mode);
BWT* initFMI(const char* BWTCodeFileName,const char* BWTOccValueFileName, const char* SAfile);
off64_t File_Size;
void File_OpenZ(const char* File_Name,const char* Mode,gzFile & Handle);
void Show_Progress(unsigned Percentage);

char Char_To_Code[256];
char Zip_Buffer[4096];
unsigned short Stats1[7];
unsigned short Stats2[7];
char NORMAL_TAGS;
char FILETYPE;
bool PRINT_DICTIONARY=true;
Mismatches_Record Mismatches;
Mismatches_Record_GIS MismatchesGIS;
Output_Record Record;
Output_Record* Hits;//[BUFFERSIZE + 1];
MMPool *mmPool;
BWT *fmi,*revfmi;
Header Head;
Offset_Record Genome_Offsets[80];
char Genome_Name_Buf[300];

int debug;
int Stat_Size;
int Mismatch_Count;
int Offset=0;
int Force_Maxhits=0;
int Genome_Count=0;
int Genome_Position;//position of offset info
int STRINGLENGTH=36;
int COLOR_JQ_CTorGA=-1;//0forCT and 1 forGA conversion

unsigned CONVERSION_FACTOR;
unsigned MAXHITS;
unsigned Total_Hits=0;
unsigned Off_Hits=0;
unsigned Total_Tags=0;
unsigned Offsets[80];
unsigned DISKBUFFERSIZE=36*1000000;

char Translated_String_Original[MAXSTRINGLENGTH+1];
char Translated_String[MAXSTRINGLENGTH+1];
char Bin_String[MAXSTRINGLENGTH+1];
char Bin_StringComplement[MAXSTRINGLENGTH+1];
char RevString[MAXSTRINGLENGTH+1];
char Org_String[MAXSTRINGLENGTH+1];
char btarray[MAXSTRINGLENGTH*4+1];
char Base_Read[MAXSTRINGLENGTH+1];
char Test_String[137];
char N[500];
char Tag_Type='*';//@ if new tag % if still in an old tag & if end...
char OUTPUT_COMPRESS=FALSE;
char VERIFY_FROM_DISK=FALSE;
char ROLLOVER,SCANBOTH;
char METH=FALSE;
int TAG_COPY_LEN;
int PLUSGAP=0;

gzFile Data_File;
//FILE* Output_File;
void* Output_File;
gzFile Output_FileG;
FILE* Binary_File;
FILE* Location_File;
FILE* Input_FileO;
FILE* Log_File;
const char* Code_To_Char="acgt";const char* Code_To_CharCAPS="ACGT";const char* Code_To_CharC="tgca";const char* Code_To_CharCAPSC="TGCA";
unsigned char* Original_Text;

char* Command_Line_Buffer;
char* Current_String;
char* INPUTFILE;
char* OUTPUTFILE;
char* BWTFILE ; 
char* OCCFILE ;
char* REVBWTINDEX;
char* REVOCCFILE;
char* REVSAFILE;
char* SAFILE;
char* BINFILE;
char* PACFILE;
char* LOCATIONFILE ; 
char* LOGFILE ; 

char MAX_MISMATCHES;
char SOLID=FALSE;
char Description[MAXDES];
char Quality[MAXTAG+1];
char cs_qual[MAXTAG+1];
char Rev_Quality[MAXTAG+1];
char Tag_Copy[MAXTAG+1];
char DP_Str[MAXTAG+1];
char Rev_Tag_Copy[MAXTAG+1];
char Tag_Copy_T[MAXTAG+1];
char New_Record[4];
char Char_To_CodeC[256];
char Char_To_CharC[256];
char Code_To_CodeC[256];
char LOADREVERSEONLY=FALSE;
char PRINT_DESC;
char FORMAT=FALSE;
char BWTFILE_DEFAULT[] = "genome.bwt";//"genome.bwt";// 
char OCCFILE_DEFAULT[] ="genome.fmv";//"genome.fmv";//
char REVBWTINDEX_DEFAULT[] ="revgenome.bwt";//"revgenome.bwt";//
char REVOCCFILE_DEFAULT[] ="revgenome.fmv";//"revgenome.fmv";//
char REVSAFILE_DEFAULT[] ="revgenome.sa";
char SAFILE_DEFAULT[] ="genome.sa";
char BINFILE_DEFAULT[] = "genome.bin";//"genome.bwt";// 
char INPUTFILE_DEFAULT[]="hits.txt";
char OUTPUTFILE_DEFAULT[]="output.txt";
char LOCATIONFILE_DEFAULT[]="location";
char LOG_DEFAULT[]="decode.log";
char USELOCATION = FALSE;
char VERIFY=FALSE;//TRUE;//FALSE;
char EXTEND=FALSE;
char PLUSSTRAND=TRUE;
char SAM=TRUE;//FALSE;
char INPUT_ZIPPED=FALSE;
char MIS_IN_PLUS_READ=FALSE;
char JQ=FALSE;
unsigned SOURCELENGTH;
char CT[4][4];

time_t Start_Time,End_Time;
map <unsigned, Ann_Info> Annotations;
map <unsigned, Ann_Info> ::iterator S,E;
//}-----------------------------  GLOBALS  -------------------------------------------------/
//{---------------------------- Command Line  -------------------------------------------------
option Long_Options[]=
{
{"help",0,NULL,'h'},
{"inputfile",1,NULL,'i'},
{"ganome",1,NULL,'g'},
{"outputfile",1,NULL,'o'},
{"buffersize",1,NULL,'b'},
{"reverseload",0,NULL,'r'},
{"verify",0,NULL,'v'},
{"extend",0,NULL,'e'},
{"verifyfromdisk",0,NULL,'V'},
{"offset",1,NULL,'O'},
{"plusstrand",0,NULL,'p'},
{"mixstrand",1,NULL,'x'},
//{"formatted",1,NULL,'f'},
{"maxhits",1,NULL,'m'},
{"zip",0,NULL,'z'},
{"location",optional_argument,NULL,'l'},
{"sam",0,NULL,'s'},
{"gis",0,NULL,'G'},
{"misplus",0,NULL,'P'},
{"logfile",1,NULL,'L'},
{0,0,0,0}
};
//}---------------------------- Command Line -------------------------------------------------
void Parse_Command_line(int argc, char* argv[]);
void Verify_Deep();
void Verify_Pairend();
void Search_Exact(struct SARange & Tag,int Start,int StringLength,BWT *fmi);
//{-----------------------------  Main  -------------------------------------------------

int main(int argc, char* argv[])
{

//{-----------------------------  COLOR TABLE  -------------------------------------------------/
	CT[0][0]=0;
	CT[0][0]=0;
	CT[1][0]=1;
	CT[1][0]=1;
	CT[2][0]=2;
	CT[2][0]=2;
	CT[3][0]=3;
	CT[3][0]=3;

	CT[0][1]=1;
	CT[0][1]=1;
	CT[1][1]=0;
	CT[1][1]=0;
	CT[2][1]=3;
	CT[2][1]=3;
	CT[3][1]=2;
	CT[3][1]=2;

	CT[0][2]=2;
	CT[0][2]=2;
	CT[2][2]=0;
	CT[2][2]=0;
	CT[1][2]=3;
	CT[1][2]=3;
	CT[3][2]=1;
	CT[3][2]=1;

	CT[0][3]=3;
	CT[0][3]=3;
	CT[3][3]=0;
	CT[3][3]=0;
	CT[2][3]=1;
	CT[2][3]=1;
	CT[1][3]=2;
	CT[1][3]=2;
//}-----------------------------  COLOR TABLE  -------------------------------------------------/
	time(&Start_Time);
	USELOCATION=TRUE;
	INPUTFILE=INPUTFILE_DEFAULT;
	OUTPUTFILE=OUTPUTFILE_DEFAULT;
	BWTFILE = BWTFILE_DEFAULT; 
	OCCFILE = OCCFILE_DEFAULT;
	REVBWTINDEX=REVBWTINDEX_DEFAULT;
	REVOCCFILE=REVOCCFILE_DEFAULT;
	REVSAFILE=REVSAFILE_DEFAULT;
	SAFILE=SAFILE_DEFAULT;
	BINFILE = BINFILE_DEFAULT;
	LOCATIONFILE=LOCATIONFILE_DEFAULT;
	LOGFILE=LOG_DEFAULT;
	Char_To_CodeC['A']=3;Char_To_CodeC['C']=2;Char_To_CodeC['G']=1;Char_To_CodeC['T']=0;Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;Char_To_CodeC['n']=0;Char_To_CodeC['N']=0;//we are using character count to store the fmicode for acgt
	Char_To_CodeC['0']=0;Char_To_CodeC['1']=1;Char_To_CodeC['2']=2;Char_To_CodeC['3']=3;
	Char_To_CharC['A']='T';Char_To_CharC['C']='G';Char_To_CharC['G']='C';Char_To_CharC['T']='A';Char_To_CharC['a']='t';Char_To_CharC['c']='g';Char_To_CharC['g']='c';Char_To_CharC['t']='a';Char_To_CharC['n']='n';Char_To_CharC['N']='N';//we are using character count to store the fmicode for acgt
	Char_To_CharC['0']='T';Char_To_CharC['1']='G';Char_To_CharC['2']='C';Char_To_CharC['3']='A';
	Code_To_CodeC[0]=3;Code_To_CodeC[1]=2;Code_To_CodeC[2]=1;Code_To_CodeC[3]=0;
	Read_INI();
	Char_To_Code['N']=0;Char_To_Code['n']=0;Char_To_Code['A']=0;Char_To_Code['C']=1;Char_To_Code['G']=2;Char_To_Code['T']=3;Char_To_Code['a']=0;Char_To_Code['c']=1;Char_To_Code['g']=2;Char_To_Code['t']=3;Char_To_Code['+']='+';Char_To_Code['-']='-';//we are using character count to store the fmicode for acgt
	Char_To_Code['0']=0;Char_To_Code['1']=1;Char_To_Code['2']=2;Char_To_Code['3']=3;
	Char_To_CodeC['N']=3;Char_To_CodeC['n']=3;Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;Char_To_CodeC['-']='-';Char_To_CodeC['+']='+';//we are using character count to store the fmicode for acgt


	Command_Line_Buffer=(char*)malloc(5000);
	Parse_Command_line(argc,argv);	
	Log_File=File_Open(LOGFILE,"w");
	printfL("\n********************************************************************************************\n*                            batdecode Version :3.00                                       *\n********************************************************************************************\n");

	File_OpenZ(INPUTFILE,"rb",Data_File);
	gz_stream *s=(gz_stream*)Data_File;
	if (s->transparent) INPUT_ZIPPED=FALSE; else INPUT_ZIPPED=TRUE;
	Input_FileO=s->file;

	fseek(Input_FileO, 0L, SEEK_END);
	File_Size = ftello64(Input_FileO);
	gzrewind(Data_File);//,0,SEEK_SET);//go top

	if (! VERIFY) 
	{
		if(OUTPUT_COMPRESS)
		{
			File_OpenZ(OUTPUTFILE,"w",Output_FileG);//Open output file...
			Output_File=Output_FileG;
		}
		else
		{
			Output_File=File_Open(OUTPUTFILE,"w");
			if(setvbuf((FILE*)Output_File,NULL,_IOFBF,DISKBUFFERSIZE*sizeof(long))) printfL("Buffer allocation failure... Disk access will be unbffered\n");
		}
	}
	//if(setvbuf(Data_File,NULL,_IOFBF,DISKBUFFERSIZE*sizeof(long))) printf("Buffer allocation failure... Disk access will be unbffered\n");
	

	if (USELOCATION)
	{
		Location_File=File_Open(LOCATIONFILE,"r");
		char* Genome_Name;unsigned Off_Cum=0,Off=0;
		while (fgets(Genome_Name_Buf,300,Location_File)!=0)
		{
			Off=atoi(Genome_Name_Buf);
			if (Genome_Count)
			{
				Annotations[Off_Cum].Size=Off;
				Annotations[Off_Cum].Name=Genome_Name;
			}
			Off_Cum+=Off;

			fgets(Genome_Name_Buf,300,Location_File);
			for(int i=0;i<40;i++) 
			{
				if (Genome_Name_Buf[i] == '\n' ||Genome_Name_Buf[i] == '\r')
				{ 
					Genome_Name_Buf[i]=0;
					Genome_Name=new char[i+1];strcpy(Genome_Name,Genome_Name_Buf);
					break;
				} 
			}
			Genome_Count++;	
		}
		if (!Genome_Count) {printfL("Error Loading genome locations..\n");exit(1);} else printfL ("%d Genomes Loaded..\n",Genome_Count);
		S=Annotations.begin();E=Annotations.end();
	}

	gzread(Data_File,&Head,sizeof(Header));
	if(!(Head.ID[0]=='B'&&Head.ID[1]=='A'&&Head.ID[2]=='T')) {printfL("Not a BAT file\n");exit(0);};
	if (Force_Maxhits) MAXHITS=Force_Maxhits; else MAXHITS=Head.MAXHITS;
	STRINGLENGTH=Head.Tag_Length;
	PRINT_DESC=Head.Print_Desc;
	FILETYPE=Head.FILETYPE;
	gzread(Data_File,&MAX_MISMATCHES,sizeof(MAX_MISMATCHES));if (MAX_MISMATCHES >5) Stat_Size=7*sizeof(unsigned short); else Stat_Size=(MAX_MISMATCHES+1)*sizeof(unsigned short);
	gzread(Data_File,&TAG_COPY_LEN,sizeof(int));
	gzread(Data_File,&ROLLOVER,sizeof(char));
	gzread(Data_File,&SCANBOTH,sizeof(char));
	if (Head.Index_Count>=100) //solid marker..
	{
		SOLID=TRUE;Head.Index_Count -= 100; printfL ("SOLiD Output ...\n");
		Char_To_CodeC['a']=0;Char_To_CodeC['c']=1;Char_To_CodeC['g']=2;Char_To_CodeC['t']=3;
		Char_To_CodeC['A']=0;Char_To_CodeC['C']=1;Char_To_CodeC['G']=2;Char_To_CodeC['T']=3;

	} 

	if(!Head.Index_Count)//!LOADREVERSEONLY)
	{
		fmi=initFMI(BWTFILE,OCCFILE,SAFILE);//Load FM indexes
		printfL("Forward index loaded ...\n");
	}

	if(SOLID)
	{
		FILE* Original_File=File_Open(PACFILE,"rb");
		Original_Text=(unsigned char*) malloc(Get_File_Size(Original_File));
		fread(Original_Text,Get_File_Size(Original_File),1,Original_File);
	}
	revfmi=initFMI(REVBWTINDEX,REVOCCFILE,REVSAFILE);//Load FM indexes
	SOURCELENGTH = revfmi->textLength;
	printfL("Reverse index loaded ...\n");
	CONVERSION_FACTOR=revfmi->textLength-STRINGLENGTH;//+1;
	if (Head.HITMODE == PAIREND) NORMAL_TAGS=FALSE; else NORMAL_TAGS=TRUE;
	printfL("Using the genome files\n %s\t %s\t %s\n %s\t %s\t %s\n", BWTFILE,OCCFILE,SAFILE,REVBWTINDEX,REVOCCFILE,REVSAFILE); 
	printfL("Input File : %s\t  Output File : %s \n",INPUTFILE, OUTPUTFILE);
	printfL("Length of a tag : %d\n", STRINGLENGTH);
	if (USELOCATION) printfL("Using location file : %s\n", LOCATIONFILE);
	if (INPUT_ZIPPED) printfL("--< Compressed input >--\n");
	if (OUTPUT_COMPRESS) printfL("=| Compressing output |=\n");
	if (Offset) printfL("Offset %d \n",Offset);
	if (PRINT_DESC==GISMODE)
	{
		if (SAM) Process_SAM();
		else if (VERIFY) Verify_GIS(); 
		else if (EXTEND) Extend_GIS(); 
		else Process_GIS();
	}
	else
	{
		switch (Head.HITMODE)
		{
			case DEFAULT:
				if (VERIFY) {printfL("File format insufficient for verification...\n");exit(0);}
				Process_Default();
				break;

			case DEEP:
				if(VERIFY) Verify_Deep();
				else if(FORMAT) Process_Deep_Formatted(); else Process_Deep();
				break;
			case PAIREND:
				if(VERIFY) Verify_Pairend();
				else if(FORMAT) Process_Pairend_Formatted(); else Process_Pairend();
				break;
		}
	}

	if (OUTPUT_COMPRESS) gzclose(Output_File);	
	printf("\r[++++++++100%%+++++++++]\n");fprintf(Log_File,"100%%\n");//progress bar....
	printfL("Total Tags/Hits[Rejected Hits] %d / %d[%d]\n",Total_Tags,Total_Hits,Off_Hits);
	time(&End_Time);printfL("\n Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
}

//}-----------------------------  Main  -------------------------------------------------

//{-----------------------------  SAM Mode  -------------------------------------------------
void Process_SAM()
{
	unsigned i, Bitsread;
	unsigned Start,End;
	unsigned Location;
	int Tail_Start;
	unsigned Previous_Tag,Hits;
	char letter;
	unsigned pos;
	int Desc_End=0;
	int Conversion_Factor,StringLength,Length_Array[3];
	unsigned Progress=0,Average_Length;
	unsigned Number_of_Tags=100;
	char Mismatches_Desc[1000];
	char Last_Orientation,Last_Half;
	char* Ins_Format;
	char* Del_Format;
	char* Mis_Format;
	char Paired_Read=FALSE;
	char* QNAME=Description+1;
	int FLAG;	
	char* RNAME;
	unsigned POS;
	int MAPQ=255; 
	char* Mismatch_Desc_Ptr;
	char M[MAXTAG];
	char* CIGAR;//[200];
	char CIGAR_Process[200];
	char MRNM='*';
	char ISIZE='0';
	char* SEQ=Tag_Copy;
	char* QUALITY=Quality;
	int IGNOREHEAD=Head.IGNOREHEAD;
	char Last_Marker;
	char Strand;
	Ann_Info A;

#define UNPAIREDF 	0
#define MINUSF 		0x10
#define UNMAPPEDF 	0x4
#define PAIREDF 	1

#define HEAD_READ	0x40
#define TAIL_READ	0x80


 
	
//--------------------------------------------------------------------------------------------
	printf("======================]\r[");//progress bar....
	for(int i=0;i<=MAXTAG;i++) {Quality[i]=Tag_Copy[i]=0;}
	if (FILETYPE==FA) Quality[0]='*';
	if (Head.HITMODE == PAIREND)
	{
		Paired_Read=TRUE;
		gzread(Data_File,&Length_Array[1],sizeof(int));
		gzread(Data_File,&Length_Array[2],sizeof(int));
	}
	else {Length_Array[1]=STRINGLENGTH;Length_Array[2]=0;}
	for (int i=0;i<MAXDES;i++) Description[i]=0;
	if (PRINT_DICTIONARY)
	{
		while (S!=E)
		{
			Ann_Info T=S->second;fprintfX(Output_File,"@SQ\tSN:%s\tLN:%d\n",T.Name,T.Size);
			S++;
		}
	}
//--------------------------------------------------------------------------------------------
	for(;;)
	{
		i=0;
		Last_Marker=Tag_Type;
		gzread(Data_File,&Tag_Type,1);
		if (Last_Marker == '@' && Tag_Type !='%') 
		{
			for(int i=0;i<TAG_COPY_LEN;i++) SEQ[i]=((SEQ[i]<'4' &&  SEQ[i]>='0') ? Code_To_CharCAPS[SEQ[i]-'0'] : SEQ[i]); 
			if (Paired_Read)
			{
				if(FILETYPE==FA) {fprintfX(Output_File,"%s\t5\t*\t0\t0\t*\t*\t0\t0\t%s\t*\n",QNAME,Tag_Copy);fprintfX(Output_File,"%s\t5\t*\t0\t0\t*\t*\t0\t0\t%s\t*\n",QNAME,Tag_Copy+Length_Array[1]+1+2*IGNOREHEAD);}
				else {fprintfX(Output_File,"%s\t5\t*\t0\t0\t*\t*\t0\t0\t%s\t*\t%s\n",QNAME,Tag_Copy,Quality);fprintfX(Output_File,"%s\t5\t*\t0\t0\t*\t*\t0\t0\t%s\t*\t%s\n",QNAME,Tag_Copy+Length_Array[1]+1+2*IGNOREHEAD,Quality+Length_Array[1]+1+2*IGNOREHEAD);}
				FLAG=1;
			}
			else 
			{
				fprintfX(Output_File,"%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",QNAME,SEQ,QUALITY);
				FLAG=0;
			}
		}
		if('&' == Tag_Type) {break;}
		if('@' == Tag_Type)//start of new tag
		{
			
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(Input_FileO);
				Average_Length=Current_Pos/Total_Tags+1;//+1 avoids divide by zero..
				Number_of_Tags=(File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/File_Size);
			}

			//fprintfX(Output_File,"@\n");
			gzgets(Data_File,Description,MAXDES);
			gzread(Data_File,Stats1,Stat_Size);
			if(!NORMAL_TAGS) gzread(Data_File,Stats2,Stat_Size);
			gzread(Data_File,Tag_Copy,TAG_COPY_LEN);if(NORMAL_TAGS && FILETYPE==FQ) gzread(Data_File,Quality,TAG_COPY_LEN);
			for(Desc_End=0;Description[Desc_End]!='\n' && Description[Desc_End]!='\r' && Description[Desc_End]!=0;Desc_End++);Description[Desc_End]=0;
			if (Paired_Read) 
			{
				for(Tail_Start=0;Tag_Copy[Tail_Start]!='\t';Tail_Start++);Tag_Copy[Tail_Start]=0;
				if(FILETYPE==FQ) {gzread(Data_File,Quality,TAG_COPY_LEN);Quality[Tail_Start]=0;}
				FLAG=1;//for(Tail_Start=Tag_Copy;*Tail_Start!='\t';Tail_Start++);*Tail_Start=0;Tail_Start++;
				for(int i=0;i<Tail_Start;i++) 
				{
					if(SOLID) Rev_Tag_Copy[i+2]=Tag_Copy[Tail_Start-i-1]; else Rev_Tag_Copy[i]=Char_To_CharC[Tag_Copy[Tail_Start-i-1]];
					if(SOLID) Rev_Quality[i+2]=Quality[Tail_Start-i-1]; else Rev_Quality[i]=Quality[Tail_Start-i-1];
				}//reverse head...
				Tail_Start++;int j=TAG_COPY_LEN-1;
				for(int i=Tail_Start;i<TAG_COPY_LEN;i++,j--) 
				{
					if(SOLID)  Rev_Tag_Copy[i+2]=Tag_Copy[j]; else Rev_Tag_Copy[i]=Char_To_CharC[Tag_Copy[j]];
					if(SOLID)  Rev_Quality[i+2]=Quality[j]; else Rev_Quality[i]=Quality[j];
				}//reverse tail...
			}//seek tail read...
			else 
			{
				for(int i=0;i<TAG_COPY_LEN;i++) 
				{
					if(SOLID) Rev_Tag_Copy[i+2]=Tag_Copy[TAG_COPY_LEN-i-1]; else Rev_Tag_Copy[i]=Char_To_CharC[Tag_Copy[TAG_COPY_LEN-i-1]]; 
					if(SOLID) Rev_Quality[i+2]=Quality[TAG_COPY_LEN-1-i]; else Rev_Quality[i]=Quality[TAG_COPY_LEN-1-i];
				}//reverse head...
				Total_Tags++;FLAG=0;
			}
			Hits=0;
			continue;
		}		
		else//Process hit....
		{
//--------------------------------------------------------------------------------------------
			//Conversion_Factor=CONVERSION_FACTOR;
			gzread(Data_File,New_Record,3);
			if (New_Record[2]) gzread(Data_File,N,New_Record[2]*2);
			gzread(Data_File,&Record,sizeof(Output_Record));
			Mismatch_Count = Record.Mismatches;
			gzread(Data_File,&MismatchesGIS,Mismatch_Count+sizeof(unsigned));
			StringLength=Length_Array[New_Record[0]];
			Conversion_Factor=revfmi->textLength-StringLength;//+1;
			if (!NORMAL_TAGS)
			{
				if(!Hits) {Last_Half=New_Record[0];Total_Tags++;}
				else if (Last_Half != New_Record[0]) {Last_Half=New_Record[0];Total_Tags++;Hits=0;}
			}
			if (SCANBOTH)
			{
				if(!Hits) Last_Orientation=New_Record[1];
				else if (Last_Orientation != New_Record[1])
				{
					Last_Orientation=New_Record[1];
					Hits=0;
				}
			}
//--------------------------- Find CIGAR ----------------------------------------------
			Mismatch_Desc_Ptr=Mismatches_Desc;
			for (int i=0;i<STRINGLENGTH;i++){M[i]=0;}//init mismatch pos..
			if (Mismatch_Count)//mismatch positions..
			{
				if (MIS_IN_PLUS_READ && New_Record[1] == '-')//give mismatches relative to plus read..
					for( int i=0;i<Mismatch_Count;i++) M[Length_Array[New_Record[0]]-MismatchesGIS.Mismatch_Pos[i]-1]=Code_To_CharCAPSC[MismatchesGIS.Mismatch_Char>>(2*i) & 3];
				else 
					for( int i=0;i<Mismatch_Count;i++) M[MismatchesGIS.Mismatch_Pos[i]]=Code_To_CharCAPS[MismatchesGIS.Mismatch_Char>>(2*i) & 3];
			}
			int j=0;
			for (int i=0;i<STRINGLENGTH;i++)//init mismatch pos..
			{
				if(M[i]){ Mismatch_Desc_Ptr+=sprintf(Mismatch_Desc_Ptr,"%d%c",j,M[i]);j=0;}
				else j++;
			}
			if (j) sprintf(Mismatch_Desc_Ptr,"%d",j);


			if(MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] == INDELMARK)//Indel...
			{
				for (int i=0;i<STRINGLENGTH;i++) CIGAR_Process[i]='M';//initialize cigar calculation...
				pos=MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-3];
				if (MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-2] == INSERTMARK)//insertion
				{
					Conversion_Factor--;StringLength++;
					CIGAR_Process[sprintf(CIGAR_Process,"%dM1I%dM",pos,StringLength-pos)]=0;
				}
				else//deletion
				{
					Conversion_Factor++;StringLength--;
					CIGAR_Process[sprintf(CIGAR_Process,"%dM1D%dM",pos,StringLength-pos-1)]=0;
				}
			}
			else
			{
				if(SOLID) CIGAR_Process[sprintf(CIGAR_Process,"%dM",StringLength+1)]=0;
				else CIGAR_Process[sprintf(CIGAR_Process,"%dM",StringLength)]=0;
			}
//--------------------------- Find CIGAR ----------------------------------------------
			for (int j=0;j<=Record.Gap && Hits< MAXHITS;j++)
			{
				if( Record.Index)//print record...
				{
					if (Record.Skip) Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
					else Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
				}
				else
				{
					if (Record.Skip) Location = fmi->saValue[Record.Start/fmi->saInterval]-Record.Skip+1;
					else Location=BWTSaValue(fmi,Record.Start);
				}
			//Orientation...	

				if ((Strand=New_Record[1])=='-') 
				{
					FLAG +=MINUSF;
					QUALITY=Rev_Quality;if(FILETYPE == FA) Rev_Quality[0]='*';
					SEQ=Rev_Tag_Copy;
				}
				else {QUALITY=Quality;SEQ=Tag_Copy;}

				if (SOLID)
				{
					//if (Strand == '+') Location -=2;
					Get_Bases(Location,STRINGLENGTH+1,Strand);
					for (int i=0;i<STRINGLENGTH+2;i++) {cs_qual[i] = QUALITY[i+2]-33;if(cs_qual[i]>60) cs_qual[i]=60;}//fix quality..
					if (!j) for(int i=0;i<STRINGLENGTH;i++) {SEQ[i+2]=Char_To_Code[SEQ[i+2]];}
					cs2nt_DP(STRINGLENGTH,Org_String,SEQ+2,cs_qual, Base_Read, btarray);
					QUALITY=cs2nt_nt_qual(STRINGLENGTH, Base_Read, SEQ+2, cs_qual, btarray);
					//for (int i=0;i<=STRINGLENGTH;i++) SEQ[i]=Code_To_CharCAPS[Base_Read[i+1]];SEQ[STRINGLENGTH+1]=0;QUALITY[STRINGLENGTH+1]=0;
					for (int i=0;i<=STRINGLENGTH;i++) SEQ[i]=Code_To_CharCAPS[Base_Read[i+1]];SEQ[STRINGLENGTH-1]=0;QUALITY[STRINGLENGTH-1]=0;
				}
				//Location_To_Genome(Location);
				Location_To_Genome(Location,A);
				if (Location+STRINGLENGTH > A.Size) 
				{
					Off_Hits++;
					FLAG=UNMAPPEDF;
				}
				CIGAR=CIGAR_Process;
				if (Paired_Read)
				{
					if (New_Record[0] == 1) FLAG += HEAD_READ; else {FLAG += TAIL_READ;QUALITY+=Length_Array[1]+1;SEQ +=Length_Array[1]+1+2*IGNOREHEAD;} //head tail.. 
					if(FILETYPE==FA) fprintfX(Output_File,"%s\t%d\t%s\t%d\t255\t%s\t*\t0\t0\t%s\t*\n",QNAME,FLAG,A.Name,Location+1,CIGAR,SEQ);
					else fprintfX(Output_File,"%s\t%d\t%s\t%d\t255\t%s\t*\t0\t0\t%s\t%s\n",QNAME,FLAG,A.Name,CIGAR,SEQ,QUALITY);
					fprintfX(Output_File,"\tNM:i:%d\tMD:Z:%s\n",Mismatch_Count,Mismatches_Desc);
					FLAG=1;
				}
				else 
				{
					fprintfX(Output_File,"%s\t%d\t%s\t%d\t255\t%s\t*\t0\t0\t%s\t%s",QNAME,FLAG,A.Name,Location+1,CIGAR,SEQ,QUALITY);
					fprintfX(Output_File,"\tNM:i:%d\tMD:Z:%s\n",Mismatch_Count,Mismatches_Desc);
					FLAG=0;
				}
				Hits++;Record.Start++;Total_Hits++;
			}
//--------------------------- Decode Mismatches ----------------------------------------------
		}

	}
}
//}-----------------------------  SAM Mode  -------------------------------------------------

//{-----------------------------  GIS Mode  -------------------------------------------------
void Process_GIS()
{
	unsigned i, Bitsread;
	unsigned Start,End;
	unsigned Location;
	int Desc_End;
	unsigned Previous_Tag,Hits;
	char letter;
	unsigned pos;
	int Conversion_Factor,StringLength,Length_Array[3];
	unsigned Progress=0,Average_Length;
	unsigned Number_of_Tags=100;
	char Mismatches_Desc[1000];
	char* Mismatch_Desc_Ptr;
	char Last_Orientation,Last_Half;
	char* Ins_Format;
	char* Del_Format;
	char* Mis_Format;
	int IGNOREHEAD=Head.IGNOREHEAD;
	char* QUALITY=Quality;
	char* SEQ=Tag_Copy;
	Ann_Info A;

	
//--------------------------------------------------------------------------------------------
	printf("======================]\r[");//progress bar....
	for(int i=0;i<=MAXTAG;i++) {Quality[i]=Tag_Copy[i]=0;}
	if (Head.HITMODE == PAIREND)
	{
		gzread(Data_File,&Length_Array[1],sizeof(int));
		gzread(Data_File,&Length_Array[2],sizeof(int));
	}
	else {Length_Array[1]=STRINGLENGTH;Length_Array[2]=0;}

	if (FORMAT) 
	{
		Ins_Format="%d<%c:";Del_Format="%d>D:";Mis_Format="%d>%c:";
		fprintfX(Output_File,"@\t%d\t%d\n",Length_Array[1]+1000,Length_Array[2]+1000);
	}
	else
	{
		Ins_Format="%d<%c\t";Del_Format="%d>D\t";Mis_Format="%d>%c\t"; 
	}
//--------------------------------------------------------------------------------------------
	for(;;)
	{
		i=0;
		gzread(Data_File,&Tag_Type,1);
		if('&' == Tag_Type) break;
		if('@' == Tag_Type)//start of new tag
		{
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(Input_FileO);
				Average_Length=Current_Pos/Total_Tags+1;//+1 avoids divide by zero..
				Number_of_Tags=(File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/File_Size);
			}

			fprintfX(Output_File,"@\n");
			gzgets(Data_File,Description,MAXDES);
			gzread(Data_File,Stats1,Stat_Size);
			if(!NORMAL_TAGS) gzread(Data_File,Stats2,Stat_Size);
			gzread(Data_File,Tag_Copy,TAG_COPY_LEN);//if(NORMAL_TAGS && FILETYPE==FQ) gzread(Data_File,Quality,TAG_COPY_LEN);
			if(FILETYPE != FQ) Quality[0]='*';
			if(FILETYPE==FQ) gzread(Data_File,Quality,TAG_COPY_LEN);
			for(Desc_End=0;Description[Desc_End]!='\n' && Description[Desc_End]!='\r' && Description[Desc_End]!=0;Desc_End++);Description[Desc_End]=0;
			if (NORMAL_TAGS){fprintfX(Output_File,"%s\t%s\t%s\t%d:%d:%d:%d:%d:%d:%d\n",Description,Tag_Copy,Quality,Stats1[0],Stats1[1],Stats1[2],Stats1[3],Stats1[4],Stats1[5],Stats1[6]);Total_Tags++;}
			else{fprintfX(Output_File,"%s\t%s\t%s\t%d:%d:%d:%d:%d:%d:%d\t%d:%d:%d:%d:%d:%d:%d\n",Description,Tag_Copy,Quality,Stats1[0],Stats1[1],Stats1[2],Stats1[3],Stats1[4],Stats1[5],Stats1[6],Stats2[0],Stats2[1],Stats2[2],Stats2[3],Stats2[4],Stats2[5],Stats2[6]);}

			Hits=0;
			if(FILETYPE != FQ) Quality[0]=0;
			continue;
		}		
		else//Process hit....
		{
//--------------------------------------------------------------------------------------------
			//Conversion_Factor=CONVERSION_FACTOR;
			gzread(Data_File,New_Record,3);
			if (New_Record[2]) gzread(Data_File,N,New_Record[2]*2);
			gzread(Data_File,&Record,sizeof(Output_Record));
			Mismatch_Count = Record.Mismatches;
			gzread(Data_File,&MismatchesGIS,Mismatch_Count+sizeof(unsigned));
			StringLength=Length_Array[New_Record[0]];
			Conversion_Factor=revfmi->textLength-StringLength;//+1;
			if (!NORMAL_TAGS)
			{
				if(!Hits) {Last_Half=New_Record[0];Total_Tags++;}
				else if (Last_Half != New_Record[0]) {Last_Half=New_Record[0];Total_Tags++;Hits=0;}
			}
			if (SCANBOTH)
			{
				if(!Hits) Last_Orientation=New_Record[1];
				else if (Last_Orientation != New_Record[1])
				{
					Last_Orientation=New_Record[1];
					Hits=0;
				}
			}
//--------------------------- Decode Mismatches ----------------------------------------------
			Mismatch_Desc_Ptr=Mismatches_Desc;
			if (New_Record[2])
			{
				for(int i=0;i<New_Record[2]*2;i+=2)
				{
					Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,"%d:%c\t",N[i],N[i+1]);
				}
			}
			if(MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] == INDELMARK)//Indel...
			{
				pos=MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-3];
				if (MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-2] == INSERTMARK)//insertion
				{
					Conversion_Factor--;StringLength++;
					Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,Ins_Format,pos,Code_To_Char[MismatchesGIS.Mismatch_Char >> 2*2]-('a'-'A'));
				}
				else//deletion
				{
					Conversion_Factor++;StringLength--;
					Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,Del_Format,pos);
				}
			}

			if (Mismatch_Count)//mismatch positions..
			{
				if (MIS_IN_PLUS_READ && New_Record[1] == '-')//give mismatches relative to plus read..
					for( int i=0;i<Mismatch_Count;i++) Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,Mis_Format,Length_Array[New_Record[0]]-MismatchesGIS.Mismatch_Pos[i]-1,Code_To_CharC[MismatchesGIS.Mismatch_Char>>(2*i) & 3]);
				else 
					for( int i=0;i<Mismatch_Count;i++) Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,Mis_Format,MismatchesGIS.Mismatch_Pos[i],Code_To_Char[MismatchesGIS.Mismatch_Char>>(2*i) & 3]);
			}
			*Mismatch_Desc_Ptr=0;
//--------------------------- Decode Mismatches ----------------------------------------------
			for (int j=0;j<=Record.Gap && Hits< MAXHITS;j++)
			{
				if( Record.Index)//print record...
				{
					if (Record.Skip) Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
					else Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
				}
				else
				{
					if (Record.Skip) Location = fmi->saValue[Record.Start/fmi->saInterval]-Record.Skip+1;
					else Location=BWTSaValue(fmi,Record.Start);
				}
				
				Location -= Offset;
				if (SOLID)
				{
					strcpy(Tag_Copy_T,Tag_Copy);
					SEQ=Tag_Copy_T;
					char First=CT[Char_To_Code[SEQ[0]]][Char_To_Code[SEQ[1]]],FirstR;
					unsigned Loc=Location;
					//if (New_Record[1] == '+') Loc --;//=2;
					Get_Bases(Loc,STRINGLENGTH+1,New_Record[1]);

					if (New_Record[1]=='-') 
					{
						QUALITY=Rev_Quality;
						SEQ=Rev_Tag_Copy;
						First=Code_To_CodeC[First];
						for(int i=0;i<TAG_COPY_LEN;i++) 
						{
							Rev_Tag_Copy[i]=Tag_Copy[TAG_COPY_LEN-i-1]; 
							Rev_Quality[i]=Quality[TAG_COPY_LEN-1-i];
						}//reverse head...
						FirstR=Org_String[STRINGLENGTH-1];
					}
					else {QUALITY=Quality;FirstR=Org_String[0];}
					if(First != FirstR && !j) 
					{
						Mismatch_Count++;
						Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,"F>F");
					}
					for (int i=0;i<STRINGLENGTH+2;i++) {cs_qual[i] = QUALITY[i]-33;if(cs_qual[i]>60) cs_qual[i]=60;}//fix quality..
					if (!j) for(int i=0;i<STRINGLENGTH;i++) {SEQ[i+2]=Char_To_Code[SEQ[i+2]];}
					cs2nt_DP(STRINGLENGTH+2,Org_String,SEQ+2,cs_qual, Base_Read, btarray);
					//QUALITY=cs2nt_nt_qual(STRINGLENGTH+2, Base_Read, SEQ, cs_qual, btarray);
					SEQ=DP_Str;
					for (int i=0;i<=STRINGLENGTH;i++) {Org_String[i]=Code_To_CharCAPS[Org_String[i]];SEQ[i]=Code_To_CharCAPS[Base_Read[i]];}
					Org_String[STRINGLENGTH+2]=SEQ[STRINGLENGTH+1]=QUALITY[STRINGLENGTH+1]=0;
				}

				if(USELOCATION) 
				{
					//Location_To_Genome(Location);
					int HT=0;
					Location_To_Genome(Location,A);
					if (Location+STRINGLENGTH > A.Size) 
					{
						Off_Hits++;
						HT= '*'-New_Record[0]-'0';
					}
					if (!PLUSSTRAND && New_Record[1]=='-') Location=Location+STRINGLENGTH-PLUSGAP;
					if (SOLID)fprintfX(Output_File,"%c\t%s\t%c\t%u\t%d\t%d\t%d\t%s\t%s\t%s\n",New_Record[0]+'0'+HT,A.Name,New_Record[1],Location,Mismatch_Count,Length_Array[New_Record[0]]+1,Record.Tag,SEQ,Org_String,Mismatches_Desc);
					else fprintfX(Output_File,"%c\t%s\t%c\t%u\t%d\t%d\t%d\t%s\n",New_Record[0]+'0'+HT,A.Name,New_Record[1],Location,Mismatch_Count,Length_Array[New_Record[0]],Record.Tag,Mismatches_Desc);
				}
				else 
				{
					fprintfX(Output_File,"%c\t%c\t%u\t%d\t%d\t%d\t%s\n",New_Record[0]+'0',New_Record[1],Location,Mismatch_Count,Length_Array[New_Record[0]],Record.Tag,Mismatches_Desc);
				}
				Hits++;Record.Start++;Total_Hits++;
			}
//--------------------------- Decode Mismatches ----------------------------------------------
		}

	}
}



//{-------------------------------------------- EXTEND GIS ---------------------------------------------

void Extend_GIS()
{
	unsigned i, Bitsread;
	unsigned Start,End;
	unsigned Location;
	int Desc_End;
	unsigned Previous_Tag,Hits;
	char letter;
	unsigned pos;
	int Conversion_Factor,StringLength,Length_Array[3];
	unsigned Progress=0,Average_Length;
	unsigned Number_of_Tags=100;
	char Mismatches_Desc[1000];
	char* Mismatch_Desc_Ptr;
	char Last_Orientation,Last_Half;
	char* Ins_Format;
	char* Del_Format;
	char* Mis_Format;
	int IGNOREHEAD=Head.IGNOREHEAD;

	char FIRSTPASS=TRUE;
	int HEAD_LENGTH,TAIL_LENGTH;
	char* String;
	char* Org_Genome;
	char* Mis_Genome=(char*)malloc(2000);

	Binary_File=fopen(BINFILE,"rb");
	unsigned Genome_File_Size=Get_File_Size(Binary_File);
	if(!(Org_Genome=(char*)malloc(Genome_File_Size))) {printf("Out of memory for genome..\n");exit(0);};
	if (!fread(Org_Genome,Genome_File_Size,1,Binary_File)) {printf ("Error reading the genome file...\n");exit(0);}
	
//--------------------------------------------------------------------------------------------
	printf("======================]\r[");//progress bar....
	for(int i=0;i<=MAXTAG;i++) {Quality[i]=Tag_Copy[i]=0;}
	if (Head.HITMODE == PAIREND)
	{
		gzread(Data_File,&Length_Array[1],sizeof(int));HEAD_LENGTH=Length_Array[1];
		gzread(Data_File,&Length_Array[2],sizeof(int));TAIL_LENGTH=Length_Array[2];
	}
	else
	{
		HEAD_LENGTH=Length_Array[1]=STRINGLENGTH;
	}

	if(FILETYPE != FQ) Quality[0]='*';

	if (FORMAT) 
	{
		Ins_Format="%d<%c:";Del_Format="%d>D:";Mis_Format="%d>%c:";
		fprintfX(Output_File,"@\t%d\t%d\n",Length_Array[1]+1000,Length_Array[2]+1000);
	}
	else
	{
		Ins_Format="%d<%c\t";Del_Format="%d>D\t";Mis_Format="%d>%c\t"; 
	}
//--------------------------------------------------------------------------------------------
	for(;;)
	{
		i=0;
		gzread(Data_File,&Tag_Type,1);
		if('&' == Tag_Type) break;
		if('@' == Tag_Type)//start of new tag
		{
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(Input_FileO);
				Average_Length=Current_Pos/Total_Tags+1;//+1 avoids divide by zero..
				Number_of_Tags=(File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/File_Size);
			}

			fprintfX(Output_File,"@\n");
			gzgets(Data_File,Description,MAXDES);
			gzread(Data_File,Stats1,Stat_Size);
			if(!NORMAL_TAGS) gzread(Data_File,Stats2,Stat_Size);
			gzread(Data_File,Tag_Copy,TAG_COPY_LEN);//if(NORMAL_TAGS && FILETYPE==FQ) gzread(Data_File,Quality,TAG_COPY_LEN);
			if (FIRSTPASS) 
			{
				STRINGLENGTH = strlen(Tag_Copy);
				if (!METH) {HEAD_LENGTH=TAIL_LENGTH=STRINGLENGTH;}
				FIRSTPASS=FALSE;
			}
			if(FILETYPE==FQ) gzread(Data_File,Quality,TAG_COPY_LEN);
			for(Desc_End=0;Description[Desc_End]!='\n' && Description[Desc_End]!='\r' && Description[Desc_End]!=0;Desc_End++);Description[Desc_End]=0;
			if (NORMAL_TAGS){fprintfX(Output_File,"%s\t%s\t%s\t%d:%d:%d:%d:%d:%d:%d\n",Description,Tag_Copy,Quality,Stats1[0],Stats1[1],Stats1[2],Stats1[3],Stats1[4],Stats1[5],Stats1[6]);Total_Tags++;}
			else{fprintfX(Output_File,"%s\t%s\t%s\t%d:%d:%d:%d:%d:%d:%d\t%d:%d:%d:%d:%d:%d:%d\n",Description,Tag_Copy,Quality,Stats1[0],Stats1[1],Stats1[2],Stats1[3],Stats1[4],Stats1[5],Stats1[6],Stats2[0],Stats2[1],Stats2[2],Stats2[3],Stats2[4],Stats2[5],Stats2[6]);}

			Hits=0;
			continue;
		}		
		else//Process hit....
		{
			memcpy(Translated_String,Tag_Copy,MAXSTRINGLENGTH);
			//for(Next_String=Translated_String;Next_String[0]!='\n' && Next_String[0]!='\r' && Next_String[0]!=0;Next_String++) if(*Next_String=='\t') *Next_String=0; else Next_String[0]=tolower(Next_String[0]);
//--------------------------------------------------------------------------------------------
			//Conversion_Factor=CONVERSION_FACTOR;
			gzread(Data_File,New_Record,3);
			if (New_Record[2]) gzread(Data_File,N,New_Record[2]*2);
			gzread(Data_File,&Record,sizeof(Output_Record));
			Mismatch_Count = Record.Mismatches;
			gzread(Data_File,&MismatchesGIS,Mismatch_Count+sizeof(unsigned));
			StringLength=Length_Array[New_Record[0]];
			Conversion_Factor=revfmi->textLength-StringLength;//+1;
			if (!NORMAL_TAGS)
			{
				if(!Hits) {Last_Half=New_Record[0];Total_Tags++;}
				else if (Last_Half != New_Record[0]) {Last_Half=New_Record[0];Total_Tags++;Hits=0;}
			}
			if (SCANBOTH)
			{
				if(!Hits) Last_Orientation=New_Record[1];
				else if (Last_Orientation != New_Record[1])
				{
					Last_Orientation=New_Record[1];
					Hits=0;
				}
			}

			if(New_Record[0]==1) {String=Translated_String+IGNOREHEAD;STRINGLENGTH=HEAD_LENGTH;}
			else {String=Translated_String+2*IGNOREHEAD+HEAD_LENGTH+1;STRINGLENGTH=TAIL_LENGTH;}

//--------------------------- Handle N's  ----------------------------------------------
			if (New_Record[2])
			{
				for(int i=0;i<New_Record[2]*2;i+=2)
				{
					String[N[i]]=N[i+1];
				}
			}
//--------------------------- Handle N's  ----------------------------------------------
			if (New_Record[1] == '-')
			{
				for (unsigned i=0;i<=STRINGLENGTH-1;i++) RevString[STRINGLENGTH-1-i]=Code_To_Char[Char_To_CodeC[String[i]]];
				String=RevString;
			}

//--------------------------- Decode Mismatches ----------------------------------------------
			/*Mismatch_Desc_Ptr=Mismatches_Desc;
			if (New_Record[2])
			{
				for(int i=0;i<New_Record[2]*2;i+=2)
				{
					Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,"%d:%c\t",N[i],N[i+1]);
				}
			}
			if(MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] == INDELMARK)//Indel...
			{
				pos=MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-3];
				if (MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-2] == INSERTMARK)//insertion
				{
					Conversion_Factor--;StringLength++;
					Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,Ins_Format,pos,Code_To_Char[MismatchesGIS.Mismatch_Char >> 2*2]-('a'-'A'));
				}
				else//deletion
				{
					Conversion_Factor++;StringLength--;
					Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,Del_Format,pos);
				}
			}

			if (Mismatch_Count)//mismatch positions..
			{
				for( int i=0;i<Mismatch_Count;i++) Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,Mis_Format,MismatchesGIS.Mismatch_Pos[i],Code_To_Char[MismatchesGIS.Mismatch_Char>>(2*i) & 3]);
			}
			*Mismatch_Desc_Ptr=0;*/
//--------------------------- Decode Mismatches ----------------------------------------------
			for (int j=0;j<=Record.Gap && Hits< MAXHITS;j++)
			{
				if( Record.Index)//print record...
				{
					if (Record.Skip) Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
					else Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
				}
				else
				{
					if (Record.Skip) Location = fmi->saValue[Record.Start/fmi->saInterval]-Record.Skip+1;
					else Location=BWTSaValue(fmi,Record.Start);
				}
				
				Location -= Offset;
				if (!PLUSSTRAND && New_Record[1]=='-') Location=Location+STRINGLENGTH-PLUSGAP;

				memcpy(Test_String,Org_Genome+Location,STRINGLENGTH);
				int Misses=0;int Mis_Count=0;
				//for (int i=0;i<StringLength;i++)
				for (int i=0;i<STRINGLENGTH;i++)
				{
					String[i]=tolower(String[i]);
					if (Test_String[i] != String[i]) 
					{
						if (METH) 
						{
							if (Test_String[i]=='c' && String[i]=='t')//Methylated..
							{
								Misses+=sprintf(Mis_Genome+Misses,"\t%d",i);
								Mis_Count++;
							}
						}
						else
						{
							Misses+=sprintf(Mis_Genome+Misses,"\t%d\t%c:%c",i,Test_String[i],String[i]);
							Mis_Count++;
						}
					}
				}
				Mis_Genome[Misses]=0;

				if(USELOCATION) 
				{
					Location_To_Genome(Location);
					fprintfX(Output_File,"%c\t%s\t%c\t%u\t%d\t%d\t%d\t%d\t%s\n",New_Record[0]+'0',Genome_Offsets[Genome_Position].Genome,New_Record[1],Location,Mismatch_Count,Length_Array[New_Record[0]],Record.Tag,Mis_Count,Mis_Genome);
				}
				else 
				{
					fprintfX(Output_File,"%c\t%c\t%u\t%d\t%d\t%d\t%d\t%s\n",New_Record[0]+'0',New_Record[1],Location,Mismatch_Count,Length_Array[New_Record[0]],Record.Tag,Mis_Count,Mis_Genome);
				}
				Hits++;Record.Start++;Total_Hits++;
			}
//--------------------------- Decode Mismatches ----------------------------------------------
		}

	}
}


//}-------------------------------------------- EXTEND GIS ---------------------------------------------


void Verify_GIS()
{
	unsigned i, Bitsread;
	unsigned Start,End;
	unsigned Location;
	int Desc_End;
	unsigned Previous_Tag,Hits;
	char letter;
	unsigned pos;
	int Conversion_Factor,StringLength,Length_Array[3];
	unsigned Progress=0,Average_Length;
	unsigned Number_of_Tags=100;
	char Mismatches_Desc[1000];
	char* Mismatch_Desc_Ptr;
	char Last_Orientation,Last_Half;
	char* Next_String;
	char* String;
	int HEAD_LENGTH,TAIL_LENGTH;
	SARange Tag;
	int IGNOREHEAD=Head.IGNOREHEAD;

	Binary_File=fopen(BINFILE,"rb");
//--------------------------------------------------------------------------------------------
	printf("======================]\r[");//progress bar....
	for(int i=0;i<=MAXTAG;i++) {Quality[i]=Tag_Copy[i]=0;}
	if (Head.HITMODE == PAIREND)
	{
		gzread(Data_File,&Length_Array[1],sizeof(int));HEAD_LENGTH=Length_Array[1];
		gzread(Data_File,&Length_Array[2],sizeof(int));TAIL_LENGTH=Length_Array[2];
	}
	else
	{
		HEAD_LENGTH=Length_Array[1]=STRINGLENGTH;
	}
//--------------------------------------------------------------------------------------------
	for(;;)
	{
		i=0;
		gzread(Data_File,&Tag_Type,1);
		if('&' == Tag_Type) break;
		if('@' == Tag_Type)//start of new tag
		{
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(Input_FileO);
				Average_Length=Current_Pos/Total_Tags+1;//+1 avoids divide by zero..
				Number_of_Tags=(File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/File_Size);
			}

			gzgets(Data_File,Description,MAXDES);
			gzread(Data_File,Stats1,Stat_Size);
			if(!NORMAL_TAGS) gzread(Data_File,Stats2,Stat_Size);
			gzread(Data_File,Tag_Copy,TAG_COPY_LEN);//if(NORMAL_TAGS && FILETYPE==FQ) gzread(Data_File,Quality,TAG_COPY_LEN);
			if(FILETYPE==FQ) gzread(Data_File,Quality,TAG_COPY_LEN);
			for(Next_String=Tag_Copy;Next_String[0]!='\n' && Next_String[0]!='\r' && Next_String[0]!=0;Next_String++) if(*Next_String=='\t') *Next_String=0; else Next_String[0]=tolower(Next_String[0]);
			Hits=0;
			if(NORMAL_TAGS) Total_Tags++;
			continue;
		}		
		else//Process hit....
		{
//--------------------------------------------------------------------------------------------
			memcpy(Translated_String,Tag_Copy,MAXSTRINGLENGTH);
			//Conversion_Factor=CONVERSION_FACTOR;
			gzread(Data_File,New_Record,3);
			if (New_Record[2]) gzread(Data_File,N,New_Record[2]*2);
			gzread(Data_File,&Record,sizeof(Output_Record));
			Mismatch_Count = Record.Mismatches;
			gzread(Data_File,&MismatchesGIS,Mismatch_Count+sizeof(unsigned));
			StringLength=Length_Array[New_Record[0]];
			Conversion_Factor=revfmi->textLength-StringLength;//+1;
			if (!NORMAL_TAGS)
			{
				if(!Hits)
				{
					Last_Half=New_Record[0];
					Total_Tags++;
				}
				else if (Last_Half != New_Record[0])
				{
					Last_Half=New_Record[0];
					Total_Tags++;
					Hits=0;
				}
			}
			if (SCANBOTH)
			{
				if(!Hits) Last_Orientation=New_Record[1];
				else if (Last_Orientation != New_Record[1]){Last_Orientation=New_Record[1];Hits=0;}
			}
//old one....
			if(New_Record[0]==1) {String=Translated_String+IGNOREHEAD;STRINGLENGTH=HEAD_LENGTH;}
			else {String=Translated_String+2*IGNOREHEAD+HEAD_LENGTH+1;STRINGLENGTH=TAIL_LENGTH;}

//--------------------------- Handle N's  ----------------------------------------------
			if (New_Record[2])
			{
				for(int i=0;i<New_Record[2]*2;i+=2)
				{
					String[N[i]]=N[i+1];
				}
			}
//--------------------------- Handle N's  ----------------------------------------------
			if (New_Record[1] == '-')
			{
				for (unsigned i=0;i<=STRINGLENGTH-1;i++) RevString[STRINGLENGTH-1-i]=Code_To_Char[Char_To_CodeC[String[i]]];
				String=RevString;
			}
//--------------------------- Decode Mismatches ----------------------------------------------
			if (Mismatch_Count)//mismatch positions..
			{
				for( int i=0;i<Mismatch_Count;i++) String[MismatchesGIS.Mismatch_Pos[i]]=Code_To_Char[MismatchesGIS.Mismatch_Char>>(2*i) & 3];
			}
			if(MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] == INDELMARK)//Indel...
			{
				pos=MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-3];
				if (MismatchesGIS.Mismatch_Pos[MAX_MISMATCHES_BOUND-2] == INSERTMARK)//insertion
				{
					Conversion_Factor--;StringLength++;
					for(int i=STRINGLENGTH+1;i>pos;i--) String[i]=String[i-1];
					String[pos]=Code_To_Char[MismatchesGIS.Mismatch_Char >> 2*2];

				}
				else//deletion
				{
					Conversion_Factor++;StringLength--;
					for(int i=pos;i<STRINGLENGTH;i++) String[i]=String[i+1];
				}
			}

			
//--------------------------- Decode Mismatches ----------------------------------------------
			for (int j=0;j<=Record.Gap && Hits< MAXHITS;j++)
			{
				if( Record.Index)//print record...
				{
					if (Record.Skip) Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
					else Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
				}
				else
				{
					if (Record.Skip) Location = fmi->saValue[Record.Start/fmi->saInterval]-Record.Skip+1;
					else Location=BWTSaValue(fmi,Record.Start);
				}
				
				Location -= Offset;
				Hits++;Record.Start++;Total_Hits++;
				if(VERIFY_FROM_DISK)
				{
					fseeko64(Binary_File,Location,SEEK_SET);
					fread(&Test_String,StringLength,1,Binary_File);
					if (strncmp(Test_String,String,StringLength))
					{
						//for (unsigned i=0;i<=StringLength-1;i++) String[i]=Code_To_Char[Char_To_CodeC[String[i]]];
						//if (strncmp(Test_String,String,StringLength)) 
						printf("Error in Tag %d\t%s \n: Location %u is\t%s\t%s\n",Record.Tag, String,Location,Test_String,Description);
					}
				}
				else
				{

					for (unsigned i=0;i<=StringLength-1;i++) Bin_String[i]=Char_To_Code[String[i]];
					Current_String=Bin_String;

					Tag.Start=0;Tag.End=SOURCELENGTH;Tag.Mismatches=0;Tag.Level=1;
					Search_Exact(Tag,-1,StringLength,revfmi);	
					if(!Tag.Start)
					{
						printf("Error in index of Tag %d\t%s\t%s\n:",Record.Tag, String,Description);
						fseeko64(Binary_File,Location,SEEK_SET);
						fread(&Test_String,StringLength,1,Binary_File);
						if (strncmp(Test_String,String,StringLength)) printf("Error in Tag %d\t%s \n: Location %u is\t%s\t%s\n",Record.Tag, String,Location,Test_String,Description);
					}
				}
			}
//--------------------------- Decode Mismatches ----------------------------------------------
		}

	}
}
//}-----------------------------  GIS Mode  -------------------------------------------------



//{-----------------------------  Verify Deep  -------------------------------------------------
void Verify_Deep()
{
	unsigned i, Bitsread;
	unsigned Start,End;
	unsigned Previous_Tag,Hits;
	unsigned Location;
	unsigned pos;
	SARange Tag;
	int Conversion_Factor,StringLength;

	unsigned Progress=0,Average_Length;
	unsigned Number_of_Tags=100;

	printf ("Verifying using %s\n", BINFILE);
	printf("======================]\r[");//progress bar....
	Binary_File=fopen(BINFILE,"rb");
	Test_String[STRINGLENGTH]=0;
	for(;;)
	{
		i=0;
		gzread(Data_File,&Tag_Type,1);
		if('&' == Tag_Type) break;
		if('@' == Tag_Type)//start of new tag
		{
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(Input_FileO);
				Average_Length=Current_Pos/Total_Tags+1;//+1 avoids divide by zero..
				Number_of_Tags=(File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/File_Size);
			}

			if(PRINT_DESC) gzgets(Data_File,Description,1000);
			gzread(Data_File,Translated_String_Original,STRINGLENGTH+1);
			Hits=0;Total_Tags++;
			continue;
		}		
		else
		{

			memcpy(Translated_String,Translated_String_Original,MAXSTRINGLENGTH);
			Conversion_Factor=CONVERSION_FACTOR;
			StringLength=STRINGLENGTH;
			gzread(Data_File,&Mismatches,sizeof(Mismatches_Record));
			gzread(Data_File,&Record,sizeof(Output_Record));

			Mismatch_Count = Record.Mismatches;

			if (Mismatch_Count)//mismatch positions..
			{
				for( int i=0;i<Mismatch_Count;i++)
				{
					Translated_String[Mismatches.Mismatch_Pos[i]]=Code_To_Char[Mismatches.Mismatch_Char>>(2*i) & 3];
				}
			}

			if(Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] == INDELMARK)//Indel...
			{
				pos=Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-3];
				if (Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-2] == INSERTMARK)//insertion
				{
					Conversion_Factor--;StringLength++;
					for(int i=STRINGLENGTH+1;i>pos;i--) 
					{
						Translated_String[i]=Translated_String[i-1];
					}
					Translated_String[pos]=Code_To_Char[Mismatches.Mismatch_Char >> 2*2];
				}
				else//deletion
				{
					Conversion_Factor++;StringLength--;
					for(int i=pos;i<STRINGLENGTH;i++) 
					{
						Translated_String[i]=Translated_String[i+1];
					}
				}
			}		

			for (int j=0;j<=Mismatches.Gap && Hits< MAXHITS;j++)
			{
				if( Record.Index)//print record...
				{
					if (Record.Skip)
					{
						Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
					}
					else
					{
						Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
					}
				}
				else
				{
					if (Record.Skip)
					{
						Location = fmi->saValue[Record.Start/fmi->saInterval]-Record.Skip+1;
					}
					else
					{
						Location=BWTSaValue(fmi,Record.Start);
					}
				}
				//Location_To_Genome(Location);
				if(VERIFY_FROM_DISK)
				{
					fseeko64(Binary_File,Location,SEEK_SET);
					fread(&Test_String,StringLength,1,Binary_File);
					if (strncmp(Test_String,Translated_String,StringLength))
					{
						for (unsigned i=0;i<=StringLength-1;i++)
						{
							Translated_String[i]=Code_To_Char[Char_To_CodeC[Translated_String[i]]];
						}

						if (strncmp(Test_String,Translated_String,StringLength))
						{
							printf("Error in Tag %d\t%s \n: Location %u is\t%s\t%s\n",Record.Tag, Translated_String,Location,Test_String,Description);
						}
					}
					Hits++;Record.Start++;Total_Hits++;
				}
				else
				{
					for (unsigned i=0;i<=StringLength-1;i++)
					{
						Bin_String[i]=Char_To_Code[Translated_String[i]];
						Bin_StringComplement[StringLength-1-i]=Char_To_CodeC[Translated_String[i]];//change later...
					}
					Tag.Start=0;Tag.End=SOURCELENGTH;Tag.Mismatches=0;Tag.Level=1;
					Current_String=Bin_String;
					Search_Exact(Tag,-1,StringLength,revfmi);	
					if(!Tag.Start)
					{
						Tag.Start=0;Tag.End=SOURCELENGTH;Tag.Mismatches=0;Tag.Level=1;
						Current_String=Bin_StringComplement;
						Search_Exact(Tag,-1,StringLength,revfmi);	

						if(!Tag.Start)
						{
							printf("Error in index of Tag %d\t%s\t%s\n:",Record.Tag, Translated_String,Description);
							fseeko64(Binary_File,Location,SEEK_SET);
							fread(&Test_String,StringLength,1,Binary_File);
							if (strncmp(Test_String,Translated_String,StringLength))
							{
								for (unsigned i=0;i<=StringLength-1;i++)
								{
									Translated_String[i]=Code_To_Char[Char_To_CodeC[Translated_String[i]]];
								}

								if (strncmp(Test_String,Translated_String,StringLength))
								{
									printf("Error in Tag %d\t%s \n: Location %u is\t%s\t%s\n",Record.Tag, Translated_String,Location,Test_String,Description);
								}
							}
						}
					}
					Hits++;Record.Start++;Total_Hits++;
				}
			}

		}

	}
}

void Search_Exact(struct SARange & Tag,int Start,int StringLength,BWT *fmi)
{
	int Level;
	unsigned long Index,Now,First,Last;
	//if (!Tag.Start) return;

	for(;;)	
	{
		Get_SARange_Fast(Current_String[Start+Tag.Level],Tag,fmi);

//printf("[%u:%u]- %d\n",Tag.Start,Tag.End,Tag.Level);
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
//}-----------------------------  Verify Deep  -------------------------------------------------

//{-----------------------------  Process pairend  -------------------------------------------------
void Verify_Pairend()
{
	unsigned i, Bitsread;
	unsigned Start,End;
	unsigned Location;
	int Desc_End;
	unsigned Previous_Tag;
	unsigned Hits=0;
	char letter;
	SARange Tag;
	unsigned pos;
	int Conversion_Factor,StringLength;
	int HEAD_LENGTH,TAIL_LENGTH;
	char Layout[2];//head/tail, orientation...
	int Genome_Length=revfmi->textLength;	
	char Pass,PM;
	char* Next_String;
	char* String;

	unsigned Progress=0,Average_Length;
	unsigned Number_of_Tags=100;

	printf ("Verifying using %s\n", BINFILE);
	printf("======================]\r[");//progress bar....
	Binary_File=fopen(BINFILE,"rb");
	gzread(Data_File,&HEAD_LENGTH,sizeof(int));
	gzread(Data_File,&TAIL_LENGTH,sizeof(int));
	for(;;)
	{
		i=0;
		gzread(Data_File,&Tag_Type,1);
		if('&' == Tag_Type) break;
		if('@' == Tag_Type)//start of new tag
		{
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(Input_FileO);
				Average_Length=Current_Pos/Total_Tags+1;//+1 avoids divide by zero..
				Number_of_Tags=(File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/File_Size);
			}

			if(PRINT_DESC) 
			{
				gzgets(Data_File,Description,1000);
				for(Desc_End=0;Description[Desc_End]!='\n' && Description[Desc_End]!='\r' && Description[Desc_End]!=0;Desc_End++);
				Description[Desc_End]=0;
			};
			//Hits=0;
			Total_Tags++;
			
			gzgets(Data_File,Translated_String_Original,MAXSTRINGLENGTH);
			for(Next_String=Translated_String_Original;Next_String[0]!='\n' && Next_String[0]!='\r' && Next_String[0]!=0;Next_String++) Next_String[0]=tolower(Next_String[0]);
			Pass=1;PM='+';
			continue;
		}		
		else
		{
			memcpy(Translated_String,Translated_String_Original,MAXSTRINGLENGTH);
			gzread(Data_File,&Layout,2);
			if(Layout[0]=='H')
			{
				String=Translated_String;STRINGLENGTH=HEAD_LENGTH;
			} 
			else
			{
				String=Translated_String+HEAD_LENGTH+1;STRINGLENGTH=TAIL_LENGTH;

				if (Pass==1) 
				{
					Hits=0;Pass=0;
				}
			}
			if (Layout[1] == '-')
			{
				for (unsigned i=0;i<=STRINGLENGTH-1;i++)
				{
					RevString[STRINGLENGTH-1-i]=Code_To_Char[Char_To_CodeC[String[i]]];
				}
				String=RevString;
				if (PM == '+')
				{
					Hits=0;PM='-';
				}
			}
			else
			{
				if (PM == '-')
				{
					Hits=0;PM='+';
				}
			}


			Conversion_Factor=Genome_Length-STRINGLENGTH;
			StringLength=STRINGLENGTH;
			gzread(Data_File,&Mismatches,sizeof(Mismatches_Record));
			gzread(Data_File,&Record,sizeof(Output_Record));
			Mismatch_Count = Record.Mismatches;

			if (Mismatch_Count)//mismatch positions..
			{
				for( int i=0;i<Mismatch_Count;i++)
				{
					String[Mismatches.Mismatch_Pos[i]]=Code_To_Char[Mismatches.Mismatch_Char>>(2*i) & 3];
				}
			}

			if(Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] == INDELMARK)//Indel...
			{
				pos=Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-3];
				if (Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-2] == INSERTMARK)//insertion
				{
					Conversion_Factor--;StringLength++;
					for(int i=STRINGLENGTH+1;i>pos;i--) 
					{
						String[i]=String[i-1];
					}
					String[pos]=Code_To_Char[Mismatches.Mismatch_Char >> 2*2];
				}
				else//deletion
				{
					Conversion_Factor++;StringLength--;
					for(int i=pos;i<STRINGLENGTH;i++) 
					{
						String[i]=String[i+1];
					}
				}
			}

			for (int j=0;j<=Mismatches.Gap && Hits< MAXHITS;j++)
			{
				if( Record.Index)//print record...
				{
					if (Record.Skip)
					{
						Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
					}
					else
					{
						Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
					}
				}
				else
				{
					if (Record.Skip)
					{
						Location = fmi->saValue[Record.Start/fmi->saInterval]-Record.Skip+1;
					}
					else
					{
						Location=BWTSaValue(fmi,Record.Start);
					}
				}

				if(VERIFY_FROM_DISK)
				{
					fseeko64(Binary_File,Location,SEEK_SET);
					fread(&Test_String,StringLength,1,Binary_File);
					if (strncmp(Test_String,Translated_String,StringLength))
					{
						for (unsigned i=0;i<=StringLength-1;i++)
						{
							Translated_String[i]=Code_To_Char[Char_To_CodeC[Translated_String[i]]];
						}

						if (strncmp(Test_String,Translated_String,StringLength))
						{
							printf("Error in Tag %d\t%s \n: Location %u is\t%s\t%s\n",Record.Tag, Translated_String,Location,Test_String,Description);
						}
					}
					Hits++;Record.Start++;Total_Hits++;
				}
				else
				{
					for (unsigned i=0;i<=StringLength-1;i++)
					{
						Bin_String[i]=Char_To_Code[String[i]];
						Bin_StringComplement[StringLength-1-i]=Char_To_CodeC[String[i]];//change later...
					}
					Tag.Start=0;Tag.End=SOURCELENGTH;Tag.Mismatches=0;Tag.Level=1;
					Current_String=Bin_String;
					Search_Exact(Tag,-1,StringLength,revfmi);	
					if(!Tag.Start)
					{
						Tag.Start=0;Tag.End=SOURCELENGTH;Tag.Mismatches=0;Tag.Level=1;
						Current_String=Bin_StringComplement;
						Search_Exact(Tag,-1,StringLength,revfmi);	

						if(!Tag.Start)
						{
							printf("Error in index of Tag %d\t%s\t%s\n:",Record.Tag, Translated_String,Description);
							fseeko64(Binary_File,Location,SEEK_SET);
							fread(&Test_String,StringLength,1,Binary_File);
							if (strncmp(Test_String,String,StringLength))
							{
								for (unsigned i=0;i<=StringLength-1;i++)
								{
									String[i]=Code_To_Char[Char_To_CodeC[String[i]]];
								}

								if (strncmp(Test_String,RevString,StringLength))
								{
									printf("Error in Tag %d\t%s \n: Location %u is\t%s\t%s\n",Record.Tag, String,Location,Test_String,Description);
								}
							}
						}
					}
				}
				Hits++;Record.Start++;Total_Hits++;
			}

		}
	}
}
//}-----------------------------  Process pairend  -------------------------------------------------

//{-----------------------------  Process Deep  -------------------------------------------------
void Process_Deep()
{
	unsigned i, Bitsread;
	unsigned Start,End;
	unsigned Location;
	int Desc_End;
	unsigned Previous_Tag,Hits;
	char letter;
	unsigned pos;
	int Conversion_Factor,StringLength;

	unsigned Progress=0,Average_Length;
	unsigned Number_of_Tags=100;
	printf("======================]\r[");//progress bar....

	for(;;)
	{
		i=0;
		gzread(Data_File,&Tag_Type,1);
		if('&' == Tag_Type) break;
		if('@' == Tag_Type)//start of new tag
		{
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(Input_FileO);
				Average_Length=Current_Pos/Total_Tags+1;//+1 avoids divide by zero..
				Number_of_Tags=(File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/File_Size);
			}

			fprintfX(Output_File,"@\n");
			//fprintf(Output_File,"@\n");
			if(PRINT_DESC) 
			{
				gzgets(Data_File,Description,1000);
				for(Desc_End=0;Description[Desc_End]!='\n' && Description[Desc_End]!='\r' && Description[Desc_End]!=0;Desc_End++);
				Description[Desc_End]=0;
				fprintfX(Output_File,"%s ",Description);
				//fprintf(Output_File,"%s ",Description);
			};
			gzread(Data_File,Translated_String,STRINGLENGTH+1);
			//fprintf(Output_File,"%s",Translated_String);Hits=0;Total_Tags++;
			fprintfX(Output_File,"%s",Translated_String);Hits=0;Total_Tags++;
			continue;
		}		
		else
		{
			//modify out...	
			Conversion_Factor=CONVERSION_FACTOR;
			StringLength=STRINGLENGTH;
			//Offset=0;
			gzread(Data_File,&Mismatches,sizeof(Mismatches_Record));
			gzread(Data_File,&Record,sizeof(Output_Record));
			Mismatch_Count = Record.Mismatches;

			//fprintf(Output_File,"\t%u %d ",Record.Tag,Mismatch_Count);
			fprintfX(Output_File,"\t%u %d ",Record.Tag,Mismatch_Count);

			if(Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] == INDELMARK)//Indel...
			{
				pos=Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-3];
				if (Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-2] == INSERTMARK)//insertion
				{
					Conversion_Factor--;StringLength++;
					//fprintf(Output_File," %d<%c ",pos,Code_To_Char[Mismatches.Mismatch_Char >> 2*2]-('a'-'A'));
					fprintfX(Output_File," %d<%c ",pos,Code_To_Char[Mismatches.Mismatch_Char >> 2*2]-('a'-'A'));
				}
				else//deletion
				{
					Conversion_Factor++;StringLength--;
					//fprintf(Output_File," %d>D ",pos);
					fprintfX(Output_File," %d>D ",pos);
				}
			}

			if (Mismatch_Count)//mismatch positions..
			{
				for( int i=0;i<Mismatch_Count;i++)
				{
					pos=Mismatches.Mismatch_Pos[i];
					letter=Mismatches.Mismatch_Char>>(2*i) & 3;
					//Translated_String[pos]=Code_To_Char[letter];
					//fprintf(Output_File," %d>%c ",pos,Code_To_Char[letter]);
					fprintfX(Output_File," %d>%c ",pos,Code_To_Char[letter]);
				}
			}

			//fprintf(Output_File," "); //seperator...
			fprintfX(Output_File," "); //seperator...
			for (int j=0;j<=Mismatches.Gap && Hits< MAXHITS;j++)
			{
				if( Record.Index)//print record...
				{
					if (Record.Skip)
					{
						Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
					}
					else
					{
						Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
					}
				}
				else
				{
					if (Record.Skip)
					{
						Location = fmi->saValue[Record.Start/fmi->saInterval]-Record.Skip+1;
					}
					else
					{
						Location=BWTSaValue(fmi,Record.Start);
					}
				}
				
				Location -= Offset;
				if (!PLUSSTRAND && Translated_String[STRINGLENGTH]=='-') Location=Location+STRINGLENGTH-PLUSGAP;
				if(USELOCATION) {Location_To_Genome(Location);fprintfX(Output_File,"%u;[%s]",Location,Genome_Offsets[Genome_Position].Genome);}
				else fprintfX(Output_File,"%u;",Location); 

				Hits++;Record.Start++;Total_Hits++;
			}
			//fprintf(Output_File,"\n");
			fprintfX(Output_File,"\n");

		}

	}
}
//}-----------------------------  Process Deep  -------------------------------------------------

//{-----------------------------  Process Deep Formatted  -------------------------------------------------
void Process_Deep_Formatted()
{
	unsigned i, Bitsread;
	unsigned Start,End;
	unsigned Location;
	int Desc_End;
	unsigned Previous_Tag,Hits;
	char letter;
	unsigned pos;
	char Init=TRUE;
	int Conversion_Factor,StringLength;

	unsigned Progress=0,Average_Length;
	unsigned Number_of_Tags=100;
	printf("======================]\r[");//progress bar....

	Total_Hits=0;
	for(;;)
	{
		i=0;
		gzread(Data_File,&Tag_Type,1);
		if('&' == Tag_Type) break;
		if('@' == Tag_Type)//start of new tag
		{
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(Input_FileO);
				Average_Length=Current_Pos/Total_Tags+1;//+1 avoids divide by zero..
				Number_of_Tags=(File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/File_Size);
			}
			//if(!Init){fprintf(Output_File,"@\n");} else {fprintf(Output_File,"@\t%d\n",STRINGLENGTH);Init=FALSE;}
			if(!Init){fprintfX(Output_File,"@\n");} else {fprintfX(Output_File,"@\t%d\n",STRINGLENGTH);Init=FALSE;}
			if(PRINT_DESC) 
			{
				gzgets(Data_File,Description,1000);
				for(Desc_End=0;Description[Desc_End]!='\n' && Description[Desc_End]!='\r' && Description[Desc_End]!=0;Desc_End++);
				Description[Desc_End]=0;
				//fprintf(Output_File,"%s",Description);
				fprintfX(Output_File,"%s",Description);
			};
			gzread(Data_File,Translated_String,STRINGLENGTH+1);
			//fprintf(Output_File,"\t%s\n",Translated_String);Hits=0;Total_Tags++;
			fprintfX(Output_File,"\t%s\n",Translated_String);Hits=0;Total_Tags++;
			continue;
		}		
		else
		{
		//modify out...	
			Conversion_Factor=CONVERSION_FACTOR;
			StringLength=STRINGLENGTH;
			//Offset=0;
			gzread(Data_File,&Mismatches,sizeof(Mismatches_Record));
			gzread(Data_File,&Record,sizeof(Output_Record));
			Mismatch_Count = Record.Mismatches;
			
			//fprintf(Output_File,"%u\t%d\t",Record.Tag,Mismatch_Count);
			fprintfX(Output_File,"%u\t%d\t",Record.Tag,Mismatch_Count);

			if(Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] == INDELMARK)//Indel...
			{
				pos=Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-3];
				if (Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-2] == INSERTMARK)//insertion
				{
					Conversion_Factor--;StringLength++;
					//fprintf(Output_File,"%d<%c\t",pos,Code_To_Char[Mismatches.Mismatch_Char >> 2*2]-('a'-'A'));
					fprintfX(Output_File,"%d<%c\t",pos,Code_To_Char[Mismatches.Mismatch_Char >> 2*2]-('a'-'A'));
				}
				else//deletion
				{
					Conversion_Factor++;StringLength--;
					//fprintf(Output_File,"%d>D\t",pos);
					fprintfX(Output_File,"%d>D\t",pos);
				}
			}

			if (Mismatch_Count)//mismatch positions..
			{
				for( int i=0;i<Mismatch_Count;i++)
				{
					pos=Mismatches.Mismatch_Pos[i];
					letter=Mismatches.Mismatch_Char>>(2*i) & 3;
					//Translated_String[pos]=Code_To_Char[letter];
					fprintfX(Output_File,"%d>%c\t",pos,Code_To_Char[letter]);
					//fprintf(Output_File,"%d>%c\t",pos,Code_To_Char[letter]);
				}
			}

			//fprintf(Output_File,"|%d\t",Mismatches.Gap+1); //seperator...
			fprintfX(Output_File,"|%d\t",Mismatches.Gap+1); //seperator...
			for (int j=0;j<=Mismatches.Gap && Hits< MAXHITS;j++)
			{
				if( Record.Index)//print record...
				{
					if (Record.Skip)
					{
						Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
					}
					else
					{
						Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
					}
				}
				else
				{
					if (Record.Skip)
					{
						Location = fmi->saValue[Record.Start/fmi->saInterval]-Record.Skip+1;
					}
					else
					{
						Location=BWTSaValue(fmi,Record.Start);
					}
				}

				Location -= Offset;
				if (!PLUSSTRAND && Translated_String[STRINGLENGTH]=='-') Location=Location+STRINGLENGTH-PLUSGAP;
				if(USELOCATION) {Location_To_Genome(Location);fprintfX(Output_File,"%s:%u;",Genome_Offsets[Genome_Position].Genome,Location);}
				else fprintfX(Output_File,"%u;",Location);

				Hits++;Record.Start++;Total_Hits++;
			}
			//fprintf(Output_File,"\n");
			fprintfX(Output_File,"\n");
		}

	}
}
//}-----------------------------  Process Deep  -------------------------------------------------

//{-----------------------------  Process pairend  -------------------------------------------------
void Process_Pairend()
{
	unsigned i, Bitsread;
	unsigned Start,End;
	unsigned Location;
	int Desc_End;
	unsigned Previous_Tag,Hits;
	char letter;
	unsigned pos;
	int Conversion_Factor,StringLength;
	int HEAD_LENGTH,TAIL_LENGTH;
	char Layout[2];//head/tail, orientation...
	int Genome_Length=revfmi->textLength;	
	char Pass;
	char PM;

	unsigned Progress=0,Average_Length;
	unsigned Number_of_Tags=100;
	//if (MAXHITS==1) MAXHITS=2;
	//MAXHITS=2*MAXHITS;
	gzread(Data_File,&HEAD_LENGTH,sizeof(int));
	gzread(Data_File,&TAIL_LENGTH,sizeof(int));

	printf("======================]\r[");//progress bar....
	fflush(stdout);
	for(;;)
	{
		i=0;
		gzread(Data_File,&Tag_Type,1);
		if('&' == Tag_Type) break;
		if('@' == Tag_Type)//start of new tag
		{
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(Input_FileO);
				Average_Length=Current_Pos/Total_Tags+1;//+1 avoids divide by zero..
				Number_of_Tags=(File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/File_Size);
			}

			//fprintf(Output_File,"@\n");
			fprintfX(Output_File,"@\n");
			if(PRINT_DESC) 
			{
				gzgets(Data_File,Description,1000);
				for(Desc_End=0;Description[Desc_End]!='\n' && Description[Desc_End]!='\r' && Description[Desc_End]!=0;Desc_End++);
				Description[Desc_End]=0;
				//fprintf(Output_File,"%s\t",Description);
				fprintfX(Output_File,"%s\t",Description);
			};
			gzgets(Data_File,Translated_String,MAXSTRINGLENGTH);
			//fprintf(Output_File,"%s",Translated_String);Hits=0;Total_Tags++;
			fprintfX(Output_File,"%s",Translated_String);Hits=0;Total_Tags++;
			Pass=1;PM='+';
			continue;
		}		
		else
		{
		//modify out...	
			gzread(Data_File,&Layout,2);

			if(Layout[0]=='H') 
			{
				STRINGLENGTH=HEAD_LENGTH;
			}
			else 
			{
				STRINGLENGTH=TAIL_LENGTH;
				if (Pass==1) 
				{
					Hits=0;
					Pass=0;
				}
			}

			if(Layout[1]=='-')
			{
				if (PM == '+')
				{
					//Hits=0;PM='-';
					PM='-';if (!ROLLOVER) Hits=0;
					if(SCANBOTH) Hits=0;
				}
			}
			else
			{
				if (PM == '-')
				{
					PM='+';if (!ROLLOVER) Hits=0;
					if(SCANBOTH) Hits=0;
				}
			}

			Conversion_Factor=Genome_Length-STRINGLENGTH;
			//printf("H %d ", STRINGLENGTH);
			StringLength=STRINGLENGTH;
			//Offset=0;
			gzread(Data_File,&Mismatches,sizeof(Mismatches_Record));
			gzread(Data_File,&Record,sizeof(Output_Record));
			Mismatch_Count = Record.Mismatches;
			
			/////fprintf(Output_File,"\t%u %d ",Record.Tag,Mismatch_Count);
			fprintfX(Output_File,"%c\t%c\t",Layout[0],Layout[1]);//,Record.Tag,Mismatch_Count);
			//fprintf(Output_File,"%c\t%c\t",Layout[0],Layout[1]);//,Record.Tag,Mismatch_Count);

			if(Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] == INDELMARK)//Indel...
			{
				pos=Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-3];
				if (Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-2] == INSERTMARK)//insertion
				{
					Conversion_Factor--;StringLength++;
					//fprintf(Output_File,"%d<%c\t",pos,Code_To_Char[Mismatches.Mismatch_Char >> 2*2]-('a'-'A'));
					fprintfX(Output_File,"%d<%c\t",pos,Code_To_Char[Mismatches.Mismatch_Char >> 2*2]-('a'-'A'));
				}
				else//deletion
				{
					Conversion_Factor++;StringLength--;
					//fprintf(Output_File,"%d>D\t",pos);
					fprintfX(Output_File,"%d>D\t",pos);
				}
			}

			fprintfX(Output_File," "); //seperator...
			//fprintf(Output_File," "); //seperator...
			for (int j=0;j<=Mismatches.Gap && Hits< MAXHITS;j++)
			{
				if( Record.Index)//print record...
				{
					if (Record.Skip)
					{
						Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
					}
					else
					{
						Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
					}
				}
				else
				{
					if (Record.Skip)
					{
						Location = fmi->saValue[Record.Start/fmi->saInterval]-Record.Skip+1;
					}
					else
					{
						Location=BWTSaValue(fmi,Record.Start);
					}
				}

				Location -= Offset;
				if (!PLUSSTRAND && Translated_String[STRINGLENGTH]=='-') Location=Location+STRINGLENGTH-PLUSGAP;
				if(USELOCATION) {Location_To_Genome(Location);fprintfX(Output_File,"%s:%u;",Genome_Offsets[Genome_Position].Genome,Location);}
				else fprintfX(Output_File,"%u;",Location);

				Hits++;Record.Start++;Total_Hits++;
			}
			//fprintf(Output_File,"\t%d\t%d\n",Mismatch_Count,STRINGLENGTH);
			fprintfX(Output_File,"\t%d\t%d\n",Mismatch_Count,STRINGLENGTH);

		}

	}
}
//}-----------------------------  Process pairend  -------------------------------------------------
//formatted
//{-----------------------------  Process pairend Formatted  -------------------------------------------------
void Process_Pairend_Formatted()
{
	unsigned i, Bitsread;
	unsigned Start,End;
	unsigned Location;
	int Desc_End;
	unsigned Previous_Tag,Hits;
	char letter;
	unsigned pos;
	int Conversion_Factor,StringLength;
	int HEAD_LENGTH,TAIL_LENGTH;
	char Layout[2];//head/tail, orientation...
	char Pass,PM;
	char Init=TRUE;
	int Genome_Length=revfmi->textLength;	
	unsigned Progress=0,Average_Length;
	unsigned Number_of_Tags=100;

	//if (MAXHITS==1) MAXHITS=2;
	//MAXHITS=2*MAXHITS;
	gzread(Data_File,&HEAD_LENGTH,sizeof(int));
	gzread(Data_File,&TAIL_LENGTH,sizeof(int));

	printf("======================]\r[");//progress bar....
	fflush(stdout);

	for(;;)
	{
		
		i=0;
		gzread(Data_File,&Tag_Type,1);
		if('&' == Tag_Type) break;
		if('@' == Tag_Type)//start of new tag
		{
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(Input_FileO);
				Average_Length=Current_Pos/Total_Tags+1;//+1 avoids divide by zero..
				Number_of_Tags=(File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/File_Size);
			}

			//if(!Init){fprintf(Output_File,"@\n");} else {fprintf(Output_File,"@\t%d\t%d\n",HEAD_LENGTH,TAIL_LENGTH);Init=FALSE;}
			if(!Init){fprintfX(Output_File,"@\n");} else {fprintfX(Output_File,"@\t%d\t%d\n",HEAD_LENGTH,TAIL_LENGTH);Init=FALSE;}
			if(PRINT_DESC) 
			{
				gzgets(Data_File,Description,1000);
				for(Desc_End=0;Description[Desc_End]!='\n' && Description[Desc_End]!='\r' && Description[Desc_End]!=0;Desc_End++);
				Description[Desc_End]=0;
				//fprintf(Output_File,"%s",Description);
				fprintfX(Output_File,"%s",Description);
			};
			gzgets(Data_File,Translated_String,MAXSTRINGLENGTH);
			//fprintf(Output_File,"\t%s\n",Translated_String);Hits=0;Total_Tags++;
			fprintfX(Output_File,"\t%s\n",Translated_String);Hits=0;Total_Tags++;
			Pass=1;PM='+';
			continue;
		}		
		else
		{
		//modify out...	
			gzread(Data_File,&Layout,2);
			if(Layout[0]=='H') 
			{
				STRINGLENGTH=HEAD_LENGTH;
			}
			else 
			{
				STRINGLENGTH=TAIL_LENGTH;
				if (Pass==1) 
				{
					Hits=0;
					Pass=0;
				}
			}

			if(Layout[1]=='-')
			{
				if (PM == '+')
				{
					PM='-';if (!ROLLOVER) Hits=0;
					if(SCANBOTH) Hits=0;
				}
			}
			else
			{
				if (PM == '-')
				{
					PM='+';if (!ROLLOVER) Hits=0;
					if(SCANBOTH) Hits=0;
				}
			}
			//if(Layout[0]=='H') STRINGLENGTH=HEAD_LENGTH; else {STRINGLENGTH=TAIL_LENGTH;if (Pass==1) {Hits=0;Pass=0;}}
			Conversion_Factor=Genome_Length-STRINGLENGTH;
			//printf("H %d ", STRINGLENGTH);
			StringLength=STRINGLENGTH;
			//Offset=0;
			gzread(Data_File,&Mismatches,sizeof(Mismatches_Record));
			gzread(Data_File,&Record,sizeof(Output_Record));
			Mismatch_Count = Record.Mismatches;
			
				fprintfX(Output_File,"%u\t%d\t",Record.Tag,Mismatch_Count);
				//fprintf(Output_File,"%u\t%d\t",Record.Tag,Mismatch_Count);
				fprintfX(Output_File,"%c\t%c\t",Layout[0],Layout[1]);//,Record.Tag,Mismatch_Count);
				//fprintf(Output_File,"%c\t%c\t",Layout[0],Layout[1]);//,Record.Tag,Mismatch_Count);

				if(Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-1] == INDELMARK)//Indel...
				{
					pos=Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-3];
					if (Mismatches.Mismatch_Pos[MAX_MISMATCHES_BOUND-2] == INSERTMARK)//insertion
					{
						Conversion_Factor--;StringLength++;
						//fprintf(Output_File,"%d<%c\t",pos,Code_To_Char[Mismatches.Mismatch_Char >> 2*2]-('a'-'A'));
						fprintfX(Output_File,"%d<%c\t",pos,Code_To_Char[Mismatches.Mismatch_Char >> 2*2]-('a'-'A'));
					}
					else//deletion
					{
						Conversion_Factor++;StringLength--;
						//fprintf(Output_File,"%d>D\t",pos);
						fprintfX(Output_File,"%d>D\t",pos);
					}
				}
				

				if (Mismatch_Count)//mismatch positions..
				{
					for( int i=0;i<Mismatch_Count;i++)
					{
						pos=Mismatches.Mismatch_Pos[i];
						letter=Mismatches.Mismatch_Char>>(2*i) & 3;
						//Translated_String[pos]=Code_To_Char[letter];
						//fprintf(Output_File,"%d>%c\t",pos,Code_To_Char[letter]);
						fprintfX(Output_File,"%d>%c\t",pos,Code_To_Char[letter]);
					}
				}
				
				/////fprintf(Output_File," "); //seperator...
				fprintfX(Output_File,"|%d\t",Mismatches.Gap+1); //seperator...
				//fprintf(Output_File,"|%d\t",Mismatches.Gap+1); //seperator...
				for (int j=0;j<=Mismatches.Gap && Hits< MAXHITS;j++)
				{
					if( Record.Index)//print record...
					{
						if (Record.Skip)
						{
							Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
						}
						else
						{
							Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
						}
					}
					else
					{
						if (Record.Skip)
						{
							Location = fmi->saValue[Record.Start/fmi->saInterval]-Record.Skip+1;
						}
						else
						{
							Location=BWTSaValue(fmi,Record.Start);
						}
					}
					Location -= Offset;
					if (!PLUSSTRAND && Translated_String[STRINGLENGTH]=='-') Location=Location+STRINGLENGTH-PLUSGAP;
					if(USELOCATION) {Location_To_Genome(Location);fprintfX(Output_File,"%s:%u;",Genome_Offsets[Genome_Position].Genome,Location);}
					else fprintfX(Output_File,"%u;",Location);

					Hits++;Record.Start++;Total_Hits++;
				}
				/////fprintf(Output_File,"\t%d\t%d\n",Mismatch_Count,STRINGLENGTH);
				fprintfX(Output_File,"\n");
				//fprintf(Output_File,"\n");

		}

	}
}
//}-----------------------------  Process pairend  -------------------------------------------------


//{-----------------------------  Process Default  -------------------------------------------------
void Process_Default()
{
	unsigned i, Bitsread;
	unsigned Start,End,Location;
	unsigned Hitcount,Last_Tag;
	char PlusMinus;
	char Indel=0;
	int Conversion_Factor;

	Hits=(Output_Record*)malloc(sizeof(Output_Record)*BUFFERSIZE+1);
	Hits[BUFFERSIZE].Start=0;
	Hitcount=0;
	for(;;)
	{
		i=0;
		Bitsread=gzread(Data_File,Hits,sizeof(Output_Record)*BUFFERSIZE);
		while(Hits[i].Start)//while not sentinel...
		{
			if (Last_Tag != Hits[i].Tag) Hitcount=0;//new tag, new count...
			Last_Tag=Hits[i].Tag;
			for (int j=0;j<=Hits[i].Gap && Hitcount< MAXHITS;j++)
			{
				Indel='M';
				Conversion_Factor=CONVERSION_FACTOR;
				if (Hits[i].Mismatches>=100) {Location=Location+STRINGLENGTH;Hits[i].Mismatches=Hits[i].Mismatches-100;PlusMinus='-';} else PlusMinus='+';
				if (Hits[i].Mismatches>=75) {Conversion_Factor--;Hits[i].Mismatches=Hits[i].Mismatches-75;Indel='I';} 
				else if (Hits[i].Mismatches>=50) {Conversion_Factor++;Hits[i].Mismatches=Hits[i].Mismatches-50;Indel='D';}

				if( Hits[i].Index)//print record...
				{
					if (Record.Skip)
					{
						Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
					}
					else
					{
						Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
					}
				}
				else
				{
					if (Record.Skip)
					{
						Location = fmi->saValue[Record.Start/fmi->saInterval]-Record.Skip+1;
					}
					else
					{
						Location=BWTSaValue(fmi,Record.Start);
					}
				}
				Location -= Offset;
				if (USELOCATION) Location_To_Genome(Location);

				//fprintf(Output_File,"%u \t %u \t %c%d \t %c[%s] \n ",Hits[i].Tag,Location,Indel,Hits[i].Mismatches,PlusMinus,Genome_Offsets[Genome_Position].Genome);
				fprintfX(Output_File,"%u \t %u \t %c%d \t %c[%s] \n ",Hits[i].Tag,Location,Indel,Hits[i].Mismatches,PlusMinus,Genome_Offsets[Genome_Position].Genome);
				Total_Hits++;

				Hitcount++;Hits[i].Start++;
			}
			i++;
		}
		if(Bitsread<BUFFERSIZE) break;
	}
}
//}-----------------------------  Process Default  -------------------------------------------------

void Location_To_Genome(unsigned & Location)
{
	Genome_Position=0;
	while ( Genome_Position< Genome_Count )
	{
		if (Location < Offsets[Genome_Position]) break;
		Genome_Position++;
	}
	Genome_Position--;
	Location=Location-Offsets[Genome_Position];
	
}

void Location_To_Genome(unsigned & Location,Ann_Info & A)
{
	map <unsigned,Ann_Info> ::iterator I;
	I=Annotations.lower_bound(Location);
	if (I->first != Location) I--;
	A=I->second;	
	Location=Location-I->first;
}
//{-----------------------------  Parse Command Line  -------------------------------------------------
void Parse_Command_line(int argc, char* argv[])
{
	int Current_Option=0;
	const char* Short_Options ="ej:szrhvVi:b:L:l::o:Gg:O:pPfMm:x:";//allowed options....
	char* This_Program = argv[0];//Current program name....
	const char* Help_String=
"Parameters:\n"
" --help | -h\t\t\t\t Print help\n"
" --inputfile | -i <filename>\t\t Name of input file\n"
" --genome | -g <filename>\t\t Name of the genome mapped against\n"
" --outputfile | -o <filename>\t\t Name of output file\n"
" --buffersize | -b <integer> \t\t Size of disk buffers\n"
" --extend | -e \t\t\t\t extend hits\n"
" --verify | -v \t\t\t\t Verify hits\n"
" --verifyfromdisk | -V \t\t\t Verify hits from the disk ...\n"
" --location | -l <filename> \t\t use this file to filter locations by region ...\n"
" --offset | -O <integer> \t\t subtract <integer> from hit location...\n"
" --plusstrand | -p \t\t\t always output + strand coordinate....\n"
" --mixstrand | -x \t\t\t output mixed strand coordinates....\n"
" --maxhits | -m <integer> \t\t Maximum number of hits to output....\n"
" --zip | -z  \t\t\t\t compress output....\n"
" --gis | -G \t\t\t\t Output in batman format...\n"
" --misplus | -P \t\t\t give mismatch info relative to read..\n"
" --logfile | -L \t\t\t give mismatch info relative to read..\n"
;
	char* Name;int Last_Dash;char* Genome_Name;
	if(argc == 1) {printf("%s \n",Help_String);exit(0);}

	for(;;)	
	{
		Current_Option=getopt_long(argc, argv, Short_Options, Long_Options, NULL);
		if (Current_Option == -1 ) break;
		switch(Current_Option)
		{
			case 'h':
				printf("%s \n",Help_String);exit(0);
			case 'j':
				JQ=TRUE;
				COLOR_JQ_CTorGA=atol(optarg);
				break;
			case 'i':
				INPUTFILE=optarg;
				break;
			case 'x':
				PLUSSTRAND=FALSE;
				PLUSGAP=atol(optarg);
				break;
			case 'P':
				MIS_IN_PLUS_READ=TRUE;
				break;
			case 'M':
				EXTEND=TRUE;
				METH=TRUE;
				break;
			case 'e':
				EXTEND=TRUE;
				break;
			/*case 'L':
				SOLID=TRUE;
				break;*/
			case 'L':
				LOGFILE=optarg;
				break;
			case 's':
				SAM=TRUE;
				break;
			case 'G':
				SAM=FALSE;
				break;
			case 'f':
				FORMAT=TRUE;
				break;
			case 'z':
				OUTPUT_COMPRESS=TRUE;
				break;
			case 'p':
				PLUSSTRAND=TRUE;
				break;
			case 'o':
				OUTPUTFILE=optarg;
				break;
			case 'm':
				Force_Maxhits=atol(optarg);
				break;
			case 'O':
				Offset=atol(optarg);
				break;
			case 'v':
				VERIFY=TRUE;
				break;
			case 'V':
				VERIFY_FROM_DISK=TRUE;
				VERIFY=TRUE;
				break;
			case 'r':
				LOADREVERSEONLY = TRUE;
				break;
			case 'b':
				DISKBUFFERSIZE=atol(optarg);
				break;
			case 'l':
				if (optarg) LOCATIONFILE=optarg;	
				USELOCATION=TRUE;
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

				REVBWTINDEX = (char*)Command_Line_Buffer;
				if(Last_Dash) Last_Dash=Genome_Name-optarg+1; else Genome_Name--;
				strncpy(REVBWTINDEX,optarg,Last_Dash);
				REVBWTINDEX[Last_Dash+0]='r';REVBWTINDEX[Last_Dash+1]='e';REVBWTINDEX[Last_Dash+2]='v';
				strcpy(REVBWTINDEX+Last_Dash+3,Genome_Name+1);
				strcat(REVBWTINDEX+Last_Dash+3,".bwt"); 

				BWTFILE=REVBWTINDEX+500;
				strncpy(BWTFILE,optarg,Last_Dash);
				strcpy(BWTFILE+Last_Dash,Genome_Name+1);
				strcat(BWTFILE+Last_Dash,".bwt"); 


				REVOCCFILE = BWTFILE+500;
				strncpy(REVOCCFILE,optarg,Last_Dash);
				REVOCCFILE[Last_Dash+0]='r';REVOCCFILE[Last_Dash+1]='e';REVOCCFILE[Last_Dash+2]='v';
				strcpy(REVOCCFILE+Last_Dash+3,Genome_Name+1);
				strcat(REVOCCFILE+Last_Dash+3,".fmv"); 


				OCCFILE=REVOCCFILE+500;			
				strncpy(OCCFILE,optarg,Last_Dash);
				strcpy(OCCFILE+Last_Dash,Genome_Name+1);
				strcat(OCCFILE+Last_Dash,".fmv"); 

				SAFILE=OCCFILE+500;			
				strncpy(SAFILE,optarg,Last_Dash);
				strcpy(SAFILE+Last_Dash,Genome_Name+1);
				strcat(SAFILE+Last_Dash,".sa");

				REVSAFILE = SAFILE+500;
				strncpy(REVSAFILE,optarg,Last_Dash);
				REVSAFILE[Last_Dash+0]='r';REVSAFILE[Last_Dash+1]='e';REVSAFILE[Last_Dash+2]='v';
				strcpy(REVSAFILE+Last_Dash+3,Genome_Name+1);
				strcat(REVSAFILE+Last_Dash+3,".sa"); 

				BINFILE=REVSAFILE+500;			
				strncpy(BINFILE,optarg,Last_Dash);
				strcpy(BINFILE+Last_Dash,Genome_Name+1);
				strcat(BINFILE+Last_Dash,".bin");

				LOCATIONFILE=BINFILE+500;			
				strncpy(LOCATIONFILE,optarg,Last_Dash);
				strcpy(LOCATIONFILE+Last_Dash,Genome_Name+1);
				strcat(LOCATIONFILE+Last_Dash,".ann.location");

				PACFILE=LOCATIONFILE+500;			
				strncpy(PACFILE,optarg,Last_Dash);
				strcpy(PACFILE+Last_Dash,Genome_Name+1);
				strcat(PACFILE+Last_Dash,".pac");




/*
				REVBWTINDEX = Command_Line_Buffer;
				REVBWTINDEX[0]='r';REVBWTINDEX[1]='e';REVBWTINDEX[2]='v';
				strcpy(REVBWTINDEX+3,optarg);
				strcat(REVBWTINDEX+3,".bwt"); 
				BWTFILE=REVBWTINDEX+3;
				REVOCCFILE = BWTFILE+500;
				
				REVOCCFILE[0]='r';REVOCCFILE[1]='e';REVOCCFILE[2]='v';
				strcpy(REVOCCFILE+3,optarg);	
				strcpy(REVOCCFILE+500,REVOCCFILE);
				REVSAFILE=REVOCCFILE+500;
				strcat(REVSAFILE,".sa");
				SAFILE=REVSAFILE+3;
				strcat(REVOCCFILE+3,".fmv"); 
				OCCFILE=REVOCCFILE+3;			

				BINFILE=REVSAFILE +500;
				strcpy(BINFILE,optarg);
				strcat(BINFILE,".bin");*/
				break;
			default:
				printf("%s \n",Help_String);
				exit(0);
		}
	}	
}

//}-----------------------------  Parse Command Line  -------------------------------------------------


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
	Handle=fopen(File_Name,Mode);
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
/*	gz_stream *s=(gz_stream*)Handle;
	if (s->transparent) INPUT_ZIPPED=FALSE; else INPUT_ZIPPED=TRUE;
	Input_FileO=s->file;

	fseek(Input_FileO, 0L, SEEK_END);
	File_Size = ftello64(Input_FileO);
	gzrewind(Handle);//,0,SEEK_SET);//go top*/

}
//}----------------------------------- FILE HANDLING ---------------------------------------------------------

//{----------------------------------- FM INDEX ROUTINES ---------------------------------------------------------
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initFMI
 *  Description:  Opens FM index fmiFile
 * =====================================================================================
 */
BWT* initFMI(const char* BWTCodeFileName,const char* BWTOccValueFileName, const char* SAFile) 

{
	BWT *fmi;
        int PoolSize = 524288;
	MMMasterInitialize(3, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);

	fmi = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, SAFile, NULL, NULL, NULL);//Load FM index

	return fmi;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_SARange
 *  Description:  gets the SA range of strings having prefix [New_Char][Range]
 * =====================================================================================
 */

//}----------------------------------- FM INDEX ROUTINES ---------------------------------------------------------

//{-----------------------------------DEBUG ROUTINE---------------------------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Print_Location
 *  Description:  Prints the hex location of start and end of Range
 * =====================================================================================

void Print_Location (unsigned Range,BWT *fmi)
{
		printf("%x :",BWTSaValue(fmi,Range));//<<":";//%d;%d>",coordinate,textPosition[Pattern_Number]);
		//printf("%x :",BWTSaValue(fmi,Range.End));//<<":";//%d;%d>",coordinate,textPosition[Pattern_Number]);
}*/
//}-----------------------------------DEBUG ROUTINE---------------------------------------------------------

void Read_INI()
{

	dictionary* Dictionary;
	Dictionary=iniparser_load("batman.ini",FALSE);
	if (!Dictionary) Dictionary=iniparser_load("~/batman.ini",FALSE);
	if (Dictionary)
	{
		PLUSSTRAND = iniparser_getint(Dictionary,"decode:plusstrand",0);
		FORMAT = iniparser_getint(Dictionary,"decode:formatted",1);
		
	}
	iniparser_freedict(Dictionary);
}

void Show_Progress(unsigned Percentage)
{
	if (Percentage >=97) return;
	printf("+%d%\b\b\b",Percentage);
	fprintf(Log_File,"%d%%->",Percentage);
	fflush(stdout);
}

void fprintfX(void* Handle, const char* Format_String, ...)
{
	va_list argptr;
	va_start(argptr, Format_String);
	if(OUTPUT_COMPRESS)
	{
		int Len;

		Zip_Buffer[sizeof(Zip_Buffer) - 1] = 0;
		Len = vsnprintf(Zip_Buffer, sizeof(Zip_Buffer), Format_String, argptr);
		gzwrite((gzFile*) Handle, Zip_Buffer, (unsigned)Len);

	}
	else
	{
		vfprintf((FILE*) Handle,Format_String,argptr);
	}
	va_end(argptr);
}


unsigned Get_File_Size(FILE* File)
{
	fseek (File , 0 , SEEK_END);
	unsigned Size = ftell (File);
	rewind (File);
	return Size;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_Bases
 *  Description:  read bases from packed file..
 * =====================================================================================
 */
void Get_Bases (unsigned Location,int StringLength,char Strand)
{
	Location--;
	for (int i=0;i<StringLength;i++)
	{
		unsigned char L= (unsigned char)(Original_Text[(Location+i)/4]<< (((Location+i) % 4) * 2)) >>6;
		Org_String[i]=L;
	}
}

static const int nst_ntnt2cs_table[] = { 4, 0, 0, 1, 0, 2, 3, 4, 0, 3, 2, 4, 1, 4, 4, 4 };

#define COLOR_MM 19
#define NUCL_MM  25

void cs2nt_JQ(int size, const char *nt_ref, char *cs_read,char* cs_qual, char *nt_read, char *btarray)
{
	int h[8], curr, last;
	int x, y, xmin, hmin, k;
	// h[0..3] and h[4..7] are the current and last best score array, depending on curr and last

	// recursion: initial value
	//for(int i=0;i<size;i++) {cs_read[i]=Char_To_Code[cs_read[i]];}
	if (nt_ref[0] >= 4) memset(h, 0, sizeof(int) << 2);
	else {
		for (x = 0; x != 4; ++x) h[x] = NUCL_MM;
		h[nt_ref[0]] = 0;
	}
	// recursion: main loop
	curr = 1; last = 0;
	for (k = 1; k <= size; ++k) {
		for (x = 0; x != 4; ++x) {
			int min = 0x7fffffff, ymin = 0;
			for (y = 0; y != 4; ++y) {
				int s = h[last<<2|y];
				if ((cs_qual[k-1]) != 63 && cs_read[k-1] != nst_ntnt2cs_table[1<<x|1<<y])
					s += ((cs_qual[k-1]) < COLOR_MM)? COLOR_MM : (cs_qual[k-1]); // color mismatch
				if ((COLOR_JQ_CTorGA==0 && nt_ref[k]==2 && x==0) || (COLOR_JQ_CTorGA==1 && nt_ref[k]==1 && x==3)) s+=NUCL_MM/2; //jq
				else if (nt_ref[k] < 4 && nt_ref[k] != x) s += NUCL_MM; // nt mismatch				
				if (s < min) {
					min = s; ymin = y;
				}
			}
			h[curr<<2|x] = min; btarray[k<<2|x] = ymin;
		}
		last = curr; curr = 1 - curr; // swap
	}
	// back trace
	hmin = 0x7fffffff; xmin = 0;
	for (x = 0; x != 4; ++x) {
		if (h[last<<2|x] < hmin) {
			hmin = h[last<<2|x]; xmin = x;
		}
	}
	nt_read[size] = xmin;
	for (k = size - 1; k >= 0; --k)
	{
		nt_read[k] = btarray[(k+1)<<2 | nt_read[k+1]];
	}
}
/*
  {A,C,G,T,N} -> {0,1,2,3,4}
  nt_ref[0..size]: nucleotide reference: 0/1/2/3/4
  cs_read[0..size-1]: color read+qual sequence: base<<6|qual; qual==63 for N
  nt_read[0..size]: nucleotide read sequence: 0/1/2/3 (returned)
  btarray[0..4*size]: backtrack array (working space)
 */

void cs2nt_DP(int size, const char *nt_ref, char *cs_read,char* cs_qual, char *nt_read, char *btarray)
{
	if (JQ) 
	{
		cs2nt_JQ(size, nt_ref, cs_read, cs_qual, nt_read, btarray);
		return;
	}
	int h[8], curr, last;
	int x, y, xmin, hmin, k;
	// h[0..3] and h[4..7] are the current and last best score array, depending on curr and last

	// recursion: initial value
	//for(int i=0;i<size;i++) {cs_read[i]=Char_To_Code[cs_read[i]];}
	if (nt_ref[0] >= 4) memset(h, 0, sizeof(int) << 2);
	else {
		for (x = 0; x != 4; ++x) h[x] = NUCL_MM;
		h[nt_ref[0]] = 0;
	}
	// recursion: main loop
	curr = 1; last = 0;
	for (k = 1; k <= size; ++k) {
		for (x = 0; x != 4; ++x) {
			int min = 0x7fffffff, ymin = 0;
			for (y = 0; y != 4; ++y) {
				int s = h[last<<2|y];
				if ((cs_qual[k-1]) != 63 && cs_read[k-1] != nst_ntnt2cs_table[1<<x|1<<y])
					s += ((cs_qual[k-1]) < COLOR_MM)? COLOR_MM : (cs_qual[k-1]); // color mismatch
				if (nt_ref[k] < 4 && nt_ref[k] != x) s += NUCL_MM; // nt mismatch
				if (s < min) {
					min = s; ymin = y;
				}
			}
			h[curr<<2|x] = min; btarray[k<<2|x] = ymin;
		}
		last = curr; curr = 1 - curr; // swap
	}
	// back trace
	hmin = 0x7fffffff; xmin = 0;
	for (x = 0; x != 4; ++x) {
		if (h[last<<2|x] < hmin) {
			hmin = h[last<<2|x]; xmin = x;
		}
	}
	nt_read[size] = xmin;
	for (k = size - 1; k >= 0; --k)
	{
		nt_read[k] = btarray[(k+1)<<2 | nt_read[k+1]];
	}
}

/*
   nt_read[0..size]: nucleotide read sequence: 0/1/2/3
   cs_read[0..size-1]: color read+qual sequence: base<<6|qual; qual==63 for N
   tarray[0..size*2-1]: temporary array
*/

char *cs2nt_nt_qual(int size, const char *nt_read, const char *cs_read, char* cs_qual, char *tarray)
{
	int k, c1, c2;
	char *t2array = tarray + size;
	// get the color sequence of nt_read
	c1 = nt_read[0];
	for (k = 1; k <= size; ++k) {
		c2 = nt_read[k]; // in principle, there is no 'N' in nt_read[]; just in case
		tarray[k-1] = (c1 >= 4 || c2 >= 4)? 4 : nst_ntnt2cs_table[1<<c1 | 1<<c2];
		c1 = c2;
	}
	for (k = 1; k != size; ++k) {
		int q = 0;
		if (tarray[k-1] == cs_read[k-1] && tarray[k] == cs_read[k]) {
			q = (int)(cs_qual[k-1]) + (int)(cs_qual[k]) + 10;
		} else if (tarray[k-1] == cs_read[k-1]) {
			q = (int)(cs_qual[k-1]) - (int)(cs_qual[k]);
		} else if (tarray[k] == cs_read[k]) {
			q = (int)(cs_qual[k]) - (int)(cs_qual[k-1]);
		} // else, q = 0
		if (q < 0) q = 0;
		if (q > 60) q = 60;
		t2array[k] =  q;
		if ((cs_qual[k-1]) == 63 || (cs_qual[k]) == 63) t2array[k] = 0;
	}
	for (k = 1; k != size; ++k) t2array[k] = (t2array[k]==63) ? 33: t2array[k]+33;
	return t2array + 1; // of size-2
}

void printfL(char* Format, ...)
{
    va_list args;
    va_start(args,Format);
    vprintf(Format,args);
    va_end(args);
    va_start(args,Format);
    vfprintf(Log_File,Format,args);
    va_end(args);
}

