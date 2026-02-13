#include <string>
#include "string.h"
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <limits.h>
#include <cassert>
#include <algorithm>
#include <cmath>
#include "ssw.h"

//chr::Bound - region under study..
std::string Chr;
unsigned Bound;
int STD=50;
int MAX_SIGLEN=10;
bool FILTER_DUP=false;

int Stat_Split_Adjust=0;
int Stat_Split_Reads=0;
std::map <std::string, unsigned> Annotations;
std::ofstream Log_File;
bool Log=false;

char Char_To_Code[255];
char Char_To_CodeC[256];

struct Cigar_Info
{
	        int M;int I;int D;int Mis;int Length;int Score;int Indel_Count;int QScore;int BQScore;
};

struct vir_hit
{
	unsigned St;
	unsigned Ed;
	char Sign;//Left/Right handedness of integration..
	char Orientation;//the mapped direction in genome..
	char Type;//A split hit or not..
	std::string Sequence;
};

struct splitread
{
	unsigned St;
	unsigned Ed;
	char Sign;
	unsigned BP;
	int Overlap;
	int Gap;
};

class Simple_Interval//basically a simple interval, but in the context of a read can store blast info..
{
	public:
		unsigned St;
		unsigned Ed;
		int H_St;
		int H_Ed;
		int V_St;
		int V_Ed;
		int V_Q_St;
		int V_Q_Ed;
		char Type;
		bool Paired;
		char Orientation;
		int Read_Length;
		float E;
		Simple_Interval();
};

Simple_Interval::Simple_Interval()
{
	St=Ed=INT_MAX;
	V_St=V_Ed=INT_MAX;
}

class Interval
{
	public:
		unsigned St;
		unsigned Ed;
		int Count;
		std::map<std::string,char> Reads;
		Interval();
		void Add_Read(std::string & S);
		void Add_Read(Interval & I);
};

Interval::Interval()
{
	St=Ed=INT_MAX;
	Count=0;
}
void Interval::Add_Read(std::string & S)
{
	/*if(Reads.count(S))
	{
		std::cerr<<"Multiple read ID"<<Chr<<"\t"<<Bound<<"\n";
		exit(EXIT_FAILURE);
	}
	else*/
	{
		Reads[S]=1;
	}
	Count++;
}

void Interval::Add_Read(Interval & Int)
{
	std::map<std::string,char>::iterator I;
	for(I=Int.Reads.begin();I!=Int.Reads.end();++I)
	{
		if(Reads.count(I->first)==0)
		{
			Reads[I->first]=1;
			Count++;
		}
		else
		{
			//std::cerr <<"DUP"<<I->first<<"\n";
		}
	}
}
struct INTPAIR
{
	int Rank;
	int Abs_Rank;
};

std::ifstream H_Fq,T_Fq;//Open input files..

double Calculate_P_Value(Interval & Int,std::map<std::string,Simple_Interval> & MapH,std::map<std::string,Simple_Interval> & MapT,double & Best_P);
void Init_Character_Codes();
bool Recover_Split(
		std::string & Des,
		vir_hit & V_Hit,
		std::ifstream & H_Fq,
		std::map<std::string,Simple_Interval> & Mapped_ReadsPT,
		std::map<std::string,Simple_Interval> & Mapped_ReadsMT,
		std::map<std::string,Simple_Interval> & Plus_Reads,
		std::map<std::string,Simple_Interval> & Minus_Reads
		);
void ssw_cigar_print(s_align* a,char* Cigar,Cigar_Info & C,char* Ref,char* Pattern,int Clip_H,int Clip_T,int StringLength);
void Read2Bin(char* Dest,const char* Source,int StringLength);
void Read2RevCBin(char* Dest,const char* Source,int StringLength);
void Split(const std::string &s, char delim, std::vector<std::string> &elems);
void Fill_Split_Info(Simple_Interval & Human,char V_Suf_Pref,int V_St,int V_Ed,int V_Q_St,int V_Q_Ed);
void Build_Possible_Clusters(std::map<unsigned,Interval> & Map,std::vector<Interval> & Clusters);
void Get_Command_Line(int argc,char* argv[],unsigned & Bound,std::string & Chr,std::ifstream & Hg_Head,std::ifstream & Hg_Tail,std::ifstream & Vir_Head,std::ifstream & Vir_Tail,std::string & Genome_File,std::ifstream & H_Fq,std::ifstream & T_Fq);
void Read_Human_Blast(std::ifstream & Blast_Map,
		std::vector<int> & Total_Rank_P,
		std::vector<int> & Total_Rank_M,
		std::map<unsigned,Interval> & Plus_Maps,
		std::map<unsigned,Interval> & Minus_Maps,
		std::map<std::string,Simple_Interval> & Plus_Reads,
		std::map<std::string,Simple_Interval> & Minus_Reads,
		std::map<std::string,INTPAIR> & Read_Status
		);//,char FR);
void Change_Human_Read(
		char & Orientation,
		char O_Sign,
		std::string & Des,
		char Pref_Suf,
		std::map<std::string,Simple_Interval> & Mapped_Reads,
		std::map<std::string,INTPAIR> & Rank_H,
		std::map<std::string,INTPAIR> & Rank_T,
		Simple_Interval **Split_Interval
		);
void Check_Correctness(
		std::map<std::string,Simple_Interval> & Plus_ReadsH,
		std::map<std::string,Simple_Interval> & Minus_ReadsH,
		std::map<std::string,Simple_Interval> & Plus_ReadsT,
		std::map<std::string,Simple_Interval> & Minus_ReadsT,
		std::map<std::string,INTPAIR> & Rank_H,
		std::map<std::string,INTPAIR> & Rank_T
		);
float Median(std::vector<int> & Rank);
unsigned Get_File_Size(FILE* File);
int Int_Median(std::vector<int> & Rank,int & O1,int & O2);
void Check_Split_Read(std::ifstream & VFile,std::map<std::string,Simple_Interval> & Mapped_Reads,std::map<std::string,splitread> & Split,std::map<std::string,vir_hit> & Vir);
void Load_Virus_Hits(std::ifstream & VFile,std::map<std::string,Simple_Interval> & Mapped_ReadsPH,std::map<std::string,Simple_Interval> & Mapped_ReadsMH,std::map<std::string,Simple_Interval> & Mapped_ReadsPT,std::map<std::string,Simple_Interval> & Mapped_ReadsMT,std::map<std::string,vir_hit> & Vir,std::map<std::string,INTPAIR> & Rank_H,std::map<std::string,INTPAIR> & Rank_T,char HT);
int Count_Split_Read(Interval & Int,std::map<std::string,splitread> & D,std::map<std::string,splitread> & A,std::map<std::string,splitread> & B,std::map<std::string,splitread> & C,std::map<std::string,INTPAIR> Rank_H,std::map<std::string,INTPAIR> Rank_T,std::vector<int> & Ranks,std::vector<int> & Abs_Ranks,int & Uniq_Count,int & Multi_Count,int & Abs_Rank_Count);
int Is_Split_Read(std::map<std::string,splitread> & A,std::map<std::string,splitread> & B,std::map<std::string,splitread> & C,std::map<std::string,splitread> & D,splitread & S,std::string & Des);
void Get_Viral_Contig(Interval & Int,std::map<std::string,vir_hit> & V,unsigned & St,unsigned & Ed,std::map<std::string,Simple_Interval> & MapH,std::map<std::string,Simple_Interval> & MapsT,char & Sign,char H_Sign,int & Split_Reads,std::map<std::string,INTPAIR> & Rank_H,std::map<std::string,INTPAIR> & Rank_T,std::string & Chr,unsigned BP);
void Switch_Contig(unsigned & S,unsigned & E,char & Sign);
inline void Switch_Orientation(char & Sign);
void Drop_Rank(std::string & Des,std::map<std::string,INTPAIR> & Rank_H,std::map<std::string,INTPAIR> & Rank_T);
bool Rank_Bad(std::string & Des,std::map<std::string,INTPAIR> & Rank_H,std::map<std::string,INTPAIR> & Rank_T);
FILE* File_Open(const char* File_Name,const char* Mode);
std::string Seek_Read(std::string & Des,std::ifstream & Fq);
void Load_Location(char* LOCATIONFILE, std::map <std::string, unsigned> & Annotations);

int HGMAPCUTOFF=10;
unsigned char *Original_Text;
std::string Genome_File;
void Get_Bases (std::string & Chr,unsigned Location,int StringLength,char* Org_String);

int main(int argc, char* argv[])
{
	std::ifstream Hg_Head,Hg_Tail,Vir_Head,Vir_Tail;//Open input files..
	Get_Command_Line(argc,argv,Bound,Chr,Hg_Head,Hg_Tail,Vir_Head,Vir_Tail,Genome_File,H_Fq,T_Fq);
	Init_Character_Codes();
	init_SSW();
	init_SSW_Clip(2 /*match*/,2 /*mismatch*/,3 /*gap_open*/,1 /*gap_extension*/);

	std::map<unsigned,Interval> Plus_Maps,Minus_Maps;
	std::map<std::string,Simple_Interval> Plus_ReadsH,Minus_ReadsH,Plus_ReadsT,Minus_ReadsT;
	float Median_Rank_P=0,Abs_Median_Rank_M;std::vector<int> Total_Rank_P;
	float Median_Rank_M=0,Abs_Median_Rank_P;std::vector<int> Total_Rank_M;
	std::map<std::string,INTPAIR> Rank_H;
	std::map<std::string,INTPAIR> Rank_T;
	std::vector<int> Ranks;
	std::vector<int> Abs_Ranks;

	FILE* Original_File=File_Open(Genome_File.c_str(),"rb");
	Original_Text=(unsigned char*) malloc(Get_File_Size(Original_File));
	if(!fread(Original_Text,Get_File_Size(Original_File),1,Original_File))fprintf (stderr,"Init(): Error loading genome...\n");

	
//Get HG reads near cluster regions..
	Read_Human_Blast(Hg_Head,Total_Rank_P,Total_Rank_M,Plus_Maps,Minus_Maps,Plus_ReadsH,Minus_ReadsH,Rank_H);//,'F');
	Read_Human_Blast(Hg_Tail,Total_Rank_P,Total_Rank_M,Plus_Maps,Minus_Maps,Plus_ReadsT,Minus_ReadsT,Rank_T);//,'R');

//Find out split reads if any..

	std::map<std::string,splitread> Split_HP,Split_HM,Split_TP,Split_TM;
	std::map<std::string,vir_hit> Viral_Hits;
	Load_Virus_Hits(Vir_Head,Plus_ReadsH,Minus_ReadsH,Plus_ReadsT,Minus_ReadsT,Viral_Hits,Rank_H,Rank_T,'H');Vir_Head.clear();Vir_Head.seekg(0,std::ios::beg);
	Load_Virus_Hits(Vir_Tail,Plus_ReadsH,Minus_ReadsH,Plus_ReadsT,Minus_ReadsT,Viral_Hits,Rank_H,Rank_T,'T');Vir_Head.clear();Vir_Tail.seekg(0,std::ios::beg);

	//Check_Correctness(Plus_ReadsT,Minus_ReadsT,Plus_ReadsH,Minus_ReadsH,Rank_H,Rank_T);


	Median_Rank_P=Median(Total_Rank_P);
	Median_Rank_M=Median(Total_Rank_M);

	std::vector<Interval> Clusters_P,Clusters_M;
	Build_Possible_Clusters(Plus_Maps,Clusters_P);
	Build_Possible_Clusters(Minus_Maps,Clusters_M);


	std::vector<Interval>::iterator I;
	std::map<unsigned,Interval> Final_Plus,Final_Minus;

	//Collapse intervals..
	for(I=Clusters_M.begin();I!=Clusters_M.end();I++)
	{
		Interval Int=*I;
		assert(Int.St<Int.Ed && Int.Ed !=INT_MAX && Int.St !=INT_MAX);
		if(!Final_Minus.count(Int.Ed))
		{
			Final_Minus[Int.Ed]=Int;
		}
		else
		{
			if(Final_Minus[Int.Ed].St>Int.St) //-------*-----(--------] where * is the current ( is the prev.
			{
				Final_Minus[Int.Ed]=Int;
			}
		}
	}

	for(I=Clusters_P.begin();I!=Clusters_P.end();I++)
	{
		Interval Int=*I;
		//std::cout << Int.St <<"-" << Int.Ed<<"\t"<<Int.Count<<"\n"; 
		assert(Int.St<Int.Ed && Int.Ed !=INT_MAX && Int.St !=INT_MAX);
		if(!Final_Plus.count(Int.Ed))
		{
			Final_Plus[Int.Ed]=Int;
		}
		else
		{
			if(Final_Plus[Int.Ed].St>Int.St) //-------[-----)-----*-- where * is the current ( is the prev.
			{
				Final_Plus[Int.Ed]=Int;
			}
		}
	}
	
	std::map<unsigned,Interval>::iterator I_Map;double Best_P;
	unsigned VSt,VEd;char VSign;int Top_Rank_Count;int Split_Reads; 
	for(I_Map=Final_Minus.begin();I_Map!=Final_Minus.end();I_Map++)
	{
		Interval Int=I_Map->second;
		double P_Value=Calculate_P_Value(Int,Minus_ReadsH,Minus_ReadsT,Best_P);
		if (Log) Log_File<<Chr<<"\t"<<Int.St<<"\t"<<Int.Ed<<std::endl;
		Get_Viral_Contig(Int,Viral_Hits,VSt,VEd,Minus_ReadsH,Minus_ReadsT,VSign,'-',Split_Reads,Rank_H,Rank_T,Chr,Int.St);//if(VSign=='-') Switch_Contig(VSt,VEd,VSign);
		assert(I_Map->first > Int.St);int Uniq_Count,Multi_Count;
		Count_Split_Read(Int,Split_HP,Split_HM,Split_TP,Split_TM,Rank_H,Rank_T,Ranks,Abs_Ranks,Uniq_Count,Multi_Count,Top_Rank_Count);
		Median_Rank_M=Median(Ranks);
		Abs_Median_Rank_M=Median(Abs_Ranks);
		std::cout <<"-\t" <<Chr<<"\t"<< Int.St     << "\t" <<Int.Ed<<"\t"<<VSign<<"\t"<<VSt<<"\t"<<VEd<<"\t"<<Median_Rank_M<<"\t"<<Int.Count<<"\t"<<Split_Reads<<"\t"<<Uniq_Count<<"\t"<<Multi_Count<<"\t"<<Top_Rank_Count<<"\t"<<P_Value<<"\t"<<Best_P<<"\n"; 
	}
	for(I_Map=Final_Plus.begin();I_Map!=Final_Plus.end();I_Map++)
	{
		Interval Int=I_Map->second;
		double P_Value=Calculate_P_Value(Int,Plus_ReadsH,Plus_ReadsT,Best_P);
		if (Log) Log_File<<Chr<<"\t"<<Int.St<<"\t"<<Int.Ed<<std::endl;
		Get_Viral_Contig(Int,Viral_Hits,VSt,VEd,Plus_ReadsH,Plus_ReadsT,VSign,'+',Split_Reads,Rank_H,Rank_T,Chr,I_Map->first);//Switch_Contig(VSt,VEd,VSign);
		assert(I_Map->first > Int.St);int Uniq_Count,Multi_Count;
		Count_Split_Read(Int,Split_HP,Split_HM,Split_TP,Split_TM,Rank_H,Rank_T,Ranks,Abs_Ranks,Uniq_Count,Multi_Count,Top_Rank_Count);
		Median_Rank_P=Median(Ranks);
		Abs_Median_Rank_P=Median(Abs_Ranks);
		std::cout <<"+\t" <<Chr<<"\t"<<I_Map->first<< "\t"<<Int.St<<"\t"<<VSign<<"\t"<<VSt<<"\t"<<VEd<<"\t"<<Median_Rank_P<<"\t"<<Int.Count<<"\t"<<Split_Reads<<"\t"<<Uniq_Count<<"\t"<<Multi_Count<<"\t"<<Top_Rank_Count<<"\t"<<P_Value<<"\t"<<Best_P<<"\n"; 
	}
	std::cerr<<"Number of split reads processed:\t"<<Stat_Split_Reads<<"\n";
	std::cerr<<"Number of split reads Re-aligned:\t"<<Stat_Split_Adjust<<"\n";

}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Read_Human_Blast
 *  Description:  Reads BLAST mapping to human..
 * =====================================================================================
 */
	
void Read_Human_Blast(
		std::ifstream & Blast_Map,
		std::vector<int> & Total_Rank_P,
		std::vector<int> & Total_Rank_M,
		std::map<unsigned,Interval> & Plus_Maps,
		std::map<unsigned,Interval> & Minus_Maps,
		std::map<std::string,Simple_Interval> & Plus_Reads,
		std::map<std::string,Simple_Interval> & Minus_Reads,
		std::map<std::string,INTPAIR> & Read_Status
		)//,char FR)
{
	std::string Line,Des;
	int Insert_Size=1000;
	int STD=50;

	while (getline(Blast_Map, Line)) 
	{
		if(Line[0]=='@')//Is it Query start?
		{
			Des=Line;//Save Description of read...
			if(!getline(Blast_Map, Line)) break;
			unsigned Read_Length=strtoul(Line.c_str(),NULL,0);//

			std::size_t Pipe_Pos = Des.find("|");//Find insert size..
			if (Pipe_Pos!=std::string::npos)
			{
				Insert_Size=strtoul((Des.substr(Pipe_Pos+1)).c_str(),NULL,0)+2*STD;
				assert (Insert_Size>=0);
			}

			bool No_Map=true;
			int Rank=0,Hit_Rank=0,Abs_Rank=0;std::string Hit_Expect;

			while(true)//Parsing hits belonging to Des..
			{
				if(!getline(Blast_Map, Line)) break;
				if(Line[0] != '>') break;
				Rank++;

				std::vector<std::string> Fields;
				Split(Line,'\t',Fields);
				if(Hit_Expect.compare(Fields[5])==0)//If tophit is same confident as next best hit, give the lowest rank..
				{
					Hit_Rank++;
				}
				/*  else
				{
					std::cout<<Des<<Hit_Expect<<"\t"<<Fields[5]<<"\n";
					assert(Abs_Rank<2);
					Abs_Rank++;
				}*/

				if(No_Map)
				{
					Fields[0].erase(0,1);//remove > from chr..
					if(Fields[0].compare(Chr)!=0)
						continue;
					unsigned St=strtoul(Fields[8].c_str(),NULL,0);//Hit Start and End..
					unsigned Ed=strtoul(Fields[9].c_str(),NULL,0);
					if(St>Ed)
						std::swap(St,Ed);
					if(St<Bound-Insert_Size || Ed>Bound+2*Insert_Size)
						continue;

					Simple_Interval Query;
					Query.Paired=false;
					Query.St=strtoul(Fields[6].c_str(),NULL,0);
					Query.Ed=strtoul(Fields[7].c_str(),NULL,0);
					Query.H_St=St;Query.H_Ed=Ed;Query.Read_Length=Read_Length;
					std::istringstream(Fields[5]) >> Query.E; 
					//Query.E=std::stof(Fields[5]);

					if(Query.St>Query.Ed)
						std::swap(Query.St,Query.Ed);
					Abs_Rank=Hit_Rank=Rank;Hit_Expect=Fields[5];

					int St_Gap=Query.St,Ed_Gap=Read_Length-Query.Ed;
					if(((Query.St!=1)||(Query.Ed!=Read_Length)) && ((St_Gap+Ed_Gap)>6))//Split read..
					{
						if(St_Gap>10 && Ed_Gap>10)//Lousy Match of hg side..
						{
							continue;//Ignore this read..
						}

						if(St_Gap>Ed_Gap && Query.St!=1) 
							Query.Type='S';
						else if(St_Gap<Ed_Gap && Query.Ed!=Read_Length)  
							Query.Type='P';

						if(Fields[11].compare("+")==0)
							Query.Orientation='+';
						else
							Query.Orientation='-';
						
						if((Query.Type=='S' && Query.Orientation=='+')||(Query.Type=='P' && Query.Orientation=='-'))//(+ suffix,- prefix)
						{
							Minus_Maps[St].Add_Read(Des);
							Minus_Maps[St].Ed=Ed;
							Total_Rank_M.push_back(Rank);
							Minus_Reads[Des]=Query;
						}
						if((Query.Type=='P' && Query.Orientation=='+')||(Query.Type=='S' && Query.Orientation=='-'))//(- suffix,+ prefix)
						{
							Plus_Maps[St].Add_Read(Des);
							Plus_Maps[St].Ed=Ed;
							Total_Rank_P.push_back(Rank);
							Plus_Reads[Des]=Query;
						}

					}
					else
					{
						Query.Type='F';
						if((Fields[11].compare("-")==0))
						{
							Minus_Maps[St].Add_Read(Des);
							Minus_Maps[St].Ed=Ed;
							Total_Rank_M.push_back(Rank);
							Minus_Reads[Des]=Query;
						}
						else
						{
							Plus_Maps[St].Add_Read(Des);
							Plus_Maps[St].Ed=Ed;
							Total_Rank_P.push_back(Rank);
							Plus_Reads[Des]=Query;
						}
					}
					No_Map=false;
				}
			}
			INTPAIR I;
			if(Rank==1)
			{
				I.Rank= -1;//Unique hit..
			}
			else if(Hit_Rank)
			{
				I.Rank=Hit_Rank;
			}
			I.Abs_Rank=Abs_Rank;
			Read_Status[Des]=I;
		}
	}
	if (Blast_Map.bad()) 
	{
		std::cerr << "File IO Error..\n";exit(EXIT_FAILURE);
	}
	else if (!Blast_Map.eof()) 
	{
		std::cerr << "Format Error..\n";exit(EXIT_FAILURE);
	}
}

void Build_Possible_Clusters(std::map<unsigned,Interval> & Map,std::vector<Interval> & Clusters)
{
	std::map<unsigned,Interval>::iterator C_St,C_Cur;
	for(C_St=Map.begin();C_St!=Map.end();++C_St)
	{
		Interval Int=C_St->second;assert(Int.Ed!=INT_MAX);
		unsigned St=C_St->first,Last_St=INT_MAX;
		unsigned Ed=Int.Ed,Last_Ed=INT_MAX;

		Interval I;
		for(C_Cur=C_St;C_Cur!=Map.end();++C_Cur)
		{
			if(FILTER_DUP)
			{
				if(Last_St != INT_MAX)
				{
					if((Last_Ed==Ed) && (Last_St==St))
					{
						continue;
					}
				}
				Last_Ed==Ed;Last_St=St;
			}
			Int=C_Cur->second;assert(Int.Ed!=INT_MAX);
			assert(Int.Ed!=INT_MAX && Int.St==INT_MAX);//C_Cur->first:Int.Ed is the interval..
			if(Int.Ed - St > 850)
			{
				break;
			}
			else
			{
				Ed=Int.Ed;
				I.Add_Read(Int);
				//I.Count+=Int.Count;
			}
		}
		I.St=St;I.Ed=Ed;
		Clusters.push_back(I);
	}
}

void Split(const std::string &s, char delim, std::vector<std::string> &elems) 
{
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) 
	{
		elems.push_back(item);
	}
}

void Get_Command_Line(int argc,char* argv[],unsigned & Bound,std::string & Chr,std::ifstream & Hg_Head,std::ifstream & Hg_Tail,std::ifstream & Vir_Head,std::ifstream & Vir_Tail,std::string & Genome_File,std::ifstream & H_Fq,std::ifstream & T_Fq)
{
	if(argc<11)
	{
		std::cerr << argv[0] << " Chr Start head.blast.hg.txt tail.blast.hg.txt head.blast.vir.txt tail.blast.vir.txt \n"; 
		std::cerr << argc << "- Need more parameters..\n";
		exit(EXIT_FAILURE);
	}
	else
	{
		if(argc>11)
		{
			Log_File.open(argv[11],std::ofstream::out|std::ofstream::app);Log=true;
			if(!Log_File.is_open())
			{
				std::cerr << "Cannot open log file";
				exit(EXIT_FAILURE);	
			}
			if(argc>12)
			{
				if(strcmp(argv[12],"TRUE")==0)
				{
					std::cerr <<"Filtering Dups..\n";
					FILTER_DUP=true;
				}
			}
		}

		Chr.assign(argv[1]);
		Bound=strtoul(argv[2],NULL,0);
		Hg_Head.open(argv[3]);Hg_Tail.open(argv[4]);Vir_Head.open(argv[5]);Vir_Tail.open(argv[6]);
		Genome_File=argv[7];
		H_Fq.open(argv[8]);T_Fq.open(argv[9]);
		if(!(H_Fq.is_open() && T_Fq.is_open() && Hg_Head.is_open() && Hg_Tail.is_open() && Vir_Head.is_open() && Vir_Tail.is_open()))
		{
			std::cerr << "Error opening files..\n";
			exit(EXIT_FAILURE);
		}
		Load_Location(argv[10],Annotations);
	}
}

float Median(std::vector<int> & Rank)
{
	float median;
	int size = Rank.size();
	if(!size)
		return INT_MAX;

	std::sort(Rank.begin(), Rank.end());

	if (size  % 2 == 0)
	{
		median = (Rank[size / 2 - 1] + Rank[size / 2]) / 2;
	}
	else 
	{
		median = Rank[size / 2];
	}

	return median;
}

int Int_Median(std::vector<int> & Rank,int & O1,int & O2)
{
	int median;
	int size = Rank.size();
	if(!size)
		return INT_MAX;

	std::sort(Rank.begin(), Rank.end());
	O1=*Rank.begin();
	O2=Rank.back();

	median = Rank[size / 2];

	return median;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Virus_Hits;
 *  Description:  Load viral hits corresponding to human mappings..
 * =====================================================================================
 */

void Load_Virus_Hits(
		std::ifstream & VFile,
		std::map<std::string,Simple_Interval> & Mapped_ReadsPH,
		std::map<std::string,Simple_Interval> & Mapped_ReadsMH,
		std::map<std::string,Simple_Interval> & Mapped_ReadsPT,
		std::map<std::string,Simple_Interval> & Mapped_ReadsMT,
		std::map<std::string,vir_hit> & Vir,
		std::map<std::string,INTPAIR> & Rank_H,
		std::map<std::string,INTPAIR> & Rank_T,
		char HT
		)
{
	std::string Line,Des;

	while (getline(VFile, Line)) 
	{
		if(Line[0]=='@')//Is it Query start?
		{
			Des=Line;//Save Description of read...

			if(!Mapped_ReadsPT.count(Des) && !Mapped_ReadsMT.count(Des) && !Mapped_ReadsPH.count(Des) && !Mapped_ReadsMH.count(Des)) continue;
			if(!getline(VFile, Line)) break;
			unsigned Read_Length=strtoul(Line.c_str(),NULL,0);//

			bool No_Map=true;
			while(true)
			{
				if(!getline(VFile, Line)) break;
				if(Line[0] != '>') break;
				if(No_Map)
				{
					std::vector<std::string> Fields;
					Split(Line,'\t',Fields);
					Fields[0].erase(0,1);//remove > from chr..
					unsigned St=strtoul(Fields[8].c_str(),NULL,0);//Hit Start and End..
					unsigned Ed=strtoul(Fields[9].c_str(),NULL,0);

					Interval Query;
					Query.St=strtoul(Fields[6].c_str(),NULL,0);
					Query.Ed=strtoul(Fields[7].c_str(),NULL,0);

					vir_hit V_Hit,Vold;
					V_Hit.St=St;V_Hit.Ed=Ed;Vold=V_Hit;
					if(V_Hit.St>V_Hit.Ed)
					{
						V_Hit.Orientation='-';
						assert(!Fields[11].compare("-"));
						//std::swap(V_Hit.St,V_Hit.Ed);
					}
					else
					{
						V_Hit.Orientation='+';
						assert(!Fields[11].compare("+"));
					}
					char Orientation=Fields[11].c_str()[0];
					Simple_Interval *Split_Interval=NULL;//=Mapped_ReadsPH[Des];

					int St_Gap=Query.St,Ed_Gap=Read_Length-Query.Ed;//TODO
					//if(St_Gap>10||(Query.St!=1)||(Query.Ed!=Read_Length))//Split read..
					int Q_St,Q_Ed;unsigned Ref_St,Ref_Ed;
					if((Query.St!=1)||(Query.Ed!=Read_Length))//Split read..
					{
						char SP;
						if(Query.St!=1)
							SP='S';
						else if(Query.Ed!=Read_Length)
							SP='P';

						if(SP=='S')
						{
							if(HT=='H')
							{	
								if((!Mapped_ReadsPH.count(Des)&&!Mapped_ReadsMH.count(Des))&&(St_Gap+Ed_Gap>HGMAPCUTOFF))//sizeable chunk does not map to HG ..
								{
									std::cerr<<Chr<<":"<<Bound<<"\t"<<Des<<"\n";
									if(!Recover_Split(Des,V_Hit,H_Fq,Mapped_ReadsPT,Mapped_ReadsMT,Mapped_ReadsPH,Mapped_ReadsMH))
									{
										Drop_Rank(Des,Rank_H,Rank_T);
									}
									else
									{
										Stat_Split_Adjust++;
									}
								}
								if(Mapped_ReadsPH.count(Des))
								{
									Change_Human_Read(Orientation,'-',Des,'P',Mapped_ReadsPH,Rank_H,Rank_T,&Split_Interval);
								}
								else if(Mapped_ReadsMH.count(Des))
								{
									Change_Human_Read(Orientation,'+',Des,'P',Mapped_ReadsMH,Rank_H,Rank_T,&Split_Interval);
								}
							}
							else 
							{
								assert(HT=='T');
								if(!Mapped_ReadsPT.count(Des)&&!Mapped_ReadsMT.count(Des)&&(St_Gap+Ed_Gap>HGMAPCUTOFF))//sizeable chunk does not map to HG ..
								{
									if(!Recover_Split(Des,V_Hit,T_Fq,Mapped_ReadsPH,Mapped_ReadsMH,Mapped_ReadsPT,Mapped_ReadsMT))
									{
										Drop_Rank(Des,Rank_H,Rank_T);
									}
									else
									{
										Stat_Split_Adjust++;
									}
									std::cerr<<Chr<<":"<<Bound<<"\t"<<Des<<"\n";
								}
								if(Mapped_ReadsPT.count(Des))
								{
									Orientation='-';
									if(Mapped_ReadsPT[Des].Type=='P') Split_Interval=&Mapped_ReadsPT[Des];
									
								}
								else if(Mapped_ReadsMT.count(Des))
								{
									Orientation='+';
									if(Mapped_ReadsMT[Des].Type=='P') Split_Interval=&Mapped_ReadsMT[Des];
								}
							}
						}
						else 
						{
							assert(SP=='P');
							if(HT=='H') 
							{
								if(!Mapped_ReadsPH.count(Des)&&Mapped_ReadsMH.count(Des)&&(St_Gap+Ed_Gap>HGMAPCUTOFF))//sizeable chunk does not map to HG ..
								{
									if(!Recover_Split(Des,V_Hit,H_Fq,Mapped_ReadsPT,Mapped_ReadsMT,Mapped_ReadsPH,Mapped_ReadsMH))
									{
										Drop_Rank(Des,Rank_H,Rank_T);
									}
									else
									{
										Stat_Split_Adjust++;
									}
									std::cerr<<Chr<<":"<<Bound<<"\t"<<Des<<"\n";
								}
								if(Mapped_ReadsMH.count(Des))
								{
									std::swap(V_Hit.St,V_Hit.Ed);
									Switch_Orientation(V_Hit.Orientation);
									Change_Human_Read(Orientation,'-',Des,'S',Mapped_ReadsMH,Rank_H,Rank_T,&Split_Interval);
								}
								else if(Mapped_ReadsPH.count(Des))
								{
									Change_Human_Read(Orientation,'+',Des,'S',Mapped_ReadsPH,Rank_H,Rank_T,&Split_Interval);
								}
							}
							else 
							{
								assert(HT=='T');
								if(!Mapped_ReadsMT.count(Des)&&!Mapped_ReadsPT.count(Des)&&(St_Gap+Ed_Gap>HGMAPCUTOFF))//sizeable chunk does not map to HG ..
								{
									if(!Recover_Split(Des,V_Hit,T_Fq,Mapped_ReadsPH,Mapped_ReadsMH,Mapped_ReadsPT,Mapped_ReadsMT))
									{
										Drop_Rank(Des,Rank_H,Rank_T);
									}
									else
									{
										Stat_Split_Adjust++;
									}
									std::cerr<<Chr<<":"<<Bound<<"\t"<<Des<<"\n";
								}
								if(Mapped_ReadsPT.count(Des))
								{
									/*Orientation='+';
									assert(!Mapped_ReadsMT.count(Des));
									if(Mapped_ReadsPT[Des].Type=='S') Split_Interval=&Mapped_ReadsPT[Des];*/
									Change_Human_Read(Orientation,'+',Des,'S',Mapped_ReadsPT,Rank_H,Rank_T,&Split_Interval);
									
								}
								else if(Mapped_ReadsMT.count(Des))
								{
									/*Orientation='-';*/
									std::swap(V_Hit.St,V_Hit.Ed);
									//if(Mapped_ReadsMT[Des].Type=='S') Split_Interval=&Mapped_ReadsMT[Des];
									Change_Human_Read(Orientation,'-',Des,'S',Mapped_ReadsMT,Rank_H,Rank_T,&Split_Interval);
								}
							}
						}
						if(Split_Interval)
						{
							if(!Rank_Bad(Des,Rank_H,Rank_T))
							{
								V_Hit.Type='S';
								Stat_Split_Reads++;
							}
							else
								V_Hit.Type='U';
							V_Hit.Sign=Orientation;Vir[Des]=V_Hit;	
							Fill_Split_Info(*Split_Interval,SP,V_Hit.St,V_Hit.Ed,Query.St,Query.Ed);//load read with split read info..
						}
					}
					else
					{
						V_Hit.Type='F';//a full hit..
						if(HT=='H') 
						{
							if(Mapped_ReadsPT.count(Des) && Mapped_ReadsMT.count(Des))
							{
								Drop_Rank(Des,Rank_H,Rank_T);

							}
							if(Mapped_ReadsPT.count(Des))
							{
								Mapped_ReadsPT[Des].Paired=true;
								Orientation='-';
								std::swap(V_Hit.St,V_Hit.Ed);
								Switch_Orientation(V_Hit.Orientation);
								V_Hit.Sign=Orientation;Vir[Des]=V_Hit;	
							}
							if(Mapped_ReadsMT.count(Des))
							{
								Mapped_ReadsMT[Des].Paired=true;
								Orientation='+';
								V_Hit.Sign=Orientation;Vir[Des]=V_Hit;	
							}
						}
						else
						{
							if(Mapped_ReadsPH.count(Des) && Mapped_ReadsMH.count(Des))
							{
								Drop_Rank(Des,Rank_H,Rank_T);

							}
							if(Mapped_ReadsPH.count(Des))
							{
								Mapped_ReadsPH[Des].Paired=true;
								Orientation='-';
								std::swap(V_Hit.St,V_Hit.Ed);
								Switch_Orientation(V_Hit.Orientation);
								V_Hit.Sign=Orientation;Vir[Des]=V_Hit;	
							}
							if(Mapped_ReadsMH.count(Des))
							{
								Mapped_ReadsMH[Des].Paired=true;
								Orientation='+';
								V_Hit.Sign=Orientation;Vir[Des]=V_Hit;	
							}
						}
					}

					No_Map=false;
				}
			}
		}
	}
	if (VFile.bad()) 
	{
		std::cerr << "File IO Error..\n";exit(EXIT_FAILURE);
	}
	else if (!VFile.eof()) 
	{
		std::cerr << "Format Error..\n";exit(EXIT_FAILURE);
	}
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Check_Split_Read
 *  Description:  Check if there are split mappings between hg and hbv..
 * =====================================================================================
 */

void Check_Split_Read(
		std::ifstream & VFile,
		std::map<std::string,Simple_Interval> & Mapped_Reads,
		std::map<std::string,splitread> & Split_Reads,
		std::map<std::string,vir_hit> & Vir
		)
{
	std::string Line,Des;

	while (getline(VFile, Line)) 
	{
		if(Line[0]=='@')//Is it Query start?
		{
			Des=Line;//Save Description of read...
			if(!Mapped_Reads.count(Des)) continue;
			if(!getline(VFile, Line)) break;

			int Q_Len=atoi(Line.c_str());//Read Length..
			Simple_Interval HG_Hit=Mapped_Reads[Des];

			bool No_Map=true;
			int Rank=0;int Best_Score=INT_MAX;
			splitread S;
			while(true)
			{
				if(!getline(VFile, Line)) break;
				if(Line[0] != '>') break;
				Rank++;
				if(No_Map)
				{
					std::vector<std::string> Fields;
					Split(Line,'\t',Fields);
					Fields[0].erase(0,1);//remove > from chr..
					unsigned St=strtoul(Fields[8].c_str(),NULL,0);//Hit Start and End..
					unsigned Ed=strtoul(Fields[9].c_str(),NULL,0);
					S.St=St;S.Ed=Ed;

					vir_hit V_Hit;
					//V_Hit.St=St;V_Hit.Ed=Ed;V_Hit.Sign=Fields[11].c_str()[0];
					V_Hit.Ed=St;V_Hit.St=Ed;V_Hit.Sign=Fields[11].c_str()[0];Switch_Orientation(V_Hit.Sign);
					Vir[Des]=V_Hit;	

					Interval Query;
					Query.St=strtoul(Fields[6].c_str(),NULL,0);
					Query.Ed=strtoul(Fields[7].c_str(),NULL,0);
					if(Query.St>Query.Ed)
						std::swap(Query.St,Query.Ed);

					int Overlap,Gap;
					if(Query.St>HG_Hit.St)//--HG-- --VIR--
					{
						if(Query.Ed<=HG_Hit.Ed)//Virus hit completely in human..
						{
							continue;
						}
						else
						{
							Overlap=Query.St-HG_Hit.Ed-1;
							Gap=Q_Len-Query.Ed;
							S.BP=St;
						}
					}
					else if(Query.St<HG_Hit.St)//--VIR-- --HG--
					{
						if(Query.Ed>=HG_Hit.Ed)//human hit completely in virus..
						{
							continue;
						}
						else
						{
							Overlap= -Query.St+HG_Hit.Ed-1;
							Gap=Q_Len-HG_Hit.Ed;
							S.BP=Ed;
						}
					}
					S.Overlap=Overlap;S.Gap=Gap;

					if(Best_Score> abs(Overlap)+Gap)
					{
						//std::cout<<Overlap<<Des<<"\n";
						Best_Score=abs(Overlap)+Gap;
						S.Sign=Fields[11][0];
						Split_Reads[Des]=S;
						//std::cerr<<"SP:"<<Des<<"\n";
					}

					No_Map=false;
				}
			}
		}
	}
	if (VFile.bad()) 
	{
		std::cerr << "File IO Error..\n";exit(EXIT_FAILURE);
	}
	else if (!VFile.eof()) 
	{
		std::cerr << "Format Error..\n";exit(EXIT_FAILURE);
	}
}

int Count_Split_Read(Interval & Int,std::map<std::string,splitread> & D,std::map<std::string,splitread> & A,std::map<std::string,splitread> & B,std::map<std::string,splitread> & C,std::map<std::string,INTPAIR> Rank_H,std::map<std::string,INTPAIR> Rank_T,std::vector<int> & Ranks,std::vector<int> & Abs_Ranks,int & Uniq_Count,int & Multi_Count,int & Abs_Rank_Count)
{
	int Split_Read_Count=0,Split_Origin;Uniq_Count=Multi_Count=0;
	splitread Split_Info;

	Abs_Rank_Count=0;
	for(std::map<std::string,char>::iterator I=Int.Reads.begin();I!=Int.Reads.end();I++)
	{
		std::string Des=I->first;
		if((Split_Origin=Is_Split_Read(A,B,C,D,Split_Info,Des)))
		{
			Split_Read_Count++;
		}

		int Rank;//Find the rank of hits and uniqness..
		if(Rank_H.count(Des))
		{
			INTPAIR I;
			I=Rank_H[Des];
			Rank=I.Rank;
			if(Rank==-1)
			{
				Uniq_Count++;
				Rank=1;
			}
			else
			{
				Multi_Count++;
			}
			Ranks.push_back(Rank);
			if(I.Abs_Rank==1) Abs_Rank_Count++;
			Abs_Ranks.push_back(I.Abs_Rank);
		}
		if(Rank_T.count(Des))
		{
			INTPAIR I;
			I=Rank_T[Des];
			Rank=I.Rank;
			if(Rank==-1)
			{
				Uniq_Count++;
				Rank=1;
			}
			else
			{
				Multi_Count++;
			}
			Ranks.push_back(Rank);
			if(I.Abs_Rank==1) Abs_Rank_Count++;
			Abs_Ranks.push_back(I.Abs_Rank);
		}
	}
	return Split_Read_Count;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_Viral_Contig
 *  Description:  Get the Viral contig belonging to interval.. 
 * =====================================================================================
 */

struct Hit
{
	int Count;
	char Sign;
	int St;
	int Ed;
};

bool Cmp_Contig(const Hit &a, const Hit &b)
{
	    return a.Count > b.Count;
}

struct V_END
{
	int Val;
	char Type;
};

void Get_Viral_Contig(Interval & Int,std::map<std::string,vir_hit> & V,unsigned & St,unsigned & Ed,std::map<std::string,Simple_Interval> & MapH,std::map<std::string,Simple_Interval> & MapsT,char & Sign,char H_Sign,int & Split_Reads,std::map<std::string,INTPAIR> & Rank_H,std::map<std::string,INTPAIR> & Rank_T,std::string & Chr,unsigned BP)
{

	Hit PP,PM,MM,MP;Split_Reads=INT_MAX;
	std::map<std::string,char>::iterator I=Int.Reads.begin();
	std::vector <V_END> StartsPP,EndsPP,StartsPM,EndsPM,StartsMP,StartsMM,EndsMP,EndsMM;
	for(;I!=Int.Reads.end();I++)//get reads involved in the integration..
	{
		std::string Des=I->first;
		if(V.count(Des))
		{
			if (Log) Log_File<<"\t"<<Des<<std::endl;
			vir_hit VHit=V[Des];
			//assert(VHit.Type!='U');
			if(VHit.Type=='U') continue;
			V_END V;
			V.Type=VHit.Type;
			if(VHit.Sign=='+')//seperate into two starts and end groups..
			{
				//assert(VHit.St<= VHit.Ed);
				if(H_Sign!='+')//Human and virus not on same side..
				{
					if(VHit.Orientation=='+')
					{
						V.Val=VHit.St;
						StartsPP.push_back(V);
						V.Val=VHit.Ed;
						EndsPP.push_back(V);
					}
					else
					{
						V.Val=VHit.St;
						StartsPM.push_back(V);
						V.Val=VHit.Ed;
						EndsPM.push_back(V);
					}
				}
			}
			else
			{
				assert(VHit.Sign=='-');//assert(VHit.St>= VHit.Ed);
				if(H_Sign!='-')//Human and virus not on same side..
				{
					if(VHit.Orientation=='+')
					{
						V.Val=VHit.St;
						StartsMP.push_back(V);
						V.Val=VHit.Ed;
						EndsMP.push_back(V);
					}
					else
					{
						V.Val=VHit.St;
						StartsMM.push_back(V);
						V.Val=VHit.Ed;
						EndsMM.push_back(V);
					}
				}
			}

		}
	}
	int PM_Count=0,PP_Count=0;//Number of hits in each array..
	int MM_Count=0,MP_Count=0;//Number of hits in each array..
	int StMP=INT_MAX,EdMP=0;//extremes of virus..
	int StPM=0,EdPM=INT_MAX;
	int StMM=0,EdMM=INT_MAX;
	int StPP=INT_MAX,EdPP=0;
	int MedSt;

	int Sp_MP=0;//extremes of virus..
	int Sp_PM=0;
	int Sp_MM=0;
	int Sp_PP=0;
	std::vector<Hit> Viral_Side;
	std::vector<Hit> Split_Side;
	std::vector<int> Sp;
	int O1=0,O2=0;

	if(!StartsMP.empty())//Find the leftmost start..
	{
		MP.Count=MP_Count=StartsMP.size();
		for(std::vector <V_END>::iterator I=StartsMP.begin();I!=StartsMP.end();++I)
		{
			if(I->Type == 'S') Sp.push_back(I->Val); 
			if(I->Val<StMP) StMP=I->Val;
		}
		MP.St=StMP;
		MP.Sign='+';
	}
	if((MedSt=Int_Median(Sp,O1,O2))!=INT_MAX)
	{
		std::cerr<<"ranki\t"<<O1<<"\t"<<MedSt<<"\t"<<O2<<"\n";
		MP.St=Sp_MP=MedSt; 
		MP.Count=Sp.size();
		Split_Side.push_back(MP);
	}
	Sp.clear();

	if(!EndsMP.empty())
	{
		MP_Count=EndsMP.size();
		for(std::vector <V_END>::iterator I=EndsMP.begin();I!=EndsMP.end();++I)
		{
			if(I->Val>EdMP) EdMP=I->Val;
		}
		MP.Ed=EdMP;
		Viral_Side.push_back(MP);
		
	}

	if(!StartsMM.empty())//Find the leftmost start..
	{
		MM.Count=MM_Count=StartsMM.size();
		for(std::vector <V_END>::iterator I=StartsMM.begin();I!=StartsMM.end();++I)
		{
			if(I->Type == 'S') Sp.push_back(I->Val); 
			if(I->Val>StMM) StMM=I->Val;
		}
		MM.St=StMM;
		MM.Sign='-';
	}
	if((MedSt=Int_Median(Sp,O1,O2))!=INT_MAX) 
	{
		std::cerr<<"ranki\t"<<O1<<"\t"<<MedSt<<"\t"<<O2<<"\n";
		MM.St=Sp_MM=MedSt; 
		MM.Count=Sp.size();
		Split_Side.push_back(MM);
	}
	Sp.clear();

	if(!EndsMM.empty())
	{
		MM_Count=EndsMM.size();
		for(std::vector <V_END>::iterator I=EndsMM.begin();I!=EndsMM.end();++I)
		{
			if(I->Val<EdMM) EdMM=I->Val;
		}
		MM.Ed=EdMM;
		Viral_Side.push_back(MM);
	}

	if(!StartsPP.empty())
	{
		PP_Count=StartsPP.size();
		for(std::vector <V_END>::iterator I=StartsPP.begin();I!=StartsPP.end();++I)
		{
			if(I->Type == 'S') Sp.push_back(I->Val); 
			if(I->Val<StPP) StPP=I->Val;
		}
		PP.St=StPP;
		PP.Sign='+';
	}
	if((MedSt=Int_Median(Sp,O1,O2))!=INT_MAX) 
	{
		std::cerr<<"ranki\t"<<O1<<"\t"<<MedSt<<"\t"<<O2<<"\n";
		PP.St=Sp_PP=MedSt; 
		PP.Count=Sp.size();
		Split_Side.push_back(PP);
	}
	Sp.clear();

	if(!EndsPP.empty())
	{
		PP.Count=PP_Count=EndsPP.size();
		for(std::vector <V_END>::iterator I=EndsPP.begin();I!=EndsPP.end();++I)
		{
			if(I->Val>EdPP) EdPP=I->Val;
		}
		PP.Ed=EdPP;
		Viral_Side.push_back(PP);
	}

	if(!StartsPM.empty())//Find the leftmost start..
	{
		PM.Count=PM_Count=StartsPM.size();
		for(std::vector <V_END>::iterator I=StartsPM.begin();I!=StartsPM.end();++I)
		{
			if(I->Type == 'S') Sp.push_back(I->Val); 
			if(I->Val>StPM) StPM=I->Val;
		}
		PM.St=StPM;
		PM.Sign='-';
	}
	if((MedSt=Int_Median(Sp,O1,O2))!=INT_MAX) 
	{
		std::cerr<<"ranki\t"<<O1<<"\t"<<MedSt<<"\t"<<O2<<"\n";
		PM.St=Sp_PM=MedSt; 
		PM.Count=Sp.size();
		Split_Side.push_back(PM);
	}
	Sp.clear();

	if(!EndsPM.empty())
	{
		PM_Count=EndsPM.size();
		for(std::vector <V_END>::iterator I=EndsPM.begin();I!=EndsPM.end();++I)
		{
			if(I->Val<EdPM) EdPM=I->Val;
		}
		PM.Ed=EdPM;
		Viral_Side.push_back(PM);
	}

	Hit Final_Hit;Final_Hit.St=0;Final_Hit.Ed=0;Final_Hit.Sign='?';
	if(!Split_Side.empty())
	{
		std::sort(Split_Side.begin(),Split_Side.end(),Cmp_Contig);
		int Size=Split_Side.size();
		assert(Size);
		Final_Hit=Split_Side[0];Final_Hit.Ed=0;
		if(Size>1)
		{
			if(Final_Hit.Count<Split_Side[1].Count +2)
			{
				Final_Hit.St=0;
				Final_Hit.Ed=0;
				Final_Hit.Sign='?';
				Final_Hit.Count=0;
			}
		}
		St=Final_Hit.St;Ed=Final_Hit.Ed;Sign=Final_Hit.Sign;
		Split_Reads=Final_Hit.Count;
	}
	else
	{
		std::sort(Viral_Side.begin(),Viral_Side.end(),Cmp_Contig);
		int Size=Viral_Side.size();
		if(Size)
		{
			Final_Hit=Viral_Side[0];
			if(Size>1)
			{
				if(Final_Hit.Count<Viral_Side[1].Count +2)
				{
					Final_Hit.St=0;
					Final_Hit.Ed=0;
					Final_Hit.Sign='?';
				}
			}
		}
		St=Final_Hit.St;Ed=Final_Hit.Ed;Sign=Final_Hit.Sign;Split_Reads=0;
	}

	assert(Split_Reads!=INT_MAX);
	/* St=0;Ed=0;Sign='?';//Return the longest contig and orientation..
	if(EdP)
	{
		if(EdM==INT_MAX)//only unique contig orientation
		{
			St=StP;Ed=EdP;Sign='+';
		}
		else
		{
			if(P_Count/M_Count >2)//More than twice reads support Plus..
			{
				St=StP;Ed=EdP;Sign='+';
			}
			else if(M_Count/P_Count >2)//More than twice reads support Minus..
			{
				St=StM;Ed=EdM;Sign='-';
			}
		}

	}
	else if(EdM)
	{
		St=StM;Ed=EdM;Sign='-';
	}*/
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Is_Split_Read
 *  Description:  Check if the Des is a split read, and if so get its details..
 * =====================================================================================
 */
int Is_Split_Read(std::map<std::string,splitread> & A,std::map<std::string,splitread> & B,std::map<std::string,splitread> & C,std::map<std::string,splitread> & D,splitread & S,std::string & Des)
{
	if(A.count(Des))
	{
		S=A[Des];
		return 1;
	}
	if(A.count(Des))
	{
		S=B[Des];
		return 1;
	}
	if(A.count(Des))
	{
		S=C[Des];
		return 1;
	}
	if(A.count(Des))
	{
		S=D[Des];
		return 1;
	}
	return 0;
}

void Switch_Contig(unsigned & S,unsigned & E,char & Sign)
{
	if(Sign=='?')
		return;
	unsigned T=S;
	S=E;E=T;
	/*if(Sign=='+')
		Sign='-';
	else
		Sign='+';*/
}

inline void Switch_Orientation(char & Sign)
{
	if(Sign=='?')
		return;
	if(Sign=='+')
		Sign='-';
	else
		Sign='+';
}

void Fill_Split_Info(Simple_Interval & Human,char V_Suf_Pref,int V_St,int V_Ed,int V_Q_St,int V_Q_Ed)
{
	if(Human.Type=='P')
	{
		assert(V_Suf_Pref!='P');
		Human.V_St=V_St;Human.V_Ed=V_Ed;
	}
	else if(Human.Type=='S')
	{
		assert(V_Suf_Pref!='S');
		Human.V_St=V_Ed;Human.V_Ed=V_St;
	}
	else
	{
		Human.V_St=INT_MAX;Human.V_Ed=INT_MAX;
	}
	Human.V_Q_St=V_Q_St;Human.V_Q_Ed=V_Q_Ed;
}

double Calculate_P_Value(Interval & Int,std::map<std::string,Simple_Interval> & MapH,std::map<std::string,Simple_Interval> & MapT,double & Best_P)
{
	std::map<std::string,char>::iterator I=Int.Reads.begin();
	double Prod=1;
	Best_P=1;
	for(;I!=Int.Reads.end();I++)//get reads involved in the integration..
	{
		std::string Des=I->first;
		if(MapH.count(Des))
		{
			double E=MapH[Des].E;
			double P=E;
			if(P>0.01)
			{
				P=exp(-E);
				P=1-P;
			}
			if(Best_P>P) Best_P=P;
			Prod=Prod*P;
			if(fabs(Prod-0)>1e-40)
			{
				Prod=1e-40;
			}
		}
		if(MapT.count(Des))
		{
			double E=MapT[Des].E;
			double P=E;
			if(P>0.01)
			{
				P=exp(-E);
				P=1-P;
			}
			if(Best_P>P) Best_P=P;
			Prod=Prod*P;
			if(fabs(Prod-0)>1e-40)
			{
				Prod=1e-40;
			}
		}
	}

	//double P_Value=1-Prod;
	double P_Value=Prod;
	return P_Value;
}
void Change_Human_Read(
		char & Orientation,
		char O_Sign,
		std::string & Des,
		char Pref_Suf,
		std::map<std::string,Simple_Interval> & Mapped_Reads,
		std::map<std::string,INTPAIR> & Rank_H,
		std::map<std::string,INTPAIR> & Rank_T,
		Simple_Interval **Split_Interval
		)
{

	Orientation=O_Sign;
	assert(Mapped_Reads.count(Des));
	if(Mapped_Reads[Des].Type==Pref_Suf) 
	{
		*Split_Interval=&Mapped_Reads[Des];
	}
	else//Human in bad pos relative to virus..
	{
//std::cout<<Mapped_Reads[Des].Type<<" Hum bad \n";
		//Drop_Rank(Des,Rank_H,Rank_T);
	}

}

void Drop_Rank(std::string & Des,std::map<std::string,INTPAIR> & Rank_H,std::map<std::string,INTPAIR> & Rank_T)
{
	if(Rank_H.count(Des))
	{
		Rank_H[Des].Rank=1000;
		Rank_H[Des].Abs_Rank=1000;
	}
	if(Rank_T.count(Des))
	{
		Rank_T[Des].Rank=1000;
		Rank_T[Des].Abs_Rank=1000;
	}
}

bool Rank_Bad(std::string & Des,std::map<std::string,INTPAIR> & Rank_H,std::map<std::string,INTPAIR> & Rank_T)
{
	bool Bad_Hit=false;
	if(Rank_H.count(Des))
	{
		if(Rank_H[Des].Rank==1000) Bad_Hit=true;
	}
	if(Rank_T.count(Des))
	{
		if(Rank_T[Des].Rank==1000) Bad_Hit=true;
	}
	return Bad_Hit;
}

/*void Check_Correctness(
		std::map<std::string,Simple_Interval> & Plus_ReadsH,
		std::map<std::string,Simple_Interval> & Minus_ReadsH,
		std::map<std::string,Simple_Interval> & Plus_ReadsT,
		std::map<std::string,Simple_Interval> & Minus_ReadsT,
		std::map<std::string,INTPAIR> & Rank_H,
		std::map<std::string,INTPAIR> & Rank_T
		)
{
	std::map<std::string,Simple_Interval>::iterator I;
	for(I=Plus_ReadsH.begin();I!=Plus_ReadsH.end();I++)
	{
		std::string Des=I->first;
		Simple_Interval S=I->second;
		if(S.Read_Length-(S.Ed-S.St)>25)
		{

		}

	}
}*/

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
		printf("File %s Cannot be opened ....\n",File_Name);
		return Handle;
	}
	else 
		return Handle;
}

unsigned Get_File_Size(FILE* File)
{
	fseek (File , 0 , SEEK_END);
	unsigned Size = ftell (File);
	rewind (File);
	return Size;
}

std::string Seek_Read(std::string & Des,std::ifstream & Fq)
{
	std::string L;
	//L.clear();
	//return L;
	while(getline(Fq, L))
	{
		if(Des.compare(1,std::string::npos,L,1,std::string::npos)==0)
		{
			getline(Fq, L);//Get Sequence..
			assert(!L.empty());
			return L;
		}
	}
	std::cerr<<"Missing read from fastq\n";exit(100);
	//return L;
}

void Load_Location(char* LOCATIONFILE, std::map <std::string, unsigned> & Annotations)
{
	std::string Line;
	std::string Chr;
	std::ifstream Location_File;
	Location_File.open(LOCATIONFILE);

	assert (Location_File.is_open());
	unsigned Off_Cum=0,Off=0;
	int Genome_Count=0;
	while (getline(Location_File,Line))
	{
		Off=atoi(Line.c_str());
		if (Genome_Count)
		{
			Annotations[Chr]=Off_Cum;
		}
		Off_Cum+=Off;

		getline(Location_File,Chr);
		Genome_Count++;	
	}
}

void Get_Bases (std::string & Chr,unsigned Location,int StringLength,char* Org_String) 
{
	assert(Location>=1);
	Location+=Annotations[Chr];//Get absolute location..
	for (int i=0;i<=StringLength;i++)
	{
		unsigned char L= (unsigned char)(Original_Text[(Location+i)/4]<< (((Location+i) % 4) * 2)) >>6;
		Org_String[i]=L;
		//std::cout<<"ACGT"[L];
	}
	//std::cout<<std::endl;
}

void Read2Bin(char* Dest,const char* Source,int StringLength)
{
	for (int i=0;i<StringLength;i++)
	{
		*Dest=Char_To_Code[Source[i]];Dest++;
	}
	*Dest=0;
}

void Read2RevCBin(char* Dest,const char* Source,int StringLength)
{
	for (int i=StringLength-1;i>=0;i--)
	{
		*Dest=Char_To_CodeC[Source[i]];Dest++;
	}
	*Dest=0;
}

void Init_Character_Codes()
{
	Char_To_Code['A']=0;Char_To_Code['C']=1;Char_To_Code['G']=2;Char_To_Code['T']=3;Char_To_Code['a']=0;Char_To_Code['c']=1;Char_To_Code['g']=2;Char_To_Code['t']=3;Char_To_Code['n']=0;Char_To_Code['N']=0;//we are using character count to store the fmicode for acgt
	Char_To_CodeC['A']=3;Char_To_CodeC['C']=2;Char_To_CodeC['G']=1;Char_To_CodeC['T']=0;Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;Char_To_CodeC['n']=0;Char_To_CodeC['N']=0;//we are using character count to store the fmicode for acgt
}

void ssw_cigar_print(s_align* a,char* Cigar,Cigar_Info & C,char* Ref,char* Pattern,int Clip_H,int Clip_T,int StringLength) 
{ //print the cigar out
	Clip_H=a->read_begin1;Clip_T=0;
	int c = 0,Tot_Length=0;char* Cigar_Ptr=Cigar;
	bool Cig_Err=false;
	if(Clip_H) Cigar_Ptr+=sprintf(Cigar_Ptr,"%dS", Clip_H); 
	for (c = 0; c < a->cigarLen; ++c) {
		int32_t letter = 0xf&*(a->cigar + c);
		int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
		Cigar_Ptr+=sprintf(Cigar_Ptr,"%d", length);
		if (letter == 0) 
		{
			Tot_Length+=length;
			*Cigar_Ptr='M';
			for(int i=0;i<length;i++)
			{
				assert(*Ref<4 && *Ref>=0);
				assert(*Pattern<4 && *Pattern>=0);
				if(*Ref!=*Pattern)C.Mis++;
				Ref++;Pattern++;
			}
			Cigar_Ptr++;C.M+=length;
		}
		else if (letter == 1)
		{
			Tot_Length+=length;
			*Cigar_Ptr='I';Cigar_Ptr++;C.I+=length;
			Pattern+=length;
		} 
		else 
		{
			*Cigar_Ptr='D';Cigar_Ptr++;C.D+=length;
			Ref+=length;
		}
		if(Cigar_Ptr-Cigar>=MAX_SIGLEN-6) 
		{
			assert(StringLength-Tot_Length-Clip_H>=0);
			Cigar_Ptr+=sprintf(Cigar_Ptr,"%dS", StringLength-Tot_Length-Clip_H); 
			Cig_Err=true;
			break;
		}
	}
	if(Clip_T && !Cig_Err) Cigar_Ptr+=sprintf(Cigar_Ptr,"%dS", Clip_T); 
	*Cigar_Ptr=0;
	//if(Cig_Err) *Cigar=0;
	C.Length=Tot_Length;
}

bool Recover_Split(
		std::string & Des,
		vir_hit & V_Hit,
		std::ifstream & H_Fq,
		std::map<std::string,Simple_Interval> & Mapped_ReadsPT,
		std::map<std::string,Simple_Interval> & Mapped_ReadsMT,
		std::map<std::string,Simple_Interval> & Plus_Reads,
		std::map<std::string,Simple_Interval> & Minus_Reads
		)
		
{
	char Suf_Pref;int Q_St;int Q_Ed;unsigned Ref_Start;unsigned Ref_End;
	bool Rev=false;
	int Insert_Size=1000;
	std::size_t Pipe_Pos = Des.find("|");//Find insert size..
	if (Pipe_Pos!=std::string::npos)
	{
		Insert_Size=strtoul((Des.substr(Pipe_Pos+1)).c_str(),NULL,0)+2*STD;
		assert (Insert_Size>=0);
	}

	char *Org_String=(char*)malloc((Insert_Size+1)*sizeof(char));
	unsigned Mate_Loc;unsigned Flank_St;
	std::cerr<<Chr<<":"<<Bound<<"\t"<<Des<<"\n";
	V_Hit.Sequence=Seek_Read(Des,H_Fq);
	if(Mapped_ReadsMT.count(Des))
	{
		Mate_Loc=Mapped_ReadsMT[Des].H_Ed;Flank_St=Mate_Loc-Insert_Size;
	}
	else if(Mapped_ReadsPT.count(Des))
	{
		Mate_Loc=Mapped_ReadsPT[Des].H_St;Flank_St=Mate_Loc;
		Rev=true;
	}
	else
	{
		return false;
	}
	Get_Bases(Chr,Flank_St,Insert_Size,Org_String);

	char *Sequence=(char*)malloc((V_Hit.Sequence.length()+1)*sizeof(char));
	if(Rev)
		Read2RevCBin(Sequence,V_Hit.Sequence.c_str(),V_Hit.Sequence.length());
	else
		Read2Bin(Sequence,V_Hit.Sequence.c_str(),V_Hit.Sequence.length());

	s_profile* p = ssw_init((int8_t*)Sequence, V_Hit.Sequence.length(), mata, n, 1);
	s_align* Aln=mengyao_ssw_core(Org_String,V_Hit.Sequence.length(),Sequence,Insert_Size,0,0/*DP*/, p);
	char Cigar[200];int Clip_H,Clip_T;Cigar_Info Cig_Info;
	ssw_cigar_print(Aln,Cigar,Cig_Info,Org_String,Sequence,Clip_H,Clip_T,V_Hit.Sequence.length());
	free (Org_String);
	free (Sequence);

	int Read_Length=V_Hit.Sequence.length();
	if(Q_St!=0 && Q_Ed!=0)
	{
		Q_St=Aln->read_begin1+1;
		Q_Ed=Aln->read_end1;
		Ref_Start=Aln->ref_begin1+Flank_St;
		Ref_End=Aln->ref_end1+Flank_St;
		if(Ref_Start>Ref_End)
			std::swap(Ref_Start,Ref_End);
		if(Ref_Start<Bound-Insert_Size || Ref_End>Bound+2*Insert_Size)
		{
			return false;
		}
		if(Q_St>Q_Ed)
			std::swap(Q_St,Q_Ed);
		int St_Gap=Q_St,Ed_Gap=Read_Length-Q_Ed;

		char Q_Type='B';
		if(St_Gap>Ed_Gap && Q_St!=1) 
			Q_Type='S';
		else if(St_Gap<Ed_Gap && Q_Ed!=Read_Length)  
			Q_Type='P';
		else return false;

		Simple_Interval Query;
		if(Rev) 
		{
			Query.Orientation='-';
			if(Q_Type=='S') 
			{
				Q_Type='P';
				int T=Q_St;
				Q_St=Read_Length-Q_Ed;
				Q_Ed=Read_Length-T;
				unsigned T1=Ref_Start;
				Ref_Start=Ref_End;
				Ref_End=T1;
			}
			else 
			{
				Q_Type='S'; 
			}

		}	
		else 
			Query.Orientation='+';

		Query.St=Q_St;Query.Ed=Q_Ed;
		Query.H_St=Ref_Start;Query.H_Ed=Ref_End;
		Query.Paired=true;
		Query.Type=Q_Type;
		Query.Read_Length=Read_Length;
		if((Query.Type=='S' && Query.Orientation=='+')||(Query.Type=='P' && Query.Orientation=='-'))//(+ suffix,- prefix)
		{
			Minus_Reads[Des]=Query;

		}
		else if((Query.Type=='P' && Query.Orientation=='+')||(Query.Type=='S' && Query.Orientation=='-'))//(- suffix,+ prefix)
		{
			Plus_Reads[Des]=Query;
		}
		else
		{
			return false;
		}

		return true;
	}
	return false;
}
