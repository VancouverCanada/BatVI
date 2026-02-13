package msa;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ariyaratnep
 * Date: 10/10/12
 * Time: 10:10 AM
 * To change this template use File | Settings | File Templates.
 */
public class BuildConsSeq
{
	/*
  final static int HASH_LENGTH=8;
  final static int SHIFT_LENGTH=3;
  final static int DIGIT_VAL=7;
  final static int MIN_OL_LENGTH=50;
  final static double MIN_MATCH_PERCENT=1.0;
  final static double OL_MATCH_RATIO=1.0;
	*/
	
  final static int HASH_LENGTH = 8;
  final static int SHIFT_LENGTH = 3;
  final static int DIGIT_VAL = 7;
  final static int MIN_OL_LENGTH = 30;
  final static double MIN_MATCH_PERCENT = 0.75; //0.96; //PY. why 0.96? 2014.02.11.
  final static double COMBINE_MATCH_PERCENT = 0.9;
  final static double OL_MATCH_RATIO = 0.5;

  static boolean init=false;
  static int c2i[]=null;
  static char i2c[]=null;

  HashEntry[] h=null;
  int nHashEntries=0;

  static int andValue=0;
  char[][] seqs=null;
  ConsSeq[] result=null;

  public static synchronized void init()
  {
    if (init) return;
    c2i=new int[256];
    c2i['A']=0; c2i['a']=0;
    c2i['C']=1; c2i['c']=1;
    c2i['G']=2; c2i['g']=2;
    c2i['T']=3; c2i['t']=3;
    c2i['N']=4; c2i['n']=4;

    i2c=new char[5];
    i2c[0]='A'; i2c[1]='C'; i2c[2]='G'; i2c[3]='T'; i2c[4]='N';

    andValue=0; for (int c1=0; c1<HASH_LENGTH; c1++) { andValue=andValue<<SHIFT_LENGTH; andValue=andValue|DIGIT_VAL; }

    init = true;
  }

  public BuildConsSeq(SWMap[] maps)
  {

    if (!init) init();

    seqs=new char[maps.length][];

    buildIndex(maps);
    HashCompCount[] comp=getComp();
    for (int c1=0; c1<h.length; c1++) h[c1]=null;
    h=null;
    buildCons(maps, comp);
    for (int c1=0; c1<seqs.length; c1++) seqs[c1]=null;
    seqs=null;
    comp=null;

  }




  private void buildCons(SWMap[] maps, HashCompCount[] comp)
  {
    int[] readLookUp=new int[maps.length]; Arrays.fill(readLookUp, -1);
    int[] offsetLookUp=new int[maps.length]; Arrays.fill(offsetLookUp, -1);
    ConsBuild[] cons=new ConsBuild[maps.length]; int consCount=0;

    int curRound=0;
    boolean cont=true;
    HashMap addedThisRound=new HashMap();
    while (cont)
    {
      cont = false;
      boolean newCons=false;
      curRound++;
      addedThisRound.clear();
      for (int c1=0; c1<comp.length; c1++)
      {
        HashCompCount cur=comp[c1];
        if (readLookUp[cur.read1]!=-1 && readLookUp[cur.read2]!=-1)  continue;
        else if ((readLookUp[cur.read1]==-1 && readLookUp[cur.read2]==-1) && !newCons)
        {
          OverlapScore s=getOverlapScore(seqs[cur.read1], seqs[cur.read2], cur.offset);
          if (scoreFail(s)) continue;
          //System.out.println(maps[cur.read1].tag+"\t"+maps[cur.read2].tag+"\t"+cur.count+"\t"+cur.offset);
          cont = true; newCons=true;
          cons[consCount]=new ConsBuild(maps.length, cur, seqs); cons[consCount].lastMod=curRound;
          readLookUp[cur.read1]=consCount; offsetLookUp[cur.read1]=0;
          readLookUp[cur.read2]=consCount; offsetLookUp[cur.read2]=cur.offset;

          addedThisRound.put(new HashReadCons(cur.read1, consCount), null);
          addedThisRound.put(new HashReadCons(cur.read2, consCount), null);

          consCount++;
        }
        else if (readLookUp[cur.read1]==-1 && readLookUp[cur.read2]!=-1)
        {
          int consIndex=readLookUp[cur.read2];

          if (cons[consIndex].lastMod>=curRound-1 && !addedThisRound.containsKey(new HashReadCons(cur.read1, consIndex)))
          {
            addedThisRound.put(new HashReadCons(cur.read1, consIndex), null);
            int relOffset=offsetLookUp[cur.read2]-cur.offset;

            boolean addSuccess=cons[consIndex].addCons(cur.read1, seqs, relOffset);
            if (addSuccess)
            {
              readLookUp[cur.read1]=consIndex;
              offsetLookUp[cur.read1]=relOffset;
              cons[consIndex].lastMod=curRound;
            }
          }
        }
        else if (readLookUp[cur.read1]!=-1 && readLookUp[cur.read2]==-1)
        {
          int consIndex=readLookUp[cur.read1];
          if (cons[consIndex].lastMod>=curRound-1 && !addedThisRound.containsKey(new HashReadCons(cur.read2, consIndex)))
          {
            addedThisRound.put(new HashReadCons(cur.read2, consIndex), null);
            int relOffset=offsetLookUp[cur.read1]+cur.offset;

            boolean addSuccess=cons[consIndex].addCons(cur.read2, seqs, relOffset);
            if (addSuccess)
            {
              readLookUp[cur.read2]=consIndex;
              offsetLookUp[cur.read2]=relOffset;
              cons[consIndex].lastMod=curRound;
            }
          }
        }
      }

    }

    buildConsSeq(maps, cons, consCount, seqs);

    readLookUp=null;
    offsetLookUp=null;
    for (int c1=0; c1<cons.length; c1++)
    {
      if (cons[c1]==null) continue;
      cons[c1].offsets=null;
      cons[c1].readids=null;
    }
    cons = null;



  }


  private void buildConsSeq(SWMap[] maps, ConsBuild[] cons, int consCount, char[][] seqs)
  {
    result=new ConsSeq[consCount];
    for (int c1=0; c1<result.length; c1++)
    {
      result[c1]=new ConsSeq(cons[c1].getCons(seqs));

      result[c1].minPos=Integer.MAX_VALUE; result[c1].maxPos=Integer.MIN_VALUE;
      for (int c2=0; c2<cons[c1].n; c2++)
      {
        int readid=cons[c1].readids[c2];
        result[c1].minPos=Math.min(result[c1].minPos, maps[readid].pos);
        result[c1].maxPos=Math.max(result[c1].maxPos, maps[readid].pos);
      }
      result[c1].mapids=new int[cons[c1].n];
      for (int c2=0; c2<result[c1].mapids.length; c2++) result[c1].mapids[c2]=cons[c1].readids[c2];
    }
  }

  public static boolean scoreFail(OverlapScore s)
  {
    if (s.matchLength<0) return true;
    if (s.olLength<MIN_OL_LENGTH) return true;
    if (s.matchLength<s.olLength*MIN_MATCH_PERCENT) return true;
    return false;
  }


  private OverlapScore getOverlapScore(char[] seq2, char[] seq1, int offset)
  {
    int count=0, match=0;

    for (int c1=0; c1<seq1.length; c1++)
    {
      int pos=offset+c1;
      if (pos<0) continue;
      if (pos>=seq2.length) break;
      count++;
      if (seq1[c1]==seq2[pos]) match++;

    }
    return new OverlapScore(count, match);
  }


  private HashCompCount[] getComp()
  {
    Hashtable<HashComp, Integer> map=new Hashtable<HashComp, Integer>();
    for (int c1=0; c1<nHashEntries; c1++)
    {
      for (int c2=c1+1; c2<nHashEntries; c2++)
      {
        if (h[c1].v!=h[c2].v) break;
        if (h[c1].readIndex==h[c2].readIndex) continue;
        HashComp comp=null;
        if (h[c1].start<=h[c2].start) comp=new HashComp(h[c2].readIndex, h[c1].readIndex, h[c2].start-h[c1].start);
        else comp=new HashComp(h[c1].readIndex, h[c2].readIndex, h[c1].start-h[c2].start);

        Integer v=map.get(comp);
        if (v==null) map.put(comp, new Integer(1));
        else map.put(comp, new Integer(v.intValue()+1));
        comp=null;
      }
    }

    Iterator<HashComp> i=map.keySet().iterator();
    HashCompCount[] result=new HashCompCount[map.size()]; int n=0;
    while (i.hasNext())
    {
      HashComp c=i.next(); Integer v=map.get(c);
      result[n++]=new HashCompCount(c.read1, c.read2, c.offset, v.intValue());
    }
    i=null;
    map.clear();
    map=null;

    Arrays.sort(result);

    return result;
  }

  private void buildIndex(SWMap[] maps)
  {
    int n=0;
    for (int c1=0; c1<maps.length; c1++)
    {
      n = n+Math.max(0,maps[c1].tag.length()-HASH_LENGTH+1);
    }
    h=new HashEntry[n];
    for (int c1=0; c1<h.length; c1++) h[c1]=new HashEntry();


    nHashEntries=0;

    for (int mapPos=0; mapPos<maps.length; mapPos++)
    {
      int curValue=0;
      seqs[mapPos]=maps[mapPos].tag.toCharArray();
      if (seqs[mapPos].length<HASH_LENGTH) continue;
      for (int seqPos=0; seqPos<HASH_LENGTH-1; seqPos++) { curValue=curValue<<SHIFT_LENGTH; curValue=curValue|c2i[seqs[mapPos][seqPos]]; }
      for (int seqPos=HASH_LENGTH-1; seqPos<seqs[mapPos].length; seqPos++)
      {
        curValue=curValue<<SHIFT_LENGTH; curValue=curValue|c2i[seqs[mapPos][seqPos]];

        if (!acceptHash(seqs[mapPos], seqPos-HASH_LENGTH+1, seqPos)) continue;

        h[nHashEntries].v=(curValue&andValue);
        h[nHashEntries].readIndex=mapPos;
        h[nHashEntries].start=seqPos-HASH_LENGTH+1;
        nHashEntries++;
      }
    }

    Arrays.sort(h, 0, nHashEntries);
  }

  private boolean acceptHash(char[] s, int from, int to)
  {
    boolean homo=true;
    for (int c1=from+1; c1<=to; c1++) if (s[c1]!=s[from]) { homo=false; break; }
    if (homo) return false;

    boolean dual=true;
    for (int c1=from+2; c1<=to; c1++) if (s[c1]!=s[c1-2]) { dual=false; break; }
    if (dual) return false;

    return true;
  }

  static String hash2String(int v)
  {
    if (!init) init();
    String s="";
    for (int c1=0; c1<HASH_LENGTH; c1++)
    {
      s = i2c[v&DIGIT_VAL]+s;
      v=v>>SHIFT_LENGTH;
    }
    return s;
  }

  public static ConsSeq[] getConsSeq(SWMap[] maps)
  {
    BuildConsSeq c=new BuildConsSeq(maps);
    ConsSeq[] result=c.result;
    c = null;
    return result;
  }
}

class HashEntry implements Comparable
{
  int v;
  int readIndex;
  int start;

  public HashEntry()
  {
    if (!BuildConsSeq.init) BuildConsSeq.init();
    v=BuildConsSeq.andValue;
  }

  public int compareTo(Object o)
  {
    HashEntry e=(HashEntry)o;
    if (this.v<e.v) return -1;
    if (this.v>e.v) return 1;
    return 0;
  }

  public String toString() { return readIndex+"\t"+start+"\t"+v+"\t"+BuildConsSeq.hash2String(v); }

}

class HashComp
{
  int read1, read2;
  int offset;

  public HashComp(int read1, int read2, int offset)
  {
    this.read1=read1;
    this.read2=read2;
    this.offset=offset;
  }

  public boolean equals(Object o)
  {
    if (o ==null) return false;
    if (o == this) return true;
    if (o.getClass() != getClass()) return false;

    HashComp h=(HashComp)o;
    if ((this.read1==h.read1 && this.read2==h.read2) && this.offset==h.offset) return true;
    else return false;
  }

  public int hashCode()
  {
    int v=0;
    v=v|read1;
    v=v<<12;
    v=v|read2;
    v=v<<12;
    v=v|offset;
    return v;
  }
}

class HashCompCount implements Comparable<HashCompCount>
{
  int read1, read2;
  int offset;
  int count=0;

  public HashCompCount(int read1, int read2, int offset, int count)
  {
    this.read1=read1;
    this.read2=read2;
    this.offset=offset;
    this.count=count;
  }

  public int compareTo(HashCompCount h)
  {
    if (this.count>h.count) return -1;
    if (this.count<h.count) return 1;
    return 0;
  }
}

class ConsBuild
{
  int n=0;
  int readids[]=null;
  int offsets[]=null;
  char[] cons=null;
  int minOffset=0;
  int lastMod=-1;

  public ConsBuild(int maxLength, HashCompCount h, char[][] seqs)
  {
    readids=new int[maxLength]; offsets=new int[maxLength];
    n=2;
    readids[0]=h.read1; offsets[0]=0;
    readids[1]=h.read2; offsets[1]=h.offset;
    cons=getCons(seqs);
  }

  public boolean addCons(int readid, char[][] seqs, int relOffset)
  {

    OverlapScore s=getOverlapScore(readid, seqs, relOffset);
    if (BuildConsSeq.scoreFail(s)) return false;
    else
    {
      readids[n]=readid;
      offsets[n]=relOffset;
      n++;
      cons=getCons(seqs);
      return true;
    }
  }

  public char[] getCons(char[][] seqs)
  {
    int minV=Integer.MAX_VALUE, maxV=Integer.MIN_VALUE;
    for (int c1=0; c1<n; c1++)
    {
      minV=Math.min(minV, offsets[c1]);
      maxV=Math.max(maxV, offsets[c1]+seqs[readids[c1]].length);
    }
    minOffset=minV;
    int cons[][]=new int[maxV-minV][];
    for (int c1=0; c1<cons.length; c1++) { cons[c1]=new int[5]; Arrays.fill(cons[c1], 0); }

    for (int c1=0; c1<n; c1++)
    {
      for (int c2=0; c2<seqs[readids[c1]].length; c2++)
      {
        int pos=c2+offsets[c1]-(minV);
        cons[pos][BuildConsSeq.c2i[seqs[readids[c1]][c2]]]++;
      }
    }
    char[] result=new char[maxV-minV]; Arrays.fill(result, '.');
    for (int c1=0; c1<cons.length; c1++)
    {
      int maxIndex=-1, maxCount=0, maxCountCount=0;
      for (int c2=0; c2<cons[c1].length; c2++)
      {
        if (cons[c1][c2]>maxCount) { maxCount=cons[c1][c2]; maxIndex=c2; maxCountCount=1; }
        else if (cons[c1][c2]==maxCount) maxCountCount++;
      }
      if (maxIndex==-1) { System.out.println("WARNING: maxIndex is -1  "+ c1+"\t"+minV+"\t"+maxV); }
      else if (maxCountCount>1) result[c1]='N';
      else result[c1]=BuildConsSeq.i2c[maxIndex];
    }
    return result;
  }

  private OverlapScore getOverlapScore(int readid, char[][] seqs, int offset)
  {
    int olLength=0, matchLength=0;
    int maxAllowedMis=(int)(((double)seqs[readid].length)*(1-BuildConsSeq.MIN_MATCH_PERCENT));
    for (int curPos=0; curPos<seqs[readid].length; curPos++)
    {
      int pos=curPos+offset-minOffset;
      if (pos<0 || pos>=cons.length) continue;
      olLength++;
      if (cons[pos]=='N' || cons[pos]==seqs[readid][curPos]) matchLength++;
      else if (olLength-matchLength>maxAllowedMis) return new OverlapScore(-1, -1);
    }
    return new OverlapScore(olLength, matchLength);

  }

  public void print(char[][] seqs)
  {

    int minOffset=Integer.MAX_VALUE; for (int c1=0; c1<n; c1++) minOffset=Math.min(minOffset, offsets[c1]);
    if (minOffset<0) minOffset=-minOffset; else minOffset=0;
    for (int c1=0; c1<n; c1++)
    {
      for (int c2=0; c2<offsets[c1]+minOffset; c2++) System.out.print(" ");
      for (int c2=0; c2<seqs[readids[c1]].length; c2++) System.out.print(seqs[readids[c1]][c2]);
      System.out.println();
    }
  }
}

class OverlapScore
{
  int olLength=0;
  int matchLength=0;

  public OverlapScore(int olLength, int matchLength) { this.olLength=olLength; this.matchLength=matchLength; }
}


class HashReadCons
{
  int readid, consid;

  public HashReadCons(int readid, int consid)
  {
    this.readid=readid;
    this.consid=consid;
  }

  public boolean equals(Object o)
  {
    if (o ==null) return false;
    if (o == this) return true;
    if (o.getClass() != getClass()) return false;

    HashReadCons h=(HashReadCons)o;
    if (this.readid==h.readid && this.consid==consid) return true;
    else return false;
  }

  public int hashCode()
  {
    int v=0;
    v=v|readid;
    v=v<<12;
    v=v|consid;
    return v;
  }
}
