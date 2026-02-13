package msa;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ariyaratnep
 * Date: 10/23/12
 * Time: 11:39 AM
 * To change this template use File | Settings | File Templates.
 */

class ConsSeq
{
  char[] seq;
  int[] forb;
  int[] againstb;
  int minPos=-1, maxPos=-1;
  int[] mapids=null;
  String seqStr=null;

  public ConsSeq(int n)
  {
    seq=new char[n]; Arrays.fill(seq, ' ');
    forb=new int[n]; Arrays.fill(forb, 0);
    againstb=new int[n]; Arrays.fill(againstb, 0);
  }

  public ConsSeq(char[] seq)
  {
    this.seq=seq;
    this.seqStr=new String(seq);
  }
}


class SWMap implements Comparable
{
  String id;
  int target;
  int pos;
  String tag;

  public SWMap() { }

  public SWMap(SWMap m)
  {
    this.id=m.id;
    this.target=m.target;
    this.pos=m.pos;
    this.tag=m.tag;
  }

  public SWMap(String line)
  {
    String fields[]=line.split("\t");
    id=fields[0];
    target=Integer.parseInt(fields[1]);
    pos=Integer.parseInt(fields[2]);
    tag=fields[3].toUpperCase();
  }

  public int compareTo(Object o)
  {
    SWMap m=(SWMap)o;
    if (this.target<m.target) return -1;
    if (this.target>m.target) return 1;
    if (this.pos<m.pos) return -1;
    if (this.pos>m.pos) return 1;
    return 0;
  }

  public static SWMap[] toArray(Vector temp)
  {
    SWMap[] result=new SWMap[temp.size()];
    for (int c1=0; c1<result.length; c1++) result[c1]=(SWMap)(temp.elementAt(c1));
    Arrays.sort(result);
    return result;
  }

  public String toString() { return id+"\t"+target+"\t"+pos+"\t"+tag; }
}