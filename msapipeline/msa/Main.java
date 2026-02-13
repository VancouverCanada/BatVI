package msa;

/**
 * Created by IntelliJ IDEA.
 * User: ariyaratnep
 * Date: 10/23/12
 * Time: 11:41 AM
 * To change this template use File | Settings | File Templates.
 */

import java.util.*;
import java.io.*;

public class Main
{
  public static void main(String args[]) throws IOException
  {
    Vector<SWMap> temp=new Vector<SWMap>();
    BufferedReader fileIn=new BufferedReader(new InputStreamReader(new FileInputStream(args[0]))); String line; String Cluster;
    PrintWriter fileOut=new PrintWriter(new FileOutputStream(args[0]+".out"));
    while ((Cluster=fileIn.readLine())!=null)
    {
	    temp.clear();
	    String id="1";
	    int Lines_Read=0;
	    while (!(line=fileIn.readLine()).equals("--"))
	    {

		    if (line.startsWith(">")) id=line.substring(1);
		    else
		    {
			    String tags[]=line.split("\t");
			    for (int c1=0; c1<tags.length; c1++)
			    {
				    Lines_Read++;
				    if (Lines_Read>1000) break;
				    SWMap m=new SWMap(); m.id=id+"_"+c1; m.target=1; m.pos=0; m.tag=tags[c1].toUpperCase();
				    temp.addElement(m);
			    }
		    }
	    }
	    SWMap[] maps=SWMap.toArray(temp);
	    System.err.println(maps.length+" reads loaded");
	    ConsSeq[] s=BuildConsSeq.getConsSeq(maps);

	    fileOut.println("*\n"+Cluster+"\n*");
	    for (int c1=0; c1<s.length; c1++)
	    {
		    fileOut.println(">"+maps[s[c1].mapids[0]].id+"_"+ s[c1].mapids.length);
		    fileOut.println(s[c1].seqStr);
	    }
    }

    fileOut.close();
  }
}
