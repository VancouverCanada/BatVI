package msa;

/**
 * 2013.10.16
 * Guan Peiyong
 * for Chandana Virus project.
 * The program ouputs <input>.log file and <input>.out file
 * which keeps the timing stats and assembly information.
 */

import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Vector;

public class BatchMain {
    private final static int MAX_NO_READS = 250;
    private final static int MIN_NO_READS = 2;

    public static void main(String args[]) throws IOException {
        String fileInName = args[0];
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(fileInName)));
        PrintWriter fileOut = new PrintWriter(new FileOutputStream(fileInName + ".out"));
        PrintWriter fileLog = new PrintWriter(new FileOutputStream(fileInName + ".log"));

        Vector<SWMap> temp = new Vector<SWMap>();
        String line;

        boolean bHeader = true;

				int counter = 0;
				
        while ((line = fileIn.readLine()) != null) {
        		
        		counter++;
        		System.out.println(counter);
        		
            SWMap m1 = getSWMap(line);
            String currId = m1.id;
            String prevId = currId;
            temp.addElement(m1);

            int sn = 1;

            while (prevId != null && prevId.equals(currId)) {
                line = fileIn.readLine();

								if(line==null) return;
								counter++;
								
        				System.out.println(counter);
        		
                SWMap m2 = getSWMap(line);
                currId = m2.id;
                if(prevId==null) return;
                
                if (prevId.equals(currId)) {
                    temp.addElement(m2);

                    if (temp.size() >= MAX_NO_READS) {
                        buildConsSeq(temp, sn, fileOut, fileLog, bHeader);
                        sn++;
                        if (bHeader) bHeader = false;
                    }

                } else {
                    buildConsSeq(temp, sn, fileOut, fileLog, bHeader);
                    sn++;
                    if (bHeader) bHeader = false;

                    temp.clear();
                    temp.addElement(m2);
                }
            }
        }
        fileIn.close();
        fileOut.close();
        fileLog.close();
    }

    private static SWMap getSWMap(String line) {
        String[] split = line.split("\t");
        String id = split[0];
        String seq = split[1];
        SWMap m = new SWMap();
        m.id = id;
        m.target = 1;
        m.pos = 0;
        m.tag = seq.toUpperCase();
        return m;
    }

    private static String getDateTime(Date date) {
        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        return dateFormat.format(date);
    }

    private static void buildConsSeq(Vector<SWMap> reads, int sn, PrintWriter fileOut, PrintWriter fileLog, boolean bHeader) {

        if (reads.size() >= MIN_NO_READS) {
            if (bHeader) {
                fileLog.println("currId" + "\t" + "sn" + "\t" + "start" + "\t" + "end" + "\t" + "contigCnt" + "\t" + "clusterSize" + "\t" + "duration");
            }

            Date start = new Date();
            long startM = System.currentTimeMillis();

            SWMap[] maps = SWMap.toArray(reads);
            ConsSeq[] s = BuildConsSeq.getConsSeq(maps);
            for (int c1 = 0; c1 < s.length; c1++) {
                //>chr1:950802_21_9_22_500_2
                //><clusterID>_<clusterSplitSN>_<contigSN>_<contigCnt>_<clusterSize>_<readsUsedForAssemblyCnt>
                fileOut.println(">" + maps[s[c1].mapids[0]].id + "_" + sn + "_" + (c1 + 1) + "_" + s.length + "_" + maps.length + "_" + s[c1].mapids.length);
                fileOut.println(s[c1].seqStr);
            }

            Date end = new Date();
            long endM = System.currentTimeMillis();
            fileLog.println(reads.elementAt(0).id + "\t" + sn + "\t" + getDateTime(start) + "\t" + getDateTime(end) + "\t" + s.length + "\t" + maps.length + "\t" + (endM - startM));

        } else {
            System.out.println("MIN_NO_READS...");
        }

        reads.clear();
        fileOut.flush();
        fileLog.flush();
    }
}

