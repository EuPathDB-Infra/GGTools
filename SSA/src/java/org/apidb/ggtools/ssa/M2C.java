package org.apidb.ggtools.ssa;

/*
Conceived and implemented by Gregory Robert Grant, University of Pennsylvania, Philadelphia, PA.
<ggrant@pcbi.upenn.edu>, <greg@grant.org>
Copyright 2009
All rights reserved.  
*/

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Random;

public class M2C {
    public static void main(String[] args) {
	System.err.println("\nWritten by Gregory R. Grant\nUniversity of Pennsylvania, 2010.");
	if(Array.getLength(args) < 3) {
	    if(Array.getLength(args) == 1) {
		if(args[0].equals("-example")) {
		    printExample();
		    System.exit(0);
		}
	    }
	    printUsage();
	    System.exit(0);
	}
	boolean FLAG = true;
	int numchunks = 1;
	BufferedWriter out_log = null;
	try {
	    long genomelength = 0;
	    int maxUniqueReadsOneSpan = 0;
	    int numUniqueReads = 0;
	    //long heapSize = Runtime.getRuntime().totalMemory();
	    //long heapMaxSize = Runtime.getRuntime().maxMemory();
	    boolean windows = false;
	    int strandIndex = -1;
	    int maxSpan = 0;
	    int normalizeToThisNumberOfReads = 0;
	    boolean nodensity = false;
	    String name = args[0];
	    int count_cutoff=1;
	    int windowSize = 1;
	    int windowDisplacement = 1;
	    int startCoordinate = 1;
	    boolean rightEndInclusion = true;
	    int startCoordinate_outfile = 1;
	    boolean rightEndInclusion_outfile = true;
	    File maskFile = new File("");
	    File repeatsFile = new File("");
	    boolean use_mask = false;
	    boolean repeat_mask = false;
	    boolean removeIdenticalReads = false;
	    boolean output_filtered_reads = false;
	    int extensionLength = 0;
	    boolean ucsc = false;
	    if(Array.getLength(args) >= 4) {
		boolean option_recognized = false;
		for(int i = 3; i<Array.getLength(args); i++) {
		    if(args[i].equals("-outputfilteredreads")) {
			output_filtered_reads = true;
			option_recognized = true;
		    }
		    if(args[i].equals("-windows")) {
			windows = true;
			option_recognized = true;
		    }
		    if(args[i].equals("-removeidenticalreads")) {
			removeIdenticalReads = true;
			option_recognized = true;
		    }
		    if(args[i].equals("-nodensity")) {
			nodensity = true;
			option_recognized = true;
		    }
		    if(args[i].equals("-ucsc")) {
			rightEndInclusion_outfile = false;
			startCoordinate_outfile = 0;
			ucsc = true;
			option_recognized = true;
		    }
		    if(args[i].equals("-openinterval_infile")) {
			rightEndInclusion = false;
			option_recognized = true;
		    }
		    if(args[i].equals("-start_coordinate_infile")) {
			try {
			    startCoordinate = Integer.parseInt(args[i+1]);
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -start_coordinate_infile must be followed by a 0 or 1, (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			if(startCoordinate != 0 && startCoordinate != 1) {
			    System.err.println("\nERROR: the argument -start_coordinate_infile must be followed by a 0 or 1, (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-openinterval_outfile")) {
			rightEndInclusion_outfile = false;
			option_recognized = true;
		    }
		    if(args[i].equals("-start_coordinate_outfile")) {
			try {
			    startCoordinate_outfile = Integer.parseInt(args[i+1]);
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -start_coordinate_outfile must be followed by a 0 or 1, (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			if(startCoordinate_outfile != 0 && startCoordinate_outfile != 1) {
			    System.err.println("\nERROR: the argument -start_coordinate_outfile must be followed by a 0 or 1, (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-ucsc")) {
			rightEndInclusion_outfile = false;
			startCoordinate_outfile = 0;
			option_recognized = true;
		    }
		    if(args[i].equals("-chunks")) {
			try {
			    numchunks = Integer.parseInt(args[i+1]);
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -numchunks must be followed by a non-negative integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			if(numchunks < 0) {
			    System.err.println("\nERROR: the argument -numchunks must be followed by a non-negative integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-countcutoff")) {
			try {
			    count_cutoff = Integer.parseInt(args[i+1]);
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -countcutoff must be followed by a non-negative integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			if(count_cutoff < 0) {
			    System.err.println("\nERROR: the argument -countcutoff must be followed by a non-negative integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-name")) {
			name = args[i+1];
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-extensionlength")) {
			try {
			    extensionLength = Integer.parseInt(args[i+1]);
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -extensionlength must be followed by a non-negative integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			if(count_cutoff < 0) {
			    System.err.println("\nERROR: the argument -extensionlength must be followed by a non-negative integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-windowsize")) {
			try {
			    windowSize = Integer.parseInt(args[i+1]);
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -windowsize must be followed by a positive integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			if(windowSize < 1) {
			    System.err.println("\nERROR: the argument -windowsize must be followed by a positive integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-windowdisplacement")) {
			try {
			    windowDisplacement = Integer.parseInt(args[i+1]);
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -windowdisplacement must be followed by a positive integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			if(windowDisplacement < 1) {
			    System.err.println("\nERROR: the argument -windowsize must be followed by a positive integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-strandindex")) {
			try {
			    strandIndex = Integer.parseInt(args[i+1]);
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -strandindex must be followed by a non-negative integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			if(strandIndex < 0) {
			    System.err.println("\nERROR: the argument -strandindex must be followed by a non-negative integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-maxspan")) {
			try {
			    maxSpan = Integer.parseInt(args[i+1]);
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -maxspan must be followed by a non-negative integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			if(maxSpan < 1) {
			    System.err.println("\nERROR: the argument -maxspan must be followed by a positive integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-norm")) {
			try {
			    normalizeToThisNumberOfReads = Integer.parseInt(args[i+1]);
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -norm must be followed by a non-negative integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			if(normalizeToThisNumberOfReads < 1) {
			    System.err.println("\nERROR: the argument -norm must be followed by a positive integer (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-mask1")) {
			try {
			    maskFile = new File(args[i+1]);
			    use_mask = true;
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -mask2 must be followed by a valid file name (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(args[i].equals("-mask2")) {
			try {
			    repeatsFile = new File(args[i+1]);
			    repeat_mask = true;
			} catch (Exception E2) {
			    System.err.println("\nERROR: the argument -mask2 must be followed by a valid file name (you have \"" + args[i+1] + "\")\n\n");
			    System.exit(1);
			}
			i++;
			option_recognized = true;
		    }
		    if(!option_recognized) {
			System.err.println("\nERROR: the argument \"" + args[i] + "\" is not recognized.\n\n");
			System.exit(1);
		    }
		    
		}
	    }
	    File sampleFile = new File("");
	    File outputFile_log = new File("");
	    
	    if(windowSize > 10000000 || windowDisplacement > 10000000) {
		System.err.println("window size and window displacement must be less or equal to 10,000,000");
		System.exit(0);
	    }
	    sampleFile = new File(args[0]);
	    String logFileName = "";
	    args[1] = args[1].replaceAll("[^a-zA-Z0-9]txt$", "");
	    args[1] = args[1].replaceAll("[^a-zA-Z0-9]Txt$", "");
	    args[1] = args[1].replaceAll("[^a-zA-Z0-9]TXT$", "");
	    logFileName = args[2];
	    outputFile_log = new File(logFileName);
	    out_log = new BufferedWriter(new FileWriter(outputFile_log));
		
	    System.err.println("");
	    System.err.println("Assessing sample file...");
	    Object[] fileSpecs = assessFile(sampleFile, rightEndInclusion, strandIndex);
	    System.err.println("Done Assessing sample file...");
	    int numChromosomes_all = (Integer) fileSpecs[0];
	    System.err.println("\n" + numChromosomes_all + " chromosomes");
	    int maxReadsOneSpan = (Integer) fileSpecs[1];
	    int numReads_sample = (Integer) fileSpecs[2];
	    strandIndex = (Integer) fileSpecs[3];
	    int skip = (Integer) fileSpecs[4];
	    @SuppressWarnings("unchecked")
	    Hashtable<String, Integer> chrName2Number = (Hashtable<String, Integer>) fileSpecs[5];
	    String[] chrNumber2Name = (String[]) fileSpecs[6];
	    
	    out_log.write("input file = " + args[0] + "\n");
	    out_log.write("number of headers lines in sample file: " + skip + "\n");
	    out_log.write("column " + strandIndex + " in sample file gives the strand (counting from 0)\n");
	    out_log.write("extension length = " + extensionLength + "\n");
	    out_log.write("start coordinate of infile = " + startCoordinate + "\n");
	    out_log.write("right end inclusion in infile = " + rightEndInclusion + "\n");
	    out_log.write("start coordinate of outfile = " + startCoordinate_outfile + "\n");
	    out_log.write("right end inclusion in outfile = " + rightEndInclusion_outfile + "\n");
	    out_log.write("num chromosomes = " + numChromosomes_all + "\n");
	    String mm = format(maxReadsOneSpan);
	    out_log.write("max reads on one chromosome in the sample file = " + mm + "\n");
	    mm = format(numReads_sample);
	    out_log.write("num reads in sample file = " + mm + "\n");
	    out_log.write("remove identical reads = " + removeIdenticalReads + "\n");
	    if(repeat_mask) {
		out_log.write("mask1 = true\n");
		out_log.write("mask1 fie = " + repeatsFile + "\n");
	    }
	    else {
		out_log.write("mask1: none\n");
	    }
	    if(use_mask) {
		out_log.write("mask2: " + maskFile + "\n");
	    }
	    else {
		out_log.write("mask2: none\n");
	    }
	    out_log.write("chromosomes:\nnum\tname\tnum reads\n");
	    out_log.flush();
	    System.err.println("");
	    int num_chromosomes_per_chunk = 1;
	    int chr_start_num = 0;
	    File outputFile_density = new File(args[1]);
	    BufferedWriter out_density = new BufferedWriter(new FileWriter(outputFile_density));
	    while(FLAG) {
		FLAG = false;
		out_density.close();
		outputFile_density = new File(args[1]);
		out_density = new BufferedWriter(new FileWriter(outputFile_density));
		try {
		    if(ucsc) {
			out_density.write("track type=bedGraph name=\"" + name  + "\" description=\"" + name + "\" visibility=full color=255,0,0 priority=20\n");
		    }
		    for(int chunk=0; chunk<numchunks; chunk++) {
			//int chunk2 = chunk + 1;
			if(numChromosomes_all % numchunks == 0)
			    num_chromosomes_per_chunk = numChromosomes_all / numchunks;
			else
			    num_chromosomes_per_chunk = numChromosomes_all / numchunks + 1;
			if(num_chromosomes_per_chunk > numChromosomes_all) {
			    num_chromosomes_per_chunk = numChromosomes_all;
			}
			chr_start_num = chunk * num_chromosomes_per_chunk;
			int chr_end_num = (chunk+1) * num_chromosomes_per_chunk - 1;
			if(chr_end_num > numChromosomes_all - 1) {
			    chr_end_num = numChromosomes_all - 1;
			}
			//int chr_start_num2 = chr_start_num+1;
			//int chr_end_num2 = chr_end_num+1;
			int numChromosomes = chr_end_num - chr_start_num + 1;
			System.err.println("Reading the sample file... trying to read " + numChromosomes + " chromosomes.");
			Object[] readfileOutput = readFile(sampleFile, numChromosomes, maxReadsOneSpan, numReads_sample, extensionLength, strandIndex, removeIdenticalReads, skip, rightEndInclusion, chrName2Number, chr_start_num, chr_end_num);
			int[][] locationArray = (int[][])readfileOutput[0];
			int[] readCnt = (int[])readfileOutput[1];
			for(int i=chr_start_num; i<=chr_end_num; i++) {
			    int ii=i+1;
			    out_log.write("  " + ii + "\t" + chrNumber2Name[i] + "\t" + readCnt[i-chr_start_num] + "\n");
			}
			long genomelength_sample = (Long)readfileOutput[2];
			long[] chrMin_sample = ((long[])readfileOutput[3]);
			long[] chrMax_sample = ((long[])readfileOutput[4]);
			int num_unmasked_sample = 0;
			int maxUniqueReadsOneSpan2 = ((Integer)readfileOutput[5]);
			if(maxUniqueReadsOneSpan2 > maxUniqueReadsOneSpan) {
			    maxUniqueReadsOneSpan = maxUniqueReadsOneSpan2;
			}
			int[][] readLength = (int[][])readfileOutput[6];
			int maxReadLength = (Integer)readfileOutput[7];
			for(int i=0; i<numChromosomes; i++) {
			    numUniqueReads = numUniqueReads + readCnt[i];
			}
			
			int totalLength = extensionLength + maxReadLength;
			genomelength = genomelength_sample + genomelength;
			int[][][] locationArray_unextended = new int[0][0][0];
			
			if(repeat_mask) {
			    //long repeat_length_sample = 0;
			    System.err.println("");
			    System.err.println("Parsing the sample file...");
			    readfileOutput = readFileReportStrand(sampleFile, numChromosomes, maxReadsOneSpan, numReads_sample, 0, strandIndex, removeIdenticalReads, skip, chrName2Number);
			    locationArray_unextended = (int[][][])readfileOutput[0];
			    System.err.println("");
			    System.err.println("Masking repeats in the sample file...");
			    Object[] repeatMaskOutput = repeatMask(locationArray_unextended, readCnt, repeatsFile, numChromosomes, readLength, maxReadsOneSpan, extensionLength, chrMin_sample, chrMax_sample, chrNumber2Name, chrName2Number);
			    locationArray = (int[][])repeatMaskOutput[0];
			    readCnt = (int[])repeatMaskOutput[1];
			    num_unmasked_sample = (Integer)repeatMaskOutput[2];
			    readLength = (int[][])repeatMaskOutput[3];
			    mm = format(num_unmasked_sample);
			    out_log.write("num reads in the sample after masking repeats: " + mm + "\n");
			    out_log.flush();
			    mm = format(genomelength);
			    out_log.flush();
			}
			
			numReads_sample = 0;
			for(int c=0; c < numChromosomes; c++) {
			    numReads_sample = numReads_sample + readCnt[c];
			}
			
			// Note, the difference between masking repeats (called mask1) and masking other regions (called mask2) of bias is that reads are repeat masked if their unextended reads overlap repeat regions, while other masking is done if the extended read overlaps the region.
			if(use_mask) {
			    Object[] returnArray = assessRepeatsFile(maskFile);
			    int maxRegionsOneSpan = (Integer)returnArray[0];
			    int numRegions = (Integer)returnArray[1];
			    int[] regionCnt = (int[])returnArray[2];
			    int numSpans_regions = (Integer)returnArray[3];
			    @SuppressWarnings("unchecked")
			    Hashtable<String, Integer> chrName2Number_mask = (Hashtable<String, Integer>) returnArray[4];
			    String[] chrNumber2Name_mask = (String[])returnArray[5];
			    
			    returnArray = readMaskFile(maskFile, numSpans_regions, maxRegionsOneSpan, numRegions, chrName2Number_mask);
			    int[][][] locations_regions = (int[][][])returnArray[0];
			    regionCnt = (int[])returnArray[1];
			    Object[] maskRegions_output = maskRegions(locationArray, readCnt, readLength, locations_regions, regionCnt, extensionLength, chrNumber2Name, chrName2Number, chrNumber2Name_mask, chrName2Number_mask);
			    locationArray = (int[][])maskRegions_output[0];
			    readCnt = (int[])maskRegions_output[1];
			    num_unmasked_sample = (Integer)maskRegions_output[2];
			    readLength = (int[][])maskRegions_output[3];
			    mm = format(num_unmasked_sample);
			    out_log.write("number reads in sample after filtering out regions in the file: " + maskFile + ": " + mm + "\n");
			    out_log.flush();
			}
			
			// normalize to command line count argument number of reads 
			if(normalizeToThisNumberOfReads > 0) {
			    int diff = numReads_sample - normalizeToThisNumberOfReads;
			    Object[] normalize_return = normalize(locationArray, readLength, readCnt, diff);
			    locationArray = (int[][])normalize_return[0];
			    readLength = (int[][])normalize_return[1];
			    readCnt = (int[])normalize_return[2];
			    numReads_sample = 0;
			    for(int c=0; c<numChromosomes; c++) {
				numReads_sample = numReads_sample + readCnt[c];
			    }
			    mm = format(numReads_sample);
			    out_log.write("number of reads in the sample after normalization: " + mm + "\n");
			    out_log.flush();
			}
			
			//String chr = "";
			if(output_filtered_reads == true) {
			    File outputFile_p = new File("filtered_reads.txt");
			    BufferedWriter out_p = new BufferedWriter(new FileWriter(outputFile_p));
			    for(int i=0; i<numChromosomes; i++) {
				for(int j=0; j<readCnt[i]; j++) {
				    long loc = locationArray[i][j];
				    long loc2 = loc + readLength[i][j];
				    out_p.write(chrNumber2Name[i] + "\t" + loc + "\t" + loc2 + "\n");
				}
			    }
			    out_p.close();
			}
			
			if(!nodensity) {
			    System.err.println("");
			    writeDensityFile(locationArray, readLength, numChromosomes, readCnt, out_density, startCoordinate, startCoordinate_outfile, rightEndInclusion_outfile, chrMin_sample, chrMax_sample, windowSize, windowDisplacement, totalLength, count_cutoff, chrNumber2Name, windows, maxSpan, chr_start_num);
			}
		    }
		} catch (OutOfMemoryError e) {
		    if(num_chromosomes_per_chunk == 1) {
			System.err.println("\nERROR: Darn, I ran completely out of memory, processing in chunks didn't work.\nTry rerunning with -Xmx500m option, or higher if you have enough memory. See the\nhelp for more details.\n");
			System.exit(0);
		    }
		    int numchrperchunkold = num_chromosomes_per_chunk;
		    while(numchrperchunkold <= num_chromosomes_per_chunk && num_chromosomes_per_chunk > 1) {
			numchunks++;
			if(numChromosomes_all % numchunks == 0)
			    num_chromosomes_per_chunk = numChromosomes_all / numchunks;
			else
			    num_chromosomes_per_chunk = numChromosomes_all / numchunks + 1;
			if(num_chromosomes_per_chunk > numChromosomes_all) {
			    num_chromosomes_per_chunk = numChromosomes_all;
			}
		    }
		    System.err.println("\nWARNING: Ran out of memory.  Going to try to reprocess in " + numchunks + " chunks, but you might\nwant to try rerunning with the -Xmx500m option, or higher if you have enough\nmemory.  See the help for more details.\n");
		    FLAG = true;
		}
	    }
	    mm = format(numUniqueReads);
	    out_log.write("num unique reads in sample = " + mm + "\n");
	    mm = format(maxUniqueReadsOneSpan);
	    out_log.write("max unique reads on one chromosome in sample = " + mm + "\n");
	    mm = format(genomelength);
	    out_log.write("genomelength = " + mm + "\n");
	    out_log.flush();
	} catch (Exception e) { 
	    printUsage();
	    e.printStackTrace(System.err);
	    System.exit(0);
	} finally {
	  closeQuietly(out_log);
	}
    }

    private static void closeQuietly(BufferedWriter writer) {
      if (writer != null)
        try { writer.close(); }
        catch (IOException e) { }
    }

    public static Object[] repeatMask(int[][][] locations_reads, int[] readCnt, File repeatsFile, int numChromosomes, int[][] readLength, int maxReadsOneSpan, int extensionLength, long[] chrMin, long[] chrMax, String[] chrNumber2Name, Hashtable<String, Integer> chrName2Number) {
	int numSpans = Array.getLength(locations_reads);
	Object[] returnArray = new Object[4];
	returnArray = assessRepeatsFile(repeatsFile);
	int maxRepeatsOneSpan = (Integer)returnArray[0];
	int numRepeats = (Integer)returnArray[1];
	int[] repeatCnt = (int[])returnArray[2];
	int numSpans_repeats = (Integer)returnArray[3];
	@SuppressWarnings("unchecked")
	Hashtable<String, Integer> chrName2Number_repeats = (Hashtable<String, Integer>) returnArray[4];
	//String[] chrNumber2Name_repeats = (String[])returnArray[5];
	returnArray = readMaskFile(repeatsFile, numSpans_repeats, maxRepeatsOneSpan, numRepeats, chrName2Number_repeats);
	int[][][] locations_repeats = (int[][][])returnArray[0];
	repeatCnt = (int[])returnArray[1];

	boolean masked = false;
	int flag = 0;
	int[][] unmasked_reads = new int[numSpans][maxReadsOneSpan];
	int[][] unmasked_readLength = new int[numSpans][maxReadsOneSpan];
	int[] num_unmasked_per_span = new int[numSpans];
	int max_unmasked_one_span = 0;
	for(int span = 0; span<numSpans; span++) {
	    num_unmasked_per_span[span] = 0;
	}
	int num_unmasked = 0;
	for(int span = 0; span<numSpans; span++) {
	    if(chrName2Number_repeats.containsKey(chrNumber2Name[span])) {
		if(readCnt[span] > 0) {
		    int[][] temparray_reads = new int[readCnt[span]][3];
		    for(int i=0; i<readCnt[span]; i++) {
			temparray_reads[i][0] = locations_reads[span][i][0];
			temparray_reads[i][1] = locations_reads[span][i][1];
			temparray_reads[i][2] = readLength[span][i];
		    }
		    Arrays.sort(temparray_reads, FirstIndexComparator);
		    int start = 0;
		    int repeatcounter = 0;
		    for(int read=0; read<readCnt[span]; read++) {
			flag=0;
			repeatcounter=start;
			masked = false;
			while(flag == 0 && repeatcounter < repeatCnt[chrName2Number_repeats.get(chrNumber2Name[span])]) {
			    if(temparray_reads[read][0] + readLength[span][read] - 1 >= locations_repeats[chrName2Number_repeats.get(chrNumber2Name[span])][repeatcounter][0] && temparray_reads[read][0] <= locations_repeats[chrName2Number_repeats.get(chrNumber2Name[span])][repeatcounter][1]) {
				masked = true;
				flag = 1;
			    }
			    if(temparray_reads[read][0] + readLength[span][read] - 1 < locations_repeats[chrName2Number_repeats.get(chrNumber2Name[span])][repeatcounter][0]) {
				flag=1;
			    }
			    if(temparray_reads[read][0] > locations_repeats[chrName2Number_repeats.get(chrNumber2Name[span])][repeatcounter][1]) {
				start++;
				if(start == numRepeats) {
				    flag = 1;
				}
			    }
			    repeatcounter++;
			}
			if(!masked) {
			    unmasked_readLength[span][num_unmasked_per_span[span]] = temparray_reads[read][2];
			    if(temparray_reads[read][1] == 1) {
				unmasked_reads[span][num_unmasked_per_span[span]] = temparray_reads[read][0];
			    }
			    if(temparray_reads[read][1] == -1) {
				unmasked_reads[span][num_unmasked_per_span[span]] = temparray_reads[read][0] - extensionLength;
			    }
			    num_unmasked++;
			    num_unmasked_per_span[span]++;
			    if(num_unmasked_per_span[span] >  max_unmasked_one_span) {
				max_unmasked_one_span++;
			    }
			}
		    }
		}
	    }
	}
	int[][] unmasked_reads_return = new int[numSpans][max_unmasked_one_span];
	int[][] unmasked_readLength_return = new int[numSpans][max_unmasked_one_span];
	for(int span=0; span<numSpans; span++) {
	    for(int i=0; i<num_unmasked_per_span[span]; i++) {
		unmasked_reads_return[span][i] = unmasked_reads[span][i];
		unmasked_readLength_return[span][i] = unmasked_readLength[span][i];
	    }
	}

	returnArray = new Object[5];
	returnArray[0] = unmasked_reads_return;
	returnArray[1] = num_unmasked_per_span;
	returnArray[2] = num_unmasked;
	returnArray[3] = unmasked_readLength_return;
	return returnArray;
    }

    public static Object[] maskRegions(int[][] locations_reads, int[] readCnt, int[][] readLength, int[][][] locations_regions, int[] regionsCnt, int extensionLength, String[] chrNumber2Name, Hashtable<String, Integer> chrName2Number, String[] chrNumber2Name_regions, Hashtable<String, Integer> chrName2Number_regions) {
	int numSpans_regions = Array.getLength(regionsCnt);
	int numRegions = 0;
	int maxRegionsOneSpan = 0;
	for(int i=0; i<numSpans_regions; i++) {
	    if(regionsCnt[i] > maxRegionsOneSpan) 
		maxRegionsOneSpan = regionsCnt[i];
	    numRegions = numRegions + regionsCnt[i];
	}
	int numSpans = Array.getLength(locations_reads);
	int maxReadsOneSpan = 0;
	for(int i=0; i<numSpans; i++) {
	    if(readCnt[i] > maxReadsOneSpan) 
		maxReadsOneSpan = readCnt[i];
	}
	boolean masked = false;
	int flag = 0;
	int[][] unmasked_reads = new int[numSpans][maxReadsOneSpan];
	int[][] unmasked_readLength = new int[numSpans][maxReadsOneSpan];
	int[] num_unmasked_per_span = new int[numSpans];
	int max_unmasked_one_span = 0;
	for(int span = 0; span<numSpans; span++) {
	    num_unmasked_per_span[span] = 0;
	}
	int num_unmasked = 0;
	for(int span = 0; span<numSpans; span++) {
	    if(chrName2Number_regions.containsKey(chrNumber2Name[span]) && readCnt[span] > 0) {
		int[][] temparray_reads = new int[readCnt[span]][2];
		for(int i=0; i<readCnt[span]; i++) {
		    temparray_reads[i][0] = locations_reads[span][i];
		    temparray_reads[i][1] = readLength[span][i];
		}
		Arrays.sort(temparray_reads, FirstIndexComparator);
		int[][] temparray_regions = new int[regionsCnt[chrName2Number_regions.get(chrNumber2Name[span])]][2];
		for(int i=0; i<regionsCnt[chrName2Number_regions.get(chrNumber2Name[span])]; i++) {
		    temparray_regions[i][0] = locations_regions[chrName2Number_regions.get(chrNumber2Name[span])][i][0];
		    temparray_regions[i][1] = locations_regions[chrName2Number_regions.get(chrNumber2Name[span])][i][1];
		}
		Arrays.sort(temparray_regions, FirstIndexComparator);
		//int len = Array.getLength(temparray_regions);
		int start = 0;
		int regioncounter = 0;
		for(int read=0; read<readCnt[span]; read++) {
		    flag=0;
		    regioncounter=start;
		    masked = false;
		    while(flag == 0 && regioncounter < regionsCnt[chrName2Number_regions.get(chrNumber2Name[span])]) {
			if(temparray_reads[read][0] + temparray_reads[read][1] + extensionLength - 1 >= temparray_regions[regioncounter][0] && temparray_reads[read][0] <= temparray_regions[regioncounter][1]) {
			    masked = true;
			    flag = 1;
			}
			if(temparray_reads[read][0] + temparray_reads[read][1] + extensionLength - 1 < temparray_regions[regioncounter][0]) {
			    flag=1;
			}
			if(temparray_reads[read][0] > temparray_regions[regioncounter][1]) {
			    start++;
			    if(start == numRegions) {
				flag = 1;
			    }
			}
			regioncounter++;
		    }
		    if(!masked) {
			unmasked_reads[span][num_unmasked_per_span[span]] = temparray_reads[read][0];
			unmasked_readLength[span][num_unmasked_per_span[span]] = temparray_reads[read][1];
			num_unmasked++;
			num_unmasked_per_span[span]++;
			if(num_unmasked_per_span[span] >  max_unmasked_one_span) {
			    max_unmasked_one_span++;
			}
		    }
		}
	    }
	}
	int[][] unmasked_reads_return = new int[numSpans][max_unmasked_one_span];
	int[][] unmasked_readLength_return = new int[numSpans][max_unmasked_one_span];
	for(int span=0; span<numSpans; span++) {
	    for(int i=0; i<num_unmasked_per_span[span]; i++) {
		unmasked_reads_return[span][i] = unmasked_reads[span][i];
		unmasked_readLength_return[span][i] = unmasked_readLength[span][i];
	    }
	}

	Object[] returnArray = new Object[4];
	returnArray[0] = unmasked_reads_return;
	returnArray[1] = num_unmasked_per_span;
	returnArray[2] = num_unmasked;
	returnArray[3] = unmasked_readLength_return;
	return returnArray;
    }

    private static Object[] readFile(File file, int numChromosomes, int maxReadsOneSpan, int numReads, int extensionLength, int strandIndex, boolean removeIdenticalReads, int skip, boolean rightEndInclusion, Hashtable<String, Integer> chrName2Number, int chr_start_num, int chr_end_num) {

        FileInputStream inFile = null;
	int[][][] locationArray = new int[2][numChromosomes][maxReadsOneSpan];
	int maxReadLength = 0;
	int[][][] readLength = new int[2][numChromosomes][maxReadsOneSpan];
	int[][] readCnt = new int[2][numChromosomes];
	int[] readCnt2 = new int[numChromosomes];
        try {
            inFile = new FileInputStream(file);
        } catch(FileNotFoundException e) {
	    System.err.println("The file " + file + " does not seem to exist...");
            System.exit(1);
        }
	long[] chrMin = new long[numChromosomes];
	long[] chrMax = new long[numChromosomes];
	for(int i=0; i<numChromosomes; i++) {
	    chrMin[i] = 1000000000;
	    chrMax[i] = 0;
	}
	for(int i=0; i<numChromosomes; i++) {
	    readCnt[0][i] = 0;
	    readCnt[1][i] = 0;
	    readCnt2[i] = 0;
	}
	//int flag = 0;
	int chrNum = 0;
        try {
            inFile = new FileInputStream(file);
            BufferedReader br = new BufferedReader(new InputStreamReader(inFile));
            String line = "";
	    for(int i=0; i<skip; i++) {
		line = br.readLine();
	    }
	    int CNTTMP = 0;
            while((line = br.readLine()) != null) {
		CNTTMP++;
		String[] temp = line.split("\\s+", 0);
		String chr = temp[0];
		chrNum = chrName2Number.get(chr) -  chr_start_num;
		if(chrName2Number.get(chr) >= chr_start_num && chrName2Number.get(chr) <= chr_end_num) {
		    if(Array.getLength(temp) < 4) {
			System.err.println("line " + CNTTMP + " is the offending one.");
		    }
		    String strand = temp[strandIndex];
		    if(strand.equals("-")) {
			if(rightEndInclusion) {
			    readLength[0][chrNum][readCnt[0][chrNum]] = Integer.parseInt(temp[2]) - Integer.parseInt(temp[1]) + 1;
			}
			else {
			    readLength[0][chrNum][readCnt[0][chrNum]] = Integer.parseInt(temp[2]) - Integer.parseInt(temp[1]);
			}
			locationArray[0][chrNum][readCnt[0][chrNum]] = Integer.parseInt(temp[1]) - extensionLength;
			if(locationArray[0][chrNum][readCnt[0][chrNum]] < chrMin[chrNum])
			    chrMin[chrNum] = locationArray[0][chrNum][readCnt[0][chrNum]];
			if(locationArray[0][chrNum][readCnt[0][chrNum]] + extensionLength + readLength[0][chrNum][readCnt[0][chrNum]] > chrMax[chrNum]) 
			    chrMax[chrNum] = locationArray[0][chrNum][readCnt[0][chrNum]] + extensionLength + readLength[0][chrNum][readCnt[0][chrNum]] - 1;
			readCnt[0][chrNum]++;
		    }
		    else {
			if(rightEndInclusion) {
			    readLength[1][chrNum][readCnt[1][chrNum]] = Integer.parseInt(temp[2]) - Integer.parseInt(temp[1]) + 1;
			}
			else {
			    readLength[1][chrNum][readCnt[1][chrNum]] = Integer.parseInt(temp[2]) - Integer.parseInt(temp[1]);
			}
			locationArray[1][chrNum][readCnt[1][chrNum]] = Integer.parseInt(temp[1]);
			if(locationArray[1][chrNum][readCnt[1][chrNum]] < chrMin[chrNum])
			    chrMin[chrNum] = locationArray[1][chrNum][readCnt[1][chrNum]];
			if(locationArray[1][chrNum][readCnt[1][chrNum]] + extensionLength + readLength[1][chrNum][readCnt[1][chrNum]] > chrMax[chrNum]) 
			    chrMax[chrNum] = locationArray[1][chrNum][readCnt[1][chrNum]] + extensionLength + readLength[1][chrNum][readCnt[1][chrNum]] - 1;
			readCnt[1][chrNum]++;
		    }
		}
	    }
            inFile.close();
        } catch(IOException e) {
            e.printStackTrace(System.err);
            System.exit(1);
        }
	for(int s=0; s<2; s++) {
	    int l1 = 0;
	    int l2 = 0;
	    int r1 = 0;
	    //int r2 = 0;
	    int j = 0;
	    for(int chr = 0; chr < numChromosomes; chr++) {
		int[][] temparray = new int[readCnt[s][chr]][2];
		for(int i=0; i<readCnt[s][chr]; i++) {
		    temparray[i][0] = locationArray[s][chr][i];
		    temparray[i][1] = readLength[s][chr][i];
		}
		Arrays.sort(temparray, FirstIndexComparator);
		j=0;
		for(int i=0; i<readCnt[s][chr]-1; i++) {
		    l1 = temparray[i][0];
		    l2 = temparray[i+1][0];
		    if(l1 != l2 || !removeIdenticalReads) {
			locationArray[s][chr][j]=l1;
			r1 = temparray[i][1];
			readLength[s][chr][j]=r1;
			j++;
		    }
		}	    
		if(readCnt[s][chr] > 0) {
		    locationArray[s][chr][j]=temparray[readCnt[s][chr]-1][0];
		    readLength[s][chr][j]=temparray[readCnt[s][chr]-1][1];
		j++;
		readCnt[s][chr] = j;
		}
	    }		
	}

	int maxReadsOneSpanFinal = 0;
	for(int chr = 0; chr < numChromosomes; chr++) {
	    if(readCnt[0][chr]+readCnt[1][chr] > maxReadsOneSpanFinal) 
		maxReadsOneSpanFinal = readCnt[0][chr]+readCnt[1][chr];
	}
	int[][] locationArray2 = new int[numChromosomes][maxReadsOneSpanFinal];
	int[][] readLength2 = new int[numChromosomes][maxReadsOneSpanFinal];
	for(int chr = 0; chr < numChromosomes; chr++) {
	    int[][] temparray = new int[readCnt[0][chr]+readCnt[1][chr]][2];
	    for(int i=0; i<readCnt[0][chr]; i++) {
		temparray[i][0] = locationArray[0][chr][i];
		temparray[i][1] = readLength[0][chr][i];
	    }
	    for(int i=0; i<readCnt[1][chr]; i++) {
		temparray[i+readCnt[0][chr]][0] = locationArray[1][chr][i];
		temparray[i+readCnt[0][chr]][1] = readLength[1][chr][i];
	    }
	    Arrays.sort(temparray, FirstIndexComparator);
	    for(int i=0; i<readCnt[0][chr] + readCnt[1][chr]; i++) {
		locationArray2[chr][i] = temparray[i][0];
		readLength2[chr][i] = temparray[i][1];
		if(readLength2[chr][i] > maxReadLength) {
		    maxReadLength = readLength2[chr][i];
		}
	    }
	    readCnt2[chr] = readCnt[0][chr] + readCnt[1][chr];
	}

	long genomelength = 0;
	for(int i=0; i<numChromosomes; i++) {
	    if(chrMax[i] >= chrMin[i])
		genomelength = genomelength + chrMax[i] - chrMin[i] + 1;
	}
	if(maxReadLength > 5000) {
	    System.err.println("WARNING: You have long reads, your longest is " + maxReadLength + " bases.\nThis program was designed for short reads. In principle it will work\nwith long reads, but it might take too long to run.\n\n");
	}

	Object[] RETURN = new Object[8];
	RETURN[0] = locationArray2;
	RETURN[1] = readCnt2;
	RETURN[2] = genomelength;
	RETURN[3] = chrMin;
	RETURN[4] = chrMax;
	RETURN[5] = maxReadsOneSpanFinal;	
	RETURN[6] = readLength2;
	RETURN[7] = maxReadLength;
	return RETURN;
    }

    private static Object[] readFileReportStrand(File file, int numChromosomes, int maxReadsOneSpan, int numReads, int extensionLength, int strandIndex, boolean removeIdenticalReads, int skip,  Hashtable<String, Integer> chrName2Number) {
        FileInputStream inFile = null;
	int[][][][] locationArray = new int[2][numChromosomes][maxReadsOneSpan][2];
	int[][] readCnt = new int[2][numChromosomes];
        try {
            inFile = new FileInputStream(file);
        } catch(FileNotFoundException e) {
	    System.err.println("The file " + file + " does not seem to exist...");
            System.exit(1);
        }
	for(int i=0; i<numChromosomes; i++) {
	    readCnt[0][i] = 0;
	    readCnt[1][i] = 0;
	}
	//int flag = 0;
	int chrNum = 0;
        try {
            inFile = new FileInputStream(file);
            BufferedReader br = new BufferedReader(new InputStreamReader(inFile));
            String line = "";
	    for(int i=0; i<skip; i++) {
		line = br.readLine();
	    }
            while((line = br.readLine()) != null) {
		String[] temp = line.split("\\s+", 0);
		String chr = temp[0];
		String strand = temp[strandIndex];
		chrNum = chrName2Number.get(chr);
		if(strand.equals("-")) {
		    locationArray[0][chrNum][readCnt[0][chrNum]][0] = Integer.parseInt(temp[1]) - extensionLength;
		    locationArray[0][chrNum][readCnt[0][chrNum]][1] = -1;
		    readCnt[0][chrNum]++;
		}
		else {
		    locationArray[1][chrNum][readCnt[1][chrNum]][0] = Integer.parseInt(temp[1]);
		    locationArray[1][chrNum][readCnt[1][chrNum]][1] = 1;
		    readCnt[1][chrNum]++;
		}
	    }
            inFile.close();
	} catch(IOException e) {
	    e.printStackTrace(System.err);
	    System.exit(1);
	}
	
	int[] maxReadsOneSpan2 = new int[2];
	for(int s=0; s<2; s++) {
	    maxReadsOneSpan2[s] = 0;
	    int l1 = 0;
	    int l2 = 0;
	    int j = 0;
	    for(int chr = 0; chr < numChromosomes; chr++) {
		int[][] temparray = new int[readCnt[s][chr]][2];
		for(int i=0; i<readCnt[s][chr]; i++) {
		    temparray[i][0] = locationArray[s][chr][i][0];
		    temparray[i][1] = locationArray[s][chr][i][1];
		}
		Arrays.sort(temparray, FirstIndexComparator);
		j=0;
		for(int i=0; i<readCnt[s][chr]-1; i++) {
		    l1 = temparray[i][0];
		    l2 = temparray[i+1][0];
		    if(l1 != l2 || !removeIdenticalReads) {
			locationArray[s][chr][j][0]=l1;
			locationArray[s][chr][j][1]=temparray[i][1];
			j++;
		    }
		}
		if(readCnt[s][chr] > 0) {
		    locationArray[s][chr][j][0]=temparray[readCnt[s][chr]-1][0];
		    locationArray[s][chr][j][1]=temparray[readCnt[s][chr]-1][1];
		}
		j++;
		readCnt[s][chr] = j;
		if(j > maxReadsOneSpan2[s])
		    maxReadsOneSpan2[s] = j;
	    }
	}

	int[][][] locationArray2 = new int[numChromosomes][maxReadsOneSpan2[0] + maxReadsOneSpan2[1]][2];
	for(int chr = 0; chr < numChromosomes; chr++) {
	    int[][] temparray = new int[readCnt[0][chr]+readCnt[1][chr]][2];
	    for(int i=0; i<readCnt[0][chr]; i++) {
		temparray[i][0] = locationArray[0][chr][i][0];
		temparray[i][1] = locationArray[0][chr][i][1];
	    }
	    for(int i=0; i<readCnt[1][chr]; i++) {
		temparray[i+readCnt[0][chr]][0] = locationArray[1][chr][i][0];
		temparray[i+readCnt[0][chr]][1] = locationArray[1][chr][i][1];
	    }
	    Arrays.sort(temparray, FirstIndexComparator);
	    for(int i=0; i<readCnt[0][chr]+readCnt[1][chr]; i++) {
		locationArray2[chr][i][0] = temparray[i][0];
		locationArray2[chr][i][1] = temparray[i][1];
	    }
	}

	Object[] RETURN = new Object[2];
	RETURN[0] = locationArray2;
	RETURN[1] = readCnt;
	return RETURN;
    }

    private static Object[] assessFile(File file, boolean rightEndInclusion, int strandIndex) {
        FileInputStream inFile = null;
        //int x = 0;
	//int chrNum = 0;
	int numChromosomes = 0;
	Object[] answer = new Object[7];
	int[] readCnt = new int[5002];
	int skip = 0;
	int maxReads = 0;
	int numReads = 0;
	String line = "";
	String firstline = "";
	Hashtable<String, Integer> chrName2Number = new Hashtable<String, Integer>();
	String[] chrNumber2Name = new String[1000];
	for(int i=0; i<1000; i++) {
	    readCnt[i] = 0;
	}
        try {
            inFile = new FileInputStream(file);
        } catch(FileNotFoundException e) {
	    System.err.println("The file " + file + " does not seem to exist...");
            System.exit(1);
        }
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(inFile));
	    firstline = br.readLine();
	    if(firstline == null) {
		System.err.println("This file is empty...");
		System.exit(1);
	    }
	    String[] temp = firstline.split("\\s+", 0);	
	    int n = Array.getLength(temp);
	    boolean testflag = true;
	    while(testflag == true) {
		testflag = false;
		try {
		    Integer.parseInt(temp[1]);
		    Integer.parseInt(temp[2]);
		} catch (Exception e3) {
		    testflag = true;
		    skip++;
		    firstline = br.readLine();
		    if(firstline == null) {
			System.err.println("\nError: Reached the end of the file without finding a proper line of data,\ncheck you have all fields in each line, there should be four ...\n");
			System.exit(1);
		    }
		    temp = firstline.split("\\s+", 0);
		    n = Array.getLength(temp);
		}
		if(n < 4) {
		    testflag = true;
		    skip++;
		    firstline = br.readLine();
		    if(firstline == null) {
			System.err.println("\nError: Reached the end of the file without finding a proper line of data,\ncheck you have all fields in each line, there should be four ...\n");
			System.exit(1);
		    }
		    temp = firstline.split("\\s+", 0);
		    n = Array.getLength(temp);
		}
	    }
            inFile.close();
	    if(strandIndex < 0) {
		int strand_cols_count = 0;
		for(int i=0; i<n; i++) {
		    if(temp[i].equals("+") || temp[i].equals("-")) {
			strandIndex = i;
			strand_cols_count++;
		    }
		}
		if(strand_cols_count > 1) {
		    System.out.print("Can't determine strand column, please enter it (starting counting from 0): ");
		    BufferedReader br2 = new BufferedReader(new InputStreamReader(System.in));
		    String userinput = null;
		    try {
			userinput = br2.readLine();
		    } catch (IOException ioe) {
			System.out.println("Sorry, there was an error in trying to input column number!");
			System.exit(1);
		    }
		    strandIndex = Integer.parseInt(userinput);
		    if(!(temp[strandIndex].equals("+") || temp[strandIndex].equals("-"))) {
			System.out.println("\nSorry, that does not look like a strand column, it must have only \"+\" and \"-\" signs.\nInstead I found the following: \"" + temp[strandIndex] + "\"");
			System.exit(0);
		    }
		}
	    }
	} catch(IOException e2) {
	    e2.printStackTrace(System.err);
	    System.exit(1);
	}
	//int linecounter = 0;
        try {
            inFile = new FileInputStream(file);
            BufferedReader br = new BufferedReader(new InputStreamReader(inFile));
	    for(int i=0; i<skip; i++) {
		line = br.readLine();
		//linecounter++;
	    }
	    //int flag = 0;
            while((line = br.readLine()) != null) {
		//linecounter++;
		//flag = 0;
		String[] temp = line.split("\\s+", 0);
		String chr = temp[0];
		if(!(chrName2Number.containsKey(chr))) {
		    chrName2Number.put(chr, numChromosomes);
		    chrNumber2Name[numChromosomes] = chr;
		    numChromosomes++;
		    if(numChromosomes > 5000) {
			System.err.println("\nERROR: Your file has more than 5000 chromosomes, that's to many...\n\n");
			System.exit(1);
		    }
		}
		numReads++;
		readCnt[chrName2Number.get(chr)]++;
		if(readCnt[chrName2Number.get(chr)] > maxReads) {
		    maxReads = readCnt[chrName2Number.get(chr)];
		}
	    }
	    inFile.close();
	} catch(IOException e) {
	    e.printStackTrace(System.err);
	    System.exit(1);
	}
	answer[0] = numChromosomes;
	answer[1] = maxReads;
	answer[2] = numReads;
	answer[3] = strandIndex;
	answer[4] = skip;
	answer[5] = chrName2Number;
	answer[6] = chrNumber2Name;

	return answer;
    }

    private static Object[] assessRepeatsFile(File file) {
	//boolean rightEndInclusion = true;
        FileInputStream inFile = null;
	int numChromosomes = 0;
        //int x = 0;
	//int chrNum = 0;
	int[] repeatCnt = new int[1000];
	Hashtable<String, Integer> chrName2Number_repeats = new Hashtable<String, Integer>();
	String[] chrNumber2Name_repeats = new String[1000];
	for(int i=0; i<1000; i++) {
	    repeatCnt[i] = 0;
	}
        try {
            inFile = new FileInputStream(file);
        } catch(FileNotFoundException e) {
	    System.err.println("The repeats file does not seem to exist...");
            System.exit(1);
        }
	int maxRepeats = 0;
	int numRepeats = 0;
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(inFile));
            String line = "";
	    //int flag = 0;
            while((line = br.readLine()) != null) {
		String[] temp = line.split("\\s+", 0);
		String chr = temp[0];
		numRepeats++;
		if(!(chrName2Number_repeats.containsKey(chr))) {
		    chrName2Number_repeats.put(chr, numChromosomes);
		    chrNumber2Name_repeats[numChromosomes] = chr;
		    numChromosomes++;
		}
		repeatCnt[chrName2Number_repeats.get(chr)]++;
		if(repeatCnt[chrName2Number_repeats.get(chr)] > maxRepeats) {
		    maxRepeats = repeatCnt[chrName2Number_repeats.get(chr)];
		}
	    }
	    inFile.close();
	} catch(IOException e2) {
            e2.printStackTrace(System.err);
            System.exit(1);
        }

	Object[] answer = new Object[6];
	answer[0] = (int)maxRepeats;
	answer[1] = numRepeats;
	answer[2] = repeatCnt;
	answer[3] = numChromosomes;
	answer[4] = chrName2Number_repeats;
	answer[5] = chrNumber2Name_repeats;
	return answer;
    }

    private static Object[] readMaskFile(File file, int numChromosomes, int maxRepeatsOneSpan, int numSpans, Hashtable<String, Integer> chrName2Number_mask) {
        FileInputStream inFile = null;
        //int x = 0;
	int[][][] locationArray = new int[numChromosomes][maxRepeatsOneSpan][2];
	int[] repeatCnt = new int[numChromosomes];
	int chrNum = 0;
	for(int i=0; i<numChromosomes; i++) {
	    repeatCnt[i] = 0;
	}
        try {
            inFile = new FileInputStream(file);
        } catch(FileNotFoundException e) {
	    System.err.println("The repeats file does not seem to exist...");
            System.exit(1);
        }
        try {
            inFile = new FileInputStream(file);
            BufferedReader br = new BufferedReader(new InputStreamReader(inFile));
            String line = "";
            while((line = br.readLine()) != null) {
		String[] temp = line.split("\\s+", 0);
		String chr = temp[0];
		chrNum = chrName2Number_mask.get(chr);
		locationArray[chrNum][repeatCnt[chrNum]][0] = Integer.parseInt(temp[1]);
		locationArray[chrNum][repeatCnt[chrNum]][1] = Integer.parseInt(temp[2]);
		repeatCnt[chrNum]++;
	    }
	    inFile.close();
	} catch(Exception e) {
	    e.printStackTrace(System.err);
	    System.exit(1);
	}
	Object[] RETURN = new Object[2];
	RETURN[0] = locationArray;
	RETURN[1] = repeatCnt;
	return RETURN;
    }

    public static Object[] normalize (int[][] locationArray, int[][] readLength, int[] readCnt, int numberToRemove) {
	// randomly removes numberToRemove things from the locationArray
	Random generator = new Random();
	int numChromosomes = Array.getLength(locationArray);
	int[] readCnt_new = new int[numChromosomes];
	int max=0;
	long genomelength = 0;
	long[] cumulative = new long[numChromosomes];
	cumulative[0] = readCnt[0];
	for(int i=0; i<numChromosomes; i++) {
	    readCnt_new[i] = 0;
	    if(readCnt[i] > max) {
		max = readCnt[i];
	    }
	    genomelength = genomelength + readCnt[i];
	    if(i > 0)
		cumulative[i] = cumulative[i-1] + readCnt[i];
	}
	int[][] locationArray_new = new int[numChromosomes][max];
	int[][] readLength_new = new int[numChromosomes][max];
	boolean[][] table = new boolean[numChromosomes][max];
	for(int c=0; c<numChromosomes; c++) {
	    for(int i=0; i<max; i++) {
		table[c][i] = true;
	    }
	}
	int cnt = 0;
	long p = 0;
	while(cnt < numberToRemove) {
	    p = nextLong(generator, genomelength);
	    int c = 0;
	    while(p > cumulative[c])
		c++;
	    int r = generator.nextInt(readCnt[c]);
	    if(table[c][r] == true) {
		table[c][r] = false;
		cnt++;
	    }
	}
	for(int i=0; i<numChromosomes; i++) {
	    for(int j=0; j<readCnt[i]; j++) {
		if(table[i][j] == true) {
		    locationArray_new[i][readCnt_new[i]] = locationArray[i][j];
		    readLength_new[i][readCnt_new[i]] = readLength[i][j];
		    readCnt_new[i]++;
		}
	    }
	}
	max = 0;
	for(int i=0; i<numChromosomes; i++) {
	    if(readCnt_new[i] > max) {
		max = readCnt_new[i];
	    }
	}
	int[][] locationArray_return = new int[numChromosomes][max];
	int[][] readLength_return = new int[numChromosomes][max];
	for(int i=0; i<numChromosomes; i++) {
	    for(int j=0; j<readCnt_new[i]; j++) {
		locationArray_return[i][j] = locationArray_new[i][j];
		readLength_return[i][j] = readLength_new[i][j];
	    }
	}

	Object[] RETURN = new Object[3];
	RETURN[0] = locationArray_return;
	RETURN[1] = readLength_return;
	RETURN[2] = readCnt_new;

	return RETURN;
    }

    public static long nextLong (Random generator, long m) {
	long d = m / 10;
	int w = (int)d;
	int r=0;
	int s=0;
	long answer;
	if(m == 0) {
	    answer = 0;
	}
	else {
	    if(w > 0) {
		r = generator.nextInt(w);
		s = generator.nextInt(10);
		answer = d * s + r;
	    }
	    else {
		answer = generator.nextInt((int)m);
	    }
	}
	return answer;
    }

    public static String format (long num) {
        String s = Long.toString(num);
        String answer = "";
        String[] temp = s.split("", 0);

        int len = Array.getLength(temp);
        int count = 0;
        for(int i=len - 1; i>0; i--) {
            answer = temp[i] + answer;
            count++;
            if(count == 3 && i>1) {
                answer = "," + answer;
                count = 0;
            }
        }
        return answer;
    }

    public static void writeDensityFile(int[][] location, int[][] readLength, int numSpans, int[] readCnt, BufferedWriter out_density, int startCoordinate, int startCoordinate_outfile, boolean rightEndInclusion_outfile, long[] chrMin, long[] chrMax, int windowSize, int windowDisplacement, int totalLength, int count_cutoff, String[] chrNumber2Name, boolean windows, int maxSpan, int chr_start_num) {
        // location is an array of start positions for the reads
	//int sizeOfArrays = 100000;
	//int max_sample = 0;
	int density_number_spans_per_window = windowSize / windowDisplacement;
	if(density_number_spans_per_window == 0) {
	    density_number_spans_per_window = 1;
	}
	int[] SPAN = new int[4];
	int first_flag = 1;
	for(int span = 0; span<numSpans; span++) {
	    int[][] localCounts = new int[density_number_spans_per_window][2];
	    for(int i = 0; i<density_number_spans_per_window; i++) {
		localCounts[i][0] = 0;
		localCounts[i][1] = 0;
	    }
	    System.err.println("Writing " + chrNumber2Name[span+chr_start_num] + ": " + readCnt[span] + " reads");
	    if(readCnt[span] > 0) {
		int[][] temparray = new int[readCnt[span]][2];
		int maxreadlength = 0;
		for(int i=0; i<readCnt[span]; i++) {
		    temparray[i][0] = location[span][i];
		    temparray[i][1] = readLength[span][i];
		    if(readLength[span][i] > maxreadlength) {
			maxreadlength = readLength[span][i];
		    }
		}
		System.err.println("The max read length on " + chrNumber2Name[span+chr_start_num] + " is " + maxreadlength);
		Arrays.sort(temparray, FirstIndexComparator);
		int spanlength = temparray[readCnt[span] - 1][0] + totalLength;
		int start = 0;
		int readcounter = 0;
		int num_over_loc = 0;
		int flag;
		for(int loc=startCoordinate; loc<spanlength+2; loc = loc + windowDisplacement) {
		    flag = 0;
		    readcounter=start;
		    num_over_loc = 0;
		    while(flag == 0 && (readcounter < readCnt[span])) {
			if(((loc + windowSize - 1) >= temparray[readcounter][0]) && (loc <= (temparray[readcounter][0] + temparray[readcounter][1] - 1))) {
			    num_over_loc++;
			}
			if((loc + windowSize - 1) < temparray[readcounter][0]) {
			    flag=1;
			}
			if(loc > (temparray[readcounter][0] + maxreadlength)) {
			    start++;
			    if(start == readCnt[span]) {
				flag = 1;
			    }
			}
			readcounter++;
		    }
		    try {
			int i = 0;
			while(i < density_number_spans_per_window && localCounts[i][0] < loc && localCounts[i][0] != 0) {
			    int ll1 = localCounts[i][0];
			    int ll2 = localCounts[i][0] + windowDisplacement;
			    if(startCoordinate_outfile == 1 && startCoordinate == 0) {
				ll1++;
				ll2++;
			    }
			    if(startCoordinate_outfile == 0 && startCoordinate == 1) {
				ll1--;
				ll2--;
			    }
			    if(windows) {
				if(localCounts[i][1] >= count_cutoff) {
				    if(rightEndInclusion_outfile) {
					int tmp = ll2 - 1;
					out_density.write(chrNumber2Name[span+chr_start_num]  + "\t" + ll1 + "\t" + tmp + "\t" + localCounts[i][1] + "\n");
				    }
				    else {
					out_density.write(chrNumber2Name[span+chr_start_num]  + "\t" + ll1 + "\t" + ll2 + "\t" + localCounts[i][1] + "\n");
				    }
				    out_density.flush();
				}
			    }
			    if(first_flag == 0 && SPAN[3] == localCounts[i][1]) { // span extends here
				SPAN[2] = ll2;
			    }
			    if(localCounts[i][1] >= count_cutoff && first_flag == 1) {  // span starts here
				first_flag = 0;
				SPAN[0] = span;
				SPAN[1] = ll1;
				SPAN[2] = ll2;
				SPAN[3] = localCounts[i][1];
			    }
			    if(((SPAN[3] != localCounts[i][1]) || ((0 < maxSpan) && (maxSpan < (SPAN[2] - SPAN[1])))) && first_flag == 0) { // span ends here
				first_flag = 1;
				if(rightEndInclusion_outfile) {
				    int tmp = SPAN[2] - 1;
				    if(!windows) {
					out_density.write(chrNumber2Name[SPAN[0]+chr_start_num]  + "\t" + SPAN[1] + "\t" + tmp + "\t" + SPAN[3] + "\n");
				    }
				}
				else {
				    if(!windows) {
					out_density.write(chrNumber2Name[SPAN[0]+chr_start_num]  + "\t" + SPAN[1] + "\t" + SPAN[2] + "\t" + SPAN[3] + "\n");
				    }
				}
				out_density.flush();
				if(localCounts[i][1] >= count_cutoff) { // new span starts
				    first_flag = 0;
				    SPAN[0] = span;
				    SPAN[1] = ll1;
				    SPAN[2] = ll2;
				    SPAN[3] = localCounts[i][1];
				}
			    }
			    i++;
			}
			for(int j=i; j<density_number_spans_per_window; j++) {
			    localCounts[j-i][0] = loc + (j-i) * windowDisplacement;
			    if(localCounts[j][1] <= num_over_loc)
				localCounts[j-i][1] = num_over_loc;
			    if(localCounts[j][1] > num_over_loc)
				localCounts[j-i][1] = localCounts[j][1];
			}
			for(int j=density_number_spans_per_window - i; j< density_number_spans_per_window; j++) {
			    localCounts[j][0] = loc + j * windowDisplacement;
			    localCounts[j][1] = num_over_loc;
			}
		    } catch (Exception e) {
			e.printStackTrace(System.err);
		    }
		}
		try {
		    int i = 0;
		    while(i < density_number_spans_per_window) {
			int ll1 = localCounts[i][0];
			int ll2 = localCounts[i][0] + windowDisplacement;
			if(startCoordinate_outfile == 1 && startCoordinate == 0) {
			    ll1++;
			    ll2++;
			}
			if(startCoordinate_outfile == 0 && startCoordinate == 1) {
			    ll1--;
			    ll2--;
			}
			if(localCounts[i][1] >= count_cutoff)
			    if(rightEndInclusion_outfile) {
				int tmp = ll2 - 1;
				out_density.write(chrNumber2Name[span+chr_start_num] + "\t" + ll1 + "\t" + tmp + "\t" + localCounts[i][1] + "\n");
			    }
			    else {
				out_density.write(chrNumber2Name[span+chr_start_num] + "\t" + ll1 + "\t" + ll2 + "\t" + localCounts[i][1] + "\n");
			    }
			i++;
		    }
		} catch (Exception e) {
		    e.printStackTrace(System.err);
		}
	    }
	}
	try {
	    out_density.flush();
	} catch (Exception e) {
	    e.printStackTrace(System.err);
	}
	System.err.println("");
    }


    public static void printUsage() {
	System.err.println("\n\n===================================================================================\n\nThis program takes short spans of genome and builds a 'depth-of-coverage' plot.  Input\nis a file of spans in genome coordinates and the output gives, for each genomic location,\nthe count of how many spans contain that location.  By default the output is given as\nspans where all bases in a span have the same count, but there are options to report\nit differntly.  There is an option (-ucsc) to output as a bedgraph file for immediate\nuploading to the UCSC browser.  Note: UCSC expects data to have start coordinate 0 and\nthe right endpoint of the span is *not* included (they call\nthis 'zero based half open').\n\nUSAGE:\n------\n\njava -jar M2C.jar <sample file> <output file> <log file> [options]\n\nPARAMETER DESCRIPTIONS:\n---------------------\n\n* <sample file> is the name of the file of reads.  See below for file format.\n\n* <output file> is the name of the file you want to save the output to.  <log file> is the name of the log file.\n\nThese are the only required arguments.\n\nOPTIONS:\n-------\n\n-name str: the name 'str' is used in the definition line of the output file, if there\nis a definition line.\n\n-countcutoff n: only spans that have count exceeding n (a non-negative integer) will\nbe output.  Default is n=1.\n\n-start_coordinate_infile n: n is 0 or 1 depending on whether the coordinate of the\nfirst base in input file is 0 or 1.  Default is n=1.\n\n-openinterval_infile: if set then the spans in the input file \"chr start end\" are\nassumed to *not* contain position \"end\".  Default is to not assume this, i.e. by\ndefault it assumes that \"end\" is included in the span.\n\n-start_coordinate_outfile n: n is 0 or 1 depending on whether the coordinate of the\nfirst base in output file is to be 0 or 1.  Default is n=1.\n\n-openinterval_outfile: if set then the spans in the output file \"chr start end\" are\nassumed to *not* contain position \"end\".  Default is to not assume this, i.e. by\ndefault it assumes that \"end\" is included in the span.\n\n-ucsc: if set, will output a format that is appropriate to load as a custom track\nin the UCSC genome browser.  this means:\n - there will be an appropriate header line\n - the coordinates in the outfile are the same as setting:\n     -start_coordinate_infile 0 -openinterval_outfile\n\n-windowsize n: compute counts in a sliding window of width n (set to 1 for resolution\nat the single base level).  Default is n=1.\n\n-windowdisplacement n: slide the window n bases to get consecutive windows.  Set this\nto 1 to output every window of length 'window length'.  Set this to be the same as\n-windowsize to get adjacent non-overlapping windows.  Default is n=1.\n\n-windows: outputs a value for each window, instead of combining windows into spans.  By\ndefault it will combine consecutive windows with the same density into a signle span.\n\n-strandIndex n: n gives the column with the strand (+ or -), starting counting from zero.\nThe strand column will be detected automatically unless it is ambiguous, in which case\nit will ask for the correct column, unless this parameter is set.\n\n-extensionlength: extension length is the total fragment length minus the read length.\nTypically used for ChIP-seq, not for RNA-seq.  Default is n = 0.\n\n-maxspan n: n is a positive integer.  By default, the program joins consecutive windows\nwith the same density into a single 'span'.   If -maxspan n is set, no span will be\nlonger than n nucleotides, so that consecutive spans might have the same density.\n\n-mask1 <filename>: An optional file of regions to mask with three tab delimited columns:\n    - chr, start, end\n    - coordinates start at 1 and right endpoint is included.\nReads are masked if their *unextended* coordinates overlap a region in the mask file.\nDownload properly formatted repeats files here: http://cbil.upenn.edu/STAR/repeats.\n\n-mask2 <filename>: An optional file of regions to mask with three tab delimited columns:\n    - chr, start, end\n    - coordinates start at 1 and right endpoint is included.\nReads are masked if their *extended* coordinates overlap a region in the mask file.\n\n-outputfilteredreads: If some filtering or masking of reads is done, set this option\nto write out the reads that were not filtered to a file called 'filtered_reads.txt'.\n\n-removeidenticalreads: if multiple reads map to the same location, keep only one.\n\n-norm n: reads will be thrown out randomly to produce a final set of n reads (if there\nare enough reads left after filtering to do so).\n\n-nodensity: do not output the density file.  Use this if you just want to output the\nfiltered reads.\n\n\nFILE FORMAT:\n------------\n\nThe sample file should have at least four tab delimited columns:\n\nchr_name   start_position   end_position   strand\n\nStrand can be any column after the end position column, it will figure it out which is\nthe strand column automatically and intervening columns will be ignored (if it cannot\ndetermine the strand column it will prompt you to input it, or you can put the column\nnumber (starting counting from 0) as an optional command line argument).  Data files\ncan have header lines, they will be ignored.\n\nFor more on file formats, run\n> M2C -example\n\nNOTES:\n------\n\n* NOTE 1: program will not work correctly if the genome length is greater than 21.4 Gb.\n\n* NOTE 2: you will probably need to increase the default memory.  After 'java' put the\noption -Xmx500m to raise it to 500 Mb, or change 500 to whatever is necessary.\nIf it hangs for a long time but doesn't crash it might run faster with more memory.\nWith 5 million reads we find it necessary to use at least 500mb.\n\nUSAGE REPEATED:\n--------------\n\njava -jar M2C.jar <sample file> <output file> <log file> [options]\n\n");
    }

    public static void printExample() {
	System.err.println("\n\n===================================================================================\n\nInput File.  This is the first 20 lines of a valid input file:\n\n-----------------------------------------------------------------\nchr17   23890145        23890251        +\nchr3    102856741       102856847       +\nchr10   74096904        74097010        -\nchr14   65198837        65198943        +\nchrM    6352            6458            +\nchr14   21136638        21136744        +\nchr7    118217611       118217717       +\nchr11   90104648        90104754        -\nchr13   49001839        49001945        -\nchr19   17474596        17474702        +\nchr1    72439224        72439330        +\nchr11   6217256         6217362         +\nchr18   62679475        62679581        +\nchr5    92487616        92487723        +\nchrX    106355619       106355726       -\nchr9    107576950       107577057       +\nchr14   70220320        70220427        +\nchr14   28860500        28860607        +\nchr1    88253626        88253733        +\nchr8    74044378        74044485        +\n...\n\nAnd if the -ucsc option is set, this is the first 20 lines of the output file:\n\n----------------------------------------------------------------\ntrack type=bedGraph name=\"example.cov\" description=\"window size = 1, window dispacement = 1, example.cov\" visibility=full color=255,0,0 priority=20\nchr1    3108511 3108619 1\nchr1    3289145 3289253 1\nchr1    3290596 3290631 1\nchr1    3292774 3292882 1\nchr1    3293203 3293276 1\nchr1    3296861 3296862 1\nchr1    3296862 3296881 2\nchr1    3296881 3296885 3\nchr1    3296885 3296902 4\nchr1    3296902 3296915 6\nchr1    3296915 3296934 7\nchr1    3296934 3296943 8\nchr1    3296943 3296946 10\nchr1    3296946 3296949 11\nchr1    3296949 3296950 12\nchr1    3296950 3296956 13\nchr1    3296956 3296969 14\nchr1    3296969 3296970 13\nchr1    3296970 3296989 12\n...\n\n* Note 1: If -ucsc is not set, then there will be no header line.\n\n* Note 2: If the -ucsc option is used, then the line:\n> chr1    3108511 3108619 1\nin the output file indicates positions 3108512 to 3108619 inclusive.\nThis is because the ucsc browser wants custom trakcs in 'zero-based half open' format.\n\n* Note 3: Unless the -extensionlength option is used, the strand makes no difference.\n\n");
    }
    
    public static Comparator<int[]> FirstIndexComparator = new Comparator<int[]>() {
	    @Override
	    public int compare(int[] array1, int[] array2) {
		return (array1[0] - array2[0]);
	    }
    };
}
