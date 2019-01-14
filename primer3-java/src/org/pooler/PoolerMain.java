package org.pooler;

import org.pooler.AllPrimers.PS_cache;

public class PoolerMain {
	//	  printf("%s [options] FASTA-file\n",arg0());
	//      puts("Options can include any of:");
	//      puts("--counts  (show score or deltaG-range pair counts for whole input)");
	boolean doCounts;

	//      puts("--self-omit (omits self-interaction and pair-interaction in --counts)");
	boolean selfOmit;

	//      puts("--print-bonds=THRESHOLD e.g. --print-bonds=1");
	double threshold;
	boolean hasThreshold;
	//      puts("--dg[=KELVIN[,mg[,cation[,dNTP]]]] to use dG instead of score");

	boolean calcDG =  false;
	float temp = Temprature.C_to_kelvin(37);
	float mg = 0;
	float cation = 50;
	float dNTP  = 0;
	// TODO 
	// table = deltaG_table(temp,mg,cation,dNTP);

	//      puts("--suggest-pools suggests a number for --pools (or use --pools=?)");
	boolean suggestPools = false;

	//      puts("--pools[=NUM[,MINS[,PREFIX]]] e.g. --pools=2,1,poolfile-");
	int numPools = 0;
	int limitMinsToRun ;
	String poolPrefix;
	//      puts("(Set prefix to a single hyphen (-) to write all to stdout)");
	//      puts("--max-count=NUM (per pool)");
	int maxCount;
	//      puts("--genome=PATH to check amplicons for overlaps in the genome (.2bit or .fa)");
	String genomeFile;
	//      puts("--amp-max=LENGTH sets max amplicon length for the overlap check (default 220)"); /* v1.35 added 0 = unlimited, not available in interactive version */
	int maxAmpliconLen;
	//      puts("--multiplx=FILE (write MultiPLX input after the --genome stage)");
	String multiplexFile;
	//      puts("--seedless (don't seed random number generator, and use the same\n           
	// one across all operating systems for reproducibility)");
	boolean seedless = false;




	public void runPooler(String fileName) {
		AllPrimers ap = new AllPrimers();
		try {

			float [] table = null;
			if(this.calcDG) {
				table = Temprature.deltaG_table(temp,mg,cation, dNTP);
			}

			ap.loadFasta(fileName);

			// FASTA file loaded OK
			if(this.doCounts) {
				if(this.selfOmit) {
					PS_cache cache = ap.PS_precalc(table,null,null,0);
					//		          if(cache.scores) {
					int[] pools = new int[ap.np]; /* all in pool 0 */
					//		            if(!memFail(pools,_memFail)) {
					if(table != null) ap.dGprintStats(pools,cache.scores);
					else ap.printStats(pools,cache.scores);
					//		              free(pools);
					//		            }
					//		            PS_free(cache);
					//		          } /* else memFail will already have printed something */
				} else {
					ap.addTags();
					if(table != null) {
						ap.dGandScoreCounts(ap,table);
					}
					else 
					{
						int[] counts = ap.counts();
						for(int i = 0 ; i < counts.length ;i++)
						{
							System.out.println(i + "\t"  + counts[i]);
						}
					}
					ap.removeTags();
				}
			}

			//		      if(numPools>0 && maxCount) {
			//		          /* preliminary: check range of maxCount */
			//		          int average=2*averagePairsRoundUp(ap.np,numPools);
			//		          if(maxCount<average) {
			//		            fputs("Can't do that: --max-count is too low!\n",stderr); suggestMax(numPools,average,ap.np);
			//		            exit(1); /* no need to free memory */
			//		          }
			//		        }
			int nAmplicons=0; boolean[] overlappingAmplicons=null; int[] primerNoToAmpliconNo=null;
			if(genomeFile != null) {
				Amplicons amplicons= new Amplicons(ap,genomeFile,maxAmpliconLen,multiplexFile,1);
				overlappingAmplicons = amplicons.getOverlappingAmplicons(); /* TODO: support ,0 here for just an all-amplicon-locations-file from command line? */
				primerNoToAmpliconNo = amplicons.primerNoToAmpliconNo; 
				nAmplicons = amplicons.nAmplicons;
			}
			if(numPools != 0 || suggestPools ) {
				PS_cache cache= ap.PS_precalc(table,overlappingAmplicons,primerNoToAmpliconNo,nAmplicons);
				if(suggestPools || numPools < 0) {
					int suggestion = PoolSplit.suggest_num_pools(ap,cache,table);
					System.err.format("Computer suggestion is %d pools.\n",suggestion);
					if(numPools < 0) {
						numPools = suggestion;
						if (maxCount != 0) {
							/* as the 'preliminary' above, but we now need to check it here (TODO: duplicate code) */
							//		                int average=2 * averagePairsRoundUp(ap.np,numPools);
							//		                if(maxCount<average) {
							//		                	 System.err.format("Can't do that: --max-count is too low!\n"); 
							//		                	 suggestMax(numPools,average,ap.np);
							//		                     throw new Exception("\"Can't do that: --max-count is too low!");
							//		                }
						}
					}
				}

				int[] pools = null;
				if (numPools != 0) 
					pools = PoolSplit.split_into_pools(ap,numPools,this.limitMinsToRun,cache,seedless,table,maxCount);
				if(pools != null) {
					if(hasThreshold) {
						ap.addTags();
						if(table != null) ap.dGprintBonds(System.out,threshold,pools,table);
						else ap.printBonds(System.out,threshold,pools);
						ap.removeTags();
					}
				}

			}
			else
			{
				if(hasThreshold) {
					ap.addTags();
					if(table != null) ap.dGprintBonds(System.out,threshold,null,table);
					else ap.printBonds(System.out,threshold,null);
					ap.removeTags();
				}
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	static public void main(String[] argv) {
		PoolerMain poolerMain = new PoolerMain();

		poolerMain.doCounts = true;
		poolerMain.calcDG = true;
		poolerMain.selfOmit = false;
		poolerMain.limitMinsToRun = 1;
//		poolerMain.numPools = -1;
		poolerMain.hasThreshold = true;
		poolerMain.threshold = -5;
//		poolerMain.genomeFile = "/data/softwares/pooler/example/hg38.2bit";
		poolerMain.runPooler("/data/softwares/pooler/example/example4.txt");


	}

}
