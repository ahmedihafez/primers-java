package org.primer3.libprimer3;

import java.util.regex.Pattern;

public class IntervalList {

	
	//	org.apache.commons.lang3.tuple.MutablePair<L, R>
	// refactor to a list array
	private int[][] pairs = new int[LibPrimer3.PR_MAX_INTERVAL_ARRAY][2];

	// 
	private int count = 0;
	public boolean oligo_overlaps_interval(int start, int len) {

		int i;
		int last = start + len - 1;
		for (i = 0; i <  count; i++)
			if (!(last < pairs[i][0]
					|| start > (pairs[i][0] + pairs[i][1] - 1)))
				return true;
		return false;
	}

	/* Add a pair of integers to an  array of intervals */
	boolean addInterval( int i1, int i2)
	{
		int c = this.count;
		if (c >= LibPrimer3.PR_MAX_INTERVAL_ARRAY) return true;
		this.pairs[c][0] = i1;
		this.pairs[c][1] = i2;
		this.count++;
		return  false;
	}


	static final String numSep = ",";
	static final String intervalSep = " ";
	
	
	public IntervalList append_interval(String datum) {
		
//		IntervalArrayT2 tar2 = new IntervalArrayT2();
		IntervalList tar2 = this;


		String[] intervalStrs = datum.split(Pattern.quote(intervalSep));
		for(String intervalStr : intervalStrs){
			String[] intervalNums = intervalStr.split(Pattern.quote(numSep));
			int i1,i2;
			i1 = i2  = -1;
			if(!intervalNums[0].isEmpty())
				i1 = Integer.parseInt(intervalNums[0]);
			if(!intervalNums[1].isEmpty())
				i2 = Integer.parseInt(intervalNums[1]);

			if(tar2.addInterval(i1, i2))
				throw new IndexOutOfBoundsException("Too many elements for tag ");
		}

		return tar2; 
	}


	/**
	 * Check intervals, and add any errors to err.
	 * Update the start of each interval to
	 * be relative to the start of the included region.
	 */
	public boolean _check_and_adjust_1_interval(String tag_name, int seq_len,
			int first_index, StringBuilder err, SeqArgs sa,
			StringBuilder warning, boolean empty_allowed) {
		return _check_and_adjust_1_interval(tag_name,this.count,this.pairs,seq_len,first_index,err,sa,warning,empty_allowed);
	}

	/**
	 * Check intervals, and add any errors to err.
	 * Update the start of each interval to
	 * be relative to the start of the included region.
	 */
	public static boolean _check_and_adjust_1_interval(String tag_name,
			int count,
			int[][] pairs,
			int seq_len,
			int first_index, StringBuilder err, SeqArgs sa,
			StringBuilder warning, boolean empty_allowed) {
		int i;
		boolean outside_warning_issued = false;

		/* Subtract the first_index from the start positions in the interval
		     array */
		for (i = 0; i < count; i++) {
			if (empty_allowed && (pairs[i][0] == -1) && (pairs[i][1] == -1))
				continue;
			if (empty_allowed && (((pairs[i][0] == -1) && (pairs[i][1] != -1)) 
					|| ((pairs[i][0] != -1) && (pairs[i][1] == -1)))) {
				err.append(tag_name + " illegal interval");
				return true;
			}
			pairs[i][0] -= first_index;
		}

		for (i=0; i < count; i++) {
			if (empty_allowed && (pairs[i][0] == -1) && (pairs[i][1] == -1))
				continue;
			if (pairs[i][0] + pairs[i][1] > seq_len) {
				err.append(tag_name + " beyond end of sequence");
				return true;
			}
			/* Cause the interval start to be relative to the included region. */
			pairs[i][0] -= sa.incl_s;
			/* Check that intervals are within the included region. */
			if (pairs[i][0] < 0
					|| pairs[i][0] + pairs[i][1] > sa.incl_l) {
				if (!outside_warning_issued) {
					warning.append(tag_name +  " outside of INCLUDED_REGION");
					outside_warning_issued = true;
				}
			}
			if (pairs[i][1] < 0) {
				err.append( "Negative " + tag_name + " length");
				return true;
			}
		}
		return false;
	}

	
	
	public int getCount() {
		return count;
	}

	public int[] getPair(int i) {
		return pairs[i];
	}

	

}