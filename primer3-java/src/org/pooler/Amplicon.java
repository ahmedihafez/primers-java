package org.pooler;

import java.util.ArrayList;
import java.util.List;

public class Amplicon {
	public int fwdIndex = -1;
	public int revIndex = -1;
	public String name;
	public int poolIndex;
	Location loc = null;
	
	ArrayList<Amplicon>  overlappingAmps = null;
	
	public class Location {
		long start;
		long end;
		String seqName;
		boolean strand;
	}
	
	public void setLocation(String seqName,long start,long end,boolean strand) {
		loc = new Location();
		loc.seqName = seqName;
		loc.start = start;
		loc.end = end;
		loc.strand = strand;
	}

	public void addOverlapingAmp(Amplicon overlappingAmp) {
		if(overlappingAmps == null) {
			overlappingAmps = new ArrayList<Amplicon>();
		}
		overlappingAmps.add(overlappingAmp);
	}

	public String getName() {
		return name;
	}

	public String getSequenceName() {
		if(loc!= null)
			return loc.seqName;
		return null;
	}
	
	public long getStart() {
		if(loc!= null)
			return loc.start;
		return -1;
	}
	
	public long getEnd() {
		if(loc!= null)
			return loc.end;
		return -1;
	}
	
	public long getLength() {
		if(loc!= null)
			return loc.end + 1 - loc.start;
		return -1;
	}
	
	public String getStrand() {
		if(loc!= null)
			return loc.strand ? "+":"-";
		return "";
	}

	public List<Amplicon> getOverlaps() {
		return this.overlappingAmps;
	}
	
}