/*
  This file is part of Java porting of Primer Pooler (https://github.com/ssb22/PrimerPooler)
  Primer Pooler (c) Silas S. Brown.  For Wen.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package org.ucam.ssb22.pooler;

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