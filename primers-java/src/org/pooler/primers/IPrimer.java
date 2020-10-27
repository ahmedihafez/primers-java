/*
  This file is part of Java porting of Primer Pooler
  Original Primer Pooler (c) Silas S. Brown.  For Wen.
  Please refer to https://github.com/ssb22/PrimerPooler
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
package org.pooler.primers;

import java.io.PrintStream;

public interface IPrimer {
	IPrimer getReverse();
	boolean isDegeneratePrimer();
	void addTag(IPrimer tag);
	
	void addTagBackward(IPrimer tag);
	
	
	void removeTag(IPrimer tag);
	void removeTagBackward(IPrimer tag);
	
	int calcScore(IPrimer primer);
	float calcDeltaG(IPrimer backwordPrimer, float[] table);
	CountResult calcCount(IPrimer backwardPrimer);
	int NumPossibilities_32bases();
	IPrimer getComplement();
	
	
	
	void printBases(PrintStream outstream);
	
	void dGprint(IPrimer backwardPrimer, float minDG, PrintStream out, float[] table);
	void print(IPrimer backwardPrimer, int maxScore, PrintStream out);

}