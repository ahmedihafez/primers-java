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
package org.ucam.ssb22.pooler.primers;

import org.ucam.ssb22.pooler.Amplicons.PrimerToFind;

public abstract class DegeneratePrimer implements IPrimer {

	@Override
	public final boolean isDegeneratePrimer() {
		return true;
	}

	public abstract boolean make2bit(PrimerToFind newPrimerPoss, int possNo, int nPoss);
}
