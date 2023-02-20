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

public class dHS {
	  float dH_kcal_per_mol;
	  float dS_cal_per_molK;
	  public dHS(float dH , float dS) {
		  this.dH_kcal_per_mol = dH;
		  this.dS_cal_per_molK = dS;
	  }
}