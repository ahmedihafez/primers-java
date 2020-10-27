/*
    This file is part of primer3 porting to java


	Original file are part of https://github.com/primer3-org/primer3
	Whitehead Institute for Biomedical Research, Steve Rozen
	(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
	All rights reserved to Primer3 authors.

    Primer3 and the libprimer3 library are free software;
    you can redistribute them and/or modify them under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software (file gpl-2.0.txt in the source
    distribution); if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
package org.primer3.libprimer3;

/** 
 * Arguments to the primer program as a whole.  Values for these arguments are
 * retained _across_ different input records.  (These are the so-called
 * "Global" arguments in the documentation.)
 */
public class OligoWeights implements Cloneable {

	private double complAny;
	private double complAnyTh;
	private double complEnd;
	private double complEndTh;
	private double endQuality;
	private double endStability;
	private double gcContentQT;
	private double gcContentLT;
	private double hairpinTh;
	private double lengthQT;
	private double lengthLT;
	private double numNs;
	private double posPenalty;
	private double repeatSimilarity;
	private double seqQuality;
	private double temperatureCutoff;
	private double temperatureQT;
	private double temperatureLT;
	private double templateMispriming;
	private double templateMisprimingTh;
	private double failureRate;

	
	
	
	
	
	
	
	
	public double getComplAny() {
		return complAny;
	}

	public void setComplAny(double complAny) {
		this.complAny = complAny;
	}

	/**
	 * @return the complAnyTh
	 */
	public double getComplAnyTh() {
		return complAnyTh;
	}

	/**
	 * @param complAnyTh the complAnyTh to set
	 */
	public void setComplAnyTh(double complAnyTh) {
		this.complAnyTh = complAnyTh;
	}

	/**
	 * @return the complEnd
	 */
	public double getComplEnd() {
		return complEnd;
	}

	/**
	 * @param complEnd the complEnd to set
	 */
	public void setComplEnd(double complEnd) {
		this.complEnd = complEnd;
	}

	/**
	 * @return the complEndTh
	 */
	public double getComplEndTh() {
		return complEndTh;
	}

	/**
	 * @param complEndTh the complEndTh to set
	 */
	public void setComplEndTh(double complEndTh) {
		this.complEndTh = complEndTh;
	}

	/**
	 * @return the endQuality
	 */
	public double getEndQuality() {
		return endQuality;
	}

	/**
	 * @param endQuality the endQuality to set
	 */
	public void setEndQuality(double endQuality) {
		this.endQuality = endQuality;
	}

	/**
	 * @return the endStability
	 */
	public double getEndStability() {
		return endStability;
	}

	/**
	 * @param endStability the endStability to set
	 */
	public void setEndStability(double endStability) {
		this.endStability = endStability;
	}

	/**
	 * @return the gcContentQT
	 */
	public double getGcContentQT() {
		return gcContentQT;
	}

	/**
	 * @param gcContentQT the gcContentQT to set
	 */
	public void setGcContentQT(double gcContentQT) {
		this.gcContentQT = gcContentQT;
	}

	/**
	 * @return the gcContentLT
	 */
	public double getGcContentLT() {
		return gcContentLT;
	}

	/**
	 * @param gcContentLT the gcContentLT to set
	 */
	public void setGcContentLT(double gcContentLT) {
		this.gcContentLT = gcContentLT;
	}

	/**
	 * @return the hairpinTh
	 */
	public double getHairpinTh() {
		return hairpinTh;
	}

	/**
	 * @param hairpinTh the hairpinTh to set
	 */
	public void setHairpinTh(double hairpinTh) {
		this.hairpinTh = hairpinTh;
	}

	/**
	 * @return the lengthQT
	 */
	public double getLengthQT() {
		return lengthQT;
	}

	/**
	 * @param lengthQT the lengthQT to set
	 */
	public void setLengthQT(double lengthQT) {
		this.lengthQT = lengthQT;
	}

	/**
	 * @return the lengthLT
	 */
	public double getLengthLT() {
		return lengthLT;
	}

	/**
	 * @param lengthLT the lengthLT to set
	 */
	public void setLengthLT(double lengthLT) {
		this.lengthLT = lengthLT;
	}

	/**
	 * @return the numNs
	 */
	public double getNumNs() {
		return numNs;
	}

	/**
	 * @param numNs the numNs to set
	 */
	public void setNumNs(double numNs) {
		this.numNs = numNs;
	}

	/**
	 * @return the posPenalty
	 */
	public double getPosPenalty() {
		return posPenalty;
	}

	/**
	 * @param posPenalty the posPenalty to set
	 */
	public void setPosPenalty(double posPenalty) {
		this.posPenalty = posPenalty;
	}

	/**
	 * @return the repeatSimilarity
	 */
	public double getRepeatSimilarity() {
		return repeatSimilarity;
	}

	/**
	 * @param repeatSimilarity the repeatSimilarity to set
	 */
	public void setRepeatSimilarity(double repeatSimilarity) {
		this.repeatSimilarity = repeatSimilarity;
	}

	/**
	 * @return the seqQuality
	 */
	public double getSeqQuality() {
		return seqQuality;
	}

	/**
	 * @param seqQuality the seqQuality to set
	 */
	public void setSeqQuality(double seqQuality) {
		this.seqQuality = seqQuality;
	}

	/**
	 * @return the temperatureCutoff
	 */
	public double getTemperatureCutoff() {
		return temperatureCutoff;
	}

	/**
	 * @param temperatureCutoff the temperatureCutoff to set
	 */
	public void setTemperatureCutoff(double temperatureCutoff) {
		this.temperatureCutoff = temperatureCutoff;
	}

	/**
	 * @return the temperatureQT
	 */
	public double getTemperatureQT() {
		return temperatureQT;
	}

	/**
	 * @param temperatureQT the temperatureQT to set
	 */
	public void setTemperatureQT(double temperatureQT) {
		this.temperatureQT = temperatureQT;
	}

	/**
	 * @return the temperatureLT
	 */
	public double getTemperatureLT() {
		return temperatureLT;
	}

	/**
	 * @param temperatureLT the temperatureLT to set
	 */
	public void setTemperatureLT(double temperatureLT) {
		this.temperatureLT = temperatureLT;
	}

	/**
	 * @return the templateMispriming
	 */
	public double getTemplateMispriming() {
		return templateMispriming;
	}

	/**
	 * @param templateMispriming the templateMispriming to set
	 */
	public void setTemplateMispriming(double templateMispriming) {
		this.templateMispriming = templateMispriming;
	}

	/**
	 * @return the templateMisprimingTh
	 */
	public double getTemplateMisprimingTh() {
		return templateMisprimingTh;
	}

	/**
	 * @param templateMisprimingTh the templateMisprimingTh to set
	 */
	public void setTemplateMisprimingTh(double templateMisprimingTh) {
		this.templateMisprimingTh = templateMisprimingTh;
	}

	/**
	 * @return the failureRate
	 */
	public double getFailureRate() {
		return failureRate;
	}

	/**
	 * @param failureRate the failureRate to set
	 */
	public void setFailureRate(double failureRate) {
		this.failureRate = failureRate;
	}
	
	
	
	
	
	
	
	// ###########################################################################################
	// Setter method from a string value (source configuration file)
	public void set_template_mispriming_th(String datum) {
		this.setTemplateMisprimingTh(Double.parseDouble(datum));
	}

	public void set_end_quality(String datum) {
		this.setEndQuality(Double.parseDouble(datum));
	}

	public void set_seq_quality(String datum) {
		this.setSeqQuality(Double.parseDouble(datum));
	}

	public void set_repeat_sim(String datum) {
		this.setRepeatSimilarity(Double.parseDouble(datum));
	}

	public void set_num_ns(String datum) {
		this.setNumNs(Double.parseDouble(datum));
	}

	public void set_hairpin_th(String datum) {
		this.setHairpinTh(Double.parseDouble(datum));	
	}

	public void set_compl_end_th(String datum) {
		this.setComplEndTh(Double.parseDouble(datum));
	}

	public void set_compl_any_th(String datum) {
		this.setComplAnyTh(Double.parseDouble(datum));
	}

	public void set_compl_end(String datum) {
		this.setComplEnd(Double.parseDouble(datum));
	}

	public void set_compl_any(String datum) {
		this.setComplAny(Double.parseDouble(datum));
	}

	public void set_length_gt(String datum) {
		this.setLengthQT(Double.parseDouble(datum));
	}

	public void set_length_lt(String datum) {
		this.setLengthLT(Double.parseDouble(datum));
	}

	public void set_gc_content_lt(String datum) {
		this.setGcContentLT(Double.parseDouble(datum));
	}

	public void set_gc_content_gt(String datum) {
		this.setGcContentQT(Double.parseDouble(datum));
	}

	public void set_temp_lt(String datum) {
		this.setTemperatureLT(Double.parseDouble(datum));
	}

	public void set_temp_gt(String datum) {
		this.setTemperatureQT(Double.parseDouble(datum));
	}
	public void set_failure_rate(String datum) {
		this.setFailureRate(Double.parseDouble(datum));
	}
	public void set_template_mispriming(String datum) {
		this.setTemplateMispriming(Double.parseDouble(datum));
	}
	public void set_end_stability(String datum) {
		this.setEndStability(Double.parseDouble(datum));
	}
	public void set_pos_penalty(String datum) {
		this.setPosPenalty(Double.parseDouble(datum));
	}

	// ###########################################################################

	@Override
	protected OligoWeights clone() throws CloneNotSupportedException {
		OligoWeights newClone = new OligoWeights();
		
		// field wise copy :: Make sure that any change in the class are reflected here
		// ####################
		newClone.setComplAny(this.getComplAny());
		newClone.setComplAnyTh(this.getComplAnyTh());
		newClone.setComplEnd(this.getComplEnd());
		newClone.setComplEndTh(this.getComplEndTh());
		newClone.setEndQuality(this.getEndQuality());
		newClone.setEndStability(this.getEndStability());
		newClone.setGcContentQT(this.getGcContentQT());
		newClone.setGcContentLT(this.getGcContentLT());
		newClone.setHairpinTh(this.getHairpinTh());
		newClone.setLengthQT(this.getLengthQT());
		newClone.setLengthLT(this.getLengthLT());
		newClone.setNumNs(this.getNumNs());
		newClone.setPosPenalty(this.getPosPenalty());
		newClone.setRepeatSimilarity(this.getRepeatSimilarity());
		newClone.setSeqQuality(this.getSeqQuality());
		newClone.setTemperatureCutoff(this.getTemperatureCutoff());
		newClone.setTemperatureQT(this.getTemperatureQT());
		newClone.setTemperatureLT(this.getTemperatureLT());
		newClone.setTemplateMispriming(this.getTemplateMispriming());
		newClone.setTemplateMisprimingTh(this.getTemplateMisprimingTh());
		newClone.setFailureRate(failureRate);
		// ##########################################
		
 		return newClone;
	}



}