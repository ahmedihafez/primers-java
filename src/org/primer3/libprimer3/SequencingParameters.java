package org.primer3.libprimer3;

public class SequencingParameters {
	
	private int lead;
	private int spacing;
	private int interval;
	private int accuracy;
	
	public void setLead(String datum) {
		this.setLead(Integer.parseInt(datum));
	}
	
	/**
	 * @return the lead
	 */
	public int getLead() {
		return lead;
	}
	
	/**
	 * @param lead the lead to set
	 */
	public void setLead(int lead) {
		this.lead = lead;
	}
	
	public void setSpacing(String datum) {
		this.setSpacing(Integer.parseInt(datum));		
	}
	
	/**
	 * @return the spacing
	 */
	public int getSpacing() {
		
		return spacing;
	}
	
	/**
	 * @param spacing the spacing to set
	 */
	public void setSpacing(int spacing) {
		this.spacing = spacing;
	}
	
	public void setInterval(String datum) {
		this.setInterval(Integer.parseInt(datum));			
	}
	
	/**
	 * @return the interval
	 */
	public int getInterval() {
		return interval;
	}
	
	/**
	 * @param interval the interval to set
	 */
	public void setInterval(int interval) {
		this.interval = interval;
	}
	
	public void setAccuracy(String datum) {
		this.setAccuracy(Integer.parseInt(datum));				
	}
	
	/**
	 * @return the accuracy
	 */
	public int getAccuracy() {
		return accuracy;
	}
	
	/**
	 * @param accuracy the accuracy to set
	 */
	public void setAccuracy(int accuracy) {
		this.accuracy = accuracy;
	}
}