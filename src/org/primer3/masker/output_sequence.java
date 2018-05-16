package org.primer3.masker;

/*
 * output_sequence: A wrapper for output sequence in case the output is given in a 
 * string variable
 */
public class output_sequence {
	
	
	// all can be converted into a StringBuffer or builder
	public char [] sequence;
	int pos;
	
	/* These will be used instead of char *sequence
	 * in case of both_separately masking direction */
	public char []sequence_fwd;
	public char []sequence_rev;
	
	
	

	public void write_char_to_output(char c, char c_other, MaskerParameters mp) {
		
		if (mp.print_sequence) {
			System.err.print(c);
		} 
		else {
			if (mp.mdir == masking_direction.both_separately) {
				this.sequence_fwd[this.pos] = c;
				this.sequence_rev[this.pos] = c_other;
			} else {
				this.sequence[this.pos] = c;
			}
			this.pos += 1;
		}
	}
	void  write_header_to_output (String header_name,MaskerParameters mp)
	{
		if (mp.print_sequence) 
			System.err.println(header_name);
		else {
			if (mp.mdir == masking_direction.both_separately) {
//				v = memcpy (output_seq.sequence_fwd + output_seq.pos, header_name, strlen(header_name));
//				if (v) v = memcpy (output_seq.sequence_rev + output_seq.pos, header_name, strlen(header_name));
			} else {
//				v = memcpy (output_seq.sequence + output_seq.pos, header_name, strlen(header_name));
			}	
//			this.pos += header_name.length();
		}
	}
	
	public static output_sequence create_output_sequence (int seq_len, masking_direction mdir)
	{
		output_sequence output_seq = new output_sequence();
		if (mdir == masking_direction.both_separately) {
			output_seq.sequence_fwd = new char [seq_len];
			output_seq.sequence_rev = new char [seq_len];
		} else {
			output_seq.sequence = new char [seq_len];
		}
		output_seq.pos = 0;
		return output_seq;
	}
	
}