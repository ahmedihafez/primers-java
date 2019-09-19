package org.pooler;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.time.Instant;

import org.primer3.masker.formula_parameters;

public class Fasta2Bit {

//	public int nSeq;
	public int byteSwap;
	public long seqPtr;
	ByteBuffer wrapped;
	public int nSeq;
	public void Fasta2Bit_(String genomeFile) {
		byte[] data;
		try {
			data = formula_parameters.readFile(genomeFile);
			int magic;
			wrapped = ByteBuffer.wrap(data);
			wrapped.order(ByteOrder.LITTLE_ENDIAN);	

			int sig = wrapped.getInt();
			if (sig !=  0x1A412743)
			{
				// TODO :: Invalid signature (is this really a .2bit genome file?)
			}
			int version =  wrapped.getInt();
			nSeq = wrapped.getInt();
			wrapped.getInt(); // reserved
			seqPtr = wrapped.position();
			// System.out.println(sig);
			// System.out.println(0x1A412743); // this LITTLE_ENDIAN
			// System.out.println(0x4327411A); // this BIG_ENDIAN



		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	BufferedInputStream bufferStream = null;
	public Fasta2Bit(String genomeFile) throws Exception {
		
		Path inputFilePath = Paths.get(genomeFile);
		try {
				InputStream input = Files.newInputStream(  inputFilePath, StandardOpenOption.READ);
				bufferStream = new BufferedInputStream(input);
				readHeader();
				parseSeqDir();
		} catch(Exception e) {

			e.printStackTrace(System.err);
			throw e;
		}
	}

	ByteBuffer readData(int startPostion , int size ) throws IOException {
		byte[] data = new byte[size];
		
		// handle startPostion here
		if(startPostion != 0 ) {
			// this need Skipping
		}
		
		bufferStream.read(data, 0, size);
		ByteBuffer wrapped = ByteBuffer.wrap(data);
		wrapped.order(ByteOrder.LITTLE_ENDIAN);
		return wrapped;
	}
	
	class SeqInfo{
		String seqname  = null;
		int offset = -1 ;
		int byteSize = -1 ;
		int index = -1;
	}
	
	SeqInfo[] seqDir;
	int currentSeqIndex = 0;
	
	class SeqData {
		SeqInfo seq;
		int currentSeqIndex;
		ByteBuffer data;
		int nBases;
		int nUnknown;
		int nBasesInByte;
		int[] unknownStart ;
		int[] unknownLen ;
		
		
		
		
		
		boolean printProgress =  true;
		
		
		
		
		int currentBaseNo= 0;
		int currentByteIndex= 0;
		byte basesLeftInByte = 0;
		byte currentByte;
		int currentUnknownRegionIndex = 0;
		long currentBuf= 0, currentValid = 0;
		// move next to new base
		
		
		Instant nextThreadUpdate = Instant.ofEpochSecond(0) ;
		
		
		public boolean next() {
			if(currentBaseNo >= nBases) return false;
			
			while(true) {
				// advance to the next valid byte and return
				if(currentBaseNo >= nBases) // if everything is correct it is not going to end here
					// as long we have advanced to a new byte we should complete 4 cycles inside here
					return true; // or false ??
				if(basesLeftInByte == 0 ) {
					currentByte = data.get();
					basesLeftInByte = 4;
					if(printProgress)
					if (Instant.now().getEpochSecond() > nextThreadUpdate.getEpochSecond() + 2)
					{
						System.err.format("\rScanning %s %2d%%",seq.seqname,(int)((float)currentBaseNo*100.0/(float)nBases));
						nextThreadUpdate = Instant.now();
					}
				}
				if(currentBaseNo < unknownStart[currentUnknownRegionIndex]) {
					  /* This function has to be FAST: seriously inner-loop.
				     For ease of binary search in amplicons.c, bases are
				     shifted into buf from LEFT (not right as in 64.h),
				     so buf contains the last few bases IN REVERSE from
				     the genome cursor (which is an 'end' cursor).  */
				  currentBuf = ((currentBuf) >>> 2) | ((((long)currentByte >>> (2*--basesLeftInByte)) & (long)3) << 62);
				  if(currentValid != ~(long)0) /* (only at start, so put the 'if' around it to save a couple of instructions) */
					  currentValid = ((currentValid) >>> 2) | ((long)3 << 62);
				 currentBaseNo++;
				 return true;
				  
				}
				else {
					--basesLeftInByte; // ignore an unknown region
					if(++currentBaseNo  == (unknownStart[currentUnknownRegionIndex] + unknownLen[currentUnknownRegionIndex] ) )
					{
						currentUnknownRegionIndex++;
						currentBuf = currentValid = 0;
					}
				}
			}
		} 
	}
	
	
	void parseSeqDir() {
//		System.out.println(seqPtr);
		// ByteBuffer wrapped = readData(0,1)	;
		try {
			seqDir = new  SeqInfo[nSeq];
			for (int i = 0 ; i < nSeq; i ++ )
			{
				int nextSeqNameLen = bufferStream.read();
				// System.out.println(nextSeqNameLen);
				byte[] seqNamechars = new byte[nextSeqNameLen];
				bufferStream.read(seqNamechars,0,nextSeqNameLen);
//				System.out.println(new String(seqNamechars));
				byte[] offsetData = new byte[4];
				bufferStream.read(offsetData,0,4);
				ByteBuffer wrapped = ByteBuffer.wrap(offsetData);
				wrapped.order(ByteOrder.LITTLE_ENDIAN);
				int offset = wrapped.getInt();
//				System.out.println(offset);

				seqDir[i] = new SeqInfo(); 
				seqDir[i].offset = offset;
				seqDir[i].index = i; // original index
				seqDir[i].seqname = new String(seqNamechars);
				
				// size should be like this not good ?? this is how much should we read in the byte array
				if(i > 0 ) {
					seqDir[i-1].byteSize = seqDir[i].offset - seqDir[i-1].offset;
//					System.out.println(allSeqs[i-1].byteSize);
				}
			}
			seqDir[seqDir.length-1].byteSize = seqDir[seqDir.length-1].offset - seqDir[seqDir.length-2].offset;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	SeqData parseNextSeq() throws IOException {
		SeqInfo curSeq = seqDir[currentSeqIndex];
		
		ByteBuffer wrapped = readData(0,curSeq.byteSize);
		
		int nBases = wrapped.getInt();
		int nUnknown = wrapped.getInt(); // len of an array holding Nnnnn blocks
		
		int[] unknownStart = new int[nUnknown +1]; // no need for + 1
		int[] unknownLen = new int[nUnknown  +1]; 
		
		for(int i=0; i<nUnknown; i++) {
			unknownStart[i] = wrapped.getInt();
		}
		unknownStart[nUnknown] = nBases;
		for(int i=0; i<nUnknown; i++) {
			unknownLen[i] = wrapped.getInt();
		}
		unknownLen[nUnknown] = 0;
		int maskedRegionLen = wrapped.getInt();
		
		wrapped.position(wrapped.position() + (maskedRegionLen*2*4)+ 4  ); // skip 'masked blocks' (TODO: are these ever relevant?) and the 4-byte 'reserved' word
		
		
		// after that is the seq bases
		// System.out.println(curSeq.seqname);
		// System.out.println(nBases);

		SeqData seqData = new SeqData();
		seqData.seq = curSeq;
		seqData.data = wrapped;
		seqData.currentSeqIndex = currentSeqIndex;
		seqData.unknownLen = unknownLen;
		seqData.unknownStart = unknownStart;
		seqData.nUnknown= nUnknown; // no need for this
		seqData.nBases = nBases;
		seqData.nBasesInByte = seqData.nBases/4;
		
		
		currentSeqIndex++;
		
		return seqData;
		
	}
	
	boolean readHeader() {
		// Header 
		// 4-byte magic-number -- Signature
		// 4-byte version
		// 4-byte -- number of seq in the file
		// 4-byte -- reserved
		// total 16-byte for now to read
		
		try {
			int headSize = 16;
			ByteBuffer wrapped = readData(0,headSize)	;		
			int sig = wrapped.getInt();
			
			if (sig !=  0x1A412743)
			{
				// TODO :: Invalid signature (is this really a .2bit genome file?)
				return false;
			}
			int version =  wrapped.getInt();
			nSeq = wrapped.getInt();
			wrapped.getInt(); // reserved
			seqPtr = headSize;
		} catch (IOException e) {
			// something is wrong 
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	
	public static void main(String[] argv ) {
		
//		Fasta2Bit testFile;
//		try {
//			Genome.go_throgh_genome(  "/data/softwares/pooler/example/hg38.2bit");
//		} catch (Exception e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
	}
	
	
}
