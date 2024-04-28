import FIFO::*;
import FIFOF::*;
import Vector::*;

import BRAM::*;
import BRAMFIFO::*;


// Sequences
typedef 32768 SeqNum;
// Motif
typedef 16 MotifLength;
typedef TMul#(MotifLength, 5) MotifSize;
typedef 16 PeNumMotif;
typedef TMul#(MotifSize, PeNumMotif) MotifRelayLength;
typedef TDiv#(SeqNum, PeNumMotif) MotifRelaySize; 	// 2048


interface ScoreCalculatorIfc;
	method Action putMotifUnchanged(Bit#(MotifRelayLength) m);
	method Action putMotifChanged(Bit#(MotifSize) m);
	method ActionValue#(Bit#(32)) get;
endinterface
(* synthesize *)
module mkScoreCalculator(ScoreCalculatorIfc); // 2048 + 2 cycles
	// I/O
	FIFO#(Bit#(MotifRelayLength)) motifUnchangedQ <- mkSizedBRAMFIFO(valueOf(MotifRelaySize));
	FIFO#(Bit#(MotifSize)) motifChangedQ <- mkFIFO;
	FIFO#(Bit#(32)) scoreQ <- mkFIFO;
	//--------------------------------------------------------------------------------------------
	// Phase 1
	//--------------------------------------------------------------------------------------------
	// Update the motif
	//--------------------------------------------------------------------------------------------
	FIFO#(Bit#(MotifRelayLength)) motifTmpQ <- mkFIFO;
	Reg#(Bool) changeMotifOn <- mkReg(True);
	rule changeMotif1( changeMotifOn );
		motifUnchangedQ.deq;
		let m = motifUnchangedQ.first;
		motifTmpQ.enq(m);
		changeMotifOn <= False;
	endrule
	rule changeMotif2; // 1 cycle
		motifTmpQ.deq;
		motifChangedQ.deq;
		let t = motifTmpQ.first;
		let m = motifChangedQ.first;
		t[79:0] = m;
		motifUnchangedQ.enq(t);
	endrule
	//--------------------------------------------------------------------------------------------
	// Pick the most frequent letter at each motif position
	//--------------------------------------------------------------------------------------------
	FIFO#(Bit#(MotifRelayLength)) motifQ <- mkSizedBRAMFIFO(valueOf(MotifRelaySize));
	FIFO#(Vector#(MotifLength, Vector#(20, Bit#(32)))) baseQ <- mkFIFO;
	Reg#(Vector#(MotifLength, Vector#(20, Bit#(32)))) baseR <- mkReg(replicate(replicate(0)));
	Reg#(Bit#(32)) pickLetterCnt <- mkReg(0);
	rule pickLetter1( !changeMotifOn ); // 2047 cycles + 1 cycle
		motifUnchangedQ.deq;
		let m = motifUnchangedQ.first;

		Vector#(MotifLength, Vector#(20, Bit#(32))) base = baseR;
		for ( Integer i = 0; i < valueOf(PeNumMotif); i = i + 1 ) begin
			Bit#(MotifSize) motif = truncate(m >> (i * valueOf(MotifSize)));
			for ( Integer j = 0; j < valueOf(MotifLength); j = j + 1 ) begin
				Bit#(5) c = truncate(motif >> (j * 5));
				if ( c == 5'b00000 ) base[j][0] = base[j][0] + 1;
				else if ( c == 5'b00001 ) base[j][1] = base[j][1] + 1;
				else if ( c == 5'b00010 ) base[j][2] = base[j][2] + 1;
				else if ( c == 5'b00011 ) base[j][3] = base[j][3] + 1;
				else if ( c == 5'b00100 ) base[j][4] = base[j][4] + 1;
				else if ( c == 5'b00101 ) base[j][5] = base[j][5] + 1;
				else if ( c == 5'b00110 ) base[j][6] = base[j][6] + 1;
				else if ( c == 5'b00111 ) base[j][7] = base[j][7] + 1;
				else if ( c == 5'b01000 ) base[j][8] = base[j][8] + 1;
				else if ( c == 5'b01001 ) base[j][9] = base[j][9] + 1;
				else if ( c == 5'b01010 ) base[j][10] = base[j][10] + 1;
				else if ( c == 5'b01011 ) base[j][11] = base[j][11] + 1;
				else if ( c == 5'b01100 ) base[j][12] = base[j][12] + 1;
				else if ( c == 5'b01101 ) base[j][13] = base[j][13] + 1;
				else if ( c == 5'b01110 ) base[j][14] = base[j][14] + 1;
				else if ( c == 5'b01111 ) base[j][15] = base[j][15] + 1;
				else if ( c == 5'b10000 ) base[j][16] = base[j][16] + 1;
				else if ( c == 5'b10001 ) base[j][17] = base[j][17] + 1;
				else if ( c == 5'b10010 ) base[j][18] = base[j][18] + 1;
				else if ( c == 5'b10011 ) base[j][19] = base[j][19] + 1;
			end
		end

		if ( pickLetterCnt == fromInteger(valueOf(MotifRelaySize)) ) begin
			baseQ.enq(base);
			pickLetterCnt <= 0;
		end else begin
			pickLetterCnt <= pickLetterCnt + 1;
		end

		motifQ.enq(m);
	endrule
	FIFO#(Vector#(MotifLength, Bit#(5))) patternQ <- mkFIFO;
	rule pickLetter2; // 1 cycle
		baseQ.deq;
		let base = baseQ.first;	

		Vector#(MotifLength, Bit#(5)) pattern = replicate(0);
		for ( Integer i = 0; i < valueOf(MotifLength); i = i + 1 ) begin
			if ( (base[i][0] >= base[i][1]) && (base[i][0] >= base[i][2]) && (base[i][0] >= base[i][3]) &&
			     (base[i][0] >= base[i][4]) && (base[i][0] >= base[i][5]) && (base[i][0] >= base[i][6]) &&
		     	     (base[i][0] >= base[i][7]) && (base[i][0] >= base[i][8]) && (base[i][0] >= base[i][9]) &&
		    	     (base[i][0] >= base[i][10]) && (base[i][0] >= base[i][11]) && (base[i][0] >= base[i][12]) && 
		     	     (base[i][0] >= base[i][13]) && (base[i][0] >= base[i][14]) && (base[i][0] >= base[i][15]) && 
		     	     (base[i][0] >= base[i][16]) && (base[i][0] >= base[i][17]) && (base[i][0] >= base[i][18]) &&
		     	     (base[i][0] >= base[i][19]) ) begin
				pattern[i] = 5'b00000;
			end else if ( (base[i][1] >= base[i][2]) && (base[i][1] >= base[i][3]) && (base[i][1] >= base[i][4]) &&
		       		      (base[i][1] >= base[i][5]) && (base[i][1] >= base[i][6]) && (base[i][1] >= base[i][7]) &&	
			      	      (base[i][1] >= base[i][8]) && (base[i][1] >= base[i][9]) && (base[i][1] >= base[i][10]) &&
			      	      (base[i][1] >= base[i][11]) && (base[i][1] >= base[i][12]) && (base[i][1] >= base[i][13]) && 
			      	      (base[i][1] >= base[i][14]) && (base[i][1] >= base[i][15]) && (base[i][1] >= base[i][16]) && 
			      	      (base[i][1] >= base[i][17]) && (base[i][1] >= base[i][18]) && (base[i][1] >= base[i][19]) ) begin
				pattern[i] = 5'b00001;
			end else if ( (base[i][2] >= base[i][3]) && (base[i][2] >= base[i][4]) && (base[i][2] >= base[i][5]) && 
				      (base[i][2] >= base[i][6]) && (base[i][2] >= base[i][7]) && (base[i][2] >= base[i][8]) &&
			     	      (base[i][2] >= base[i][9]) && (base[i][2] >= base[i][10]) && (base[i][2] >= base[i][11]) && 
			      	      (base[i][2] >= base[i][12]) && (base[i][2] >= base[i][13]) && (base[i][2] >= base[i][14]) && 
			      	      (base[i][2] >= base[i][15]) && (base[i][2] >= base[i][16]) && (base[i][2] >= base[i][17]) && 
			      	      (base[i][2] >= base[i][18]) && (base[i][2] >= base[i][19]) ) begin
				pattern[i] = 5'b00010;
			end else if ( (base[i][3] >= base[i][4]) && (base[i][3] >= base[i][5]) && (base[i][3] >= base[i][6]) && 
				      (base[i][3] >= base[i][7]) && (base[i][3] >= base[i][8]) && (base[i][3] >= base[i][9]) && 
			      	      (base[i][3] >= base[i][10]) && (base[i][3] >= base[i][11]) && (base[i][3] >= base[i][12]) && 
			      	      (base[i][3] >= base[i][13]) && (base[i][3] >= base[i][14]) && (base[i][3] >= base[i][15]) &&
			     	      (base[i][3] >= base[i][16]) && (base[i][3] >= base[i][17]) && (base[i][3] >= base[i][18]) && 
			      	      (base[i][3] >= base[i][19]) ) begin
				pattern[i] = 5'b00011;
			end else if ( (base[i][4] >= base[i][5]) && (base[i][4] >= base[i][6]) && (base[i][4] >= base[i][7]) &&
		        	      (base[i][4] >= base[i][8]) && (base[i][4] >= base[i][9]) && (base[i][4] >= base[i][10]) &&
		      		      (base[i][4] >= base[i][11]) && (base[i][4] >= base[i][12]) && (base[i][4] >= base[i][13]) &&
				      (base[i][4] >= base[i][14]) && (base[i][4] >= base[i][15]) && (base[i][4] >= base[i][16]) &&	
			      	      (base[i][4] >= base[i][17]) && (base[i][4] >= base[i][18]) && (base[i][4] >= base[i][19]) ) begin
				pattern[i] = 5'b00100;
			end else if ( (base[i][5] >= base[i][6]) && (base[i][5] >= base[i][7]) && (base[i][5] >= base[i][8]) &&
		       		      (base[i][5] >= base[i][9]) && (base[i][5] >= base[i][10]) && (base[i][5] >= base[i][11]) &&
		      		      (base[i][5] >= base[i][12]) && (base[i][5] >= base[i][13]) && (base[i][5] >= base[i][14]) &&
			      	      (base[i][5] >= base[i][15]) && (base[i][5] >= base[i][16]) && (base[i][5] >= base[i][17]) &&
			      	      (base[i][5] >= base[i][18]) && (base[i][5] >= base[i][19]) ) begin
				pattern[i] = 5'b00101;
			end else if ( (base[i][6] >= base[i][7]) && (base[i][6] >= base[i][8]) && (base[i][6] >= base[i][9]) && 
				      (base[i][6] >= base[i][10]) && (base[i][6] >= base[i][11]) && (base[i][6] >= base[i][12]) &&
			      	      (base[i][6] >= base[i][13]) && (base[i][6] >= base[i][14]) && (base[i][6] >= base[i][15]) &&
			      	      (base[i][6] >= base[i][16]) && (base[i][6] >= base[i][17]) && (base[i][6] >= base[i][18]) &&
			      	      (base[i][6] >= base[i][19]) ) begin
				pattern[i] = 5'b00110;
			end else if ( (base[i][7] >= base[i][8]) && (base[i][7] >= base[i][9]) && (base[i][7] >= base[i][10]) && 
				      (base[i][7] >= base[i][11]) && (base[i][7] >= base[i][12]) && (base[i][7] >= base[i][13]) &&
			      	      (base[i][7] >= base[i][14]) && (base[i][7] >= base[i][15]) && (base[i][7] >= base[i][16]) &&
			      	      (base[i][7] >= base[i][17]) && (base[i][7] >= base[i][18]) && (base[i][7] >= base[i][19]) ) begin
				pattern[i] = 5'b00111;
			end else if ( (base[i][8] >= base[i][9]) && (base[i][8] >= base[i][10]) && (base[i][8] >= base[i][11]) &&
				      (base[i][8] >= base[i][12]) && (base[i][8] >= base[i][13]) && (base[i][8] >= base[i][14]) &&
			      	      (base[i][8] >= base[i][15]) && (base[i][8] >= base[i][16]) && (base[i][8] >= base[i][17]) &&
			      	      (base[i][8] >= base[i][18]) && (base[i][8] >= base[i][19]) ) begin
				pattern[i] = 5'b01000;
			end else if ( (base[i][9] >= base[i][10]) && (base[i][9] >= base[i][11]) && (base[i][9] >= base[i][12]) &&
		       		      (base[i][9] >= base[i][13]) && (base[i][9] >= base[i][14]) && (base[i][9] >= base[i][15]) &&
			      	      (base[i][9] >= base[i][16]) && (base[i][9] >= base[i][17]) && (base[i][9] >= base[i][18]) &&
			      	      (base[i][9] >= base[i][19]) ) begin
				pattern[i] = 5'b01001;
			end else if ( (base[i][10] >= base[i][11]) && (base[i][10] >= base[i][12]) && (base[i][10] >= base[i][13]) &&
		       		      (base[i][10] >= base[i][14]) && (base[i][10] >= base[i][15]) && (base[i][10] >= base[i][16]) &&
			      	      (base[i][10] >= base[i][17]) && (base[i][10] >= base[i][18]) && (base[i][10] >= base[i][19]) ) begin
				pattern[i] = 5'b01010;
			end else if ( (base[i][11] >= base[i][12]) && (base[i][11] >= base[i][13]) && (base[i][11] >= base[i][14]) && 
				      (base[i][11] >= base[i][15]) && (base[i][11] >= base[i][16]) && (base[i][11] >= base[i][17]) &&
			      	      (base[i][11] >= base[i][18]) && (base[i][11] >= base[i][19]) ) begin
				pattern[i] = 5'b01011;
			end else if ( (base[i][12] >= base[i][13]) && (base[i][12] >= base[i][14]) && (base[i][12] >= base[i][15]) && 
				      (base[i][12] >= base[i][16]) && (base[i][12] >= base[i][17]) && (base[i][12] >= base[i][18]) &&
			      	      (base[i][12] >= base[i][19]) ) begin
				pattern[i] = 5'b01100;
			end else if ( (base[i][13] >= base[i][14]) && (base[i][13] >= base[i][15]) && (base[i][13] >= base[i][16]) && 
				      (base[i][13] >= base[i][17]) && (base[i][13] >= base[i][18]) && (base[i][13] >= base[i][19]) ) begin
				pattern[i] = 5'b01101;
			end else if ( (base[i][14] >= base[i][15]) && (base[i][14] >= base[i][16]) && (base[i][14] >= base[i][17]) &&
		       		      (base[i][14] >= base[i][18]) && (base[i][14] >= base[i][19]) ) begin
				pattern[i] = 5'b01110;
			end else if ( (base[i][15] >= base[i][16]) && (base[i][15] >= base[i][17]) && (base[i][15] >= base[i][18]) &&
			      	      (base[i][15] >= base[i][19]) ) begin
				pattern[i] = 5'b01111;
			end else if ( (base[i][16] >= base[i][17]) && (base[i][16] >= base[i][18]) && (base[i][16] >= base[i][19]) ) begin
				pattern[i] = 5'b10000;
			end else if ( (base[i][17] >= base[i][18]) && (base[i][17] >= base[i][19]) ) begin
				pattern[i] = 5'b10001;
			end else if ( (base[i][18] >= base[i][19]) ) begin
				pattern[i] = 5'b10010;
			end else begin
				pattern[i] = 5'b10011;
			end
		end

		patternQ.enq(pattern);
	endrule
	//--------------------------------------------------------------------------------------------
	// Phase 2
	//--------------------------------------------------------------------------------------------
	// Compare between each motif and the picked string
	// Get the score via Hamming Distance
	//--------------------------------------------------------------------------------------------
	Reg#(Vector#(MotifLength, Bit#(5))) patternR <- mkReg(replicate(0));
	Reg#(Bit#(32)) scoreR <- mkReg(0);
	Reg#(Bit#(1)) getScoreCnt1 <- mkReg(0);
	Reg#(Bit#(32)) getScoreCnt2 <- mkReg(0);
	rule getScore; // 2048 cycles
		motifQ.deq;
		let m = motifQ.first;

		Vector#(MotifLength, Bit#(5)) pattern = replicate(0);
		if ( getScoreCnt1 == 0 ) begin
			patternQ.deq;
			let p = patternQ.first;
			pattern = p;
			patternR <= p;
		end else begin
			pattern = patternR;
		end
		
		Bit#(32) score = scoreR;
		for ( Integer i = 0; i < valueOf(PeNumMotif); i = i + 1 ) begin
			Bit#(MotifSize) motif = truncate(m >> (i * valueOf(MotifSize)));
			for ( Integer j = 0; j < valueOf(MotifLength); j = j + 1 ) begin
				Bit#(5) c = truncate(motif >> (j * 5));
				if ( c != pattern[j] ) score = score + 1;
			end
		end

		if ( getScoreCnt2 + 1 == fromInteger(valueOf(MotifRelaySize)) ) begin
			scoreQ.enq(score);
			getScoreCnt1 <= 0;
			getScoreCnt2 <= 0;
		end else begin
			scoreR <= score;
			getScoreCnt1 <= 1;
			getScoreCnt2 <= getScoreCnt2 + 1;
		end
	endrule


	method Action putMotifUnchanged(Bit#(MotifRelayLength) m);
		motifUnchangedQ.enq(m);
	endmethod
	method Action putMotifChanged(Bit#(MotifSize) m);
		motifChangedQ.enq(m);
	endmethod
	method ActionValue#(Bit#(32)) get;
		scoreQ.deq;
		return scoreQ.first;
	endmethod
endmodule
