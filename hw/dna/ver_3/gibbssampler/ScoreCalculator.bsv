import FIFO::*;
import FIFOF::*;
import Vector::*;

import BRAM::*;
import BRAMFIFO::*;


// Sequences
typedef 131072 SeqNum;
typedef 1000 SeqLength;
typedef 2048 SeqStoredSize;
typedef TMul#(SeqLength, 2) SeqSize;
// Motif
typedef 16 MotifLength;
typedef TMul#(MotifLength, 2) MotifSize;
typedef 16 PeNumMotif;
typedef TMul#(MotifSize, PeNumMotif) MotifRelayLength;
typedef TDiv#(SeqNum, PeNumMotif) MotifRelaySize; // 8192


interface ScoreCalculatorIfc;
	method Action putMotifUnchanged(Bit#(MotifRelayLength) m);
	method Action putMotifChanged(Bit#(32) m);
	method ActionValue#(Bit#(32)) get;
endinterface
(* synthesize *)
module mkScoreCalculator(ScoreCalculatorIfc); // 8192 + 2 cycles
	// I/O
	FIFO#(Bit#(MotifRelayLength)) motifUnchangedQ <- mkSizedBRAMFIFO(valueOf(MotifRelaySize));
	FIFO#(Bit#(32)) motifChangedQ <- mkFIFO;
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
		t[31:0] = m;
		motifUnchangedQ.enq(t);
	endrule
	//--------------------------------------------------------------------------------------------
	// Pick the most frequent letter at each motif position
	//--------------------------------------------------------------------------------------------
	FIFO#(Bit#(MotifRelayLength)) motifQ <- mkSizedBRAMFIFO(valueOf(MotifRelaySize));
	FIFO#(Vector#(MotifLength, Vector#(4, Bit#(32)))) baseQ <- mkFIFO;
	Reg#(Vector#(MotifLength, Vector#(4, Bit#(32)))) baseR <- mkReg(replicate(replicate(0)));
	Reg#(Bit#(32)) pickLetterCnt <- mkReg(0);
	rule pickLetter1( !changeMotifOn ); // 8191 cycles + 1 cycle
		motifUnchangedQ.deq;
		let m = motifUnchangedQ.first;

		Vector#(MotifLength, Vector#(4, Bit#(32))) base = baseR;
		for ( Integer i = 0; i < valueOf(PeNumMotif); i = i + 1 ) begin
			Bit#(32) motif = truncate(m >> (i * 32));
			for ( Integer j = 0; j < valueOf(MotifLength); j = j + 1 ) begin
				Bit#(2) c = truncate(motif >> (j * 2));
				if ( c == 2'b00 ) base[j][0] = base[j][0] + 1;
				else if ( c == 2'b01 ) base[j][1] = base[j][1] + 1;
				else if ( c == 2'b10 ) base[j][2] = base[j][2] + 1;
				else if ( c == 2'b11 ) base[j][3] = base[j][3] + 1;
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
	FIFO#(Vector#(MotifLength, Bit#(2))) patternQ <- mkFIFO;
	rule pickLetter2; // 1 cycle
		baseQ.deq;
		let base = baseQ.first;	

		Vector#(MotifLength, Bit#(2)) pattern = replicate(0);
		for ( Integer i = 0; i < valueOf(MotifLength); i = i + 1 ) begin
			if ( (base[i][0] >= base[i][1]) && (base[i][0] >= base[i][2]) && (base[i][0] >= base[i][3]) ) begin
				pattern[i] = 2'b00;
			end else if ( (base[i][1] >= base[i][2]) && (base[i][1] >= base[i][3]) ) begin
				pattern[i] = 2'b01;
			end else if ( base[i][2] >= base[i][3] ) begin
				pattern[i] = 2'b10;
			end else begin
				pattern[i] = 2'b11;
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
	Reg#(Vector#(MotifLength, Bit#(2))) patternR <- mkReg(replicate(0));
	Reg#(Bit#(32)) scoreR <- mkReg(0);
	Reg#(Bit#(1)) getScoreCnt1 <- mkReg(0);
	Reg#(Bit#(32)) getScoreCnt2 <- mkReg(0);
	rule getScore; // 8192 cycles
		motifQ.deq;
		let m = motifQ.first;

		Vector#(MotifLength, Bit#(2)) pattern = replicate(0);
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
			Bit#(32) motif = truncate(m >> (i * 32));
			for ( Integer j = 0; j < valueOf(MotifLength); j = j + 1 ) begin
				Bit#(2) c = truncate(motif >> (j * 2));
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
	method Action putMotifChanged(Bit#(32) m);
		motifChangedQ.enq(m);
	endmethod
	method ActionValue#(Bit#(32)) get;
		scoreQ.deq;
		return scoreQ.first;
	endmethod
endmodule
