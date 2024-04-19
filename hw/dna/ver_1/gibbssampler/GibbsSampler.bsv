import FIFO::*;
import FIFOF::*;
import Vector::*;

import BRAM::*;
import BRAMFIFO::*;

import ScoreCalculator::*;
import PssmMaker::*;
import Profiler::*;


// Sequences
typedef 32768 SeqNum;
typedef 1000 SeqLength;
typedef 2048 SeqStoredSize;
typedef TMul#(SeqLength, 2) SeqSize;
typedef TMul#(SeqNum, SeqStoredSize) DataSize;
// Motifs
typedef 16 MotifLength;
typedef TMul#(MotifLength, 2) MotifSize;
typedef 64 PeNumMotif;
typedef TMul#(MotifSize, PeNumMotif) MotifRelayLength;
typedef TDiv#(SeqNum, PeNumMotif) MotifRelaySize;


interface GibbsSamplerIfc;
	method Action putSequence(Bit#(SeqSize) s);
	method Action putMotif(Bit#(MotifRelayLength) m);
endinterface
(* synthesize *)
module mkGibbsSampler(GibbsSamplerIfc); // 1574 cycles
	// Cycle Counter
	Reg#(Bit#(32)) cycleCount <- mkReg(0);
	rule incCycleCounter;
		cycleCount <= cycleCount + 1;
	endrule

	// Score Calculator
	ScoreCalculatorIfc scoreCalculator <- mkScoreCalculator;

	// PSSM Maker
	PssmMakerIfc pssmMaker <- mkPssmMaker;

	// Profiler
	ProfilerIfc profiler <- mkProfiler;

	// I/O
	FIFO#(Bit#(SeqSize)) sequenceQ <- mkFIFO;
	FIFO#(Bit#(MotifRelayLength)) motifQ <- mkSizedBRAMFIFO(valueOf(MotifRelaySize));
	//--------------------------------------------------------------------------------------------
	// Get PSSM First
	//--------------------------------------------------------------------------------------------
	Reg#(Bit#(32)) makePssmCnt <- mkReg(0);
	rule makePssm;
		motifQ.deq;
		let m = motifQ.first;
		pssmMaker.putMotif(m);
		scoreCalculator.putMotifUnchanged(m);
		if ( makePssmCnt + 1 == fromInteger(valueOf(MotifRelaySize)) ) begin
			pssmMaker.putStartPe(0);
			makePssmCnt <= 0;
		end else begin
			if ( makePssmCnt == 0 ) begin
				pssmMaker.putStartPe(1);
			end else begin
				pssmMaker.putStartPe(0);
			end
			makePssmCnt <= makePssmCnt + 1;
		end
	endrule

	rule relaySequenceToProfiler;
		sequenceQ.deq;
		let s = sequenceQ.first;
		profiler.putSequence(s);
	endrule
	rule relayPssmToProfiler; // 570 + 1 cyckes
		let p <- pssmMaker.get;
		profiler.putPssm(p);
	endrule

	rule calScore; // 487 + 1 cycles
		let m <- profiler.get;
		scoreCalculator.putMotifChanged(m);
	endrule

	rule getResult; // 514 + 1 cycles 
		let s <- scoreCalculator.get;
		$write("\033[1;33mCycle %1d -> \033[1;33m[GibbsSampler]: \033[0m: GibbsSampler finished!\n", cycleCount);
	endrule


	method Action putSequence(Bit#(SeqSize) s);
		sequenceQ.enq(s);
	endmethod
	method Action putMotif(Bit#(MotifRelayLength) m);
		motifQ.enq(m);
	endmethod
endmodule
