import FIFO::*;
import FIFOF::*;
import Vector::*;

import BRAM::*;
import BRAMFIFO::*;

import FloatingPoint::*;
import Float32::*;

import RandomGenerator::*;

import ProfilerPhase1::*;
import ProfilerPhase2::*;


// Sequences
typedef 32768 SeqNum;
typedef 1000 SeqLength;
typedef 2048 SeqStoredSize;
typedef TMul#(SeqLength, 2) SeqSize;
// Motifs
typedef 16 MotifLength;
typedef 32 PeNumProfiler;
// Probabilities
typedef TSub#(SeqLength, MotifLength) ProbSizeTmp;
typedef TAdd#(ProbSizeTmp, 1) ProbSize; // 985
typedef 31 IterCnt; // (ProbSize / PeNum) + 1
typedef 25 RmndCnt; // ProbSize % PeNum


interface ProfilerIfc;
        method Action putSequence(Bit#(SeqSize) s);
	method Action putPssm(Vector#(MotifLength, Vector#(4, Bit#(32))) p);
        method ActionValue#(Bit#(32)) get;
endinterface
(* synthesize *)
module mkProfiler(ProfilerIfc);
	// Cycle Counter
	Reg#(Bit#(32)) cycleCount <- mkReg(0);
	rule incCycleCounter;
		cycleCount <= cycleCount + 1;
	endrule
	
	ProfilerPhase1Ifc profilerPhase1 <- mkProfilerPhase1;
	ProfilerPhase2Ifc profilerPhase2 <- mkProfilerPhase2;

        // Profiler I/O
	FIFO#(Bit#(SeqSize)) sequenceQ <- mkFIFO;
	Reg#(Bit#(SeqSize)) sequenceR <- mkReg(0);
	FIFO#(Vector#(MotifLength, Vector#(4, Bit#(32)))) pssmQ <- mkFIFO;
	Reg#(Vector#(MotifLength, Vector#(4, Bit#(32)))) pssmR <- mkReg(replicate(replicate(0)));
	FIFO#(Bit#(32)) resultQ <- mkFIFO;

	Reg#(Bit#(32)) relaySeqCnt <- mkReg(0);
	rule relaySequence;
		if ( relaySeqCnt == 0 ) begin
			sequenceQ.deq;
			let s = sequenceQ.first;
			profilerPhase1.putSequence(s);
			profilerPhase2.putSequence(s);
			sequenceR <= s;
			relaySeqCnt <= relaySeqCnt + 1;
		end else begin
			let s = sequenceR;
			profilerPhase1.putSequence(s);
			if ( relaySeqCnt + 1 == fromInteger(valueOf(IterCnt)) ) begin
				relaySeqCnt <= 0;
			end else begin
				relaySeqCnt <= relaySeqCnt + 1;
			end
		end
	endrule

	Reg#(Bit#(32)) relayPssmCnt <- mkReg(0);
	rule relayPssmNPeNum;
		if ( relayPssmCnt == 0 ) begin
			pssmQ.deq;
			let p = pssmQ.first;
			profilerPhase1.putPssm(p);
			profilerPhase1.putPeNum(fromInteger(valueOf(PeNumProfiler)));
			pssmR <= p;
			relayPssmCnt <= relayPssmCnt + 1;
		end else begin
			let p = pssmR;
			profilerPhase1.putPssm(p);
			if ( relayPssmCnt + 1 == fromInteger(valueOf(IterCnt)) ) begin
				profilerPhase1.putPeNum(fromInteger(valueOf(RmndCnt)));
				relayPssmCnt <= 0;
			end else begin
				profilerPhase1.putPeNum(fromInteger(valueOf(PeNumProfiler)));
				relayPssmCnt <= relayPssmCnt + 1;
			end
		end
	endrule

	rule getProbabilites;
		let p <- profilerPhase1.getProb;
		profilerPhase2.putProb(p);
	endrule

	rule getSum; // 443 + 1 cycles
		let s <- profilerPhase1.getSum;
		profilerPhase2.putSum(s);
	endrule

	rule getResult; // 42 + 1 cycles
		let r <- profilerPhase2.get;
		resultQ.enq(r);
	endrule


	method Action putSequence(Bit#(SeqSize) s);
                sequenceQ.enq(s);
	endmethod
	method Action putPssm(Vector#(MotifLength, Vector#(4, Bit#(32))) p);
		pssmQ.enq(p);
	endmethod
        method ActionValue#(Bit#(32)) get;
		resultQ.deq;
              	return resultQ.first;
        endmethod
endmodule

