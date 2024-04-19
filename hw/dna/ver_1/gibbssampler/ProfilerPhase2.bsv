import FIFO::*;
import FIFOF::*;
import Vector::*;

import BRAM::*;
import BRAMFIFO::*;

import FloatingPoint::*;
import Float32::*;

import RandomGenerator::*;


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
typedef 31 IterCnt; // (ProbSize / PeNumProfiler) + 1
typedef 25 RmndCnt; // ProbSize % PeNumProfiler


interface ProfilerPhase2Ifc;
        method Action putSequence(Bit#(SeqSize) s);
	method Action putProb(Vector#(PeNumProfiler, Bit#(32)) p);
        method Action putSum(Bit#(32) s);
	method ActionValue#(Bit#(32)) get;
endinterface
(* synthesize *)
module mkProfilerPhase2(ProfilerPhase2Ifc);
	// Random Generator
	RandomGeneratorFpIfc fpRndNumGenerator <- mkRandomGeneratorFp;

	// FpOp
	FpPairIfc#(32) fpDiv <- mkFpDiv32;
	FpPairIfc#(32) fpSub <- mkFpSub32;

        // Profiler Phase2 I/O
	FIFO#(Bit#(2000)) sequenceQ <- mkFIFO;
	FIFO#(Vector#(PeNumProfiler, Bit#(32))) probQ <- mkSizedBRAMFIFO(valueOf(IterCnt));
	FIFO#(Bit#(32)) sumQ <- mkFIFO;
	FIFO#(Bit#(32)) resultQ <- mkFIFO;
       
	// Generate a random value between 0.0 and 1.0 
	Reg#(Bool) genRandNumOn <- mkReg(True);
	FIFO#(Bit#(32)) randNumQ <- mkFIFO;
	rule genRandNum1( genRandNumOn );
		fpRndNumGenerator.req;
	endrule
	rule genRandNum2;
		let r <- fpRndNumGenerator.get;
		randNumQ.enq(r);
		genRandNumOn <= False;
	endrule

	rule getPartialSum; // 29 cycle
		probQ.deq;
		sumQ.deq;
		let p = probQ.first;
		let s = sumQ.first;
		fpDiv.enq(p[0], s);
	endrule

	rule pickRandPosition; // 12 cycle
		randNumQ.deq;
		fpDiv.deq;
		let d = fpDiv.first;
		let r = randNumQ.first;
		fpSub.enq(d, r);
	endrule

	rule phase25; // 1 cycle
		fpSub.deq;
		sequenceQ.deq;
		let s = sequenceQ.first;
		Bit#(32) updatedMotif = truncate(s);
		resultQ.enq(updatedMotif);
	endrule


        method Action putSequence(Bit#(SeqSize) s);
                sequenceQ.enq(s);
	endmethod
	method Action putProb(Vector#(PeNumProfiler, Bit#(32)) p);
		probQ.enq(p);
	endmethod
	method Action putSum(Bit#(32) s);
		sumQ.enq(s);
	endmethod
        method ActionValue#(Bit#(32)) get;
		resultQ.deq;
              	return resultQ.first;
        endmethod
endmodule

