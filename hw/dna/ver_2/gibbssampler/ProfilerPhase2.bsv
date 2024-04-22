import FIFO::*;
import FIFOF::*;
import Vector::*;

import BRAM::*;
import BRAMFIFO::*;

import FloatingPoint::*;
import Float32::*;

import RandomGenerator::*;


// Sequences
typedef 65536 SeqNum;
typedef 1000 SeqLength;
typedef 2048 SeqStoredSize;
typedef TMul#(SeqLength, 2) SeqSize;
// Motifs
typedef 16 MotifLength;
typedef 8 PeNumProfiler;
// Probabilities
typedef TSub#(SeqLength, MotifLength) ProbSizeTmp;
typedef TAdd#(ProbSizeTmp, 1) ProbSize; // 985
typedef 124 IterCnt; // (ProbSize / PeNumProfiler) + 1
typedef 7 RmndCnt; // ProbSize % PeNumProfiler


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
	FpPairIfc#(32) fpAdd <- mkFpAdd32;
	FpPairIfc#(32) fpSub <- mkFpSub32;

        // Profiler Phase2 I/O
	FIFO#(Bit#(SeqSize)) sequenceQ <- mkFIFO;
	Reg#(Bit#(SeqSize)) sequenceR <- mkReg(0);
	FIFO#(Vector#(PeNumProfiler, Bit#(32))) probQ <- mkSizedBRAMFIFO(valueOf(IterCnt));
	Reg#(Vector#(PeNumProfiler, Bit#(32))) probR <- mkReg(replicate(0));
	FIFO#(Bit#(32)) sumQ <- mkFIFO;
	Reg#(Bit#(32)) sumR <- mkReg(0);
	FIFO#(Bit#(32)) resultQ <- mkFIFO;
       
	// Generate a random value between 0.0 and 1.0 
	Reg#(Bool) genRandNumOn <- mkReg(True);
	FIFO#(Bit#(32)) randNumQ <- mkFIFO;
	Reg#(Bit#(32)) randNumR <- mkReg(0);
	rule genRandNum1( genRandNumOn );
		fpRndNumGenerator.req;
	endrule
	rule genRandNum2;
		let r <- fpRndNumGenerator.get;
		randNumQ.enq(r);
		genRandNumOn <= False;
	endrule

	FIFO#(Bit#(32)) partialSumTmpQ <- mkFIFO;
	FIFO#(Bit#(32)) partialSumQ <- mkSizedBRAMFIFO(valueOf(ProbSize));
	Reg#(Bool) getPartialSum1On <- mkReg(True);
	Reg#(Bool) getPartialSum2On <- mkReg(False);
	Reg#(Bit#(32)) getPartialSumCnt1 <- mkReg(0);
	Reg#(Bit#(32)) getPartialSumCnt2 <- mkReg(0);
	rule getPartialSum1( getPartialSum1On );
		if ( getPartialSumCnt1 == 0 ) begin
			probQ.deq;
			let p = probQ.first;
			fpAdd.enq(p[0], p[1]);
			partialSumQ.enq(p[0]);
			probR <= p;
			getPartialSumCnt1 <= getPartialSumCnt1 + 1;
		end else begin
			partialSumTmpQ.deq;
			let ps = partialSumTmpQ.first;
			let p = probR;
			let idx = getPartialSumCnt1 + 1;
			fpAdd.enq(ps, p[idx]);
			if ( getPartialSumCnt2 + 1 == fromInteger(valueOf(IterCnt)) ) begin
				if ( getPartialSumCnt1 + 2 == fromInteger(valueOf(RmndCnt)) ) begin
					getPartialSumCnt1 <= 0;
					getPartialSumCnt2 <= 0;
				end else begin
					getPartialSumCnt1 <= getPartialSumCnt1 + 1;
				end
			end else begin
				if ( getPartialSumCnt1 + 2 == fromInteger(valueOf(PeNumProfiler)) ) begin
					getPartialSumCnt1 <= 0;
					getPartialSumCnt2 <= getPartialSumCnt2 + 1;
				end else begin
					getPartialSumCnt1 <= getPartialSumCnt1 + 1;
				end
			end
		end
		getPartialSum1On <= False;
		getPartialSum2On <= True;
	endrule
	rule getPartialSum2( getPartialSum2On );
		fpAdd.deq;
		let ps = fpAdd.first;
		partialSumTmpQ.enq(ps);
		partialSumQ.enq(ps);
		getPartialSum1On <= True;
		getPartialSum2On <= False;
	endrule

	Reg#(Bit#(32)) divCnt <- mkReg(0);
	rule divPartialByTotal; // 29 cycle
		if ( divCnt == 0 ) begin
			partialSumQ.deq;
			sumQ.deq;
			let p = partialSumQ.first;
			let s = sumQ.first;
			fpDiv.enq(p, s);
			sumR <= s;
			divCnt <= divCnt + 1;
		end else begin
			partialSumQ.deq;
			let p = partialSumQ.first;
			let s = sumR;
			fpDiv.enq(p, s);
			if ( divCnt + 1 == fromInteger(valueOf(ProbSize)) ) begin
				divCnt <= 0;
			end else begin
				divCnt <= divCnt + 1;
			end
		end	
	endrule

	Reg#(Bit#(32)) cmprCnt <- mkReg(0);
	rule cmprWithRandNum; // 12 cycle
		if ( cmprCnt == 0 ) begin
			randNumQ.deq;
			fpDiv.deq;
			let r = randNumQ.first;
			let d = fpDiv.first;
			fpSub.enq(d, r);
			randNumR <= r;
			cmprCnt <= cmprCnt + 1;
		end else begin
			fpDiv.deq;
			let r = randNumR;
			let d = fpDiv.first;
			fpSub.enq(d, r);
			if ( cmprCnt + 1 == fromInteger(valueOf(ProbSize)) ) begin
				cmprCnt <= 0;
			end else begin
				cmprCnt <= cmprCnt + 1;
			end
		end
	endrule

	Reg#(Bit#(32)) pickCnt <- mkReg(0);
	rule pickMotif; // 1 cycle
		if ( pickCnt == 0 ) begin
			sequenceQ.deq;
			fpSub.deq;
			let s = sequenceQ.first;
			let f = fpSub.first;
			if ( f[31] == 0 ) begin
				Bit#(32) motif = truncate(s >> pickCnt);
				resultQ.enq(motif);
			end else begin
				sequenceR <= s;
				pickCnt <= pickCnt + 1;
			end
		end else begin
			fpSub.deq;
			let s = sequenceR;
			let f = fpSub.first;
			if ( f[31] == 0 ) begin
				Bit#(32) motif = truncate(s >> pickCnt);
				resultQ.enq(motif);
			end else begin
				pickCnt <= pickCnt + 1;
			end
		end
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
