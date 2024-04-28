import FIFO::*;
import FIFOF::*;
import Vector::*;

import BRAM::*;
import BRAMFIFO::*;

import FloatingPoint::*;
import Float32::*;


// Sequences
typedef 32768 SeqNum;
typedef 300 SeqLength;
typedef TMul#(SeqLength, 5) SeqSize;
// Motifs
typedef 16 MotifLength;
typedef TMul#(MotifLength, 5) MotifSize;
typedef 8 PeNumProfiler;
// Pipeline to get probabilities
typedef TDiv#(16, 2) CalProbPipe1; 	     		// 8
typedef TDiv#(CalProbPipe1, 2) CalProbPipe2; 		// 4
typedef TDiv#(CalProbPipe2, 2) CalProbPipe3; 		// 2
// Pipeline to get sum of probabilities
typedef TDiv#(PeNumProfiler, 2) SumProbPipe1;		// 4
typedef TDiv#(SumProbPipe1, 2) SumProbPipe2; 		// 2
// Probabilities
typedef TSub#(SeqLength, MotifLength) ProbSizeTmp;
typedef TAdd#(ProbSizeTmp, 1) ProbSize; 		// 285
typedef 36 IterCnt; 					// (ProbSize / PeNumProfiler) + 1
typedef 3 RmndCnt; 					// ProbSize % PeNumProfiler


interface ProfilerPhase1Ifc;
        method Action putSequence(Bit#(SeqSize) s);
	method Action putPssm(Vector#(MotifLength, Vector#(20, Bit#(32))) p);
	method Action putPeNum(Bit#(32) n);
        method ActionValue#(Vector#(PeNumProfiler, Bit#(32))) getProb;
	method ActionValue#(Bit#(32)) getSum;
endinterface
(* synthesize *)
module mkProfilerPhase1(ProfilerPhase1Ifc);
	// Cycle Counter
	Reg#(Bit#(32)) cycleCount <- mkReg(0);
	rule incCycleCounter;
		cycleCount <= cycleCount + 1;
	endrule

	// Pipeline to get the probabilities
	Vector#(PeNumProfiler, Vector#(CalProbPipe1, FpPairIfc#(32))) fpMultPipe1 <- replicateM(replicateM(mkFpMult32));
	Vector#(PeNumProfiler, Vector#(CalProbPipe2, FpPairIfc#(32))) fpMultPipe2 <- replicateM(replicateM(mkFpMult32));
	Vector#(PeNumProfiler, Vector#(CalProbPipe3, FpPairIfc#(32))) fpMultPipe3 <- replicateM(replicateM(mkFpMult32));
	Vector#(PeNumProfiler, FpPairIfc#(32)) fpMultPipe4 <- replicateM(mkFpMult32);
	
	// Pipeline to get a sum of probabilities
	Vector#(PeNumProfiler, FpPairIfc#(32)) fpAdd <- replicateM(mkFpAdd32);
	Vector#(SumProbPipe1, FpPairIfc#(32)) fpAddStage1 <- replicateM(mkFpAdd32);
	Vector#(SumProbPipe2, FpPairIfc#(32)) fpAddStage2 <- replicateM(mkFpAdd32);
	FpPairIfc#(32) fpAddStage3 <- mkFpAdd32;

        // Profiler Phase1 I/O
	FIFO#(Bit#(SeqSize)) sequenceQ <- mkSizedBRAMFIFO(valueOf(IterCnt));
	FIFO#(Vector#(MotifLength, Vector#(20, Bit#(32)))) pssmQ <- mkSizedBRAMFIFO(valueOf(IterCnt));
	FIFO#(Bit#(32)) peCntQ <- mkSizedBRAMFIFO(valueOf(IterCnt));
	FIFO#(Vector#(PeNumProfiler, Bit#(32))) probQ <- mkSizedBRAMFIFO(valueOf(IterCnt));
	FIFO#(Bit#(32)) sumQ <- mkFIFO;
	//--------------------------------------------------------------------------------------------
	// Get probabilities
	//--------------------------------------------------------------------------------------------	
	Reg#(Bit#(32)) getProb1Cnt <- mkReg(0);
	rule getProb1;
		sequenceQ.deq;
		pssmQ.deq;
		peCntQ.deq;
		let s = sequenceQ.first;
		let p = pssmQ.first;
		let n = peCntQ.first;

		Vector#(PeNumProfiler, Vector#(MotifLength, Bit#(32))) prob = replicate(replicate(0));
		for ( Bit#(32) i = 0; i < fromInteger(valueOf(PeNumProfiler)); i = i + 1 ) begin 
			Bit#(SeqSize) sqnc = s >> (fromInteger(valueOf(PeNumProfiler)) * getProb1Cnt);
			Bit#(MotifSize) m = truncate(sqnc >> i);
			if ( i < n ) begin
				for ( Bit#(32) j = 0; j < fromInteger(valueOf(MotifLength)); j = j + 1 ) begin
					Bit#(5) c = truncate(m >> (j * 5));
					if ( c == 5'b00000 ) prob[i][j] = p[j][0];
					else if ( c == 5'b00001 ) prob[i][j] = p[j][1];
					else if ( c == 5'b00010 ) prob[i][j] = p[j][2]; 
					else if ( c == 5'b00011 ) prob[i][j] = p[j][3];
					else if ( c == 5'b00100 ) prob[i][j] = p[j][4];
					else if ( c == 5'b00101 ) prob[i][j] = p[j][5];
					else if ( c == 5'b00110 ) prob[i][j] = p[j][6];
					else if ( c == 5'b00111 ) prob[i][j] = p[j][7];
					else if ( c == 5'b01000 ) prob[i][j] = p[j][8];
					else if ( c == 5'b01001 ) prob[i][j] = p[j][9];
					else if ( c == 5'b01010 ) prob[i][j] = p[j][10];
					else if ( c == 5'b01011 ) prob[i][j] = p[j][11];
					else if ( c == 5'b01100 ) prob[i][j] = p[j][12];
					else if ( c == 5'b01101 ) prob[i][j] = p[j][13];
					else if ( c == 5'b01110 ) prob[i][j] = p[j][14];
					else if ( c == 5'b01111 ) prob[i][j] = p[j][15];
					else if ( c == 5'b10000 ) prob[i][j] = p[j][16];
					else if ( c == 5'b10001 ) prob[i][j] = p[j][17];
					else if ( c == 5'b10010 ) prob[i][j] = p[j][18];
					else if ( c == 5'b10011 ) prob[i][j] = p[j][19];
				end
			end
			
			for ( Integer j = 0; j < valueOf(CalProbPipe1); j = j + 1 ) begin
				fpMultPipe1[i][j].enq(prob[i][j*2], prob[i][j*2+1]);
			end
		end

		if ( getProb1Cnt == 0 ) begin
			getProb1Cnt <= getProb1Cnt + 1;
 			$write("\033[1;33mCycle %1d -> \033[1;33m[ProfilerPhase1]: \033[0m: Started!\n", cycleCount);
		end else begin
			if ( getProb1Cnt + 1 == fromInteger(valueOf(IterCnt)) ) begin
				getProb1Cnt <= 0;
			end else begin
				getProb1Cnt <= getProb1Cnt + 1;
			end
		end
	endrule
	rule getProb2;
		Vector#(PeNumProfiler, Vector#(CalProbPipe1, Bit#(32))) r = replicate(replicate(0));
		for ( Integer i = 0; i < valueOf(PeNumProfiler); i = i + 1 ) begin 
			for ( Integer j = 0; j < valueOf(CalProbPipe1); j = j + 1 ) begin
				fpMultPipe1[i][j].deq;
				r[i][j] = fpMultPipe1[i][j].first;
			end

			for ( Integer j = 0; j < valueOf(CalProbPipe2); j = j + 1 ) begin
				fpMultPipe2[i][j].enq(r[i][j*2], r[i][j*2+1]);
			end
		end
	endrule
	rule getProb3;
		Vector#(PeNumProfiler, Vector#(CalProbPipe2, Bit#(32))) r = replicate(replicate(0));
		for ( Integer i = 0; i < valueOf(PeNumProfiler); i = i + 1 ) begin 
			for ( Integer j = 0; j < valueOf(CalProbPipe2); j = j + 1 ) begin
				fpMultPipe2[i][j].deq;
				r[i][j] = fpMultPipe2[i][j].first;
			end

			for ( Integer j = 0; j < valueOf(CalProbPipe3); j = j + 1 ) begin
				fpMultPipe3[i][j].enq(r[i][j*2], r[i][j*2+1]);
			end
		end
	endrule
	rule getProb4;
		Vector#(PeNumProfiler, Vector#(CalProbPipe3, Bit#(32))) r = replicate(replicate(0));
		for ( Integer i = 0; i < valueOf(PeNumProfiler); i = i + 1 ) begin 
			for ( Integer j = 0; j < valueOf(CalProbPipe3); j = j + 1 ) begin
				fpMultPipe3[i][j].deq;
				r[i][j] = fpMultPipe3[i][j].first;
			end

			fpMultPipe4[i].enq(r[i][0], r[i][1]);
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Add #PeNumProfiler probabilities at the same time first
	//--------------------------------------------------------------------------------------------	
	Reg#(Bool) getSumTmp1On <- mkReg(True);
	Reg#(Bool) getSumTmp2On <- mkReg(False);
	FIFO#(Vector#(PeNumProfiler, Bit#(32))) sumTmp1Q <- mkFIFO;
	FIFO#(Vector#(PeNumProfiler, Bit#(32))) sumTmp2Q <- mkFIFO;
	Reg#(Bit#(32)) getSumTmp1Cnt <- mkReg(0);
	Reg#(Bit#(32)) getSumTmp2Cnt <- mkReg(0);
	rule getSumTmp1( getSumTmp1On );
		Vector#(PeNumProfiler, Bit#(32)) r = replicate(0);
		for ( Integer i = 0; i < valueOf(PeNumProfiler); i = i + 1 ) begin
			fpMultPipe4[i].deq;
			r[i] = fpMultPipe4[i].first;
		end

		if ( getSumTmp1Cnt == 0 ) begin
			sumTmp1Q.enq(r);
			getSumTmp1Cnt <= getSumTmp1Cnt + 1;
		end else begin
			sumTmp1Q.deq;
			let s = sumTmp1Q.first;
			for ( Integer i = 0; i < valueOf(PeNumProfiler); i = i + 1 ) begin
				fpAdd[i].enq(r[i], s[i]);
			end
			if ( getSumTmp1Cnt + 1 == fromInteger(valueOf(IterCnt)) ) begin
				getSumTmp1Cnt <= 0;
			end else begin
				getSumTmp1Cnt <= getSumTmp1Cnt + 1;
			end
			getSumTmp1On <= False;
			getSumTmp2On <= True;
		end
		probQ.enq(r);
	endrule
	rule getSumTmp2( getSumTmp2On );
		Vector#(PeNumProfiler, Bit#(32)) s = replicate(0);
		for ( Integer i = 0; i < valueOf(PeNumProfiler); i = i + 1 ) begin
			fpAdd[i].deq;
			s[i] = fpAdd[i].first;
		end
		if ( getSumTmp2Cnt + 2 == fromInteger(valueOf(IterCnt)) ) begin
			sumTmp2Q.enq(s);
			getSumTmp2Cnt <= 0;
		end else begin
			sumTmp1Q.enq(s);
			getSumTmp2Cnt <= getSumTmp2Cnt + 1;
		end
		getSumTmp1On <= True;
		getSumTmp2On <= False;
	endrule
	//--------------------------------------------------------------------------------------------
	// Add result #PeNumProfiler intermediate probabilities then
	//--------------------------------------------------------------------------------------------	
	rule getSum1;
		sumTmp2Q.deq;
		let p = sumTmp2Q.first;
		for ( Integer i = 0; i < valueOf(SumProbPipe1); i = i + 1 ) begin
			fpAddStage1[i].enq(p[i*2], p[i*2+1]);
		end
	endrule
	rule getSum2;
		Vector#(SumProbPipe1, Bit#(32)) s = replicate(0);
		for ( Integer i = 0; i < valueOf(SumProbPipe1); i = i + 1 ) begin
			fpAddStage1[i].deq;
			s[i] = fpAddStage1[i].first;
		end

		for ( Integer i = 0; i < valueOf(SumProbPipe2); i = i + 1 ) begin
			fpAddStage2[i].enq(s[i*2], s[i*2+1]);
		end
	endrule
	rule getSum3;
		Vector#(SumProbPipe2, Bit#(32)) s = replicate(0);
		for ( Integer i = 0; i < valueOf(SumProbPipe2); i = i + 1 ) begin
			fpAddStage2[i].deq;
			s[i] = fpAddStage2[i].first;
		end
		fpAddStage3.enq(s[0], s[1]);
	endrule
	rule getSum4;
		fpAddStage3.deq;
		sumQ.enq(fpAddStage3.first);
		$write("\033[1;33mCycle %1d -> \033[1;33m[ProfilerPhase1]: \033[0m: Finished!\n", cycleCount);
	endrule


        method Action putSequence(Bit#(SeqSize) s);
                sequenceQ.enq(s);
	endmethod
	method Action putPssm(Vector#(MotifLength, Vector#(20, Bit#(32))) p);
		pssmQ.enq(p);
	endmethod
	method Action putPeNum(Bit#(32) n);
		peCntQ.enq(n);
	endmethod
	method ActionValue#(Vector#(PeNumProfiler, Bit#(32))) getProb;
		probQ.deq;
		return probQ.first;
	endmethod
        method ActionValue#(Bit#(32)) getSum;
		sumQ.deq;
              	return sumQ.first;
        endmethod
endmodule
