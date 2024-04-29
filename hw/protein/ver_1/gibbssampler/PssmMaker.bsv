import FIFO::*;
import FIFOF::*;
import Vector::*;

import BRAM::*;
import BRAMFIFO::*;

import FloatingPoint::*;
import Float32::*;
import TypeConverter::*;


// Sequences
typedef 32768 SeqNum;
// Motifs
typedef 16 MotifLength;
typedef TMul#(MotifLength, 5) MotifSize;
typedef 16 PeNumMotif;
typedef TMul#(MotifSize, PeNumMotif) MotifRelayLength;
typedef TDiv#(SeqNum, PeNumMotif) MotifRelaySize; // 2048
// Bases
typedef 20 BaseNum;
typedef 4 PeNumBase;
typedef TDiv#(BaseNum, PeNumBase) BaseRelaySize;


interface PssmMakerIfc;
	method Action putMotif(Bit#(MotifRelayLength) m);
	method Action putStartPe(Bit#(32) s);
	method ActionValue#(Vector#(MotifLength, Vector#(BaseNum, Bit#(32)))) get;
endinterface
(* synthesize *)
module mkPssmMaker(PssmMakerIfc); // 2106 cycles
	// Cycle Counter
	Reg#(Bit#(32)) cycleCount <- mkReg(0);
	rule incCycleCounter;
		cycleCount <= cycleCount + 1;
	endrule

	// FpOp
	Vector#(MotifLength, Vector#(4, FpPairIfc#(32))) fpDiv <- replicateM(replicateM(mkFpDiv32));

	// TypeConverter
	Vector#(MotifLength, Vector#(4, UINTtoFLOATIfc)) typeConverterBase <- replicateM(replicateM(mkUINTtoFLOAT));
	Vector#(MotifLength, UINTtoFLOATIfc) typeConverterTotal <- replicateM(mkUINTtoFLOAT);

	// I/O
	FIFO#(Bit#(MotifRelayLength)) motifQ <- mkFIFO;
	FIFO#(Bit#(32)) startPeQ <- mkFIFO;
	FIFO#(Vector#(MotifLength, Vector#(BaseNum, Bit#(32)))) pssmQ <- mkFIFO;
	//--------------------------------------------------------------------------------------------
	// Make PSSM (Position Specific Score Matrix)
	//--------------------------------------------------------------------------------------------
	// Stage 1
	//--------------------------------------------------------------------------------------------
	FIFO#(Vector#(MotifLength, Vector#(BaseNum, Bit#(32)))) baseForTotalCalQ <- mkFIFO;
	FIFO#(Vector#(MotifLength, Vector#(BaseNum, Bit#(32)))) baseForEachCalQ <- mkFIFO;
	Reg#(Vector#(MotifLength, Vector#(BaseNum, Bit#(32)))) baseR <- mkReg(replicate(replicate(1)));
	Reg#(Bit#(32)) makePssmCnt <- mkReg(0);
	rule makePssm1; // 2048 cycles
		motifQ.deq;
		startPeQ.deq;
		let motif = motifQ.first;
		let startPe = startPeQ.first;

		Vector#(MotifLength, Vector#(BaseNum, Bit#(32))) base = baseR;
		for ( Bit#(32) i = 0; i < fromInteger(valueOf(PeNumMotif)); i = i + 1 ) begin
			if ( i >= startPe ) begin
				Bit#(MotifSize) m = truncate(motif >> (i * fromInteger(valueOf(MotifSize))));
				for ( Bit#(32) j = 0; j < fromInteger(valueOf(MotifLength)); j = j + 1 ) begin
					Bit#(5) c = truncate(m >> (j * 5));
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
		end

		if ( makePssmCnt + 1 == fromInteger(valueOf(MotifRelaySize)) ) begin
			baseForTotalCalQ.enq(base);
			baseForEachCalQ.enq(base);
			makePssmCnt <= 0;
		end else begin
			if ( makePssmCnt == 0 ) begin
				$write("\033[1;33mCycle %1d -> \033[1;33m[PssmMaker]: \033[0m: PssmMaker started!\n", cycleCount);
			end
			baseR <= base;
			makePssmCnt <= makePssmCnt + 1;
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Stage 2
	//--------------------------------------------------------------------------------------------
	Reg#(Vector#(MotifLength, Vector#(BaseNum, Bit#(32)))) baseForEachCalR <- mkReg(replicate(replicate(0)));
	Reg#(Bit#(32)) makePssm2EachCnt <- mkReg(0);
	rule makePssm2_Total; // 28 cycle
		baseForTotalCalQ.deq;
		let base = baseForTotalCalQ.first;

		for ( Integer i = 0; i < valueOf(MotifLength); i = i + 1 ) begin
			Bit#(32) total = base[i][0] + base[i][1] + base[i][2] + base[i][3] + base[i][4] +
					 base[i][5] + base[i][6] + base[i][7] + base[i][8] + base[i][9] +
					 base[i][10] + base[i][11] + base[i][12] + base[i][13] + base[i][14] +
					 base[i][15] + base[i][16] + base[i][17] + base[i][18] + base[i][19];
			typeConverterTotal[i].enq(total);
		end
	endrule
	rule makePssm2_Each;
		if ( makePssm2EachCnt == 0 ) begin
			baseForEachCalQ.deq;
			let base = baseForEachCalQ.first;
			for ( Bit#(32) i = 0; i < fromInteger(valueOf(MotifLength)); i = i + 1 ) begin
				for ( Bit#(32) j = 0; j < 4; j = j + 1 ) begin
					typeConverterBase[i][j].enq(base[i][j + makePssm2EachCnt]);
				end
			end
			baseForEachCalR <= base;
			makePssm2EachCnt <= makePssm2EachCnt + 1;
		end else begin
			let base = baseForEachCalR;
			for ( Bit#(32) i = 0; i < fromInteger(valueOf(MotifLength)); i = i + 1 ) begin
				for ( Bit#(32) j = 0; j < 4; j = j + 1 ) begin
					typeConverterBase[i][j].enq(base[i][j + makePssm2EachCnt]);
				end
			end
			if ( makePssm2EachCnt + 1 == fromInteger(valueOf(BaseRelaySize)) ) begin
				makePssm2EachCnt <= 0;
			end else begin
				makePssm2EachCnt <= makePssm2EachCnt + 1;
			end
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Stage 3
	//--------------------------------------------------------------------------------------------
	Reg#(Vector#(MotifLength, Bit#(32))) totalR <- mkReg(replicate(0));
	Reg#(Bit#(32)) makePssm3Cnt <- mkReg(0);
	rule makePssm3; // 29 cycle
		if ( makePssm3Cnt == 0 ) begin
			Vector#(MotifLength, Bit#(32)) total = replicate(0);
			for ( Integer i = 0; i < valueOf(MotifLength); i = i + 1 ) begin
				total[i] <- typeConverterTotal[i].get;
				Vector#(4, Bit#(32)) base = replicate(0);
				for ( Integer j = 0; j < 4; j = j + 1 ) begin
					base[j] <- typeConverterBase[i][j].get;
					fpDiv[i][j].enq(base[j], total);
				end
			end
			totalR <= total;
			makePssm3Cnt <= makePssm3Cnt + 1;
		end else begin
			let total = totalR;
			for ( Integer i = 0; i < valueOf(MotifLength); i = i + 1 ) begin
				Vector#(4, Bit#(32)) base = replicate(0);
				for ( Integer j = 0; j < 4; j = j + 1 ) begin
					base[j] <- typeConverterBase[i][j].get;
					fpDiv[i][j].enq(base[j], total);
				end
			end
			if ( makePssm3Cnt + 1 == fromInteger(valueOf(BaseRelaySize)) ) begin
				makePssm3Cnt <= 0;
			end else begin
				makePssm3Cnt <= makePssm3Cnt + 1;
			end
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Stage 4
	//--------------------------------------------------------------------------------------------
	Reg#(Vector#(MotifLength, Vector#(BaseNum, Bit#(32)))) pssmR <- mkReg(replicate(replicate(0)));
	Reg#(Bit#(32)) makePssm4Cnt <- mkReg(0);
	rule makePssm4; // 1 cycle
		Vector#(MotifLength, Vector#(BaseNum, Bit#(32))) pssm = pssmR;
		for ( Bit#(32) i = 0; i < valueOf(MotifLength); i = i + 1 ) begin
			for ( Bit#(32) j = 0; j < 4; j = j + 1 ) begin
				fpDiv[i][j].deq;
				pssm[i][j + makePssm4Cnt] = fpDiv[i][j].first;
			end
		end
		if ( makePssm4Cnt + 1 == fromInteger(valueOf(BaseRelaySize)) ) begin
			pssmQ.enq(pssm);
			makePssm4Cnt <= 0;
		end else begin
			pssmR <= pssm;
			makePssm4Cnt <= makePssm4Cnt + 1;
			$write("\033[1;33mCycle %1d -> \033[1;33m[PssmMaker]: \033[0m: PssmMaker finished!\n", cycleCount);
		end
	endrule


	method Action putMotif(Bit#(MotifRelayLength) m);
		motifQ.enq(m);
	endmethod
	method Action putStartPe(Bit#(32) s);
		startPeQ.enq(s);
	endmethod
	method ActionValue#(Vector#(MotifLength, Vector#(BaseNum, Bit#(32)))) get;
		pssmQ.deq;
		return pssmQ.first;
	endmethod
endmodule
