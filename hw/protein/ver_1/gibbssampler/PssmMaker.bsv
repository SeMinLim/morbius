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
typedef TDiv#(SeqNum, PeNumMotif) MotifRelaySize; 	// 2048


interface PssmMakerIfc;
	method Action putMotif(Bit#(MotifRelayLength) m);
	method Action putStartPe(Bit#(32) s);
	method ActionValue#(Vector#(MotifLength, Vector#(20, Bit#(32)))) get;
endinterface
(* synthesize *)
module mkPssmMaker(PssmMakerIfc); // 2106 cycles
	// Cycle Counter
	Reg#(Bit#(32)) cycleCount <- mkReg(0);
	rule incCycleCounter;
		cycleCount <= cycleCount + 1;
	endrule

	// FpOp
	Vector#(MotifLength, Vector#(20, FpPairIfc#(32))) fpDiv <- replicateM(replicateM(mkFpDiv32));

	// TypeConverter
	Vector#(MotifLength, Vector#(20, UINTtoFLOATIfc)) typeConverterBase <- replicateM(replicateM(mkUINTtoFLOAT));
	Vector#(MotifLength, UINTtoFLOATIfc) typeConverterTotal <- replicateM(mkUINTtoFLOAT);

	// I/O
	FIFO#(Bit#(MotifRelayLength)) motifQ <- mkSizedBRAMFIFO(valueOf(MotifRelaySize));
	FIFO#(Bit#(32)) startPeQ <- mkSizedBRAMFIFO(valueOf(MotifRelaySize));
	FIFO#(Vector#(MotifLength, Vector#(20, Bit#(32)))) pssmQ <- mkFIFO;
	//--------------------------------------------------------------------------------------------
	// Make PSSM (Position Specific Score Matrix)
	//--------------------------------------------------------------------------------------------
	FIFO#(Vector#(MotifLength, Vector#(20, Bit#(32)))) baseQ <- mkFIFO;
	Reg#(Vector#(MotifLength, Vector#(20, Bit#(32)))) baseR <- mkReg(replicate(replicate(1)));
	Reg#(Bit#(32)) makePssmCnt <- mkReg(0);
	rule makePssm1; // 2048 cycles
		motifQ.deq;
		startPeQ.deq;
		let motif = motifQ.first;
		let startPe = startPeQ.first;

		Vector#(MotifLength, Vector#(20, Bit#(32))) base = baseR;
		for ( Bit#(32) i = 0; i < fromInteger(valueOf(PeNumMotif)); i = i + 1 ) begin
			if ( i >= startPe ) begin
				Bit#(MotifSize) m = truncate(motif >> (i * valueOf(MotifSize)));
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
			baseQ.enq(base);
			makePssmCnt <= 0;
		end else begin
			if ( makePssmCnt == 0 ) begin
				$write("\033[1;33mCycle %1d -> \033[1;33m[PssmMaker]: \033[0m: PssmMaker started!\n", cycleCount);
			end
			baseR <= base;
			makePssmCnt <= makePssmCnt + 1;
		end
	endrule
	rule makePssm2; // 28 cycle
		baseQ.deq;
		let base = baseQ.first;

		Vector#(MotifLength, Bit#(32)) total = replicate(0);
		for ( Integer i = 0; i < valueOf(MotifLength); i = i + 1 ) begin
			total[i] = base[i][0] + base[i][1] + base[i][2] + base[i][3] + base[i][4] +
				   base[i][5] + base[i][6] + base[i][7] + base[i][8] + base[i][9] +
				   base[i][10] + base[i][11] + base[i][12] + base[i][13] + base[i][14] +
				   base[i][15] + base[i][16] + base[i][17] + base[i][18] + base[i][19];
			typeConverterTotal[i].enq(total[i]);
			for ( Integer j = 0; j < 20; j = j + 1 ) begin
				typeConverterBase[i][j].enq(base[i][j]);
			end
		end
	endrule
	rule makePssm3; // 29 cycle
		for ( Integer i = 0; i < valueOf(MotifLength); i = i + 1 ) begin
			let total <- typeConverterTotal[i].get;
			Vector#(20, Bit#(32)) base = replicate(0);
			for ( Integer j = 0; j < 20; j = j + 1 ) begin
				base[j] <- typeConverterBase[i][j].get;
				fpDiv[i][j].enq(base[j], total);
			end
		end
	endrule
	rule makePssm4; // 1 cycle
		Vector#(MotifLength, Vector#(20, Bit#(32))) pssm = replicate(replicate(0));
		for ( Integer i = 0; i < valueOf(MotifLength); i = i + 1 ) begin
			for ( Integer j = 0; j < 20; j = j + 1 ) begin
				fpDiv[i][j].deq;
				pssm[i][j] = fpDiv[i][j].first;
			end
		end
		pssmQ.enq(pssm);
		$write("\033[1;33mCycle %1d -> \033[1;33m[PssmMaker]: \033[0m: PssmMaker finished!\n", cycleCount);
	endrule


	method Action putMotif(Bit#(MotifRelayLength) m);
		motifQ.enq(m);
	endmethod
	method Action putStartPe(Bit#(32) s);
		startPeQ.enq(s);
	endmethod
	method ActionValue#(Vector#(MotifLength, Vector#(20, Bit#(32)))) get;
		pssmQ.deq;
		return pssmQ.first;
	endmethod
endmodule
