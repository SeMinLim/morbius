import FIFO::*;
import FIFOF::*;
import BRAM::*;
import BRAMFIFO::*;
import Clocks::*;
import Vector::*;

import PcieCtrl::*;

import DRAMController::*;
import DRAMArbiter::*;

import Serializer::*;

import GibbsSampler::*;


// Dataset 1
// DMA
typedef 12 DmaReadSeqUnit;
typedef TMul#(DmaReadSeqUnit, 5) DmaReadSeqUnitSize;
typedef 8 DmaReadMtfUnit;
typedef TMul#(DmaReadMtfUnit, 5) DmaReadMtfUnitSize;
typedef 800 DmaReadSeqSize;
typedef 64 DmaReadMtfSize;
typedef TAdd#(DmaReadSeqSize, DmaReadMtfSize) DmaReadSize;
// Sequences
typedef 32768 SeqNum;
typedef 300 SeqLength;
typedef TMul#(SeqLength, 5) SeqSize;
typedef TMul#(SeqNum, 1536) DataSeqSize;
typedef TDiv#(SeqLength, DmaReadSeqUnit) ForOneSeq;
// Motifs
typedef 16 MotifLength;
typedef TMul#(MotifLength, 5) MotifSize;
typedef 16 PeNumMotif;
typedef TMul#(MotifSize, PeNumMotif) MotifRelayLength;
typedef TDiv#(SeqNum, PeNumMotif) MotifRelaySize;
typedef TMul#(MotifRelaySize, 1536) DataMotifSize;
typedef TDiv#(MotifRelayLength, DmaReadMtfUnitSize) ForOneMotifRelay;
typedef TMul#(MotifRelaySize, 3) ForDramFitterMtf;
// DRAM Write
typedef TDiv#(DataSeqSize, 512) DramWriteSeqSize;
typedef TDiv#(DataMotifSize, 512) DramWriteMtfSize;
// DRAM Read
typedef 3 DramReadSeqSize;
typedef TDiv#(DataMotifSize, 512) DramReadMtfSize;
typedef 0 DramSeqAddrStart;
typedef TDiv#(DataSeqSize, 512) DramMtfAddrStart;


interface HwMainIfc;
endinterface
module mkHwMain#(PcieUserIfc pcie, DRAMUserIfc dram) 
	(HwMainIfc);
	Clock curClk <- exposeCurrentClock;
	Reset curRst <- exposeCurrentReset;

	Clock pcieclk = pcie.user_clk;
	Reset pcierst = pcie.user_rst;

	// Cycle Counter
	Reg#(Bit#(32)) cycleCount <- mkReg(0);
	rule incCycleCounter;
		cycleCount <= cycleCount + 1;
	endrule
	
	// Serializer & DeSerializer
	// For both
	SerializerIfc#(128, 2) serializer128b64b <- mkSerializer;
	// For sequences
	DeSerializerIfc#(DmaReadSeqUnitSize, ForOneSeq) deserializer60b1500b <- mkDeSerializer;
	SerializerIfc#(1500, 3) serializer1500b500b <- mkSerializer;
	DeSerializerIfc#(500, 3) deserializer500b1500b <- mkDeSerializer;
	// For motifs
	DeSerializerIfc#(DmaReadMtfUnitSize, ForOneMotifRelay) deserializer40b1280b <- mkDeSerializer;
	DeSerializerIfc#(512, 3) deserializer512b1536b <- mkDeSerializer;

	// DRAMArbiter
	DRAMArbiterIfc#(2) dramArbiter <- mkDRAMArbiter(dram);

	// Rule On/Off
	Reg#(Bool) preprocessorSeqOn <- mkReg(True);
	Reg#(Bool) preprocessorMtfOn <- mkReg(False);
	Reg#(Bool) dramWriteSeqOn <- mkReg(False);
	Reg#(Bool) dramWriteMtfOn <- mkReg(False);
	Reg#(Bool) dramReadSeqOn <- mkReg(True);
	Reg#(Bool) dramReadMtfOn <- mkReg(False);

	// GibbsSampler
	GibbsSamplerIfc gibbsSampler <- mkGibbsSampler;

	// I/O
	FIFO#(Bit#(SeqSize)) sequenceQ <- mkFIFO;
	FIFO#(Bit#(MotifRelayLength)) motifQ <- mkSizedBRAMFIFO(valueOf(MotifRelaySize));
	//--------------------------------------------------------------------------------------------
	// Get commands from the host via PCIe
	//--------------------------------------------------------------------------------------------
	FIFOF#(Bit#(32)) dmaReadReqQ <- mkFIFOF;
	rule pcieDataWriter;
		let w <- pcie.dataReceive;

		let d = w.data;
		let a = w.addr;
		let off = (a >> 2);

		if ( off == 0 ) begin
			dmaReadReqQ.enq(d);
		end
	endrule
	FIFOF#(Bit#(32)) statusQ <- mkFIFOF;
	rule pcieDataReader;
		let r <- pcie.dataReq;

		Bit#(4) a = truncate(r.addr>>2);
		if ( a == 0 ) begin
			if ( statusQ.notEmpty ) begin
				let s = statusQ.first;
				statusQ.deq;
				pcie.dataSend(r, s);
			end else begin 
				pcie.dataSend(r, 7777);
			end
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Read DMA buffer
	//--------------------------------------------------------------------------------------------
	Reg#(Bit#(32)) dmaReadCnt1 <- mkReg(0);
	Reg#(Bit#(32)) dmaReadCnt2 <- mkReg(0);
	rule dmaReadRequester( dmaReadReqQ.notEmpty && (dmaReadCnt1 == 0) );
		dmaReadReqQ.deq;
		let d = dmaReadReqQ.first;
		pcie.dmaReadReq(0, truncate(d));
		dmaReadCnt1 <= d;
		if ( dmaReadCnt2 == 0 ) begin
			$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DMA read started!\n", cycleCount);
		end
	endrule
	rule dmaReader( dmaReadCnt1 != 0 );
		let d <- pcie.dmaReadWord;
		serializer128b64b.put(d);
		if ( dmaReadCnt1 == 1 ) begin
			if ( dmaReadCnt2 + 1 == fromInteger(valueOf(DmaReadSize)) ) begin
				statusQ.enq(1);
				dmaReadCnt2 <= 0;
				$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DMA read finished!\n", cycleCount);
			end else begin
				dmaReadCnt2 <= dmaReadCnt2 + 1;
			end
		end
		dmaReadCnt1 <= dmaReadCnt1 - 1;
	endrule
	//--------------------------------------------------------------------------------------------
	// Truncate zero-padding & Put values to DeSerializer to store to DRAM
	//--------------------------------------------------------------------------------------------
	Reg#(Bit#(32)) preprocessorSeqCnt1 <- mkReg(0);
	Reg#(Bit#(32)) preprocessorSeqCnt2 <- mkReg(0);
	rule preprocessorSeq( preprocessorSeqOn );
		let v <- serializer128b64b.get;
		deserializer60b1500b.put(truncate(v));
		if ( preprocessorSeqCnt1 + 1 == fromInteger(valueOf(ForOneSeq)) ) begin
			if ( preprocessorSeqCnt2 + 1 == fromInteger(valueOf(SeqNum)) ) begin
				preprocessorSeqCnt1 <= 0;
				preprocessorSeqCnt2 <= 0;
				preprocessorSeqOn <= False;
				preprocessorMtfOn <= True;
			end else begin
				preprocessorSeqCnt1 <= 0;
				preprocessorSeqCnt2 <= preprocessorSeqCnt2 + 1;
			end
		end else begin
			preprocessorSeqCnt1 <= preprocessorSeqCnt1 + 1;
		end
	endrule
	Reg#(Bit#(32)) dramFitterSeqCnt <- mkReg(0);
	rule dramFitterSeq;
		let v <- deserializer60b1500b.get;
		serializer1500b500b.put(v);
		if ( dramFitterSeqCnt == 0 ) begin
			dramWriteSeqOn <= True;
			dramFitterSeqCnt <= dramFitterSeqCnt + 1;
		end else begin
			if ( dramFitterSeqCnt + 1 == fromInteger(valueOf(SeqNum)) ) begin
				dramFitterSeqCnt <= 0;
			end else begin
				dramFitterSeqCnt <= dramFitterSeqCnt + 1;
			end
		end
	endrule
	Reg#(Bit#(32)) preprocessorMtfCnt1 <- mkReg(0);
	Reg#(Bit#(32)) preprocessorMtfCnt2 <- mkReg(0);
	rule preprocessorMtf( preprocessorMtfOn );
		let v <- serializer128b64b.get;
		deserializer40b1280b.put(truncate(v));
		if ( preprocessorMtfCnt1 + 1 == fromInteger(valueOf(ForOneMotifRelay)) ) begin
			if ( preprocessorMtfCnt2 + 1 == fromInteger(valueOf(MotifRelaySize)) ) begin
				preprocessorMtfCnt1 <= 0;
				preprocessorMtfCnt2 <= 0;
				preprocessorSeqOn <= True;
				preprocessorMtfOn <= False;
			end else begin
				preprocessorMtfCnt1 <= 0;
				preprocessorMtfCnt2 <= preprocessorMtfCnt2 + 1;
			end
		end else begin
			preprocessorMtfCnt1 <= preprocessorMtfCnt1 + 1;
		end
	endrule
	FIFO#(Bit#(512)) dramFitterMtfQ <- mkSizedBRAMFIFO(valueOf(ForDramFitterMtf));
	Reg#(Bit#(MotifRelayLength)) dramFitterMtfR <- mkReg(0);
	Reg#(Bit#(32)) dramFitterMtfCnt <- mkReg(0);
	rule dramFitterMtf;
		if ( dramFitterMtfCnt == 0 ) begin
			let v <- deserializer40b1280b.get;
			dramFitterMtfQ.enq(truncate(v));
			dramFitterMtfR <= (v >> 512);
			dramFitterMtfCnt <= dramFitterMtfCnt + 1;
		end else begin
			let v = dramFitterMtfR;
			dramFitterMtfQ.enq(truncate(v));
			if ( dramFitterMtfCnt + 1 == 3 ) begin
				dramFitterMtfR <= 0;
				dramFitterMtfCnt <= 0;
			end else begin
				dramFitterMtfR <= (v >> 512);
				dramFitterMtfCnt <= dramFitterMtfCnt + 1;
			end
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Store the sequences and motifs to DRAM
	//--------------------------------------------------------------------------------------------
	Reg#(Bit#(32)) dramWriteSeqCnt <- mkReg(0);
	rule dramWriterSeq( dramWriteSeqOn );
		if ( dramWriteSeqCnt == 0 ) begin
			dramArbiter.users[0].cmd(fromInteger(valueOf(DramSeqAddrStart)), fromInteger(valueOf(DramWriteSeqSize)), True);
			dramWriteSeqCnt <= dramWriteSeqCnt + 1;
			$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM write started! [SEQ]\n", cycleCount);
		end else begin
			let d <- serializer1500b500b.get;
			dramArbiter.users[0].write(zeroExtend(d));
			if ( dramWriteSeqCnt == fromInteger(valueOf(DramWriteSeqSize)) ) begin
				dramWriteSeqCnt <= 0;
				dramWriteSeqOn <= False;
				dramWriteMtfOn <= True;
				$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM write finished! [SEQ]\n", cycleCount);
			end else begin
				dramWriteSeqCnt <= dramWriteSeqCnt + 1;
			end
		end
	endrule
	Reg#(Bit#(32)) dramWriteMtfCnt <- mkReg(0);
	rule dramWriterMtf( dramWriteMtfOn );
		if ( dramWriteMtfCnt == 0 ) begin
			dramArbiter.users[0].cmd(fromInteger(valueOf(DramMtfAddrStart)), fromInteger(valueOf(DramWriteMtfSize)), True);
			dramWriteMtfCnt <= dramWriteMtfCnt + 1;
			$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM write started! [MTF]\n", cycleCount);
		end else begin
			dramFitterMtfQ.deq;
			let d = dramFitterMtfQ.first;
			dramArbiter.users[0].write(d);
			if ( dramWriteMtfCnt == fromInteger(valueOf(DramWriteMtfSize)) ) begin
				dramWriteMtfCnt <= 0;
				dramWriteMtfOn <= False;
				dramReadSeqOn <= True;
				$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM write finished! [MTF]\n", cycleCount);
			end else begin
				dramWriteMtfCnt <= dramWriteMtfCnt + 1;
			end
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Read one sequence and the motifs from DRAM 
	// Send 500bits values to DeSerializer
	//--------------------------------------------------------------------------------------------
	Reg#(Bit#(32)) dramReadSeqCnt <- mkReg(0);
	rule dramReaderSeq( dramReadSeqOn );
		if ( dramReadSeqCnt == 0 ) begin
			dramArbiter.users[0].cmd(fromInteger(valueOf(DramSeqAddrStart)), fromInteger(valueOf(DramReadSeqSize)), False);
			dramReadSeqCnt <= dramReadSeqCnt + 1;
			$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM read started! [SEQ]\n", cycleCount);
		end else begin
			let d <- dramArbiter.users[0].read;
			deserializer500b1500b.put(truncate(d));
			if ( dramReadSeqCnt == fromInteger(valueOf(DramReadSeqSize)) ) begin
				dramReadSeqCnt <= 0;
				dramReadSeqOn <= False;
				dramReadMtfOn <= True;
				$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM read finished! [SEQ]\n", cycleCount);
			end else begin
				dramReadSeqCnt <= dramReadSeqCnt + 1;
			end
		end
	endrule
	Reg#(Bit#(32)) dramReadMtfCnt <- mkReg(0);
	rule dramReaderMtf( dramReadMtfOn );
		if ( dramReadMtfCnt == 0 ) begin
			dramArbiter.users[0].cmd(fromInteger(valueOf(DramMtfAddrStart)), fromInteger(valueOf(DramReadMtfSize)), False);
			dramReadMtfCnt <= dramReadMtfCnt + 1;
			$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM read started! [MTF]\n", cycleCount);
		end else begin
			let d <- dramArbiter.users[0].read;
			deserializer512b1536b.put(d);
			if ( dramReadMtfCnt == fromInteger(valueOf(DramReadMtfSize)) ) begin
				dramReadMtfCnt <= 0;
				dramReadMtfOn <= False;
				dramReadSeqOn <= True;
				$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM read finished! [Mtf]\n", cycleCount);
			end else begin
				dramReadMtfCnt <= dramReadMtfCnt + 1;
			end
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Relay a sequence and motifs to GibbsSampler
	//--------------------------------------------------------------------------------------------
	rule relaySequence;
		let s <- deserializer500b1500b.get;
		gibbsSampler.putSequence(s);
	endrule
	rule relayMotifs;
		let m <- deserializer512b1536b.get;
		gibbsSampler.putMotif(truncate(m));
	endrule
endmodule
