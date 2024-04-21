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
// Sequences
typedef 32768 SeqNum;
typedef 1000 SeqLength;
typedef 2048 SeqStoredSize;
typedef TMul#(SeqLength, 2) SeqSize;
typedef TMul#(SeqNum, SeqStoredSize) DataSize;
// Motifs
typedef 16 MotifLength;
typedef TMul#(MotifLength, 2) MotifSize;
typedef 16 PeNumMotif;
typedef TMul#(SeqNum, MotifSize) DataMotifSize;
typedef TMul#(MotifSize, PeNumMotif) MotifRelayLength;
typedef TDiv#(SeqNum, PeNumMotif) MotifRelaySize;
// DMA
typedef 1024 DmaReadSeqSize;
typedef 16 DmaReadMtfSize;
typedef TAdd#(DmaReadSeqSize, DmaReadMtfSize) DmaReadSize;
// DRAM Write
typedef TDiv#(DataSize, 512) DramWriteSeqSize;
typedef TDiv#(DataMotifSize, 512) DramWriteMtfSize;
typedef TAdd#(DramWriteSeqSize, DramWriteMtfSize) DramWriteSize;
// DRAM Read
typedef TDiv#(SeqStoredSize, 512) DramReadSeqSize;
typedef TDiv#(DataMotifSize, 512) DramReadMtfSize;
typedef 0 DramSeqAddrStart;
typedef TDiv#(DataSize, 512) DramMtfAddrStart;


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
	DeSerializerIfc#(128, 4) deserializer128b512b <- mkDeSerializer;
	DeSerializerIfc#(512, 4) deserializer512b2048b <- mkDeSerializer;

	// DRAMArbiter
	DRAMArbiterIfc#(2) dramArbiter <- mkDRAMArbiter(dram);

	// Rule On/Off
	Reg#(Bool) dramWriteOn <- mkReg(False);
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
		deserializer128b512b.put(d);
		dmaReadCnt1 <= dmaReadCnt1 - 1;
		if ( dmaReadCnt1 == 1 ) begin
			if ( dmaReadCnt2 + 1 == fromInteger(valueOf(DmaReadSize)) ) begin
				statusQ.enq(1);
				dmaReadCnt2 <= 0;
				$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DMA read finished!\n", cycleCount);
			end else begin
				dmaReadCnt2 <= dmaReadCnt2 + 1;
			end
		end
		dramWriteOn <= True;
	endrule
	//--------------------------------------------------------------------------------------------
	// Store the sequences and motifs to DRAM
	//--------------------------------------------------------------------------------------------
	Reg#(Bit#(32)) dramWriteCnt <- mkReg(0);
	rule dramWriter( dramWriteOn );
		if ( dramWriteCnt == 0 ) begin
			dramArbiter.users[0].cmd(0, fromInteger(valueOf(DramWriteSize)), True);
			dramWriteCnt <= dramWriteCnt + 1;
			$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: Burst DRAM write started!\n", cycleCount);
		end else begin
			let d <- deserializer128b512b.get;
			dramArbiter.users[0].write(d);
			if ( dramWriteCnt == fromInteger(valueOf(DramWriteSize)) ) begin
				dramWriteCnt <= 0;
				dramWriteOn <= False;
				dramReadSeqOn <= True;
				$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: Burst DRAM write finished!\n", cycleCount);
			end else begin
				dramWriteCnt <= dramWriteCnt + 1;
			end
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Read one sequence and the motifs from DRAM 
	// Send 512bits values to DeSerializer
	//--------------------------------------------------------------------------------------------
	Reg#(Bit#(32)) dramReadSeqCnt <- mkReg(0);
	rule dramReaderSeq( dramReadSeqOn );
		if ( dramReadSeqCnt == 0 ) begin
			dramArbiter.users[0].cmd(fromInteger(valueOf(DramSeqAddrStart)), fromInteger(valueOf(DramReadSeqSize)), False);
			dramReadSeqCnt <= dramReadSeqCnt + 1;
			$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: Burst DRAM read started! [SEQ]\n", cycleCount);
		end else begin
			let d <- dramArbiter.users[0].read;
			deserializer512b2048b.put(d);
			if ( dramReadSeqCnt == fromInteger(valueOf(DramReadSeqSize)) ) begin
				dramReadSeqCnt <= 0;
				dramReadSeqOn <= False;
				dramReadMtfOn <= True;
				$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: Burst DRAM read finished! [SEQ]\n", cycleCount);
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
			$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: Burst DRAM read started! [MTF]\n", cycleCount);
		end else begin
			let d <- dramArbiter.users[0].read;
			gibbsSampler.putMotif(d);
			if ( dramReadMtfCnt == fromInteger(valueOf(DramReadMtfSize)) ) begin
				dramReadMtfCnt <= 0;
				dramReadMtfOn <= False;
				dramReadSeqOn <= True;
				$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: Burst DRAM read finished! [Mtf]\n", cycleCount);
			end else begin
				dramReadMtfCnt <= dramReadMtfCnt + 1;
			end
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Relay a sequence to GibbsSampler
	//--------------------------------------------------------------------------------------------
	rule relaySequence;
		let s <- deserializer512b2048b.get;
		gibbsSampler.putSequence(truncate(s));
	endrule
endmodule
