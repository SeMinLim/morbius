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
typedef 256 SeqLength;
typedef TMul#(SeqLength, 5) SeqSize;
typedef	TMul#(512, 3) SeqStoredSize;
typedef TMul#(SeqNum, SeqStoredSize) DataSeqSize;
// Motifs
typedef 16 MotifLength;
typedef TMul#(MotifLength, 5) MotifSize;
typedef 16 PeNumMotif;
typedef TMul#(MotifSize, PeNumMotif) MotifRelayLength;
typedef TDiv#(SeqNum, PeNumMotif) MotifRelaySize;
typedef TMul#(512, 3) MotifRelayStoredSize;
typedef TMul#(MotifRelaySize, MotifRelayStoredSize) DataMotifSize;
// DMA
typedef 640 DmaReadSeqSize;
typedef 40 DmaReadMtfSize;
typedef TAdd#(DmaReadSeqSize, DmaReadMtfSize) DmaReadSize;
// DRAM Write
typedef TDiv#(DataSeqSize, 512) DramWriteSeqSize;
typedef TDiv#(DataMotifSize, 512) DramWriteMtfSize;
// DRAM Read
typedef TDiv#(SeqStoredSize, 512) DramReadSeqSize;
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
	DeSerializerIfc#(128, 10) deserializer128b1280b <- mkDeSerializer;
	DeSerializerIfc#(512, 3) deserializer512b1536b <- mkDeSerializer;

	// DRAMArbiter
	DRAMArbiterIfc#(2) dramArbiter <- mkDRAMArbiter(dram);

	// Rule On/Off
	Reg#(Bool) dramWriteOn <- mkReg(False);
	Reg#(Bool) dramReadSeqOn <- mkReg(False);
	Reg#(Bool) dramReadMtfOn <- mkReg(False);
	Reg#(Bool) relaySeqOn <- mkReg(True);
	Reg#(Bool) relayMtfOn <- mkReg(False);

	// GibbsSampler
	GibbsSamplerIfc gibbsSampler <- mkGibbsSampler;
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
		deserializer128b1280b.put(d);
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
		dramWriteOn <= True;
	endrule
	//--------------------------------------------------------------------------------------------
	// Store the sequences and motifs to DRAM
	//--------------------------------------------------------------------------------------------
	Reg#(Bit#(SeqSize)) dramWriteR <- mkReg(0);
	Reg#(Bit#(32)) dramWriteCnt1 <- mkReg(0);
	Reg#(Bit#(32)) dramWriteCnt2 <- mkReg(0);
	rule dramWriter( dramWriteOn );
		if ( dramWriteCnt1 == 0 ) begin
			dramArbiter.users[0].cmd(0, fromInteger(valueOf(DramWriteSize)), True);
			dramWriteCnt1 <= dramWriteCnt1 + 1;
			$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM write started!\n", cycleCount);
		end else begin
			if ( dramWriteCnt2 == 0 ) begin
				let d <- deserializer128b1280b.get;
				dramArbiter.users[0].write(truncate(d));
				dramWriteR <= (d >> 512);
				dramWriteCnt1 <= dramWriteCnt1 + 1;
				dramWriteCnt2 <= dramWriteCnt2 + 1;
			end else if ( dramWriteCnt2 == 1 ) begin
				let d = dramWriteR;
				dramArbiter.users[0].write(truncate(d));
				dramWriteR <= (d >> 512);
				dramWriteCnt1 <= dramWriteCnt1 + 1;
				dramWriteCnt2 <= dramWriteCnt2 + 1;
			end else if ( dramWriteCnt2 == 2 ) begin
				let d = dramWriteR;
				dramArbiter.users[0].write(truncate(d));
				if ( dramWriteCnt1 == fromInteger(valueOf(DramWriteSize)) ) begin
					dramWriteCnt1 <= 0;
					dramWriteOn <= False;
					dramReadSeqOn <= True;
					$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM write finished!\n", cycleCount);
				end else begin
					dramWriteCnt1 <= dramWriteCnt1 + 1;
				end
				dramWriteCnt2 <= 0;
			end
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Read one sequence and the motifs from DRAM 
	//--------------------------------------------------------------------------------------------
	Reg#(Bit#(32)) dramReadSeqCnt <- mkReg(0);
	rule dramReaderSeq( dramReadSeqOn );
		if ( dramReadSeqCnt == 0 ) begin
			dramArbiter.users[0].cmd(fromInteger(valueOf(DramSeqAddrStart)), fromInteger(valueOf(DramReadSeqSize)), False);
			dramReadSeqCnt <= dramReadSeqCnt + 1;
			$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM read started! [SEQ]\n", cycleCount);
		end else begin
			let d <- dramArbiter.users[0].read;
			deserializer512b1536b.put(d);
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
				$write("\033[1;33mCycle %1d -> \033[1;33m[HwMain]: \033[0m: DRAM read finished! [Mtf]\n", cycleCount);
			end else begin
				dramReadMtfCnt <= dramReadMtfCnt + 1;
			end
		end
	endrule
	//--------------------------------------------------------------------------------------------
	// Relay a sequence and motifs to GibbsSampler
	//--------------------------------------------------------------------------------------------
	rule relaySeq( relaySeqOn );
		let s <- deserializer512b1536b.get;
		gibbsSampler.putSequence(truncate(s));
		relaySeqOn <= False;
		relayMtfOn <= True;
	endrule
	Reg#(Bit#(32)) relayMtfCnt <- mkReg(0);
	rule relayMtf( relayMtfOn );
		let m <- deserializer512b1536b.get;
		gibbsSampler.putMotif(truncate(m));
		if ( relayMtfCnt + 1 == fromInteger(valueOf(MotifRelaySize)) ) begin
			relayMtfCnt <= 0;
			relayMtfOn <= False;
		end else begin
			relayMtfCnt <= relayMtfCnt + 1;
		end
	endrule
endmodule
