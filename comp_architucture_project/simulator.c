#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#define CMD_LEN (50)
#define Memory_Size (1048576)
#define num_of_cores 4


int Memory[Memory_Size]; ///data_memory///

//// ****   Linked List to help us to make bus transactions done in round robin order  **** ////

typedef struct node {
	int val;
	struct node* next;
} node_t;

void push(node_t* head, int val) {  // adding the core who asked for bus transaction while bus is busy to the end of queue
	if (head == NULL) {
		head->val = val;
		head->next = NULL;
	}
	else {
		node_t* current = head;
		while (current->next != NULL) {
			current = current->next;
		}
		current->next = (node_t*)malloc(sizeof(node_t));
		current->next->val = val;
		current->next->next = NULL;
	}
}

void pop(node_t* head) { // popping out the head of queue witch is the core thats just finish his bus transaction
	if (head == NULL) {
		return;
	}
	else if (head->next == NULL) {
		free(head);
	}
	node_t* prev_head = head;
	head = prev_head->next;
	head->val = prev_head->next->val;
	head->next = prev_head->next->next;
	free(prev_head);
}

node_t* queue = NULL;  // global queue so all cores can get to it



////////**   These structures we need to simulate the multi-core proccesor   **////////
typedef struct IF_ID { ///register IF/ID
	int IR; //the insruction
	int NXT_PC; //next pc
	int stall;
}IF_ID;

typedef struct ID_EX { // ID/EX register //
	int NXT_PC; // next pc
	int IR;
	int rs;
	int rt;
	int rd;
	int dataRs;
	int dataRt;
	int IMM; // immediate
	int stall;
}ID_EX;

typedef struct EX_MEM { // EX/MEM register //
	int ALU;
	int rt;
	int rd;
	int IR;
	int dataRt;
	int stall;
	int NXT_PC;
}EX_MEM;

typedef struct MEM_WB { // MEM/WB register //
	int ALU;
	int MD;
	int rd;
	int IR;
	int stall;
	int NXT_PC;
}MEM_WB;

typedef struct Core {     // The Core structure //
	int DSRAM_Cache[64][4];
	int TSRAM_Cache[64];
	int CMD_Memory[1024];
	int regs[16];
	int pc;
	int fetch;
	int decode;
	int WaitingForMemory;
	int index;
	int num_cycles; 		  // number of cycles
	int num_instructions;     // number of instructions
	int num_read_hit;         // number of read hit
	int num_write_hit;	      // number of write hit
	int num_read_miss;	      // number of read miss
	int num_write_miss;	      // number of write miss
	int num_decode_stall;	  // number of decode stalls
	int num_mem_stall;	      // number of memory stalls
	int finish;
	IF_ID* if_id;
	ID_EX* id_ex;
	EX_MEM* ex_mem;
	MEM_WB* mem_wb;
} Core;

typedef struct reg {    // This structure include the informtion about the four Cores and the BUS //
	Core* core[4];
	int bus_origid;
	int bus_cmd;
	int bus_addr;
	int bus_data;
	int bus_shared;
	int Tran_clk; //Transaction time
	int clk;
	int bus_is_busy; // busy==1 ---> there is a core that use the bus
}Reg;

//// **** Here you can find the copy finctions and they are used to check step by step if the code works properly each clock cycle **** ////

void If_Id_copy(IF_ID* old, IF_ID* new) {  // copy IF/ID Register //
	int i = 0;
	i = new->IR;
	old->IR = i;
	i = new->NXT_PC;
	old->NXT_PC = i;
	i = new->stall;
	old->stall = i;
}

void Id_Ex_copy(ID_EX* old, ID_EX* new) {  // copy ID/EX Register //
	int i = 0;
	i = new->IR;
	old->IR = i;
	i = new->NXT_PC;
	old->NXT_PC = i;
	i = new->rs;
	old->rs = i;
	i = new->rt;
	old->rt = i;
	i = new->rd;
	old->rd = i;
	i = new->dataRs;
	old->dataRs = i;
	i = new->dataRt;
	old->dataRs = i;
	i = new->IMM;
	old->IMM = i;
	i = new->stall;
	old->stall = i;
}

void Ex_Mem_copy(EX_MEM* old, EX_MEM* new) {  // copy EX/MEM Register //
	int i = 0;
	i = new->IR;
	old->IR = i;
	i = new->NXT_PC;
	old->NXT_PC = i;
	i = new->rt;
	old->rt = i;
	i = new->rd;
	old->rd = i;
	i = new->dataRt;
	old->dataRt = i;
	i = new->ALU;
	old->ALU = i;
	i = new->stall;
	old->stall = i;
}

void Mem_Wb_copy(MEM_WB* old, MEM_WB* new) {  // copy MEM/WB Register //
	int i = 0;
	i = new->IR;
	old->IR = i;
	i = new->NXT_PC;
	old->NXT_PC = i;
	i = new->rd;
	old->rd = i;
	i = new->ALU;
	old->ALU = i;
	i = new->MD;
	old->MD = i;
	i = new->stall;
	old->stall = i;
}

void copyCore(Core* oldCore, Core* newCore) {  // Core copy //
	int i, j, k = 0;
	for (j = 0; j < 16; j++) {
		i = newCore->regs[j];
		oldCore->regs[j] = i;
	}
	for (j = 0; j < 64; j++) {
		for (k = 0; k < 4; k++) {
			i = oldCore->DSRAM_Cache[j][k];
			newCore->DSRAM_Cache[j][k] = i;
		}
	}
	for (j = 0; j < 64; j++) {
		i = newCore->TSRAM_Cache[j];
		oldCore->TSRAM_Cache[j] = i;
	}
	for (j = 0; j < 1024; j++) {
		i = newCore->CMD_Memory[j];
		oldCore->CMD_Memory[j] = i;
	}
	i = newCore->pc;
	oldCore->pc = i;
	i = newCore->fetch;
	oldCore->fetch = i;
	i = newCore->decode;
	oldCore->decode = i;
	i = newCore->finish;
	oldCore->finish = i;
	If_Id_copy(oldCore->if_id, newCore->if_id);
	Id_Ex_copy(oldCore->id_ex, newCore->id_ex);
	Ex_Mem_copy(oldCore->ex_mem, newCore->ex_mem);
	Mem_Wb_copy(oldCore->mem_wb, newCore->mem_wb);
}

void Regs_copy(Reg* oldReg, Reg* newReg) {  // Regs copy  //
	int i = 0;
	i = newReg->bus_addr;
	oldReg->bus_addr = i;
	i = newReg->bus_cmd;
	oldReg->bus_cmd = i;
	i = newReg->bus_data;
	oldReg->bus_data = i;
	i = newReg->bus_is_busy;
	oldReg->bus_is_busy;
	i = newReg->bus_origid;
	oldReg->bus_origid = i;
	i = newReg->clk;
	oldReg->clk = i;
	i = newReg->Tran_clk;
	oldReg->Tran_clk = i;
	for (i = 0; i < num_of_cores; i++) {
		copyCore(oldReg->core[i], newReg->core[i]);
	}
}

//// **** In this part we have the create and initialize functions **** ////

IF_ID* create_IF_ID() {
	IF_ID *if_id = (IF_ID*)malloc(sizeof(IF_ID));
	if_id->IR = 0;
	if_id->NXT_PC = 0;
	if_id->stall = 1;
	return if_id;
}

ID_EX* create_ID_EX() {
	ID_EX *id_ex = (ID_EX*)malloc(sizeof(ID_EX));
	id_ex->IR = 0;
	id_ex->NXT_PC = 0;
	id_ex->IMM = 0;
	id_ex->rd = -1;
	id_ex->rs = -1;
	id_ex->rt = -1;
	id_ex->dataRs = 0;
	id_ex->dataRt = 0;
	id_ex->stall = 1;
	return id_ex;
}

EX_MEM* create_EX_MEM() {
	EX_MEM *ex_mem = (EX_MEM*)malloc(sizeof(EX_MEM));
	ex_mem->IR = 0;
	ex_mem->ALU = 0;
	ex_mem->rt = -2;
	ex_mem->rd = -2;
	ex_mem->dataRt = 0;
	ex_mem->stall = 1;
	ex_mem->NXT_PC = 0;
	return ex_mem;
}

MEM_WB* create_MEM_WB() {
	MEM_WB *mem_wb = (MEM_WB*)malloc(sizeof(MEM_WB));
	mem_wb->IR = 0;
	mem_wb->ALU = 0;
	mem_wb->MD = 0;
	mem_wb->rd = -3;
	mem_wb->stall = 1;
	mem_wb->NXT_PC = 0;
	return mem_wb;
}

void get_imem(Core* core, FILE* imem) {   //copy the file imem to cmd_memory array  //
	char cmdLine[CMD_LEN];
	int i = 0;
	while (!feof(imem)) {
		fgets(cmdLine, CMD_LEN, imem);
		core->CMD_Memory[i++] = strtol(cmdLine, NULL, 16);
	}
}

Core* create_core(FILE *file){
	MEM_WB *mem_wb = create_MEM_WB();
	IF_ID *if_id = create_IF_ID();
	ID_EX *id_ex = create_ID_EX();
	EX_MEM *ex_mem = create_EX_MEM();
	int i, j = 0;
	Core *core = (Core*)malloc(sizeof(Core));
	for (i = 0; i < 16; i++)
		core->regs[i] = 0;
	core->pc = 0;
	for (i = 0; i < 64; i++) {
		for(j = 0; j < 4; j++)
		core->DSRAM_Cache[i][j] = 0;
	}
	for (i = 0; i < 64; i++)
		core->TSRAM_Cache[i] = 0;
	for (i = 0; i < 1024; i++)
		core->CMD_Memory[i] = 0;
	get_imem(core, file);      // init cmd_memory //
	core->fetch = 0;
	core->decode = 0;
	core->finish = 0;
	core->WaitingForMemory = 0;
	core->num_cycles = 0;
	core->num_decode_stall = 0;
	core->num_instructions = 0;
	core->num_mem_stall = 0;
	core->num_read_hit = 0;
	core->num_read_miss = 0;
	core->num_write_hit = 0;
	core->num_write_miss = 0;
	core->if_id = if_id;
	core->id_ex = id_ex;
	core->ex_mem = ex_mem;
	core->mem_wb = mem_wb;
	return core;
}

Reg* create_reg(Core** core) {
	int i;
	Reg *reg = (Reg*)malloc(sizeof(Reg));
	reg->bus_origid = 0;
	reg->bus_cmd = 0;
	reg->bus_addr = 0;
	reg->bus_data = 0;
	reg->Tran_clk = 0;
	reg->bus_is_busy = 0;
	reg->clk = 0;
	for (i = 0; i < num_of_cores; i++) {
		reg->core[i] = core[i];
	}
	return reg;
}

void get_memory(FILE* memin) {    //  copy the memin file to memory array  //
	char cmdLine[CMD_LEN];
	int i = 0;
	while (!feof(memin)) {
		fgets(cmdLine, CMD_LEN, memin);
		Memory[i++] = strtol(cmdLine, NULL, 16);
	}
}

void openFiles(FILE** files, char *args[]) {
	int i = 0;
	fopen_s(files, args[1], "r");
	fopen_s(files + 1, args[2], "r");
	fopen_s(files + 2, args[3], "r");
	fopen_s(files + 3, args[4], "r");
	fopen_s(files + 4, args[5], "r");
	for (i = 5; i < 27; i++)
		fopen_s(files + i, args[i + 1], "w");
}
void closeFiles(FILE** files) {
	int i;
	for (i = 0; i < 27; i++) {
		if (files[i]) {
			fclose(files[i]);
		}
	}
}

int signExt(int num) {      // sign extension for immediate
	int mask = 0x800;
	num = num & 0xFFF;
	if (mask & num) {
		num += 0xFFFFF000;
	}
	return num;
}

//// ** The following functions implement the operations that are supported by the Cores and should be simulated ** ////

void add(Core* oldCore, Core* newCore, int rs_data, int rt_data) {
	newCore->ex_mem->ALU = rs_data + rt_data;
}

void sub(Core* oldCore, Core* newCore, int rs_data, int rt_data) {
	newCore->ex_mem->ALU = rs_data - rt_data;
}

void and(Core* oldCore, Core* newCore, int rs_data, int rt_data) {
	newCore->ex_mem->ALU = rs_data & rt_data;
}

void or(Core* oldCore, Core* newCore, int rs_data, int rt_data) {
	newCore->ex_mem->ALU = rs_data | rt_data;
}

void xor(Core* oldCore, Core* newCore, int rs_data, int rt_data) {
	newCore->ex_mem->ALU = rs_data ^ rt_data;
}

void mul(Core* oldCore, Core* newCore, int rs_data, int rt_data) {
	newCore->ex_mem->ALU = rs_data * rt_data;
}

void sll(Core* oldCore, Core* newCore, int rs_data, int rt_data) {
	newCore->ex_mem->ALU = rs_data << rt_data;
}

void sra(Core* oldCore, Core* newCore, int rs_data, int rt_data) {
	newCore->ex_mem->ALU = rs_data >> rt_data;
}

void srl(Core* oldCore, Core* newCore, int rs_data, int rt_data) {
	newCore->ex_mem->ALU = (int)((unsigned)rs_data >> rt_data);
}

void beq(Core* oldCore, Core* newCore, int rd, int rs_data, int rt_data) {
	if (rs_data == rt_data) {
		newCore->pc = oldCore->regs[rd] & 0x3FF;
	}	
}

void bne(Core* oldCore, Core* newCore, int rd, int rs_data, int rt_data) {
	if (rs_data != rt_data) {
		newCore->pc = oldCore->regs[rd] & 0x3FF;
	}
}

void blt(Core* oldCore, Core* newCore, int rd, int rs_data, int rt_data) {
	if (rs_data < rt_data) {
		newCore->pc = oldCore->regs[rd] & 0x3FF;
	}
}
void bgt(Core* oldCore, Core* newCore, int rd, int rs_data, int rt_data) {
	if (rs_data > rt_data) {
		newCore->pc = oldCore->regs[rd] & 0x3FF;
	}
}

void ble(Core* oldCore, Core* newCore, int rd, int rs_data, int rt_data) {
	if (rs_data <= rt_data) {
		newCore->pc = oldCore->regs[rd] & 0x3FF;
	}
}

void bge(Core* oldCore, Core* newCore, int rd, int rs_data, int rt_data) {
	if (rs_data >= rt_data) {
		newCore->pc = oldCore->regs[rd] & 0x3FF;
	}
}

void jal(Core* oldCore, Core* newCore, int rd) {
	newCore->regs[15] = oldCore->pc;
	newCore->pc = oldCore->regs[rd] & 0x3FF;
}

void lw(Core* oldCore, Core* newCore, int rs_data, int rt_data) {
	newCore->ex_mem->ALU = rs_data + rt_data; 
}

void sw(Core* oldCore, Core* newCore, int rs_data, int rt_data) {
	newCore->ex_mem->ALU = rs_data + rt_data;
}
//// *** ----------------- Bus and Memory transactions functions ---------------- *** ////

void BusRd_transaction(Reg* oldReg, Reg* newReg, int CoreIndex, int ALU) {
	newReg->bus_cmd = 1;
	newReg->bus_origid = CoreIndex;
	newReg->bus_addr = ALU & 0xFFFFFFFC; // address of the first word in the block
	newReg->bus_data = 0;
	newReg->Tran_clk = oldReg->clk;
}

/// recieving data from another Core-Cache by doing flush on the BUS
void Flush_transaction_BusRd(Reg* oldReg, Reg* newReg, int CoreIndex, int modified_core, int tag, int block_index, int offset_in_block) {
	newReg->bus_cmd = 3;
	newReg->bus_origid = modified_core;
	newReg->bus_addr = oldReg->bus_addr;
	int tmp_data[4] = { 0 };
	for (int i = 0; i < 4; i++)
		tmp_data[i] = oldReg->core[modified_core]->DSRAM_Cache[block_index][i];// copying data from modified core cache
	newReg->bus_data = tmp_data; // putting data on the bus
	for (int i = 0; i < 4; i++) {
		newReg->core[CoreIndex]->DSRAM_Cache[block_index][i] = tmp_data[i]; // updating the data in cache
		if (i == offset_in_block)
			newReg->core[CoreIndex]->mem_wb->MD = tmp_data[offset_in_block]; // updating data in mem_wb register
	}
	newReg->core[CoreIndex]->TSRAM_Cache[block_index] = tag | (1 << 12); // change the state invalid --> shared
	newReg->core[modified_core]->TSRAM_Cache[block_index] = tag | (1 << 12); // change the state modified --> shared
	newReg->bus_is_busy = 0; // bus is free
	newReg->Tran_clk = oldReg->clk;
	Memory[oldReg->bus_addr] = tmp_data; // update data in the main memory 
}

void check_BusRd(Reg* oldReg, Reg* newReg, int CoreIndex, int tag, int block_index, int offset_in_block) {  // when cores see BusRd transaction
	int block_tag = 0;
	int i = 0;
	int state = 0;
	for (i = 0; i < num_of_cores; i++) {
		if (i != CoreIndex){
			block_tag = oldReg->core[i]->TSRAM_Cache[block_index] & 0xFFF;
			state = (oldReg->core[i]->TSRAM_Cache[block_index] & 0x3FFF) >> 12;
			if (state == 3 && block_tag == tag) {
				newReg->bus_shared = 1;
				Flush_transaction_BusRd(oldReg, newReg, CoreIndex, i, tag, block_index, offset_in_block); // the core do Flush on the Bus if condition is True
				// stall is canceld
				newReg->core[CoreIndex]->WaitingForMemory = 0;
			}
			if (state == 2 && block_tag == tag) {
				newReg->core[CoreIndex]->TSRAM_Cache[block_index] = tag | (2 << 12);  // change state Exclusive --> Shared
			}
		}
		
	}	 
}

void Memory_Flush(Reg* oldReg, Reg* newReg, int CoreIndex, int tag, int block_index, int offset_in_block) {
	newReg->bus_cmd = 3;
	newReg->bus_origid = 4;
	newReg->bus_addr = oldReg->bus_addr;
	newReg->bus_data = Memory[oldReg->bus_addr]; // bring the data from the main memory
	newReg->Tran_clk = oldReg->Tran_clk;
	for (int i = 0; i < 4; i++) {
		newReg->core[CoreIndex]->DSRAM_Cache[block_index][i] = Memory[oldReg->bus_addr + i]; // update cache
		if (i == offset_in_block)
			newReg->core[CoreIndex]->mem_wb->MD = Memory[oldReg->bus_addr + i]; /// update the mem-wb register while updating cache
	}
	newReg->bus_is_busy = 0; // bus is free
}

void read_miss_tran(Reg* oldReg, Reg* newReg, Core* newCore, int CoreIndex, int MESI, int** Dsram, int block_tag, int cache_tag, int ALU, int block_index, int offset_in_block) {
	if (oldReg->bus_is_busy == 0 && newReg->bus_is_busy == 0) {
		newCore->num_read_miss++;
		newReg->bus_is_busy = 1;
		if (MESI == 0 || (MESI == 1 || MESI == 2) && block_tag != cache_tag) {
			BusRd_transaction(oldReg, newReg, CoreIndex, ALU);
		}
		else {
			if (MESI == 3 && block_tag != cache_tag) {
				int address = (cache_tag << 8) | (block_index << 2);
				Memory[address] = Dsram[block_index];
				BusRd_transaction(oldReg, newReg, CoreIndex, ALU);
			}
		}
	}
	else {
		if (oldReg->bus_origid == CoreIndex && oldReg->Tran_clk == oldReg->clk - 1)
			check_BusRd(oldReg, newReg, CoreIndex, cache_tag, block_index, offset_in_block); // loking for data in other cores, checking bus shared
		else {
			/* after 16 clock cycles are passed the first word of the block should get to cache */
			if (oldReg->bus_is_busy == 1 && oldReg->bus_origid == CoreIndex && oldReg->clk == oldReg->Tran_clk + 15) {
				Memory_Flush(oldReg, newReg, CoreIndex, block_tag, block_index, offset_in_block); // Memory_flush transaction if data not found in other cores
				newReg->core[CoreIndex]->TSRAM_Cache[block_index] = block_tag | (2 << 12); // convert state ---> Exclusive
				// stall is canceld
				newReg->core[CoreIndex]->WaitingForMemory = 0;
			}
				
		}
	}
}

void BusRdX_transaction(Reg* oldReg, Reg* newReg, int CoreIndex, int ALU) {
	newReg->bus_origid = CoreIndex;
	newReg->bus_cmd = 2;
	newReg->bus_addr = ALU & 0xFFFFFFC; // address of the first word in the block
	newReg->bus_data = 0;
	newReg->bus_is_busy = 1;
	newReg->Tran_clk = oldReg->clk;
}

void Flush_transaction_BusRdX(Reg* oldReg, Reg* newReg, int CoreIndex, int modified_core, int tag, int block_index, int offset_in_block) {
	newReg->bus_cmd = 3;
	newReg->bus_origid = modified_core;
	newReg->bus_addr = oldReg->bus_addr;
	int tmp_data[4] = { 0 };
	for (int i = 0; i < 4; i++)
		tmp_data[i] = oldReg->core[modified_core]->DSRAM_Cache[block_index][i];// copying data from modified core cache
	newReg->bus_data = tmp_data; // putting data on the bus
	for (int i = 0; i < 4; i++) {
		newReg->core[CoreIndex]->DSRAM_Cache[block_index][i] = tmp_data[i]; // updating the data in cache
		if (i == offset_in_block)
			newReg->core[CoreIndex]->mem_wb->MD = tmp_data[offset_in_block]; // updating data in mem_wb register
	}
	newReg->core[modified_core]->TSRAM_Cache[block_index] = tag | (0 << 13); // convert state Modified --> Invalid
	newReg->core[CoreIndex]->TSRAM_Cache[block_index] = tag | (3 << 12); // convert state Invalid --> Modified
	newReg->bus_is_busy = 0; // bus is free
	newReg->Tran_clk = oldReg->clk;
	Memory[oldReg->bus_addr] = tmp_data; // update data in the main memory 
}

void check_BusRdX(Reg* oldReg, Reg* newReg, int CoreIndex, int tag, int block_index, int offset_in_block) {
	int block_tag = 0;
	int i = 0;
	int state = 0;
	for (i = 0; i < num_of_cores; i++) {
		if (i != CoreIndex) {
			block_tag = oldReg->core[i]->TSRAM_Cache[block_index] & 0xFFF;
			state = (oldReg->core[i]->TSRAM_Cache[block_index] & 0x3FFF) >> 12;
			if (state == 3 && block_tag == tag) {
				newReg->bus_shared = 1;
				Flush_transaction_BusRdX(oldReg, newReg, CoreIndex, i, tag, block_index, offset_in_block); // the core do Flush on the Bus if condition is True
				// stall is canceld
				newReg->bus_shared = 0;
				newReg->core[CoreIndex]->WaitingForMemory = 0;
			}
			if (state == 2 && block_tag == tag) {
				newReg->core[CoreIndex]->TSRAM_Cache[block_index] = tag | (1 << 12);  // change state Shared --> Invalid
			}
			if (state == 2 && block_tag == tag) {
				newReg->core[CoreIndex]->TSRAM_Cache[block_index] = tag | (0 << 12);  // change state Exclusive --> Invalid
			}
		}

	}
}

void write_miss_trans(Reg* oldReg, Reg* newReg, Core* newCore, int CoreIndex, int MESI, int** Dsram, int rd, int block_tag, int cache_tag, int ALU, int block_index, int offset_in_block) {
	if (oldReg->bus_is_busy == 0 && newReg->bus_is_busy == 0) {
		if (MESI == 1 && block_tag == cache_tag) { // Shared block is a write miss so we write and do BusRdX on the BUS  
			newCore->num_write_miss++;
			BusRdX_transaction(oldReg, newReg, CoreIndex, ALU);
		}
		else {   // write miss
			newCore->num_write_miss++;
			if (MESI == 3 && block_tag != cache_tag) {
				int address = (cache_tag << 8) | (block_index << 2);
				Memory[address] = Dsram[block_index];
			}
			BusRdX_transaction(oldReg, newReg, CoreIndex, ALU);
		}
	}
	else {
		if (oldReg->bus_origid == CoreIndex && oldReg->Tran_clk == oldReg->clk - 1)
			check_BusRdX(oldReg, newReg, CoreIndex, cache_tag, block_index, offset_in_block);
		else {
			/* after 16 clock cycles are passed the first word of the block should get to cache */
			if (oldReg->bus_is_busy == 1 && oldReg->bus_origid == CoreIndex && oldReg->clk == oldReg->Tran_clk + 15) {
				Memory_Flush(oldReg, newReg, CoreIndex, block_tag, block_index, offset_in_block); // Memory_flush transaction if data not found in other cores
				newReg->core[CoreIndex]->DSRAM_Cache[block_index][offset_in_block] = oldReg->core[CoreIndex]->regs[rd]; // modifieng the value as needed
				newReg->core[CoreIndex]->TSRAM_Cache[block_index] = block_tag | (3 << 12); // convert state ---> Modified
				// stall is canceld
				newReg->core[CoreIndex]->WaitingForMemory = 0;
			}
		}
	}
}

//// *** -----------------------------------Pipline---------------------------------- *** ////

void fetch(Core* oldCore, Core* newCore) {
	if (oldCore->fetch == -1)
		return;
	newCore->if_id->IR = newCore->CMD_Memory[oldCore->pc];
	newCore->if_id->NXT_PC = oldCore->pc + 1;
	if (newCore->if_id->stall == -1)
		newCore->if_id->stall = 1;
	newCore->if_id->stall = 0;
}

void decode(Core* oldCore, Core* newCore) {
	if (oldCore->if_id->stall == 1)
		return;
	if (oldCore->decode == -1) {
		oldCore->if_id->stall = 1;
		return;
	}
	newCore->id_ex->stall = 0;
	newCore->if_id->stall = 1;
	int cmd = oldCore->if_id->IR;
	int imm = cmd & 0xFFF;
	int rt = (cmd >> 12) & 0xF;
	int rs = (cmd >> 16) & 0xF;
	int rd = (cmd >> 20) & 0xF;
	int OpCode = (cmd >> 24) & 0xFF;
	oldCore->regs[1] = imm;
	newCore->regs[0] = 0;
	newCore->id_ex->rt = rt;
	newCore->id_ex->rs = rs;
	newCore->id_ex->rd = rd;
	newCore->id_ex->dataRt = oldCore->regs[rt];
	newCore->id_ex->dataRs = oldCore->regs[rs];
	newCore->id_ex->IR = oldCore->if_id->IR;
	newCore->id_ex->IMM = imm;
	newCore->id_ex->NXT_PC = oldCore->if_id->NXT_PC;
	// branch resolution //
	int rs_data = oldCore->regs[rs];
	int rt_data = oldCore->regs[rt];
	switch (OpCode) {
	case 9:
		beq(oldCore, newCore, rd, rs_data, rt_data);
		break;
	case 10:
		bne(oldCore, newCore, rd, rs_data, rt_data);
		break;
	case 11:
		blt(oldCore, newCore, rd, rs_data, rt_data);
		break;
	case 12:
		bgt(oldCore, newCore, rd, rs_data, rt_data);
		break;
	case 13:
		ble(oldCore, newCore, rd, rs_data, rt_data);
		break;
	case 14:
		bge(oldCore, newCore, rd, rs_data, rt_data);
		break;
	case 15:
		jal(oldCore, newCore, rd);
		break;
	case 20:  // halt
		newCore->fetch = -1;
		break;
	default:
		break;
	}
}

void exec(Core* oldCore, Core* newCore) {
	if (oldCore->id_ex->stall == 1) {
		newCore->id_ex->rd = -2;
		return;
	}
	int cmd = oldCore->id_ex->IR;
	int OpCode = (cmd >> 24) & 0xFF;
	int rs_data = oldCore->id_ex->dataRs;
	int rt_data = oldCore->id_ex->dataRt;
	int rd = oldCore->id_ex->rd;
	int imm = oldCore->id_ex->IMM;
	newCore->regs[1] = imm;
	switch (OpCode)
	{
	case 0:
		add(oldCore, newCore, rs_data, rt_data);
		break;
	case 1:
		sub(oldCore, newCore, rs_data, rt_data);
		break;
	case 2:
		and(oldCore, newCore, rs_data, rt_data);
		break;
	case 3:
		or(oldCore, newCore, rs_data, rt_data);
		break;
	case 4:
		xor(oldCore, newCore, rs_data, rt_data);
		break;
	case 5:
		mul(oldCore, newCore, rs_data, rt_data);
		break;
	case 6:
		sll(oldCore, newCore, rs_data, rt_data);
		break;
	case 7:
		sra(oldCore, newCore, rs_data, rt_data);
		break;
	case 8:
		srl(oldCore, newCore, rs_data, rt_data);
		break;
	case 16:
		lw(oldCore, newCore, rs_data, rt_data);
		break;
	case 17:
		sw(oldCore, newCore, rs_data, rt_data);
		break;
	default:
		break;
	}
	newCore->ex_mem->rt = oldCore->id_ex->rt;
	newCore->ex_mem->rd = rd;
	newCore->ex_mem->IR = cmd;
	newCore->ex_mem->NXT_PC = oldCore->id_ex->NXT_PC;
	newCore->ex_mem->dataRt = oldCore->id_ex->dataRt;
	newCore->ex_mem->stall = 0;
	newCore->id_ex->stall = 1;
}

void mem(Core* oldCore, Core* newCore, int coreIndex, Reg* oldReg, Reg* newReg) {
	if (oldCore->ex_mem->stall == 1) {
		newCore->mem_wb->rd = -3;  // trashing the pipline //
		return;
	}
	newCore->mem_wb->ALU = oldCore->ex_mem->ALU;
	newCore->mem_wb->IR = oldCore->ex_mem->IR;
	newCore->mem_wb->rd = oldCore->ex_mem->rd;
	newCore->mem_wb->NXT_PC = oldCore->mem_wb->NXT_PC;
	int cmd = oldCore->ex_mem->IR;
	int OpCode = (cmd >> 24) & 0xFF;
	int Dsram[64][4];
	int Tsram[64];
	int rd = oldCore->ex_mem->rd;
	int ALU = oldCore->ex_mem->ALU;
	int block_tag = 0;    // block_tag help us to know if a specific block is in the cache or not  //
	int block_index = 0;  // block index help us to know in wich index at the cache the block saved if he is the cache  //
	int offset_in_block = 0;  // each block have 4 words,the offset tells us wich word the block we need, offset = 2 LSB of the address // 
	memcpy(Dsram, oldCore->DSRAM_Cache, 256 * sizeof(int));
	memcpy(Tsram, oldCore->TSRAM_Cache, 64 * sizeof(int));
	if (OpCode == 16 || OpCode == 17) {
		block_tag = ALU >> 8;
		block_index = (ALU & 0xFF) >> 2;
		offset_in_block = ALU & 0x3;
		int cache_tag = Tsram[block_index] & 0xFFF;
		int MESI = (int)((unsigned)Tsram[block_index] >> 12) & 0x3;  // Modified = 3, Exclusive = 2, Shraed = 1, Invalid = 0 //
		if (OpCode == 16) {  // lw  //
			if (MESI != 0 && block_tag == cache_tag) { // --- read hit //
				newCore->mem_wb->MD = Dsram[block_index][offset_in_block];
				newCore->mem_wb->stall = 0;
				newCore->num_read_hit++;
			}
			else {  // read miss --> stall
				/* push coreindex to queue*/
				push(queue, coreIndex);
				while (1) {
					/* if the coreindex is the queue head then do what need to do */
					if (queue->val == coreIndex) {
						newCore->WaitingForMemory = 1;
						read_miss_tran(oldReg, newReg, newCore, coreIndex, MESI, Dsram, block_tag, cache_tag, ALU, block_index, offset_in_block);
						/* pop head out of the queue */
						break;
					}
					pop(queue);
				}
				
				
			}
		}
		else {  // OpCode = 17 --> sw  //
			if ((MESI == 3 || MESI == 2) && block_tag == cache_tag) { // --- write hit, state: modified or exclusive //
				newCore->DSRAM_Cache[block_index][offset_in_block] = oldCore->regs[rd];
				newCore->mem_wb->stall = 0;
				newCore->num_write_hit++;
			}
			else {  // write miss ---> stall
				/* should perform according to the queue */
				push(queue, coreIndex);
				while (1) {
					if (queue->val == coreIndex) {
						newCore->WaitingForMemory = 1;
						write_miss_trans(oldReg, newReg, newCore, coreIndex, MESI, Dsram, rd, block_tag, cache_tag, ALU, block_index, offset_in_block);
						break;
					}
				}
				pop(queue);
			}
		}
	}
}

void wb(Core* oldCore, Core* newCore) {
	if (oldCore->mem_wb->stall == 1)
		return;
	newCore->num_instructions++;
	int cmd = oldCore->mem_wb->IR;
	int ALU = oldCore->mem_wb->ALU;
	int rd = oldCore->mem_wb->rd;
	int MD = oldCore->mem_wb->MD;
	int OpCode = (cmd >> 24) & 0xFF;
	if (OpCode >= 0 && OpCode < 9)
		newCore->regs[rd] = ALU;
	else {
		if (OpCode == 16)
			newCore->regs[rd] = MD;

	}
	newCore->mem_wb->stall = 1;
}

//// **** -------------------------------------- Memory and Decode Stalls ------------------------------------- **** ////

void mem_stall(Core* oldCore, Core* newCore) { // Memory stall if to_stall = 1, to_stall = 0 no need to stall
	newCore->num_mem_stall++;
	newCore->mem_wb->stall = 1; // wb is stalled for the next clock cycle
	// saving the previous state in regs before mem
	If_Id_copy(oldCore->if_id, newCore->if_id);
	Id_Ex_copy(oldCore->id_ex, newCore->id_ex);
	Ex_Mem_copy(oldCore->ex_mem, newCore->ex_mem);
	newCore->pc = oldCore->pc;
}

void perform_decode_stall(Reg* oldReg, Reg* newReg, int CoreIndex) {
	If_Id_copy(oldReg->core[CoreIndex]->if_id, newReg->core[CoreIndex]->if_id); // recover previous if_id state
	newReg->core[CoreIndex]->pc = oldReg->core[CoreIndex]->pc; // take back old pc
	newReg->core[CoreIndex]->id_ex->stall = 1; // perform stall in execute next clock cycle
	if (newReg->core[CoreIndex]->WaitingForMemory == 1)
		newReg->core[CoreIndex]->num_decode_stall++;
}

void check_decode_stall(Reg* oldReg, Reg* newReg, int CoreIndex) {
	int id_opcode = newReg->core[CoreIndex]->id_ex->IR >> 24 & 0xFF;
	int ex_opcode = newReg->core[CoreIndex]->ex_mem->IR >> 24 & 0xFF;
	int mem_opcode = newReg->core[CoreIndex]->mem_wb->IR >> 24 & 0xFF;
	int wb_opcode = oldReg->core[CoreIndex]->mem_wb->IR >> 24 & 0xFF;

	if ((newReg->core[CoreIndex]->ex_mem->rd == newReg->core[CoreIndex]->id_ex->rs ||  // pipeline ex_mem.rd == (id_ex.rs or id_ex.rt) Read after Write hazrad
		newReg->core[CoreIndex]->ex_mem->rd == newReg->core[CoreIndex]->id_ex->rt) && newReg->core[CoreIndex]->ex_mem->rd != 0 && newReg->core[CoreIndex]->ex_mem->rd != 1
		&& ex_opcode != 17) //execute opcode != sw , these stalls doesnt apply at sw 
	{
		perform_decode_stall(oldReg, newReg, CoreIndex);
	}
	else if ((newReg->core[CoreIndex]->mem_wb->rd == newReg->core[CoreIndex]->id_ex->rs ||  // pipeline mem_wb.rd == (id_ex.rs or id_ex.rt) Read after Write hazard
		newReg->core[CoreIndex]->mem_wb->rd == newReg->core[CoreIndex]->id_ex->rt) && newReg->core[CoreIndex]->mem_wb->rd != 0 && newReg->core[CoreIndex]->ex_mem->rd != 1
		&& mem_opcode != 17)//memory opcode != sw 
	{
		perform_decode_stall(oldReg, newReg, CoreIndex);
	}
	else if ((id_opcode == 17) && // decode ==> sw 
		(wb_opcode != 17 && wb_opcode != 19) // wb != (sw) (if sw executing at wb, such stalls doesn't apply)
		&& newReg->core[CoreIndex]->id_ex->rd == oldReg->core[CoreIndex]->mem_wb->rd)//if wb pipe rd  == decode_exec pipe rd
	{
		perform_decode_stall(oldReg, newReg, CoreIndex);
	}
	else
		if ((id_opcode == 17) && //decode ==> sw 
			(mem_opcode != 17) // wb != (sw) (if sw executing at wb, such stalls doesn't apply)
			&& newReg->core[CoreIndex]->id_ex->rd == newReg->core[CoreIndex]->mem_wb->rd)//if memory pipe rd  == decode_exec pipe rd
		{
			perform_decode_stall(oldReg, newReg, CoreIndex);
		}

		else
			if ((id_opcode == 17) && //decode ==> sw or sc
				(ex_opcode != 17) // wb != (sw) (if sw executing at wb, such stalls doesn't apply)
				&& newReg->core[CoreIndex]->id_ex->rd == newReg->core[CoreIndex]->ex_mem->rd) //if memory pipe rd  == decode_exec pipe rd
			{
				perform_decode_stall(oldReg, newReg, CoreIndex);
			}
}

//// **** ----------------------------------- Writing To Output Files ---------------------------------- **** ////

Write_to_Coretrace(FILE* core_trace, Reg* newReg, Core* oldCore, Core* newCore) {
	char pc_fetch[50] = "---";
	if (oldCore->fetch != -1)
		sprintf(pc_fetch, "%03x", oldCore->pc);
	char pc_decode[50] = "---";
	if (oldCore->if_id->stall == 0)
		sprintf(pc_decode, "%03x", oldCore->if_id->NXT_PC - 1);
	char pc_exec[50] = "---";
	if (oldCore->id_ex->stall == 0)
		sprintf(pc_exec, "%03x", oldCore->id_ex->NXT_PC - 1);
	char pc_mem[50] = "---";
	if (oldCore->ex_mem->stall == 0)
		sprintf(pc_mem, "%03x", oldCore->ex_mem->NXT_PC - 1);
	char pc_wb[50] = "---";
	if (oldCore->mem_wb->stall == 0)
		sprintf(pc_wb, "%03x", oldCore->mem_wb->NXT_PC - 1);
	fprintf(core_trace, "%d %03s %03s %03s %03s %03s ", newReg->clk, pc_fetch, pc_decode, pc_exec, pc_mem, pc_wb);
	int i = 0;
	for (; i < 13; i++) {
		fprintf(core_trace, "%08x ", oldCore->regs[i + 2]); 
	}
	fprintf(core_trace, "%08x\n", oldCore->regs[15]);
}

void Write_to_BUStrace(Reg* newReg, FILE* bustrace) {
	if (newReg->bus_cmd == 0)
		return;
	fprintf(bustrace, "%d %01X %01X %05X %08X\n", newReg->clk, newReg->bus_origid, newReg->bus_cmd, newReg->bus_addr, newReg->bus_data);
}
void Write_to_Regout(Reg* oldReg, FILE** files) {
	int i, j;
	for (i = 0; i < num_of_cores; i++) {
		for (j = 0; j < 14; j++) {
			fprintf(files[6 + i], "%08x\n", oldReg->core[i]->regs[j + 2]);
		}
	}
}

void Write_to_Memout(FILE* memOut) {
	int i;
	for (i = 0; i < Memory_Size; i++)
		fprintf(memOut, "%08X\n", Memory[i]);
}

void Write_to_Stats(Reg* newReg, FILE** files) {
	int i, j;
	for (i = 0; i < num_of_cores; i++) {
		fprintf(files[23 + i], "cycles %d\ninstructions %d\nread_hit% d\nwrite_hit %d\nread_miss %d\nwrite_miss %d\ndecode_stall %d\nmem_stall %d\n", newReg->core[i]->num_cycles + 1, newReg->core[i]->num_instructions
			, newReg->core[i]->num_read_hit, newReg->core[i]->num_write_hit, newReg->core[i]->num_read_miss, newReg->core[i]->num_write_miss
			, newReg->core[i]->num_decode_stall, newReg->core[i]->num_mem_stall);
	}
}

void Write_to_Dsram(Reg* newReg, FILE** files) {
	int i, j, k;
	for (j = 0; j < num_of_cores; j++)
		for (i = 0; i < 64; i++) {
			for (k = 0; k < 4; k++)
				fprintf(files[15 + j], "%08X\n", newReg->core[j]->DSRAM_Cache[i][k]);
		}			
}

void write_to_Tsram(Reg* newReg, FILE** files) {
	int i, j;
	for (j = 0; j < num_of_cores; j++)
		for (i = 0; i < 64; i++)
			fprintf(files[19 + j], "%08X\n", newReg->core[j]->TSRAM_Cache[i]);
}

//// **** ------------------------------------- Operating System --------------------------------------- **** ////

void Clock_Cycle(Reg* oldReg, Reg* newReg, FILE** files) {
	for (int i = 0; i < num_of_cores; i++) {
		if (oldReg->core[i]->finish)
			continue;
		newReg->core[i]->pc = oldReg->core[i]->pc + 1;
		wb(oldReg->core[i], newReg->core[i]);
		mem(oldReg->core[i], newReg->core[i], i, oldReg, newReg);
		exec(oldReg->core[i], newReg->core[i]);
		decode(oldReg->core[i], newReg->core[i]);
		fetch(oldReg->core[i], newReg->core[i]);

		if (newReg->core[i]->WaitingForMemory) // if core is waiting for memory and memory not done yet
			mem_stall(oldReg->core[i], newReg->core[i]);
		Write_to_Coretrace(files[10 + i], newReg, oldReg->core[i], newReg->core[i]);

		if (newReg->core[i]->fetch == -1 && newReg->core[i]->if_id->stall &&  // if all conditions are met the has finished
			newReg->core[i]->id_ex->stall && newReg->core[i]->ex_mem->stall && // fetch, decode, exec, mem and wb are done
			newReg->core[i]->mem_wb->stall)
		{
			newReg->core[i]->finish = 1;
			newReg->core[i]->num_cycles = oldReg->clk;
		}
	}
}

int At_Least_OneCore_NotFinished(Reg* oldReg) {
	return !(oldReg->core[0]->finish) || !(oldReg->core[1]->finish) || !(oldReg->core[2]->finish) || !(oldReg->core[3]->finish);
}

int main(int argc, char* argv[]) {
	int i;
	Core* oldCore[4];
	Core* newCore[4];
	FILE* files[27];
	openFiles(files, argv);
	get_memory(files[4]); // copy memin file to Memory
	for (i = 0; i < num_of_cores; i++)
		newCore[i] = create_core(files[i]);
	Reg* newReg = create_reg(newCore); /// newReg ---> D
	for (i = 0; i < num_of_cores; i++)
		oldCore[i] = create_core(files[i]); /// oldReg ---> Q
	Reg* oldReg = create_reg(oldCore);
	while (At_Least_OneCore_NotFinished(oldReg)) {
		Clock_Cycle(oldReg, newReg, files);
		Write_to_BUStrace(newReg, files[14]);
		newReg->bus_cmd = 0;
		newReg->clk++;
		Regs_copy(oldReg, newReg);
	}
	Write_to_Regout(oldReg, files);
	Write_to_Memout(files[5]);
	Write_to_Stats(newReg, files);
	Write_to_Dsram(newReg, files);
	write_to_Tsram(newReg, files);
	closeFiles(files);
	return 0;
}