# `make all` OR `make verilate` -> `make build` -> `make sim` -> `make waves`

TOP_MODULE=attentionSystolicArray
RTL_PATH=../rtl
PARAM_HEADER=params.h

COMMON_SRC = ../common/common.cpp


# Parameters to pass (default values; override via command line)
N ?= 4
d ?= 4
K ?= 10

.PHONY: all
all: clean gen_params_header verilate build sim waves

.PHONY: gen_params_header
gen_params_header: $(PARAM_HEADER)

$(PARAM_HEADER):
	@echo "// Automatically generated parameters for C++ testbench" > $@
	@echo "#ifndef PARAMS_H" >> $@
	@echo "#define PARAMS_H" >> $@
	@echo "#define N $(N)" >> $@
	@echo "#define d $(d)" >> $@
	@echo "#define K $(K)" >> $@
	@echo "#endif // PARAMS_H" >> $@

.PHONY: verilate
verilate: gen_params_header .stamp.verilate

.PHONY: build
build: obj_dir/VattentionSystolicArray

.PHONY: sim
sim: waveform.vcd

.PHONY: waves
waves: waveform.vcd
	@echo
	@echo "### WAVES ###"
	gtkwave waveform.vcd

waveform.vcd: ./obj_dir/V$(TOP_MODULE)
	@echo
	@echo "### SIMULATING ###"
	@./obj_dir/V$(TOP_MODULE) +verilator+rand+reset+2

./obj_dir/V$(TOP_MODULE): .stamp.verilate
	@echo
	@echo "### BUILDING SIM ###"
	make -C obj_dir -f V$(TOP_MODULE).mk V$(TOP_MODULE)

.stamp.verilate: $(RTL_PATH)/$(TOP_MODULE).sv tb_$(TOP_MODULE).cpp
	@echo
	@echo "### VERILATING ###"
	verilator -Wall --trace --x-assign unique --x-initial unique -cc -I../rtl/ \
		$(RTL_PATH)/$(TOP_MODULE).sv --exe tb_$(TOP_MODULE).cpp $(COMMON_SRC) \
		-GN=$(N) -Gd=$(d) -GK=$(K) \
		-CFLAGS "-I../common"
	@touch .stamp.verilate

.PHONY: lint
lint: $(RTL_PATH)/$(TOP_MODULE).sv
	verilator --lint-only -I../rtl/ $(RTL_PATH)/$(TOP_MODULE).sv

.PHONY: clean
clean:
	rm -rf .stamp.*;
	rm -rf ./obj_dir
	rm -rf waveform.vcd
	rm -f $(PARAM_HEADER)