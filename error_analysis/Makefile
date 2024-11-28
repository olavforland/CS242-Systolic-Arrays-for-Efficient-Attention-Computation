# Makefile for error_analysis

TOP_MODULE=attentionSystolicArray
RTL_PATH=../rtl
PARAM_HEADER=params.h

# Parameters (default values; can be overridden via command line)
N ?= 4
d ?= 4
K ?= 10
MAGNITUDE_SCALE ?= 1.0

# Common source file
COMMON_SRC = ../common/common.cpp

.PHONY: all
all: clean gen_params_header verilate build sim

.PHONY: gen_params_header
gen_params_header: $(PARAM_HEADER)

$(PARAM_HEADER):
	@echo "// Automatically generated parameters for C++ testbench" > $@
	@echo "#ifndef PARAMS_H" >> $@
	@echo "#define PARAMS_H" >> $@
	@echo "#define N $(N)" >> $@
	@echo "#define d $(d)" >> $@
	@echo "#define K $(K)" >> $@
	@echo "#define MAGNITUDE_SCALE $(MAGNITUDE_SCALE)" >> $@
	@echo "#endif // PARAMS_H" >> $@

.PHONY: verilate
verilate: gen_params_header .stamp.verilate

.stamp.verilate: $(RTL_PATH)/$(TOP_MODULE).sv attention_error.cpp $(COMMON_SRC)
	@echo
	@echo "### VERILATING ###"
	verilator -Wall --trace --x-assign unique --x-initial unique -cc -I$(RTL_PATH)/ \
		$(RTL_PATH)/$(TOP_MODULE).sv --exe attention_error.cpp $(COMMON_SRC) \
		-GN=$(N) -Gd=$(d) -GK=$(K) \
		-CFLAGS "-I../common"
	@touch .stamp.verilate

.PHONY: build
build: obj_dir/V$(TOP_MODULE)

obj_dir/V$(TOP_MODULE): .stamp.verilate
	@echo
	@echo "### BUILDING SIMULATION ###"
	make -C obj_dir -f V$(TOP_MODULE).mk V$(TOP_MODULE)

.PHONY: sim
sim: obj_dir/V$(TOP_MODULE)
	@echo
	@echo "### RUNNING SIMULATION ###"
	./obj_dir/V$(TOP_MODULE) +verilator+rand+reset+2

.PHONY: clean
clean:
	rm -rf .stamp.*
	rm -rf ./obj_dir
	rm -f waveform.vcd
	rm -f $(PARAM_HEADER)
	rm -f results/error_results.csv

.PHONY: lint
lint:
	verilator --lint-only -I$(RTL_PATH)/ $(RTL_PATH)/$(TOP_MODULE).sv
