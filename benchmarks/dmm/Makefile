# Compilers
CC_X86=g++

# Other knobs
CODE=mkl_dense.cpp
MATRIX_CODE=mkl_dense_matrix.cpp
VECTOR_CODE=mkl_dense_vector.cpp
BIN=mkl_dense
MATRIX_BIN=mkl_dense_matrix
VECTOR_BIN=mkl_dense_vector

STDFLAGS=c++11
# MKLROOT=/opt/intel/compilers_and_libraries_2018.3.222/linux/mkl
MKLROOT=/home/fengsy/intel/compilers_and_libraries_2019.0.117/linux/mkl
IDIR=${MKLROOT}/include
CFLAGS=-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl  -DMKL_ILP64 -m64 -I${MKLROOT}/include
#CFLAGS=-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_tbb_thread -lmkl_core -ltbb -lstdc++ -lpthread -lm -ldl -DMKL_ILP64 -m64 -I${MKLROOT}/include
CFLAGS_STATIC=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -DMKL_ILP64 -m64 -I${MKLROOT}/include
#LIB=/usr/local/lib/libpapi.a

# No need to touch the following
define colorecho
      @tput setaf 6
      @tput bold
      @echo $1
      @tput sgr0
endef
DONE_STRING=$(call colorecho, " [DONE]")


local:	$(CODE)
	$(CC_X86) -std=$(STDFLAGS)  $(CFLAGS) $(CODE) -o $(BIN) -g
	@echo -n "Creating binary ./$(BIN) for $@"
	@$(DONE_STRING)

spmm:	$(CODE)
	$(CC_X86) $(CFLAGS) inner_mkl_spmm.cpp -o $(BIN) -g
	@echo -n "Creating binary ./$(BIN) for $@"
	@$(DONE_STRING)

papi:	$(CODE)
	$(CC_X86) $(CFLAGS) inner_mkl_papi.cpp -o inner_mkl_papi ~/papi-5.5.1/src/libpapi.a -g
	@echo -n "Creating binary ./inner_mkl_papi for $@"
	@$(DONE_STRING)

old:	$(CODE)
	$(CC_X86) $(CFLAGS) inner_mkl_spmm1.cpp -o $(BIN) -g
	@echo -n "Creating binary ./$(BIN) for $@"
	@$(DONE_STRING)

static:	$(CODE)
	$(CC_X86) $(CFLAGS_STATIC) $(CODE) $(LIB) -static -o $(BIN) -g
	@echo -n "Creating binary ./$(BIN) for $@"
	@$(DONE_STRING)

matrix:	$(MATRIX_CODE)
	$(CC_X86) -std=$(STDFLAGS)  $(CFLAGS) $(MATRIX_CODE) -o $(MATRIX_BIN) -g
	@echo -n "Creating binary ./$(MATRIX_BIN) for $@"
	@$(DONE_STRING)

matrix_power:	$(MATRIX_CODE)
	$(CC_X86) -std=$(STDFLAGS)  $(CFLAGS) $(MATRIX_CODE) -DPOWER=1 -o $(MATRIX_BIN) -g
	@echo -n "Creating binary ./$(MATRIX_BIN) for $@"
	@$(DONE_STRING)

vector:	$(VECTOR_CODE)
	$(CC_X86) -std=$(STDFLAGS)  $(CFLAGS) $(VECTOR_CODE) -o $(VECTOR_BIN) -g
	@echo -n "Creating binary ./$(VECTOR_BIN) for $@"
	@$(DONE_STRING)

vector_power:	$(VECTOR_CODE)
	$(CC_X86) -std=$(STDFLAGS)  $(CFLAGS) $(VECTOR_CODE) -DPOWER=1 -o $(VECTOR_BIN) -g
	@echo -n "Creating binary ./$(VECTOR_BIN) for $@"
	@$(DONE_STRING)

clean:
	@rm -f $(BIN) $(MATRIX_BIN) out.txt out_mod.txt
	@echo -n "Cleaning up"
	@$(DONE_STRING)

.PHONY:	clean
