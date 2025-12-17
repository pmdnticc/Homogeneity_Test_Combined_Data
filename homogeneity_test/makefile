COMPILER=gfortran
FC=$(COMPILER)
OPT=-O2
#DEBUG= -ggdb -pg
ifdef DEBUG
OPT=-O0
#FPE=-ffpe-trap=invalid,zero,overflow
endif
FFLAGS=$(OPT) $(FPE)
FF=$(FC) -std=legacy -ffixed-line-length-none $(FFLAGS) -J$(OBJ) -I$(OBJ) -I$(INC) -I$(LAPACK)
OBJ=obj
CODEHOME=$(PWD)
SOURCEDIR=$(CODEHOME)/src
LAPACK=$(CODEHOME)/lapack-3.12.1
INC=$(CODEHOME)/Inc
OBJDIR=$(CODEHOME)/$(OBJ)
BINDIR=$(CODEHOME)/Bin
VPATH=$(DIRS)
DIRS=$(CODEHOME):$(SOURCEDIR):$(OBJDIR)


%.o: %.f $(OBJDIR)
	$(FF) -c -o $(OBJDIR)/$@ $<


%.o: %.f90 $(OBJDIR)
	$(FF) -c -o $(OBJDIR)/$@ $<


FILES = cernroutines.o
FILES+= read_input.o simu_data.o quartic_rt.o cdflib.o mle.o
FILES+= matrix_inv.o Fisher_Info.o test_statistics.o
FILES+= write2sas.o
FILES+= emp_TIE_power.o




main.exe: $(FILES)
	$(FF) $(patsubst %,$(OBJDIR)/%,$(FILES)) -L$(LAPACK) -o $@ -llapack -lrefblas
	mv main.exe $(BINDIR)
	@echo '    ----> Executable is main.exe <----'

.PHONY: clean

clean:
	rm -f $(patsubst %,$(OBJDIR)/%,$(FILES))


