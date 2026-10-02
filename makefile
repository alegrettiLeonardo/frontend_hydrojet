#=======================================================================
# Makefile for HYDROJET (Fortran)
# baseado no modelo de RL Sartori
#=======================================================================

# dirs (ajuste se o nome da pasta de fontes for "src" em vez de "scr")
DIR_SRC = ./scr
DIR_BIN = ./bin
DIR_OBJ = ./tmp
DIR_MOD = $(DIR_OBJ)
DIR_DOC = ./doc

# folders a criar
BUILD_DIRS = $(DIR_BIN) $(DIR_OBJ) $(DIR_DOC)

# "versão" do código
SVN_REV = "2025.07"

# compilador Fortran
FC = gfortran

# flags de compilação
FFLAGS = -O2 -ffpe-trap=invalid,zero,overflow -fbacktrace -cpp -Wall \
         -std=legacy -fallow-argument-mismatch \
         -J$(DIR_MOD) -D'SVN_REV=$(SVN_REV)'


# opções de link
LDFLAGS =
LIBS    =

# sufixos
.SUFFIXES: .f .for .f90 .o

# onde procurar arquivos por extensão
vpath %.f   $(DIR_SRC)
vpath %.for $(DIR_SRC)
vpath %.f90 $(DIR_SRC)
vpath %.o   $(DIR_OBJ)

# nome do executável
PROJ = $(DIR_BIN)/hydrojet

# primeira regra
all: $(BUILD_DIRS) $(PROJ)

# cria pastas (se não existirem)
$(BUILD_DIRS):
	mkdir -p $(BUILD_DIRS)

#=======================================================================
# Regras de compilação genéricas
#=======================================================================

# .f -> .o (Fortran 77 fixo)
.f.o:
	@echo "====> compiling $@"
	$(FC) -c $(FFLAGS) -ffixed-line-length-none $< -o $(DIR_OBJ)/$@

# .f90 -> .o
.f90.o:
	@echo "====> compiling $@"
	$(FC) -c $(FFLAGS) -ffixed-line-length-none $< -o $(DIR_OBJ)/$@

# .for -> .o
.for.o:
	@echo "====> compiling $@"
	$(FC) -c $(FFLAGS) -ffixed-line-length-none $< -o $(DIR_OBJ)/$@

.PHONY: clean cleanall

# limpa objetos
clean:
	@echo "====> removing all inside $(DIR_OBJ)"
	rm -Rf $(DIR_OBJ)/*

cleanall:
	@echo "====> removing all inside $(DIR_OBJ)"
	rm -Rf $(DIR_OBJ)/*
	@echo "====> removing all inside $(DIR_BIN)"
	rm -Rf $(DIR_BIN)/*
	@echo "====> removing all inside $(DIR_DOC)"
	rm -Rf $(DIR_DOC)/*

#=======================================================================
# SEÇÃO ESPECÍFICA DO PROJETO HYDROJET
#=======================================================================

# Objetos usados no projeto (sem caminho)
FOBJ = \
	fileplot.o \
	foilsubs.o \
	initopst.o \
	findprop.o \
	optionst.o \
	compe2t.o \
	tiltsubs.o \
	prntrest.o \
	calcpadt.o \
	calcther.o \
	readwrit.o \
	mipropst.o \
	firstt.o \
	calcmesh.o \
	addcoefs.o \
	hydrojet.o \
	supportp.o \
	guespp.o \
	turbcoef.o \
	inputp.o \
	pjet.o \
	setpexit.o \
	compe1t.o \
	calcsoln.o \
	filmwt.o \
	spline.o

# Regra de link do executável
$(PROJ): $(FOBJ)
	@echo "====> building $@"
	$(FC) -o $@ $(FFLAGS) \
	      $(patsubst %, $(DIR_OBJ)/%, $(FOBJ)) \
	      $(LDFLAGS) $(LIBS)

# (Opcional) Dependências mais detalhadas podem ser colocadas aqui, se quiser.
# Exemplo:
# hydrojet.o : hydrojet.f params.o optionst.o
