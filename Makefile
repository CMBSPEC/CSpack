#===============================================================================
# set to C++ compiler that should be used
#===============================================================================
OPENMPVAR = -fopenmp
LIBOPENMP = -lgomp

CC= g++-mp-7
#CC= g++

GSL=gsl
GSL_INC_PATH = /opt/local/include/
GSL_LIB_PATH = /opt/local/lib/
#GSL_INC_PATH = /usr/local/include/
#GSL_LIB_PATH = /usr/local/lib/

GSL_LIB = $(GSL_LIB_PATH)/lib$(GSL).a
LIB_GSL = $(GSL)
LIBSGSL= -L$(GSL_LIB_PATH) -l$(LIB_GSL) -lgslcblas

ifeq ($(RUNSCINET),yes)
GSL          = gsl
GSL_INC_PATH = $(SCINET_GSL_INC)
GSL_LIB_PATH = $(SCINET_GSL_LIB)
LIBSGSL = -L$(GSL_LIB_PATH) -l$(GSL) -lgslcblas
endif

#===============================================================================
# OPENMP
#===============================================================================
USE_OPENMP   = -D OPENMP_ACTIVATED
ifeq ($(RUNSCINET),yes)
LIBOPENMP = -lgomp
endif

#===============================================================================
# grace lib (for outputs with xmgrace from the diffusion modules)
#===============================================================================
# If grace support should be used uncomment the next 3 definitions
# and set them accordingly.
#===============================================================================
USE_GRACE   = -D GRACE_DEFINED

ifeq ($(USE_GRACE),-D GRACE_DEFINED)
GR          = grace_np
GR_INC_PATH = /opt/local/include/
GR_LIB_PATH = /opt/local/lib/
LIBSGR  = -L$(GR_LIB_PATH)  -l$(GR)
endif

#===============================================================================
#===============================================================================

LIBS = $(LIBOPENMP) $(LIBSGSL) $(LIBSGR) -lm


#===============================================================================
# compiler and linker flags
#===============================================================================
CXXFLAGS =  -Wall -pedantic -O3 $(OPENMPVAR)
CXXFLAGSLOC = $(CXXFLAGS) $(USE_GRACE) $(USE_OPENMP) -D SOLVERPATH=\"$(PWD)/\"
LXXFLAGS =

#===============================================================================
# object-files
#===============================================================================
OBJSLIB = ./tools/tools_for_cs.o \
          ./tools/Patterson.o \
          ./kesolver_tools/Compton_Kernel.o \
          ./tools/routines.o \
          ./tools/define_PDE_Kernel.o \
          ./tools/Integration_routines.GSL.o \
          ./tools/CSpack_functions.o \
          ./kesolver_tools/CSpack.o \
          ./tools/ODE_solver_Rec.o \
          ./tools/PDE_solver.o

OBJS = ./tools/parser.o

#===============================================================================
# program	 
#===============================================================================
all: thermalization_solverlib thermalization_solver

thermalization_solver: libthermalization_solver.a $(OBJS) main_cs.o
			 @echo "Linking..."
			 $(CC) $(LXXFLAGS) -L. -lthermalization_solver $(OBJS)  main_cs.o $(LIBS) -o thermalization_solver


lib: thermalization_solverlib

thermalization_solverlib: $(OBJSLIB)
			 @echo "\n Creating thermalization_solver lib library\n"
			 ar rvs ./libthermalization_solver.a $?

clean:
	rm -f ./tools/*.o ./kesolver_tools/*.o *.o

tidy:
	rm -f ./tools/*.o ./kesolver_tools/*.o ./tools/*.o~
	rm -f *.o *~ thermalization_solver libthermalization_solver.a

wipeDS:
	find . -type f -name \.DS_Store -print | xargs rm

./main_cs.o: ./tools/Grace_functions.cpp

#===============================================================================
# rules
#===============================================================================
INC_PATH = -I./include -I. \
           -I$(GSL_INC_PATH) -I.
		   
.cpp.o:
	@echo "Producing object-file $@"
	$(CC) $(CXXFLAGSLOC) $(INC_PATH) -c $< -o $@
	
#===============================================================================
#===============================================================================
