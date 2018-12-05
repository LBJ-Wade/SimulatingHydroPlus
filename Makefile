PWD = `pwd`

DEBUG =
OPTIMIZATION = -O3
FLOWTRACE =
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE)
COMPILER = g++
LIBS = -lm
INCLUDES = -I.

.SUFFIXES: .o .cpp .h

.cpp.o :
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c $*.cpp
S1 =\
paramreader.cpp\
VH1+1.cpp
OBJ =\
paramreader.o\
VH1+1.o
EXE =\
vh

$(EXE) : $(OBJ) 
	$(COMPILER) $(OBJ) -o $(EXE)  $(LIBS)

# clean up misc files
clean :
	rm -f $(EXE) $(EXE).exe *\.o *~ #* core 
