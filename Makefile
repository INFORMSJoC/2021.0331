#-------------------------------------------------------------------
# Compiler selection and flags
#--------------------------------------------------------------------
CC  = gcc
CCC = g++
CCCFLAGS = -std=c++11

#--------------------------------------------------------------------
#Library flags
#--------------------------------------------------------------------
CCLNFLAGS = -lgmp
CCCLNFLAGS = -lgmpxx

#--------------------------------------------------------------------
# Applications
#---------------------------------------------------------------------
APPS = ROU
all: $(APPS)

ROU: src/ROU.cpp 
	$(CCC) -o ./src/ROU src/ROU.cpp $(CCLNFLAGS) $(CCCLNFLAGS) $(CCCFLAGS) 

#--------------------------------------------------------------------
# Additional commands
#---------------------------------------------------------------------
clean:
	rm -f $(APPS)
