
#Compiler and compiler flags
CXX := mpic++
CXXFLAGS := -Wall -std=c++11 -g
LIB := -lconfig++
INC := -I include -I /usr/local/include/eigen3 -I/usr/local/Cellar/boost/1.74.0/include

#Directories
SRCDIR := src
BUILDDIR := build
TARGETDIR := bin

#Targets
TARGETS = $(TARGETDIR)/DaMaSCUS-Simulator $(TARGETDIR)/DaMaSCUS-Analyzer $(TARGETDIR)/DaMaSCUS-SD

#Source files
SRCEXT := cpp
COMMONSRC := $(shell find $(SRCDIR) -maxdepth 1 -type f -name '*.$(SRCEXT)')
SIMSRC :=$(COMMONSRC) $(shell find $(SRCDIR)/simulation -type f -name *.$(SRCEXT))
ANASRC :=$(COMMONSRC) $(shell find $(SRCDIR)/analysis -type f -name *.$(SRCEXT))
SDSRC :=$(COMMONSRC) $(shell find $(SRCDIR)/simulation_SD -type f -name *.$(SRCEXT))

#Object files
SIMOBJ := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SIMSRC:.$(SRCEXT)=.o))
ANAOBJ := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(ANASRC:.$(SRCEXT)=.o))
SDOBJ := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SDSRC:.$(SRCEXT)=.o))

.PHONY: all simulator analyzer clean codecov sdatm 

all: CXXFLAGS += -O2 
all: $(TARGETS)

test: LIB+= --coverage
test: $(TARGETS)

sdatm: CXXFLAGS += -O2
sdatm: $(TARGETDIR)/DaMaSCUS-SD


simulator: CXXFLAGS += -O2
simulator: $(TARGETDIR)/DaMaSCUS-Simulator

analyzer: CXXFLAGS += -O2
analyzer: $(TARGETDIR)/DaMaSCUS-Analyzer

$(TARGETDIR)/DaMaSCUS-Simulator: $(SIMOBJ)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LIB)

$(TARGETDIR)/DaMaSCUS-Analyzer: $(ANAOBJ)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LIB)

$(TARGETDIR)/DaMaSCUS-SD: $(SDOBJ)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)/simulation/
	@mkdir -p $(BUILDDIR)/analysis/
	@mkdir -p $(BUILDDIR)/simulation_SD/
	$(CXX) $(CXXFLAGS) $(INC) $(LIB) -o $@ -c $<

codecov:
	

clean:
	find . -type f -name '*.o' -delete
	find . -type f -name '*.gcno' -delete
	find . -type f -name '*.gcda' -delete
	find . -type f -name '*.gcov' -delete
	rm -f $(TARGETS)
