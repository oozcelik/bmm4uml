
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --glibs) -lRooFit -lRooFitCore -lMinuit -lSPlot 

all: bmm4fitter bmm4toystudy bmm4validation  

bmm4fitter: bmm4fitter.cc bmm4common.h 
	g++ -Wall -Wextra -O3 -o $@ bmm4fitter.cc $(ROOTFLAGS) $(ROOTLIBS)

bmm4toystudy: bmm4toystudy.cc bmm4common.h
	g++ -Wall -Wextra -O3 -o $@ bmm4toystudy.cc $(ROOTFLAGS) $(ROOTLIBS)

bmm4validation: bmm4validation.cc bmm4common.h
	g++ -Wall -Wextra -O3 -o $@ bmm4validation.cc $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm -f bmm4fitter bmm4toystudy bmm4validation *~ 

