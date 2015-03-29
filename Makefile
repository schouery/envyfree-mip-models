MAKEFLAGS+=" -j 16"
LIBFORMAT     = static_pic
CPLEXDIR      = /opt/cplex/cplex125/cplex
CONCERTDIR    = /opt/cplex/cplex125/concert
ifeq ($(shell uname), Darwin)
	SYSTEM        = x86-64_osx
	CCC = clang++ -O0
	CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -stdlib=libstdc++
	CCLNFLAGS = -lconcert -lilocplex -lcplex -m64 -lm -lpthread -framework CoreFoundation -framework IOKit
else
	SYSTEM     = x86-64_sles10_4.1
	CCC = g++
	CCOPT = -m64 -O3 -fPIC -fexceptions -DNDEBUG -DIL_STD -DTHREADS=16  -DWORKMEM=65536
	# CCOPT = -m64 -O3 -fPIC -fexceptions -DNDEBUG -DIL_STD -DTHREADS=2 -DWORKMEM=32768
	# CCOPT = -m64 -O3 -fPIC -fexceptions -DNDEBUG -DIL_STD -DTHREADS=1 -DWORKMEM=4096
	CCLNFLAGS = -lilocplex -lcplex -lconcert -lm -m64 -pthread
endif
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 


# Mantenha uma linha em branco no comeco desse arquivo
all: envyfree	 pricedown-modified

envyfree: envyfree.o graph.o utils.o utility.o stm.o profit.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) envyfree.o graph.o utils.o utility.o stm.o profit.o -o envyfree $(CCLNFLAGS)

envyfree.o: envyfree.cpp
	$(CCC) -c $(CCFLAGS) envyfree.cpp -o envyfree.o

graph.o: graph.cpp graph.h
	$(CCC) -c $(CCFLAGS) graph.cpp -o graph.o

utils.o: utils.cpp utils.h
	$(CCC) -c $(CCFLAGS) utils.cpp -o utils.o

utility.o: utility.cpp utility.h
	$(CCC) -c $(CCFLAGS) utility.cpp -o utility.o

profit.o: profit.cpp profit.h
	$(CCC) -c $(CCFLAGS) profit.cpp -o profit.o

stm.o: stm.cpp
	$(CCC) -c $(CCFLAGS) stm.cpp -o stm.o

pricedown-modified: pricedown-modified.cc
	g++ -std=gnu++0x pricedown-modified.cc -o pricedown-modified


clean:
	/bin/rm -rf *.o *~ envyfree pricedown-modified
