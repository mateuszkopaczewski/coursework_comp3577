ROOTDIR=$(shell pwd)
OUTPUTDIR=$(ROOTDIR)/paraview-output/

.PHONY: all cleanall clean clean_paraview
all: step-1-gcc step-2-gcc step-3-gcc step-4-gcc step-1-icpc step-2-icpc step-3-icpc step-4-icpc
step-%: step-%-gcc step-%-icpc

# Target to be used with the GNU Compiler Collection.
step-%-gcc step-%-gcc.o NBodySimulation-gcc.o: CXX?=g++
step-%-gcc step-%-gcc.o NBodySimulation-gcc.o: CXXFLAGS?=-fopenmp -O3 -march=native -std=c++0x -fno-math-errno
NBodySimulation-gcc.o: NBodySimulation.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
step-%-gcc.o: step-%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
step-%-gcc: NBodySimulation-gcc.o step-%-gcc.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Target to be used with the Intel C++ compiler.
# In order to use this compiler on Hamilton, you should first add the
# corresponding module with
#     $ module add intel/2021.4
step-%-icpc step-%-icpc.o NBodySimulation-icpc.o: CXX?=icpc
step-%-icpc step-%-icpc.o NBodySimulation-icpc.o: CXXFLAGS?=-qopenmp -O3 -xHost -std=c++0x
NBodySimulation-icpc.o: NBodySimulation.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
step-%-icpc.o: step-%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
step-%-icpc: NBodySimulation-icpc.o step-%-icpc.o
	$(CXX) $(CXXFLAGS) -o $@ $^

.silent: cleanall clean clean_paraview

cleanall: clean clean_paraview

clean:
	rm -rf $(ROOTDIR)/step-*-gcc $(ROOTDIR)/step-*-icpc $(ROOTDIR)/*.o

clean_paraview:
	if test -d "$(OUTPUTDIR)"; then \
		find $(OUTPUTDIR) -iname "result-*.vtp" -delete; \
		rm -rf $(OUTPUTDIR)/result.pvd; \
	fi

test:
	./validate.sh
	python3 validate.py
