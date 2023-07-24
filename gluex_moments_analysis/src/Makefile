CONFINCLUDEPY = $(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('CONFINCLUDEPY'))")
NUMPYINCLUDE = $(shell python3 -c "import numpy; print(numpy.get_include())")

# ---- Link ---------------------------
C_buildMomentsMatrix.so: C_buildMomentsMatrix.o
	gcc -mtune=native -O4 -mavx2 -shared -o C_buildMomentsMatrix.so C_buildMomentsMatrix.o

# ---- gcc C compile ------------------
C_buildMomentsMatrix.o: C_buildMomentsMatrix.c C_buildMomentsMatrix.h
	gcc -c -O4 -fPIC -mtune=native -mavx2 -I$(CONFINCLUDEPY) -I$(NUMPYINCLUDE) C_buildMomentsMatrix.c

trial_model_C.so: trial_model.C trial_model.h
	root -l -b -q Load_trial_model.C

DSelector_etapi0_moments_C.so: DSelector_etapi0_moments.C DSelector_etapi0_moments.h trial_model_C.so
	root -l -b -q Load_DSelector.C

generated_moments: generated_moments.C generated_moments.h trial_model_C.so
	g++ -std=c++11 -O4 -o $@ -I$(shell root-config --incdir) -I$(HALLD_RECON_HOME)/$(BMS_OSNAME)/include -I$(ROOT_ANALYSIS_HOME)/$(BMS_OSNAME)/include -I$(HDF5ROOT)/include $< $(shell root-config --libs) -lMathMore -L $(HALLD_RECON_HOME)/$(BMS_OSNAME)/lib -L $(ROOT_ANALYSIS_HOME)/$(BMS_OSNAME)/lib trial_model_C.so -lHDDM -lxstream -lz -lbz2 -L$(HDF5ROOT/lib) -lhdf5_hl_cpp -lhdf5_hl -lhdf5_cpp -lhdf5 

clean:
	rm -rf *.o *.so *.pcm *.d
