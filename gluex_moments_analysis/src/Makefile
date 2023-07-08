# ---- Link ---------------------------
C_buildMomentsMatrix.so: C_buildMomentsMatrix.o
	gcc -mtune=native -O4 -mavx2 -shared -o C_buildMomentsMatrix.so C_buildMomentsMatrix.o

# ---- gcc C compile ------------------
C_buildMomentsMatrix.o: C_buildMomentsMatrix.c C_buildMomentsMatrix.h
	gcc -c -O4 -fPIC -mtune=native -mavx2 -I/usr/include/python3.6m -I/usr/local/lib64/python3.6/site-packages/numpy/core/include C_buildMomentsMatrix.c

generated_moments: generated_moments.C generated_moments.h
	g++ -std=c++11 -O4 -o $@ -I$(shell root-config --incdir) -I$(HALLD_RECON_HOME)/$(BMS_OSNAME)/include -I$(HDF5ROOT)/include $< $(shell root-config --libs) -lMathMore -L $(HALLD_RECON_HOME)/$(BMS_OSNAME)/lib -lHDDM -lxstream -lz -lbz2 -L$(HDF5ROOT/lib) -lhdf5_hl_cpp -lhdf5_hl -lhdf5_cpp -lhdf5 
