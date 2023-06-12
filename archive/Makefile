# ---- Link ---------------------------
C_buildMomentsMatrix.so: C_buildMomentsMatrix.o
	gcc -mtune=native -O4 -mavx2 -shared -o C_buildMomentsMatrix.so C_buildMomentsMatrix.o

# ---- gcc C compile ------------------
C_buildMomentsMatrix.o: C_buildMomentsMatrix.c C_buildMomentsMatrix.h
	gcc -c -O4 -fPIC -mtune=native -mavx2 -I/usr/include/python3.6m -I/usr/lib64/python3.6/site-packages/numpy/core/include C_buildMomentsMatrix.c

