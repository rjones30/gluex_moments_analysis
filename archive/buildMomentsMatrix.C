void buildMomentsMatrix(double *M, int m[3], int n[3], double *Ym[3], double *Yn[3])
{
  int mdim = m[0] * m[1] * m[2];
  for (int im0=0; im0 < m[0]; ++im0) {
    for (int in0=0; in0 < n[0]; ++in0) {
      double Y00 = Ym[0][im0] * Yn[0][in0];
      for (int im1=0; im1 < m[1]; ++im1) {
        int mm = (im0 * m[1] + im1) * m[2];
        for (int in1=0; in1 < n[1]; ++in1) {
          double Y0011 = Y00 * Ym[1][im1] * Yn[1][in1];
          for (int im2=0; im2 < m[2]; ++im2, ++mm) {
            double Y00112 = Y0011 * Ym[2][im2];
            int nn = (in0 * n[1] + in1) * n[2];
	    int mn = mm * ncols + nn;
            for (int in2=0; in2 < n[2]; ++in2, ++mn) {
              M[mn] += Y012 * Yn[2][in2];
	    }
	  }
	}
      }
    }
  }
}
