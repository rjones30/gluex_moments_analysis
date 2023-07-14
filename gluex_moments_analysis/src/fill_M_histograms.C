#include <TH1D.h>
#include <iostream>

void fill_M_histograms(const double *M, int m, int mdim, TH1D *hMabovediag, TH1D *hMbelowdiag)
{
  int mm(m * mdim + m);
  int mn(m * mdim);
  int n;
  for (n=0; n < m; ++n, ++mn) {
    int nn(n * mdim + n);
    hMbelowdiag->Fill(M[mn] / sqrt(M[mm] * M[nn]));
  }
  ++mn;
  for (n=m+1; n < mdim; ++n, ++mn) {
    int nn(n * mdim + n);
    hMabovediag->Fill(M[mn] / sqrt(M[mm] * M[nn]));
  }
}
