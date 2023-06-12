#include <TH1D.h>
#include <iostream>

void fill_M_histograms(const double *M, const double *Mvar, int m, int mdim, TH1D *hMabovediag, TH1D *hMbelowdiag, TH1D *hEabovediag, TH1D *hEbelowdiag)
{
  int mm(m * mdim + m);
  int mn(m * mdim);
  int n;
  for (n=0; n < m; ++n, ++mn) {
    int nn(n * mdim + n);
    hMbelowdiag->Fill(M[mn] / sqrt(M[mm] * M[nn]));
    hEbelowdiag->Fill(sqrt(Mvar[mn] / (M[mm] * M[nn])));
  }
  ++mn;
  for (n=m+1; n < mdim; ++n, ++mn) {
    int nn(n * mdim + n);
    hMabovediag->Fill(M[mn] / sqrt(M[mm] * M[nn]));
    hEabovediag->Fill(sqrt(Mvar[mn] / (M[mm] * M[nn])));
  }
}
