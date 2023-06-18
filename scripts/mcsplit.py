#!/usr/bin/env python3

import hddm_s

n = 0
nout = 0
for rec in hddm_s.istream("root://cn440.storrs.hpc.uconn.edu/Gluex/resilient/simulation/moments-6-2023/eta_pi0_p_x10.hddm"):
   if n % 10000 == 0:
      fout = hddm_s.ostream(f"eta_pi0_p_x10_{nout}.hddm")
      nout += 1
   fout.write(rec)
   n += 1
