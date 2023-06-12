#!/usr/bin/env python3

import hddm_s

events_per_file = 10000
nfout = 0
nevents = 0
for rec in hddm_s.istream("root://nod25.phys.uconn.edu/Gluex/simulation/moments-6-2023/eta_pi0_p.hddm"):
   if nevents % events_per_file == 0:
      fout = hddm_s.ostream("eta_pi0_p_{0:03d}.hddm".format(nfout))
      print(f"writing output file {nfout}")
      nfout += 1
   fout.write(rec)
   nevents += 1
