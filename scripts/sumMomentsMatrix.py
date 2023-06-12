#!/usr/bin/env python3
#
# sumMomentsMatrix.py - reads from a set of saved acceptance moments computed
#                       on subsets of a Monte Carlo signal sample and saves
#                       the summed M and Mvar matrices in the same h5 format
#                       under the name acceptance_sum.h5 in the workdir.
#
# author: richard.t.jones at uconn.edu
# version: june 10, 2023

import numpy as np
import h5py
import sys

def usage():
  print("usage: sumMomentsMatrix.py <file1.h5> [...]")
  sys.exit(1)

if len(sys.argv) < 2 or sys.argv[1][0] == '-':
  usage()

Msum = 0
Mvarsum = 0
subset_start={0:0}
subset_end={0:0}
for finh5 in sys.argv[1:]:
  with h5py.File(finh5, 'r') as h5f:
    #if Msum:
    #  Msum += h5f["Moments"][:,:]
    #  Mvarsum += h5f["Moments_variance"][:,:]
    #else:
    #  Msum = h5f["Moments"][:,:]
    #  Mvarsum = h5f["Moments_variance"][:,:]
    count = h5f["accepted_subset"]
    skip = h5f["skipped_accepted"] 
    stop = skip + count
    if stop in subset_start and skip in subset_end:
      block = subset_start.pop(stop)
      start = subset_end.pop(skip)
      subset_start[start] += count + block
      subset_end[stop + block] = start
    elif stop in subset_start:
      block = subset_start.pop(stop)
      subset_start[skip] = block + count
      subset_end[stop + block] = skip
    elif skip in subset_end:
      start = subset_end.pop(skip)
      subset_end[stop] = start
      subset_start[start] += count
    else:
      subset_start[skip] = count
      subset_stop[stop] = skip

print("event ranges found:")
for start in subset_start:
  print("  (", start, ":", start + subset_start[start], ")")

""" 
if True:
  h5out = h5py.File("acceptance_moments_{0}_{1}.h5".format(event_count, sequence_no), "w")
  h5out.create_dataset("Moments", data=M)
  h5out.create_dataset("Moments_variance", data=Mvar)
  h5out.create_dataset("accepted_subset", data=event_count)
  h5out.create_dataset("skipped_accepted", data=skip_count)
  h5out.close()
"""
