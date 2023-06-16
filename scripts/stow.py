#!/usr/bin/env python3
#
# stow.py - save any stashed output files to the output storage area.
#
# author: richard.t.jones at uconn.edu
# version: march 16, 2023

import subprocess
import sys
import re

store = "https://cn440.storrs.hpc.uconn.edu:2843/Gluex/resilient/simulation/moments-6-2023"
stash = "/home/www/docs/halld/moments-6-2023/"

def usage():
   print("Usage: stow.py")
   sys.exit(1)

if len(sys.argv) > 1:
   usage()

templates = ["eta_pi0_p_{}_geant4.hddm",
             "eta_pi0_p_{}_smeared.hddm",
             "eta_pi0_p_{}_rest.hddm",
             "eta_pi0_p_{}_rest.root",
             "eta_pi0_p_x10_{}_geant4.hddm",
             "eta_pi0_p_x10_{}_smeared.hddm",
             "eta_pi0_p_x10_{}_rest.hddm",
             "eta_pi0_p_x10_{}_rest.root",
            ]
regexps = [re.compile(template.format("([0-9]*)")) for template in templates]

lsproc = subprocess.Popen(["find", stash, "-maxdepth", "1", "-type", "f", "-mmin", "+10", "-exec", "ls", "-l", "{}", ";"], stdout=subprocess.PIPE)
for line in lsproc.communicate()[0].split(b'\n'):
   sline = line.decode('utf-8').split()
   if len(sline) < 9:
      continue
   for i in range(len(regexps)):
      fname = sline[8].split('/')[-1]
      m = regexps[i].match(fname)
      if m:
         print(f"gfal-copy -f --copy-mode streamed file://{stash}/{fname} {store}/{fname} && /bin/rm {stash}/{fname}")
