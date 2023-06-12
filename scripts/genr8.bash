#!/bin/bash
#genr8 -Agenr8_output.txt -M10000000 -r10000 < omegadelta2.input
#genr8 -Agenr8_output.txt -M10000 -r10000 < omega3pi.input
genr8 -Agenr8_output.txt -M10000000 -r71000 < eta_pi0_p.input
genr8_2_hddm -V"0 0 50 80" genr8_output.txt
