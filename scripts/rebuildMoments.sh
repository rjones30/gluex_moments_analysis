#!/bin/bash
#
# rebuildMoments.sh - script to automate the parallel generation of a mock
#                     data sample based on an arbitrary PWA model, for the
#                     purpose of demonstrating the use of support vectors
#                     to extract acceptance-corrected spherical moments
#                     for multi-particle decay reactions using the GlueX
#                     detector simulation.
#
# author: richard.t.jones at uconn.edu
# version: june 27, 2023

maxprocesses=1

function usage {
    echo "Usage: (inside gluex container) [/srv]> rebuildMoments.sh [<mX_min>,<Mx_max> <|t_min|,|t_max|>]"
    echo " where kinematic bounds on the resonance invariant mass <mX_min>,<mX_max> and Mandelstam t,"
    echo " <t_min>,<t_max> may be specified on the commandline or generated automatically if left blank."
    exit 1
}

function clean_exit {
    retcode=$?
    echo "Error in rebuildMoments.sh - $1"
    echo "quitting with exit code $retcode"
    exit $retcode
}

if [ $# = 2 ]; then
    mspec=$1
    tspec=$2
elif [ $# != 0 ]; then
    usage
elif [ $(pwd) != "/srv" ]; then
    usage
fi

source scripts/setup_root.sh
source scripts/setup_scripts.sh

function build_moments {
    mspec=$1
    tspec=$2
    workdir="rebuildMoments_${mspec}_${tspec}"
    #rm -rf $workdir || clean_exit "unable to clear space for workdir $workdir"
    #mkdir $workdir || clean_exit "unable to create workdir $workdir"
    #cp src/C_buildMomentsMatrix* $workdir || clean_exit "unable to populate workdir $workdir"
    cd $workdir || clean_exit "unable to descend into workdir $workdir"
    cat >rebuild.py <<EOI
import sys
import analyzeMomentsMatrix as ana
h=ana.analyze_moments(massEtaPi0_limits=($mspec), abst_limits=($tspec), model=1)
sys.exit(0)
EOI
    python3 rebuild.py || clean_exit "script rebuild.py failed"
}

pidlist=""
if [ "$mspec" != "" ]; then
    build_moments $mspec $tspec >rebuildMoments_${mspec}_${tspec}.log 2>&1 & pidlist="$pidlist $!"
else
    for tspec in 0.0,0.3 0.3,0.6 0.6,1.2; do
        for mspec in 0.9,1.2 1.2,1.5 1.5,1.8 1.8,2.1; do
            build_moments $mspec $tspec >rebuildMoments_${mspec}_${tspec}.log 2>&1 & pidlist="$pidlist $!"
	    if [ $(echo $pidlist | wc -w) -ge $maxprocesses ]; then
		set $pidlist
		wait $1 || clean_exit "child proces $pid crashed"
                pidlist="$2 $3 $4 $5 $6 $7 $8 $9"
            fi
        done
    done
fi

for pid in $pidlist; do
    wait $pid || clean_exit "child process $pid crashed"
done
