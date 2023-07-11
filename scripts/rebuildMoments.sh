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
    rm -rf $workdir || clean_exit "unable to clear space for workdir $workdir"
    if [ ! -d $workdir ]; then
        mkdir -p $workdir || clean_exit "unable to create workdir $workdir"
        cp src/C_buildMomentsMatrix* $workdir || clean_exit "unable to populate workdir $workdir"
    fi
    cd $workdir || clean_exit "unable to descend into workdir $workdir"
    cat >rebuild.py <<EOI
import sys
import gluex_moments_analysis.analyzeMomentsMatrix as ana
h=ana.analyze_moments(massEtaPi0_limits=($mspec), abst_limits=($tspec), model=1)
sys.exit(0)
EOI
    python3 rebuild.py || clean_exit "script rebuild.py failed"
}

pidlist=""
if [ "$mspec" != "" ]; then
    build_moments $mspec $tspec >rebuildMoments_${mspec}_${tspec}.log 2>&1 & pidlist="$pidlist $!"
else
    seqno=0
    for tspec in 0.0,0.2 0.2,0.3 0.3,0.45 0.45,0.6 0.6,0.9 0.9,1.2 1.2,1.6 1.6,2.4; do
    #for tspec in 0.0,0.3 0.3,0.6 0.6,1.2 1.2,2.5; do
        for mspec in 0.6,0.7 0.7,0.8 0.8,0.9 0.9,1.0 1.0,1.1 1.1,1.2 1.2,1.3 1.3,1.4 1.4,1.5 1.5,1.6 1.6,1.7 1.7,1.8 1.8,1.9 1.9,2.0 2.0,2.1 2.1,2.2 2.2,2.3 2.3,2.4; do
        #for mspec in 0.6,0.75 0.75,0.9 0.9,1.05 1.05,1.2 1.2,1.35 1.35,1.5 1.5,1.65 1.65,1.8 1.8,1.95 1.95,2.1 2.1,2.25 2.25,2.4; do
        #for mspec in 0.6,0.9 0.9,1.2 1.2,1.5 1.5,1.8 1.8,2.1 2.1,2.5; do
            seqno=$(expr $seqno + 1)
            inode=$(expr \( $seqno / 4 \) % 38 + 410)
            if [ $(hostname) != "cn$inode" ]; then
                continue
            fi
            logfile="rebuildMoments_${mspec}_${tspec}.log"
	    rm -f $logfile
	    if [ ! -f $logfile ]; then
                build_moments $mspec $tspec >$logfile 2>&1 & pidlist="$pidlist $!"
	        if [ $(echo $pidlist | wc -w) -ge $maxprocesses ]; then
		    set $pidlist
		    wait $1 || clean_exit "child proces $pid crashed"
                    pidlist="$2 $3 $4 $5 $6 $7 $8 $9"
                fi
            fi
        done
    done
fi

for pid in $pidlist; do
    wait $pid || clean_exit "child process $pid crashed"
done
