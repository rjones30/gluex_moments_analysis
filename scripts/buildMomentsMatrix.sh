#!/bin/bash
#
# buildMomentsMatrix.sh - construct a Monte Carlo estimate of the acceptance of GlueX
#                         for a specific multi-particle reaction in the form a matrix
#                         of spherical harmonic moments, formed as a sum over a Monte
#                         Carlo reconstructed sample, generated flat in CM angles.
#
# author: richard.t.jones at uconn.edu
# version: june 9, 2023
#
# usage: [ run within a GlueX singularity container ]
#        $ ./buildMomentsMatrix.sh <seqNo> [<offset>]

nevents=21000
nthreads=28
memoryGB=32
diskGB=32

treename="etapi0_moments"
treefile="etapi0_moments.root"
xrootdurl="root://cn440.storrs.hpc.uconn.edu"
httpsURL="https://cn441.storrs.hpc.uconn.edu:2843"
remotepath="Gluex/resilient/simulation/moments-6-2023"
inputURL="https://gryphn.phys.uconn.edu/halld/moments-6-2023"
wget="wget --ca-directory=/etc/grid-security/certificates"

function usage() {
    echo "Usage: buildMomentsMatrix.sh <seqNo> [<offset]"
    echo "  where <seqNo> = job sequence number, 0,1,2,..99"
    echo "        <offset> = job sequence number offset, default 0"
    exit 1
}

function clean_exit() {
    ls -l 
    rm -f *.o *.so *.h5
    if [ "$1" = "" -o "$1" = "0" ]; then
    	cd $jobdir
		rm -rf $workdir
        echo "Successful exit from workscript."
        exit 0
    fi
    echo "Error $1 in workscript, $2"
    while true; do
        msg=$(echo "Error $1 in workscript, $2" | sed 's/ /_/g')
        eval $($wget -O- "$inputURL/scripts/onerror?msg=$msg" 2>/dev/null)
        sleep 10
    done
    exit $1
}

function save_output() {
    maxretry=5
    retry=0
    while [[ $retry -le $maxretry ]]; do
        gfal-copy -f --copy-mode streamed file://`pwd`/$1 $httpsURL/$remotepath/$2 2>gfal-copy.err
        retcode=$?
        if [[ -s gfal-copy.err ]]; then
            cat gfal-copy.err
            retcode=$(expr $retcode + 256)
        fi
        rm gfal-copy.err
        if [[ $retcode = 0 ]]; then
            rm $1
            return 0
        elif [[ $retry -lt $maxretry ]]; then
            retry=$(expr $retry + 1)
            echo "gfal-copy returned error code $retcode, waiting $retry minutes before retrying"
            sleep $(expr $retry \* 60)
        else
            retry=$(expr $retry + 1)
            echo "gfal-copy returned error code $retcode, giving up"
        fi
    done
    return $retry
}

if [ $# -lt 1 -o $# -gt 2 ]; then
    usage
    exit 
fi

if [[ $# = 1 ]]; then
   seqNo=$(expr $1 + 0)
elif [[ $# = 2 ]]; then
   seqNo=$(expr $1 + $2 + 0)
else
   usage
fi

echo "job $simType $seqNo running on" $(hostname)
export PYTHONPATH="/home/rtj02001/.local/lib/python3.8/site-packages"
python --version
python -m site

jobdir=$(pwd)
workdir=/local/slurm-job-$SLURM_JOBID
cp * $workdir
cd $workdir

make
echo "buildMomentsMatrix.py $xrootdurl/$remotepath/$treefile:$treename $nevents $seqNo"
./buildMomentsMatrix.py $xrootdurl/$remotepath/$treefile:$treename $nevents $seqNo \
 || clean_exit $? "buildMomentsMatrix.py crashed"

ls -l
cp x509up_u7896 /tmp/x509up_u$UID
outfile=$(echo acceptance_moments_*.h5)
save_output $outfile $outfile || clean_exit $? "save of $outfile failed"
clean_exit
