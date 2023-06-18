#!/bin/bash
#
# workscript.sh - GlueX simulation samples script for production
#                 of benchmark Monte Carlo data on the osg.
#
# author: richard.t.jones at uconn.edu
# version: february 27, 2021
#
# usage: [ run within a GlueX singularity container ]
#        $ ./workscript.sh <simType> <seqNo> [<offset>]

runNo=71000
nthreads=1
nevents=10000

xrootdserver="cn440.storrs.hpc.uconn.edu"
xrootdURL="root://$xrootdserver"
remotepath="/Gluex/resilient/simulation/moments-6-2023"
gsiftpURL="gsiftp://cn440.storrs.hpc.uconn.edu"
httpsURL="https://cn441.storrs.hpc.uconn.edu:2843"
inputURL="https://gryphn.phys.uconn.edu/halld/moments-6-2023"
wget="wget --ca-directory=/etc/grid-security/certificates"

function usage() {
    echo "Usage: workscript.sh <simType> <seqNo> [<offset]"
    echo "  where <simType> = simulation type, eg. eta_pi0_p"
    echo "        <seqNo> = job sequence number, 1,2,..."
    echo "        <offset> = job sequence number offset, default 0"
    exit 1
}

function clean_exit() {
    ls -l 
    rm -f setup sim.sh control.in randoms worklist smear.root sample*.hddm dana_rest.hddm hd_root.root tree*.root *.hbook *.rz *.dat tmp.* *.astate hd_recon.config genr8.hddm
    if [ "$1" = "" -o "$1" = "0" ]; then
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
            break
        elif [[ $retry -lt $maxretry ]]; then
            retry=$(expr $retry + 1)
            echo "gfal-copy returned error code $retcode, waiting $retry minutes before retrying"
            sleep $(expr $retry \* 60)
        else
            retry=$(expr $retry + 1)
            echo "gfal-copy returned error code $retcode, giving up"
        fi
    done
    # fall through to allow job file transfer return results, failure not fatal
    mv $1 $(basename $2)
    return 0
}

if [ $# -lt 2 -o $# -gt 3 ]; then
    usage
    exit 
fi

if [[ $1 = "g4" ]]; then
   controlin=control.in
   simType=g4
elif [[ $1 = "g31" ]]; then
   controlin=control.in
   simType=g31
elif [[ $1 = "g34" ]]; then
   controlin=control.in
   simType=g34
elif [[ $1 = "eta_pi0_p" ]]; then
   controlin=control.in
   simType="$1"
elif [[ $1 = "eta_pi0_p_x10" ]]; then
   controlin=control.in
   simType="$1"
else
   usage
fi
if [[ $# = 2 ]]; then
   seqNo=$(expr $2 + 0)
elif [[ $# = 3 ]]; then
   seqNo=$(expr $2 + $3 + 0)
else
   usage
fi
seqNo3=$(echo $seqNo | awk '{printf("%03d",$1)}')

echo "job $simType $seqNo running on" $(hostname)

$wget -O setup $inputURL/scripts/setup.sh 2>/dev/null || clean_exit $? "cannot fetch setup.sh from web server"
$wget -O control.in $inputURL/config/$controlin 2>/dev/null || clean_exit $? "cannot fetch control.in from web server"
$wget -O randoms $inputURL/config/randoms.in 2>/dev/null || clean_exit $? "cannot fetch randoms.in from web server"
$wget -O hd_recon.config $inputURL/config/hd_recon.config 2>/dev/null || clean_exit $? "cannot fetch hd_recon.config from web server"

source ./setup || clean_exit $? "cannot execute setup.sh"

if false; then
    awk '/^TRIG/{print "TRIG  '$nevents'   number of events to simulate";next}\
         /^WROUT/{print "WROUT    1   0   0   hddm output only";next}\
         /^RUNNO/{print "RUNNO    '$runNo'   run number for output events";next}\
         /^NPRIEV/{print "NPRIEV  5       only print first few events";next}\
         /^EPHLIM/{print "EPHLIM  8.0 9.0 energy range in GeV";next}\
         /^RNDMSEQ/{print "RNDMSEQ    '$RSEQ'     random number sequence";next}\
                   {print}' \
         run.ffr.template > fort.15
    bggen || clean_exit $? "bggen failed"
else
    xrdcp $xrootdURL/$remotepath/${simType}_$seqNo.hddm genr8.hddm || clean_exit $? "cannot fetch $xrootdURL/$remotepath/${simType}_$seqNo.hddm"
fi

rseeds=$(head -n $seqNo randoms | tail -n1)
sed -i "s/^RNDM.*/RNDM $rseeds/" control.in || clean_exit $? "RNDM card not found in control.in"
sed -i "s/^RUN[NG].*/RUNNO $runNo/" control.in || clean_exit $? "RUN card not found in control.in"
sed -i "s/^OUTF.*/OUTFILE 'sample.hddm'/" control.in || clean_exit $? "OUTF card not found in control.in"
sed -i "s/^TRIG.*/TRIG $nevents/" control.in || clean_exit $? "TRIG card not found in control.in"

export CCDB_CONNECTION="sqlite:////cvmfs/oasis.opensciencegrid.org/gluex/private/ccdb_2023-04-12.sqlite"
export JANA_CALIB_URL=$CCDB_CONNECTION
export JANA_GEOMETRY_URL="ccdb://GEOMETRY/main_HDDS.xml"
export JANA_CALIB_CONTEXT="variation=mc"
#export JANA_CALIB_CONTEXT="variation=fdcwires_test"

if [ "$simType" = "g4" -o "$simType" = "g4t" -o "$simType" = "eta_pi0_p" -o "$simType" = "eta_pi0_p_x10" ]; then
    simApp=hdgeant4
    cat setup > sim.sh
    echo "which hdgeant4" >> sim.sh
    echo "hdgeant4 -t $nthreads" >> sim.sh
else
    simApp=hdgeant
    cat setup > sim.sh
    echo "export TMPDIR=$(pwd)" >> sim.sh
    echo "which hdgeant" >> sim.sh
    echo "hdgeant -xml=ccdb://GEOMETRY/main_HDDS.xml,run=$runNo" >> sim.sh
fi

bash sim.sh || clean_exit $? "$simApp crashed"
source ./setup
mcsmear -Pprint -PJANA:BATCH_MODE=1 \
        -PNTHREADS=$nthreads \
        -PTHREAD_TIMEOUT_FIRST_EVENT=3600 \
        -PTHREAD_TIMEOUT=600 \
        sample.hddm || clean_exit $? "mcsmear crashed"

hd_root --config=hd_recon.config \
        --nthreads=$nthreads \
        -PJANA:BATCH_MODE=1 \
        -PNTHREADS=$nthreads \
        -PTHREAD_TIMEOUT_FIRST_EVENT=3600 \
        -PTHREAD_TIMEOUT=600 \
        -PTRK:SAVE_TRUNCATED_DEDX=1 \
       sample_smeared.hddm || clean_exit $? "hd_ana crashed"

save_output sample.hddm ${simType}_${seqNo}_geant4.hddm || clean_exit $? "save of sample.hddm failed"
save_output sample_smeared.hddm ${simType}_${seqNo}_smeared.hddm || clean_exit $? "save of sample_smeared.hddm failed"
save_output dana_rest.hddm ${simType}_${seqNo}_rest.hddm || clean_exit $? "save of dana_rest.hddm failed"
save_output hd_root.root ${simType}_${seqNo}_rest.root || clean_exit $? "save of hd_root.hddm failed"
clean_exit
