#!/bin/bash
export ROOTSYS=/cvmfs/oasis.opensciencegrid.org/gluex/root-6.22.06/x86_64
if [[ -d $ROOTSYS ]]; then
    export PATH=$ROOTSYS/bin:$PATH
    export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH
fi
