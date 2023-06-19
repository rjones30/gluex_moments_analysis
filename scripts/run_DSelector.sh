#!/bin/bash

function usage {
    echo "usage: run_DSelector.sh <input.root> <output.root>"
}

if [[ $# != 2 ]]; then
    usage
    exit 1
fi

workdir=DSelector_work_$$
mkdir $workdir
cp src/DSelector_etapi0_moments.[Ch] $workdir
cd $workdir
xrdcp $1 input.root
cat <<EOI >doit.C
{
    TFile *finput = TFile::Open("input.root");
    gROOT->ProcessLine(".x ../src/Load_DSelector.C");
    gROOT->ProcessLine(".L DSelector_etapi0_moments.C++");
    gROOT->ProcessLine("etapi0__B4_T1_S1_M7_M17_F4_Tree->Process(\"DSelector_etapi0_moments\")");
    gROOT->ProcessLine(".q)");
}
EOI
root -l -b doit.C
gfal-copy -f etapi0_moments.root https://cn440.storrs.hpc.uconn.edu:2843/Gluex/resilient/simulation/moments-6-2023/$2
