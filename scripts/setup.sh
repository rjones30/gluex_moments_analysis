# add pointers to Geant4 libraries
export G4ROOT=/cvmfs/oasis.opensciencegrid.org/gluex/geant4.10.04.p02/x86_64
export GEANT4PY=$G4ROOT/src/environments/g4py
export G4WORKDIR=/cvmfs/oasis.opensciencegrid.org/gluex/HDGeant4/jlab
if [[ -x $G4ROOT/bin/geant4-site.sh ]]; then
    . $G4ROOT/bin/geant4-site.sh >/dev/null
else
    . $G4ROOT/bin/geant4.sh >/dev/null
fi
if [[ -x $G4ROOT/share/Geant4*/geant4make/geant4make-site.sh ]]; then
    . $G4ROOT/share/Geant4*/geant4make/geant4make-site.sh
else
    . $G4ROOT/share/Geant4*/geant4make/geant4make.sh
fi
if [[ -d $GEANT4PY/lib64 ]]; then
    PYTHONPATH=$PYTHONPATH:$GEANT4PY/lib64:$G4WORKDIR/g4py
else
    PYTHONPATH=$PYTHONPATH:$GEANT4PY/lib:$G4WORKDIR/g4py
fi
export DIRACXX=/cvmfs/oasis.opensciencegrid.org/gluex/Diracxx
export HDF5ROOT=/cvmfs/oasis.opensciencegrid.org/gluex/hdf5-1.12.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5ROOT/lib
export LD_LIBRARY_PATH=$DIRACXX:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$G4ROOT/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$G4WORKDIR/tmp/Linux-g++/hdgeant4:$LD_LIBRARY_PATH
export PATH=$G4WORKDIR/bin/Linux-g++:$PATH

export XERCESCROOT=/group/halld/Software/builds/Linux_CentOS7-x86_64-gcc4.8.5-cntr/xerces-c/xerces-c-3.1.4
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$XERCESCROOT/lib

export OSNAME=Linux_CentOS7-x86_64-gcc4.8.5

# add pointers to halld_sim
export HALLD_SIM_HOME=/cvmfs/oasis.opensciencegrid.org/gluex/halld_sim
export JANA_PLUGIN_PATH=$HALLD_SIM_HOME/$OSNAME/plugins:$JANA_PLUGIN_PATH
export LD_LIBRARY_PATH=$HALLD_SIM_HOME/$OSNAME/lib:$LD_LIBRARY_PATH
export PATH=$HALLD_SIM_HOME/$OSNAME/bin:$PATH

# add pointers to halld_recon
export HALLD_RECON_HOME=/cvmfs/oasis.opensciencegrid.org/gluex/halld_recon
export JANA_PLUGIN_PATH=$HALLD_RECON_HOME/$OSNAME/plugins
export LD_LIBRARY_PATH=$HALLD_RECON_HOME/$OSNAME/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$HALLD_RECON_HOME/$OSNAME/python3:$HALLD_RECON_HOME/$OSNAME/python2:$PYTHONPATH
export PATH=$HALLD_RECON_HOME/$OSNAME/bin:$PATH

# override the sqlite file if needed
export RCDB_CONNECTION=sqlite:////group/halld/www/halldweb/html/dist/rcdb.sqlite
#export CCDB_CONNECTION=sqlite:////cvmfs/oasis.opensciencegrid.org/gluex/private/ccdb_2021-3-11.sqlite
export JANA_CALIB_URL=$CCDB_CONNECTION

# configure the version of ROOT
export ROOTSYS=/cvmfs/oasis.opensciencegrid.org/gluex/root-6.22.06/x86_64-debug
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH

# configure ccdb
export CCDB_ROOT=/cvmfs/oasis.opensciencegrid.org/gluex/group/halld/Software/builds/Linux_CentOS7-x86_64-gcc4.8.5-cntr/ccdb/ccdb_1.06.07
export LD_LIBRARY_PATH=$CCDB_ROOT/lib:$LD_LIBRARY_PATH

# configure evio
EVIOROOT=/cvmfs/oasis.opensciencegrid.org/gluex/group/halld/Software/builds/Linux_CentOS7-x86_64-gcc4.8.5-cntr/evio/evio-4.4.6/Linux-x86_64
export LD_LIBRARY_PATH=$EVIOROOT/lib:$LD_LIBRARY_PATH

# configure jana
export JANA_ROOT=/cvmfs/oasis.opensciencegrid.org/gluex/JANA
export LD_LIBRARY_PATH=$JANA_ROOT/lib:$LD_LIBRARY_PATH
export JANA_PLUGIN_PATH=$JANA_PLUGIN_PATH:$JANA_ROOT/$OSNAME/plugins

if [ "$jobType" = "mergehists" ]; then
    export XROOTD_HOME=/cvmfs/oasis.opensciencegrid.org/gluex/xrootd/3.3.2
    export LD_LIBRARY_PATH=$ROOTSYS/lib:$XROOTD_HOME/lib64:$LD_LIBRARY_PATH
    export LD_PRELOAD=$XROOTD_HOME/lib64/libXrdPosixPreload.so:$XROOTD_HOME/lib64/libXrdPosix.so:$XROOTD_HOME/lib64/libXrdClient.so:$XROOTD_HOME/lib64/libXrdUtils.so
fi
