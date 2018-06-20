
echo "Setting up the OPT data analysis environment ..."

export OPTDIR="/micehome/jwalls/OPT_8.04"
export PYTHONPATH="$PYTHONPATH:$OPTDIR/opt_tools/:/$OPTDIR/pyCTSim/"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$OPTDIR/libFDR/"
export PATH="$PATH:$OPTDIR/opt_tools/:$OPTDIR/libFDR/"