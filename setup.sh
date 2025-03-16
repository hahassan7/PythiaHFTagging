# need to have correct $ROOTSYS, $FASTJET and $PYTHIA8 are needed in Makefile
# In case you have installed all from Aliroot installation including fastjet packages based on https://dberzano.github.io/alice/install-aliroot/
# no need to setup $ROOTSYS, $FASTJET need to have only $PYTHIA8
export GCC_TOOLCHAIN_ROOT=/scratch/project_2003583/focal-full-sim/alicesw/sw/slc8_x86-64/GCC-Toolchain/latest/
export LD_LIBRARY_PATH=$GCC_TOOLCHAIN_ROOT/lib:$GCC_TOOLCHAIN_ROOT/lib64${LD_LIBRARY_PATH+:$LD_LIBRARY_PATH}
export PATH=$GCC_TOOLCHAIN_ROOT/bin:$GCC_TOOLCHAIN_ROOT/libexec/bin${PATH+:$PATH}

source /scratch/project_2003583/focal-full-sim/alicesw/sw/slc8_x86-64/ROOT/latest/bin/thisroot.sh
#source /usr/local/bin/thisroot.sh
export CGAL_ROOT=/scratch/project_2003583/focal-full-sim/alicesw/sw/slc8_x86-64/cgal/latest
export PATH=${PATH}:$CGAL_ROOT/bin
export GMP_ROOT=/scratch/project_2003583/focal-full-sim/alicesw/sw/slc8_x86-64/GMP/latest
export PATH=${PATH}:$GMP_ROOT/bin
export PYTHIA8=/scratch/project_2003583/focal-full-sim/alicesw/sw/slc8_x86-64/pythia/latest
export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc
export PATH=${PATH}:$PYTHIA8/bin
export FASTJET=/scratch/project_2003583/focal-full-sim/alicesw/sw/slc8_x86-64/fastjet/latest
export PATH=$PATH:${FASTJET}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$PYTHIA8/lib:$FASTJET/lib:$CGAL_ROOT/lib:$GMP_ROOT/lib
