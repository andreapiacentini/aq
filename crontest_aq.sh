#!/bin/bash -l
if [ -z $AQ_JEDI_CRO ] ; then
  if [ -f $HOME/loadaq_crontest.sh ] ; then
    source $HOME/loadaq_crontest.sh
    echo "AQ CRONTAB TEST: loaded execution environment"
  else
    echo "AQ CRONTAB TEST: execution environment not loaded"
    exit 1
  fi
fi

if [ -z $AQ_JEDI_CRO ] ; then
  echo "AQ CRONTAB TEST: execution environment not set"
  exit 1
fi
\rm -fr $AQ_JEDI_CRO/*
\rm -fr $AQ_JEDI_BLD/*
\rm -fr $AQ_JEDI_DIR/*
cd $AQ_JEDI_SRC
git pull  > $AQ_JEDI_CRO/git_output.out 2>&1
\rm -fr bundle/*
git checkout -- bundle/CMakeLists.txt >> $AQ_JEDI_CRO/git_output.out 2>&1
git status >> $AQ_JEDI_CRO/git_output.out 2>&1
./install_aq.ksh > $AQ_JEDI_CRO/install_output.out 2>&1
install_ret=$?
if [ $install_ret = 0 ] ; then
   cd $AQ_JEDI_BLD/aq
   if [ $AQ_JEDI_SUBMIT = slurm ] ; then
       rm -f HPC_JOB HPC_OUT
       cat > HPC_JOB <<EOF
#!/bin/bash
#SBATCH -J aq-crontest
#Number of nodes
#SBATCH -N 1
#SBATCH -p debug
#SBATCH --ntasks-per-node=36
#SBATCH -o HPC_OUT
#SBATCH -e HPC_OUT

source $HOME/loadaq_crontest.sh
export KMP_DETERMINISTIC_REDUCTION=true
export I_MPI_PIN_DOMAIN=omp
export ROMIO_HINTS=/softs/Modules/modulefiles/mpi/intelmpi/romio-hints
export I_MPI_HYDRA_BOOTSTRAP=ssh
cd $AQ_JEDI_BLD/aq
ctest --output-on-failure > $AQ_JEDI_CRO/ctest_output.out 2>&1
EOF
       sbatch HPC_JOB
       g_ret=-99
       while [ $g_ret -ne 0 ] ; do
           grep 'tests passed' $AQ_JEDI_CRO/ctest_output.out > /dev/null 2>&1
           g_ret=$?
           sleep 5
       done
       ret_st=`grep 'tests passed' $AQ_JEDI_CRO/ctest_output.out | awk '{print $4}'`
   else
       ctest --output-on-failure > $AQ_JEDI_CRO/ctest_output.out 2>&1
       ret_st=$?
   fi
   if [ $ret_st -ne 0 ] ; then
       status=Failed
       message="Failed on ctest"
   else
       status=Succeeded
       message="OK"
   fi
else
   status=Failed
   message="Failed on install"
fi
mailx -s "AQ automatic testing $status at `date`" $1 <<EOF
$message
EOF
