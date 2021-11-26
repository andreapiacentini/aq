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
   ctest --output-on-failure > $AQ_JEDI_CRO/ctest_output.out 2>&1
   if [ $? -ne 0 ] ; then
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

