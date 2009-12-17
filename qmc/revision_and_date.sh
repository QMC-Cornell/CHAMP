     svn_revision_number=`svn info |grep "Revision"`
     svn_date=`svn info |grep "Last Changed Date"`
     compilation_date=`date`
     compilation_hostname=`hostname -f`
     echo " write(6,'(a)') 'SVN $svn_revision_number, $svn_date'" >revision_and_date.inc
     echo " write(6,'(a)') 'Compiled by $USER on $compilation_date on host $compilation_hostname'" >>revision_and_date.inc
