     git_commit_number=`git log -1|grep commit`
     git_date=`git log -1|grep Date`
     compilation_date=`date`
     compilation_hostname=`hostname -f`
     echo " write(6,'(a)') 'GIT $git_commit_number, $git_date'" >revision_and_date.inc
     echo " write(6,'(a)') 'Compiled by $USER on $compilation_date on host $compilation_hostname'" >>revision_and_date.inc
