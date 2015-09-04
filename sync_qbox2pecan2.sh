#!/bin/sh


logdir='/Work/Tech/rsync/logs/'
dryrun=''
backup='b'
flags=""
exclude=""
src="/Work/Research/Macrosystems/ed2-clark_compare/"
dst="pecan2:/fs/data2/rykelly/ed2-clark_compare"
backupdir='/fs/data2/rykelly/rsync/backups/'

# Exclude flag for edoutputs dirs (will be applied to either sync direction)
exclude_outputs="--exclude */0.Junk/ --exclude */1.Saved_nosync --exclude */tempdelete/ --exclude */tempsave/"



#------ Generate timestamp
timestamp=`date +%Y%m%d_%H%M%S`
printf "Syncing \n"
printf "  $src to \n"
printf "  $dst \n"
printf "Timestamp is: \n"
printf "  $timestamp"
printf "\n************************************************************\n\n"
printf "Rsync transcript:\n"


#------ Run rsync 
rsync -azvP${dryrun}${backup} --delete \
  --backup-dir=${backupdir}${timestamp} \
  --log-file=${logdir}${timestamp}.log \
  ${exclude} \
  "${src}" "${dst}"