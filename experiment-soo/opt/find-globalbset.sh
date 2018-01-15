#!/bin/bash
log=`find ./ -name debug.log`
for i in $log;
do
et=`grep 'Updated global best:'   $i |tail -1 |awk '{printf "%s\n", $8}'`
echo $i   $et
done

