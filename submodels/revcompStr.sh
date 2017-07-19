#!/bin/bash
#assume argument is one line of ATGC sequence.

#fwd=`cat $1 | tr -d '\r'`
fwd=$1
rvs=`echo $fwd | rev`
tmp1=`echo $rvs | sed 's/[Aa]/X/g'`
tmp2=`echo $tmp1 | sed 's/[Gg]/Y/g'`
rvcmp_4=`echo $tmp2 | sed 's/[Tt]/A/g'`
rvcmp_3=`echo $rvcmp_4 | sed 's/[Cc]/G/g'`
rvcmp_2=`echo $rvcmp_3 | sed 's/X/T/g'`
rvcmp=`echo $rvcmp_2 | sed 's/Y/C/g'`
echo $rvcmp
