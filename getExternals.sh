#!/bin/bash

[ -x ThirdParty ] || mkdir ThirdParty 
[ -x Data ]       || mkdir Data 

awk '{printf("echo; echo Downloading %s; svn co %s %s\n", $1, $2, $1);}' <Externals > svn_externals.sh

chmod 700 svn_externals.sh && \
sh svn_externals.sh && \
rm svn_externals.sh
