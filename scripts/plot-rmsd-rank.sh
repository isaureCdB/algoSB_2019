#!/bin/bash -i
d=`dirname $0`

name=${1%%.*}

gracebat -hdevice PNG -autoscale none -par $d/plot-rmsd-rank.par \
 $1 -world 0.9 0.1 50000 100 -printfile plot-$name.png
