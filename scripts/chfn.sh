#!/bin/bash
# Changes the extensions of a bunch of files in a directory
dir=$1
startext=$2
newext=$3
for file in ${dir}/*.${startext}; do
	mv -- "${file}" "${file%.${startext}}.${newext}"
done
