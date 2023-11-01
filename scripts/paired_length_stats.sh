#!/bin/bash

zcat "$1" | awk 'BEGIN{x=""; y=0; }{if (NR%4==1) {x=$1}; if (NR%4==2) {y=length($0); print x, y;}}' | sort -t " " -k1,1 > "$1"".txt"

zcat "$2" | awk 'BEGIN{x=""; y=0; }{if (NR%4==1) {x=$1}; if (NR%4==2) {y=length($0); print x, y;}}' | sort -t " " -k1,1 > "$2"".txt"

join "$1"".txt" "$2"".txt" | cut -d " " -f 2-3
