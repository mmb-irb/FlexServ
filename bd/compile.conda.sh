#!/bin/bash
"${FC}" -c -O -ffixed-line-length-none bdrandom.f
"${FC}" -O -o bd -ffixed-line-length-none bd2.f bdrandom.o
