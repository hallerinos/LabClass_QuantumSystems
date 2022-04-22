#!/bin/bash
pandoc --metadata-file=./pandoc_latex.yml ./$1.md -o $1.pdf