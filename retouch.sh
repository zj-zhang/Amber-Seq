#!/bin/bash
# re-touch files by generating order
# so to avoid snakemake re-run from start
# ZZ
# 2020.3.23
project=$1
touch outputs/$project/config/*
find outputs/$project/nas_search/ -exec touch {} \;
find outputs/$project/nas_final/ -exec touch {} \;
find outputs/$project/vep/ -exec touch {} \;
find outputs/$project/asb/* -exec touch {} \;
