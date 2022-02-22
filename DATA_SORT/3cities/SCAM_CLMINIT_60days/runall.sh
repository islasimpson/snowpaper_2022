#!/bin/bash
for ifile in `ls *_FLDS_*.py` ; do
  python $ifile
done
