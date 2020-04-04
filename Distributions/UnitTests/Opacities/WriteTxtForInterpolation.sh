#!/bin/bash
# --- check if tables exists and write out dataList
if [ -e "dataList.txt" ]; then
rm dataList.txt
fi
for f in *.profile; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
for f in wl-EOS-*.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
for f in *AbEm.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
for f in *Iso.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
for f in *NES.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
for f in *Pair.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
