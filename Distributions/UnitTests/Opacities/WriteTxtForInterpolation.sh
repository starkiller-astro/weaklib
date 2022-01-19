#!/bin/bash
# --- check if tables exists and write out dataList
if [ -e "dataList.txt" ]; then
rm dataList.txt
fi
if [ -e *.profile ]; then
  for f in *.profile; do
    [ -e "$f" ] && echo $f >> dataList.txt
  done
else
  echo 'No profile!'>> dataList.txt
fi
for f in wl-EOS-*.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
for f in wl-*Ab*.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
for f in wl-*Iso.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
for f in wl-*NES.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
for f in wl-*Pair.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
for f in wl-*Brem.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
