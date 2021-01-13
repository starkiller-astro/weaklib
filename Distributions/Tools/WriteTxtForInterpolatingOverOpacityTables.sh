#!/bin/bash
# --- check if tables exists and write out dataList
if [ -e "dataList.txt" ]; then
rm dataList.txt
fi

if [ -e wl-EOS-*.h5 ]; then
for f in wl-EOS-*.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
else
  echo 'EOS' >> dataList.txt
fi

if [ -e *AbEm.h5 ]; then
for f in *AbEm.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
else
  echo 'AbEm' >> dataList.txt
fi

if [ -e *Iso.h5 ]; then
for f in *Iso.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
else
  echo 'Iso' >> dataList.txt
fi

if [ -e *NES.h5 ]; then
for f in *NES.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
else
  echo 'NES' >> dataList.txt
fi

if [ -e *Pair.h5 ]; then
for f in *Pair.h5; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
else
  echo 'Pair' >> dataList.txt
fi

if [ -e "profilelist.txt" ]; then
   echo "profilelist.txt" >> dataList.txt
elif [ -e *.profile ]; then
for f in *.profile; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
elif [ -e *_hdf5_chk* ]; then
# --- check if flash checkout file exists and write out dataList
for f in *_hdf5_chk*; do
  [ -e "$f" ] && echo $f >> dataList.txt
done
else
  echo 'pro' >> dataList.txt
fi
