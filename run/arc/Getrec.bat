echo off
echo Getting recruitment estimates...
del rec.dat
for /L %%i in (1,1,%2) do awk "NR>=480&&NR<520{print $0}" %1%%i.rep >> rec.dat

@echo deleting unzipped rep files
rem del %1*.rep
