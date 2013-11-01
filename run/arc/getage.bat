echo off
@echo Getting Age comps in first and last year
del Age1.dat
for /L %%i in (1,1,%2) do awk "NR==287{print $0}"  %1%%i.rep >>Age1.dat

del Age2.dat
for /L %%i in (1,1,%2) do awk "NR==326{print $0}"  %1%%i.rep >> Age2.dat

