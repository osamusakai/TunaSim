echo off
echo Getting Fishing mortality rates
del Fmort1.dat
for /L %%i in (1,1,%2) do awk "NR==163{x=$0} NR==245{y=$0} END{print x,y} "  %1%%i.rep >>Fmort1.dat

del Fmort2.dat
for /L %%i in (1,1,%2) do awk "NR==241{x=$0} NR==283{y=$0} END{print x,y} "  %1%%i.rep >>Fmort2.dat

