echo off
echo Selectivity for fishery 1
del sel.dat
for /L %%i in (1,1,%2) do awk "NR==12{print $0}"  %1%%i.rep >>  sel.dat
del sel2.dat
for /L %%i in (1,1,%2) do awk "NR==94{print $0}"  %1%%i.rep >>  sel2.dat
