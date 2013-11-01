echo off
echo Extracting relative rbiomass from rep files
del rbiom.dat
for /L %%i in (1,1,%2) do awk "NR>=371&&NR<411{print NR-370, $0}" %1%%i.rep >> rbiom.dat
