echo off
echo Extracting catchabilities...
del f1q.dat
for /L %%i in (1,1,%2) do awk "NR>=666&&NR<706{print NR-665, $4}" %1%%i.rep >> f1q.dat
del f2q.dat
for /L %%i in (1,1,%2) do awk "NR>=707&&NR<727{print NR-706+20, $4}" %1%%i.rep >> f2q.dat
rem FOR %%i IN (1r*.rep) DO @echo %%i
REM awk "NR>=330&&NR<370{print NR-329, $0}" %i >> biom.dat
