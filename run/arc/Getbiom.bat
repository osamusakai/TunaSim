echo off
echo Extracting biomass from rep files
del biom.dat
rem for %%i in (1r*.rep) do awk "NR>=330&&NR<370{print NR-329, $0}" %%i >> biom.dat
for /L %%i in (1,1,%2) do awk "NR>=330&&NR<370{print NR-329, $0}" %1%%i.rep >> biom.dat
rem FOR %%i IN (1r*.rep) DO @echo %%i
REM awk "NR>=330&&NR<370{print NR-329, $0}" %i >> biom.dat
