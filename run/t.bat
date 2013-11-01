
rem for %%i in (1r*.rep) do awk "NR>=330&&NR<370{print NR-329, $0}" %%i >> biom.dat
for /L %%i in (1,1,100) do awk "NR>=330&&NR<370{print NR-329, $0}" %1%%i.rep >> biom.dat
rem FOR %%i IN (1r*.rep) DO @echo %%i
REM awk "NR>=330&&NR<370{print NR-329, $0}" %i >> biom.dat
rem catsim -iseed 1232 -ind %1%2catsim.dat

rem call arctrue %1%2
