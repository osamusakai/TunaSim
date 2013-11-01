echo off 
REM echo # Sample size >t
REM echo  400 >>t
REM sed -f sedcmd %1catsim.dat >>t
REM del %1catsim.dat
REM move t %1catsim.dat
rem for /L %%i in (1,1,30) do awk "NR>=666&&NR<706{print NR-665, $4}" %1%%i.rep >> f1q.dat
FOR %%i IN (*ccatsim.dat) DO call sedit %%i
REM FOR %%i IN (*rcatsim.dat) DO dir %%i
REM FOR %%i IN (*mcatsim.dat) DO call sedit %%i
REM FOR %%i IN (*ccatsim.dat) DO call sedit %%i
REM awk "NR>=330&&NR<370{print NR-329, $0}" %i >> biom.dat
