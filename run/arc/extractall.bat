echo off
@echo Args: run number, type, and n simulations (e.g., extract 1 r 50)
pause
rem unzip -j %2results.zip *%2*.rep

del biom.dat
for /L %%i in (1,1,13) do for /L %%j in (1,1,50) do awk "NR>=330&&NR<370{print "%%i%2",NR-329, $0}" %%ir%%j.rep >> biom.dat
unzip -j rresults.zip *trbiom.dat
del tbiom.dat
for /L %%i in (1,1,%2) do awk "{print %%i,$1,$2}" %%i%2tbiom.dat >> tbiom.dat
pause 
copy %1rtbiom.dat tbiom.dat

rem call getrbiom %1%2 50
rem call getsel %1%2 50
rem call getf %1%2 50
rem call getq %1%2 50
Rem call getage %1%2 50
rem call getrec %1%2 50
Rem del %1%2*.rep

Rem unzip -j %2results.zip %1%2*.par
Rem call getgrwth %1%2 50
Rem del %1%2*.par

Rem call gettrue %1 %2 

Rem echo off
Rem echo Extracting biomass from rep files
Rem del biom.dat
rem for %%i in (1r*.rep) do awk "NR>=330&&NR<370{print NR-329, $0}" %%i >> biom.dat
Rem for /L %%i in (1,1,%2) do awk "NR>=330&&NR<370{print NR-329, $0}" %1%%i.rep >> biom.dat
rem FOR %%i IN (1r*.rep) DO @echo %%i
REM awk "NR>=330&&NR<370{print NR-329, $0}" %i >> biom.dat




REM write true data files...
unzip -j rresults.zip %1rt*.dat
copy %1rtend.dat tend.dat
copy %1rtinit.dat tinit.dat
copy %1rtrec.dat trec.dat
copy %1rtbiom.dat tbiom.dat
copy %1rtrbiom.dat trbiom.dat
copy %1rtq1.dat tq1.dat
copy %1rtq2.dat tq2.dat
copy %1rtf1.dat tf1.dat
copy %1rtf2.dat tf2.dat
del %1rt*.dat
