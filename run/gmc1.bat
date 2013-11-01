REM gmult catsim.frq init.fit 01.fit
REM  estimate the age -dependent variance
gmult1 catsim.frq 01.fit 02.fit -switch 1 1 16 1 >nul
REM estimate K
gmult1 catsim.frq 02.fit 03.fit -switch 1 1 14 1 >nul
gmult1 catsim.frq 03.fit 04.fit -switch 1 -999 10 1 >nul
gmult1 catsim.frq 04.fit final.fit -switch 1 1 190 1 >nul
echo Done.... 
REM estimate M
REM gmult1 catsim.frq 03.fit 04.fit -switch 1 2 33 1
