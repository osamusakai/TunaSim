REM gmult catsim.frq init.fit 01.fit
REM  estimate the age -dependent variance
echo Stage 1
gmult1 catsim.frq 01.fit 02.fit -switch 1 1 16 1 >nul
echo Stage 1
REM estimate K
gmult1 catsim.frq 02.fit 03.fit -switch 1 1 14 1 >nul
echo Stage 2
gmult1 catsim.frq 03.fit final.fit -switch 1 1 190 1 >nul
echo Done.... 
REM estimate M
REM gmult1 catsim.frq 03.fit 04.fit -switch 1 2 33 1
