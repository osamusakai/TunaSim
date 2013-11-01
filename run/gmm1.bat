REM  estimate the age -dependent variance
REM echo Stage 1
gmult1 catsim.frq 01.fit 02.fit -switch 1 1 16 1 >nul
REM echo Stage 2
REM estimate K
gmult1 catsim.frq 02.fit 03.fit -switch 1 1 14 1 >nul
REM estimate selectivities as size-based
gmult1 catsim.frq 03.fit 04.fit -switch 1 -999 26 2 >nul
REM estimate M
REM echo Stage 4
gmult1 catsim.frq 04.fit 05.fit -switch 1 2 33 1 >nul
REM echo Stage 5
gmult1 catsim.frq 05.fit final.fit -switch 1  1 190 1 >nul
echo Done.... 
