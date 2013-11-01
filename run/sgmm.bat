REM gmult catsim.frq init.fit 01.fit
REM  estimate the age -dependent variance
REM echo Stage 1
REM gmult1 catsim.frq 01.fit 02.fit -switch 1 1 16 1 >nul
REM echo Stage 2
estimate selectivities as size-based
gmult1 catsim.frq 04.fit 05.fit -switch 1 -999 26 2 >nul
REM estimate M
REM echo Stage 4
gmult1 catsim.frq 05.fit 06.fit -switch 1 2 33 1 >nul
REM echo Stage 5
gmult1 catsim.frq 06.fit final.fit -switch 1  1 190 1 >nul
echo Done.... 
