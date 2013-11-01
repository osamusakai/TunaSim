rem gmult1 catsim.out t2.fit init.fit rem uncomment this to get startup
rem gmult1 catsim.out init.fit 01.fit
REM  estimate the age -dependent variance
rem gmult1 catsim.frq 01.fit 02.fit -switch 1 1 16 1
REM estimate K
rem gmult1 catsim.frq 02.fit 03.fit -switch 1 1 14 1
rem Make some output
gmult1 catsim.frq 03.fit 04.fit -switch 2 1 189 1 1 190 1
REM estimate M
rem gmult1 catsim.frq 03.fit 04.fit -switch 1 2 33 1


