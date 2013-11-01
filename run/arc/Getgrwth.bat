echo getting growth parameters 
del growth.dat
for /L %%i in (1,1,%2) do awk "NR==50{print 4, $1} NR==284{print 1, $1} NR==285{print 2, $1} NR==286{print 3,$1}" %1%%i.par >> growth.dat

echo deleting unzipped par files
rem del %1*.par
