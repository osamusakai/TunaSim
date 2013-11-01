rem invoke simulator
catsim %2 %3
rem invoke analyser
copy catsim.out catsim.frq
rem invoke estimator with M part....
call gmd1  
rem invoke archiver
call arcit %1

