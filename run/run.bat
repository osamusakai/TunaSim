REM run simulation saver stuff for mfcl
rem invoke simulator
catsim %2 %3
rem invoke analyser
copy catsim.out catsim.frq
call gm  
rem invoke archiver
call arcit %1

