REM echo # Catch_sd    >t
REM echo  0 >>t
echo # Rho_sd    >t
echo  0 >>t
REM sed -f sedcmd %1catsim.dat >>t
cat %1 >>t
del %1
move t %1
