echo off
@echo Args: run number, type, and n simulations (e.g., extract 1 r 50)
rem pause
unzip -j rresults.zip %1%2*.rep
call getbiom %1%2 50
call getrbiom %1%2 50
call getsel %1%2 50
call getf %1%2 50
call getq %1%2 50
call getage %1%2 50
call getrec %1%2 50
del %1%2*.rep

unzip -j rresults.zip %1%2*.par
call getgrwth %1%2 50
del %1%2*.par

call gettrue %1 %2 
