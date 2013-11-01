catsim.dat - simulator scenario (input data for simulation)
catsim.frq - simulated data set
graph.frq  - fit to catsim.frq data set
04.fit     - fit file

main.bat calls 8 different model scenarios 
{
   copy scenario file (e.g., 1rcatsim.dat) to catsim.dat for basis of simulator
   runallm.bat calls 100 simulations passes names, random number seeds as args.
   {
     run.bat does the simulation/estimation loop 
     {
       call simulator (calls catsim.exe with rn seed etc)
       copy new data file to catsim.frq
       call analyze.bat (using catsim.frq)
       {
         do mfcl steps for consecutive fits
         archive results (arcit.bat; note need subdirectory "arc")
       }
   }
   zip up this scenario (packit.bat)
}
