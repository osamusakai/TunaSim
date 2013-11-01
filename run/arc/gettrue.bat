
REM write true data files...
unzip -j rresults.zip %1%2t*.dat
copy %1%2tend.dat tend.dat
copy %1%2tinit.dat tinit.dat
copy %1%2trec.dat trec.dat
copy %1%2tbiom.dat tbiom.dat
copy %1%2trbiom.dat trbiom.dat
copy %1%2tq1.dat tq1.dat
copy %1%2tq2.dat tq2.dat
copy %1%2tf1.dat tf1.dat
copy %1%2tf2.dat tf2.dat
del %1%2t*.dat
