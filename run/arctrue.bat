copy tbiom.dat arc\%1tbiom.dat
copy trbiom.dat arc\%1trbiom.dat
copy trec.dat arc\%1trec.dat
copy ttot.dat arc\%1ttot.dat
copy tinit.dat arc\%1tinit.dat
copy tend.dat arc\%1tend.dat
copy tF1.dat arc\%1tF1.dat
copy tF2.dat arc\%1tF2.dat

zip -j arc\rresults.zip arc\%1t*.dat


