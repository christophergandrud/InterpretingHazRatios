/////////////////
// Preprocess Simmons et al. (2006) data set for replication
// Christopher Gandrud 
// 28 April 2013
/////////////////

///// Load Data /////
use "~/Dropbox/Hazard_Ratio_Research_Note/bits_io2006_rev5_2008/bits_io2006_rev5_08.dta"

//// stset the data using instructions from bit_io2006.readme.txt ////
stset year, id(dyad)  failure(bit) origin(time atrisk)

//// Rename stata Surv variables ////
rename _t0 begin
rename _t end
rename _d event
rename _st NoEvent

//// Save to csv file ////
outsheet using "SimmonIO2006.csv", comma replace