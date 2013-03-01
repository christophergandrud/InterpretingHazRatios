/////////////////////
// Export Replication Data Set From Gandrud (2012) that can be used with simPH
// Christopher Gandrud
// 1 March 2013
// Based on original source code available at: https://raw.github.com/christophergandrud/FinRegGov/master/SourceCode/PredictedHazardsCIF.do
// Output: FinSurvData.csv
/////////////////////


//// Download data ////
use "http://dl.dropbox.com/u/12581470/code/Replicability_code/Financial_Supervision_Governance_Replication/public_fin_trans_data.dta", clear


///// Remove 1987 due to SE lag missing /////
drop if year == 1987
   
//// Convert spatial effects into percents to ease interpretation ////      
gen percent_se_cbss_ocbu = se_cbss_ocbu*100
gen percent_se_eu_ocbu = se_eu_ocbu*100
gen percent_se_basel_ocbu = se_basel_ocbu*100
    
//// Extract 1st imputed data set ////
mi extract 1, clear

//// stset the data ////
stset year, id(country) failure(reg_4state == 3) enter(reg_4state == 2) exit(reg_4state == 3 4) origin(min)

//// Rename Surv variables ////
rename _t0 begin
rename _t end
rename _d event
rename _st NoEvent

//// Save to FinSurvData.csv /////  
outsheet using "FinSurvData.csv", comma replace
