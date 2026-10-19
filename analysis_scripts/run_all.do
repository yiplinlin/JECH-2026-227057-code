version 19.5
args input out
if `"`input'"'=="" local input "data/data.dta"
if `"`out'"'=="" local out "results"
confirm file "`input'"
capture mkdir "`out'"
foreach package in ftools reghdfe ppmlhdfe boottest {
    which `package'
}
do "analysis_scripts/Table_2_and_Table_S13.do" "`input'" "`out'"
foreach script in primary_and_event_studies cutoffs_and_pretrends country_inclusion cohort_trends_and_age callaway_santanna_overall poisson_heterogeneity poisson_linear_comparison poisson_bounded_variation {
    do "analysis_scripts/`script'.do" "`input'" "`out'"
}
