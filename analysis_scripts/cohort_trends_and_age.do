version 19.5
set more off
set linesize 180
set varabbrev off
set processors 1
args input folder
if `"`input'"'=="" local input "data/data.dta"
if `"`folder'"'=="" local folder "results"
local reps 5000
capture mkdir "`folder'"
local folder "`folder'/cohort_trends_and_age"
capture mkdir "`folder'"
global r2out "`folder'"
global r2reps `reps'
capture log close _all
log using "`folder'/linear_diagnostics.txt", text replace
display "R2_LINEAR_START REPS=`reps' SEED=20260708 ORIGINAL_TABLE2_OUTCOMES"
tempname R C E A
postfile `R' str12 outcome str30 spec double N clusters b se ci_lo ci_hi p_cluster p_wcb reps feasible ci_parts using "`folder'/linear_results.dta", replace
postfile `C' str12 outcome str30 spec int part double lo hi using "`folder'/wild_intervals.dta", replace
postfile `E' str12 outcome str8 window int event double b se byte omitted using "`folder'/event_estimability.dta", replace
postfile `A' str12 outcome double age_residual_sd age_residual_maxabs linear_age_b_delta rank_base rank_with_linear_age using "`folder'/age_identification.dta", replace
global r2res `R'
global r2ci `C'
capture program drop report_model
program define report_model
args outcome spec
tempname n g b se df p lo hi ci wp nr nf dist
scalar `n'=e(N)
scalar `g'=e(N_clust)
scalar `b'=_b[exposed_before15]
scalar `se'=_se[exposed_before15]
scalar `df'=e(df_r)
assert `se'>0 & `se'<.
scalar `p'=2*ttail(`df',abs(`b'/`se'))
scalar `lo'=`b'-invttail(`df',.025)*`se'
scalar `hi'=`b'+invttail(`df',.025)*`se'
estimates save "$r2out/fit_`outcome'_`spec'.ster", replace
display "MODEL_START OUTCOME=`outcome' SPEC=`spec' N=" `n' " CLUSTERS=" `g'
boottest exposed_before15, cluster(country_id) reps($r2reps) seed(20260708) weighttype(rademacher) nograph svmat
return list
scalar `wp'=r(p)
scalar `nr'=r(reps)
matrix `dist'=r(dist)
mata: st_numscalar("`nf'", sum(st_matrix("`dist'"):<.))
matrix `ci'=r(CI)
local parts=rowsof(`ci')
assert `nr'==$r2reps
assert `nf'==$r2reps
forvalues j=1/`parts' {
post $r2ci ("`outcome'") ("`spec'") (`j') (`ci'[`j',1]) (`ci'[`j',2])
}
post $r2res ("`outcome'") ("`spec'") (`n') (`g') (`b') (`se') (`lo') (`hi') (`p') (`wp') (`nr') (`nf') (`parts')
display "MODEL_COMPLETE OUTCOME=`outcome' SPEC=`spec' B=" %16.10f `b' " P_WCB=" %12.8f `wp'
end
use "`input'", clear
keep if primary_sample==1 & inlist(data_part,1,2)
gen double age_calendar=survey_year-female_birth_cohort
assert age==age_calendar if data_part==1
gen double cohort_time=female_birth_cohort-1990
gen int evcat=event15_bin+11
tempfile source
save `source'
foreach outcome in union18 birth18 birth20 births infant under5 {
use `source', clear
local y union_before18
local sample analytic_union
local part 1
local wt weight
local wtype pw
if "`outcome'"=="birth18" {
local y birth_before18
local sample analytic_birth18
}
if "`outcome'"=="birth20" {
local y birth_before20
local sample analytic_birth20
}
if inlist("`outcome'","births","infant","under5") {
local part 2
local wt weighted_women_n
local wtype aw
local sample primary_sample
local y births_15_19_w_per1000
if "`outcome'"=="infant" local y infant_deaths_fr_15_19_W_per1000
if "`outcome'"=="under5" local y under5_deaths_fr_15_19_W_per1000
}
keep if data_part==`part' & `sample'==1 & !missing(`y',exposed_before15,`wt',country_id,survey_id,cohort_id) & `wt'>0
quietly regress `y' exposed_before15 i.survey_id i.cohort_id [`wtype'=`wt'], vce(cluster country_id)
scalar base_b=_b[exposed_before15]
scalar rank0=e(rank)
if "`outcome'"=="union18" assert abs(base_b-(-.064163725))<.000001
if "`outcome'"=="birth18" assert abs(base_b-(-.052670960))<.000001
if "`outcome'"=="birth20" assert abs(base_b-(-.079728124))<.000001
if "`outcome'"=="births" assert abs(base_b-(-50.081071))<.00001
if "`outcome'"=="infant" assert abs(base_b-(-6.978656526))<.000001
if "`outcome'"=="under5" assert abs(base_b-(-10.723971581))<.000001
report_model `outcome' baseline
quietly regress age_calendar i.survey_id i.cohort_id [`wtype'=`wt'], vce(cluster country_id)
predict double age_res if e(sample), resid
quietly summarize age_res
scalar age_sd=r(sd)
gen double abs_age_res=abs(age_res)
quietly summarize abs_age_res
scalar age_max=r(max)
quietly regress `y' exposed_before15 i.survey_id i.cohort_id c.age_calendar [`wtype'=`wt'], vce(cluster country_id)
scalar bdelta=_b[exposed_before15]-base_b
scalar rank1=e(rank)
assert abs(bdelta)<.000001
post `A' ("`outcome'") (age_sd) (age_max) (bdelta) (rank0) (rank1)
drop age_res abs_age_res
quietly regress `y' exposed_before15 i.survey_id i.cohort_id i.country_id#c.cohort_time [`wtype'=`wt'], vce(cluster country_id)
report_model `outcome' country_linear_trends
quietly regress `y' exposed_before15 i.survey_id i.cohort_id i.age_calendar [`wtype'=`wt'], vce(cluster country_id)
report_model `outcome' age_indicators
quietly regress `y' exposed_before15 i.survey_id i.cohort_id i.country_id#c.cohort_time i.age_calendar [`wtype'=`wt'], vce(cluster country_id)
report_model `outcome' trends_and_age
foreach window in full narrow {
local condition "evcat<."
local first 1
local last 21
if "`window'"=="narrow" {
local condition "inrange(event15_bin,-6,8)"
local first 5
local last 19
}
quietly regress `y' ib10.evcat i.survey_id i.cohort_id [`wtype'=`wt'] if `condition', vce(cluster country_id)
estimates save "`folder'/event_`outcome'_`window'.ster", replace
forvalues j=`first'/`last' {
local ev=`j'-11
local b=_b[`j'.evcat]
local se=_se[`j'.evcat]
local omitted=(`se'==0 & `j'!=10)
post `E' ("`outcome'") ("`window'") (`ev') (`b') (`se') (`omitted')
}
}
if `part'==1 {
keep if inrange(age,20,24)
quietly regress `y' exposed_before15 i.country_id i.cohort_id [pw=weight], vce(cluster country_id)
report_model `outcome' s7_young_countryFE
quietly regress `y' exposed_before15 i.survey_id i.cohort_id [pw=weight], vce(cluster country_id)
report_model `outcome' s7_young_surveyFE
keep if !missing(rural,wealth_quintile)
quietly regress `y' exposed_before15 i.survey_id i.cohort_id [pw=weight], vce(cluster country_id)
report_model `outcome' s7_complete_unadjusted
quietly regress `y' exposed_before15 i.rural i.wealth_quintile i.survey_id i.cohort_id [pw=weight], vce(cluster country_id)
report_model `outcome' s7_complete_adjusted
}
}
postclose `R'
postclose `C'
postclose `E'
postclose `A'
use "`folder'/linear_results.dta", clear
assert _N==36
assert reps==`reps' & feasible==`reps'
list outcome spec N b p_wcb, noobs clean
display "R2_LINEAR_COMPLETE MODELS=36 REPS=`reps'"
log close
