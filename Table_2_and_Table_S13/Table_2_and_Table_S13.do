version 19.5
clear all
set more off
set linesize 255
set varabbrev off
set processors 1
args input out
if `"`input'"'=="" local input "data/data.dta"
if `"`out'"'=="" local out "results"
confirm file "`input'"
capture mkdir "`out'"
local joint "`out'"
local output "`joint'/Table_2"
capture mkdir "`output'"
capture log close _all
log using "`joint'/Table_2_and_Table_S13.log", text replace name(combined)
display "COMBINED_RUN_START DATE=`c(current_date)' TIME=`c(current_time)'"
display "STATA_VERSION=`c(stata_version)' REPS=5000 SEED=20260708"
display "INPUT=`input'"
which boottest
display "TABLE_2_START"
global r2out "`output'"
global r2reps 5000
tempname R C
postfile `R' str12 outcome str30 spec double N clusters b se ci_lo ci_hi p_cluster p_wcb reps feasible ci_parts using "`output'/Table_2_results.dta", replace
postfile `C' str12 outcome str30 spec int part double lo hi using "`output'/Table_2_intervals.dta", replace
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
regress `y' exposed_before15 i.survey_id i.cohort_id [`wtype'=`wt'], vce(cluster country_id)
scalar base_b=_b[exposed_before15]
scalar rank0=e(rank)
if "`outcome'"=="union18" assert abs(base_b-(-.064163725))<.000001
if "`outcome'"=="birth18" assert abs(base_b-(-.052670960))<.000001
if "`outcome'"=="birth20" assert abs(base_b-(-.079728124))<.000001
if "`outcome'"=="births" assert abs(base_b-(-50.081071))<.00001
if "`outcome'"=="infant" assert abs(base_b-(-6.978656526))<.000001
if "`outcome'"=="under5" assert abs(base_b-(-10.723971581))<.000001
report_model `outcome' baseline
}
postclose `R'
postclose `C'
use "`output'/Table_2_results.dta", clear
assert _N==6
assert clusters==16 & reps==5000 & feasible==5000
generate double display_scale=cond(inlist(outcome,"union18","birth18","birth20"),100,1)
foreach variable in b se ci_lo ci_hi {
    replace `variable'=`variable'*display_scale
}
list, noobs
export delimited using "`output'/Table_2_results.csv", replace
use "`output'/Table_2_intervals.dta", clear
generate double display_scale=cond(inlist(outcome,"union18","birth18","birth20"),100,1)
replace lo=lo*display_scale
replace hi=hi*display_scale
list, noobs
export delimited using "`output'/Table_2_intervals.csv", replace
macro drop r2out r2reps r2res r2ci
display "TABLE_2_COMPLETE MODELS=6 REPS_PER_MODEL=5000"
local out "`joint'/Table_S13"
capture mkdir "`out'"
display "TABLE_S13_START"
display "S13_WCB_START DATE=`c(current_date)' TIME=`c(current_time)'"
display "DATA=`input'"
display "MODEL_SCOPE=Table S13 panels A and B"
display "INFERENCE=Null-imposed wild bootstrap-t, Rademacher weights, country clusters, 5000 replications, seed 20260708"
display "CONFIDENCE_INTERVALS=95 percent, inversion of two-sided bootstrap tests"
display "PANEL_A_DISPLAY_UNIT=percentage points; PANEL_B_DISPLAY_UNIT=per 1000 women"
display "NO_DATA_OR_EXISTING_ANALYSIS_FILES_ARE_MODIFIED"
which boottest
tempname R C
postfile `R' str1 panel str40 outcome str20 unit long N int clusters double estimate double se double ci_cluster_low double ci_cluster_high double p_cluster double p_wild long reps long feasible int ci_parts double ci_wild_low double ci_wild_high double delta_previous double delta_fe_b double delta_fe_se double ci_endpoint_max_deviation using "`out'/s13_wild_results.dta", replace
postfile `C' str1 panel str40 outcome int component double lower double upper using "`out'/s13_wild_intervals.dta", replace
global s13_post_results `R'
global s13_post_intervals `C'
capture program drop report_s13
program define report_s13
    args panel outcome wtype wt scale expectedn expectedb
    tempname n g b se df pc lo hi delta bb ss ci wp nr nf dd null ci_low ci_high orig endpoint_error endpoint_p
    tempvar sample
    scalar `n'=e(N)
    scalar `g'=e(N_clust)
    scalar `b'=_b[exposed_before15]
    scalar `se'=_se[exposed_before15]
    scalar `df'=e(df_r)
    scalar `pc'=2*ttail(`df',abs(`b'/`se'))
    scalar `lo'=`b'-invttail(`df',.025)*`se'
    scalar `hi'=`b'+invttail(`df',.025)*`se'
    scalar `delta'=`b'-`expectedb'
    assert `n'==`expectedn'
    assert `g'==16
    assert abs(`delta')<.0000005
    gen byte `sample'=e(sample)
    estimates store `orig'
    quietly regress `outcome' exposed_before15 i.survey_id i.cohort_id [`wtype'=`wt'] if `sample', vce(cluster country_id)
    assert e(sample)==`sample'
    scalar `bb'=_b[exposed_before15]-`b'
    scalar `ss'=_se[exposed_before15]-`se'
    assert abs(`bb')<.0000001
    assert abs(`ss')<.0000001
    estimates restore `orig'
    estimates drop `orig'
    display "S13_MODEL_START PANEL=`panel' OUTCOME=`outcome' N=" `n' " CLUSTERS=" `g'
    display "PREVIOUS_POINT_DIFFERENCE=" %20.12g `delta'
    display "SAME_SAMPLE_FE_EQUIVALENCE_B=" %20.12g `bb' " SE=" %20.12g `ss'
    boottest exposed_before15, cluster(country_id) reps(5000) seed(20260708) weighttype(rademacher) level(95) nograph svmat
    return list
    scalar `wp'=r(p)
    scalar `nr'=r(reps)
    scalar `null'=r(null)
    matrix `dd'=r(dist)
    matrix `ci'=r(CI)
    mata: st_numscalar("`nf'",sum(st_matrix("`dd'"):<.))
    assert `nr'==5000
    assert `nf'==5000
    assert `null'==1
    assert inrange(`wp',0,1)
    local parts=rowsof(`ci')
    assert `parts'>=1
    scalar `endpoint_error'=0
    forvalues j=1/`parts' {
        post $s13_post_intervals ("`panel'") ("`outcome'") (`j') (`scale'*`ci'[`j',1]) (`scale'*`ci'[`j',2])
        forvalues side=1/2 {
            if `ci'[`j',`side']<. {
                local boundary: display %24.17g `ci'[`j',`side']
                quietly boottest exposed_before15=`boundary', cluster(country_id) reps(5000) seed(20260708) weighttype(rademacher) noci nograph
                scalar `endpoint_p'=r(p)
                scalar `endpoint_error'=max(`endpoint_error',abs(`endpoint_p'-.05))
                assert abs(`endpoint_p'-.05)<.003
                display "S13_CI_ENDPOINT_CHECK OUTCOME=`outcome' PART=`j' SIDE=`side' NULL=" %20.12g `ci'[`j',`side'] " P=" %10.6f `endpoint_p'
            }
        }
    }
    scalar `ci_low'=.
    scalar `ci_high'=.
    if `parts'==1 {
        scalar `ci_low'=`scale'*`ci'[1,1]
        scalar `ci_high'=`scale'*`ci'[1,2]
    }
    local unit "per 1000 women"
    if "`panel'"=="A" local unit "percentage points"
    post $s13_post_results ("`panel'") ("`outcome'") ("`unit'") (`n') (`g') (`scale'*`b') (`scale'*`se') (`scale'*`lo') (`scale'*`hi') (`pc') (`wp') (`nr') (`nf') (`parts') (`ci_low') (`ci_high') (`delta') (`bb') (`ss') (`endpoint_error')
    display "S13_MODEL_COMPLETE PANEL=`panel' OUTCOME=`outcome' ESTIMATE=" %16.9f (`scale'*`b') " WILD_P=" %10.6f `wp' " WILD_CI_LOW=" %16.9f `ci_low' " WILD_CI_HIGH=" %16.9f `ci_high' " CI_PARTS=`parts' VALID_REPS=" `nf'
end
use "`input'", clear
keep if primary_sample==1
capture drop overlap_main
gen byte overlap_main=1 if data_part<.
assert overlap_main==1 if inlist(data_part,1,2)
local macros gdp_pc_constant_survey health_exp_gdp_survey urban_pop_pct_survey under5_mortality_rate_survey rule_of_law_survey ucdp_total_deaths_survey dtp3_coverage_survey
preserve
keep if data_part==1
foreach v of local macros {
    bysort country_id survey_id: assert `v'==`v'[1]
}
keep country_id survey_id `macros'
bysort country_id survey_id: keep if _n==1
isid country_id survey_id
tempfile macro_crosswalk
save `macro_crosswalk'
restore
foreach y in union_before18 birth_before18 birth_before20 {
    local n 241506
    local b -.058244
    if "`y'"=="union_before18" {
        local n 238175
        local b -.070745
    }
    if "`y'"=="birth_before20" local b -.086463
    regress `y' exposed_before15 `macros' i.survey_id i.cohort_id [pw=weight] if data_part==1 & overlap_main==1 & `y'<. & exposed_before15<. & weight>0, vce(cluster country_id)
    report_s13 A `y' pw weight 100 `n' `b'
}
keep if data_part==2 & overlap_main==1
drop `macros'
merge m:1 country_id survey_id using `macro_crosswalk', assert(match) nogen
isid country_id survey_id cohort_id
foreach y in births_15_19_w_per1000 births_15_24_w_per1000 infant_deaths_fr_15_19_W_per1000 under5_deaths_fr_15_19_W_per1000 {
    local n 500
    local b -57.659048343
    if "`y'"=="births_15_24_w_per1000" {
        local n 250
        local b -53.684217234
    }
    if "`y'"=="infant_deaths_fr_15_19_W_per1000" {
        local n 450
        local b -7.507276830
    }
    if "`y'"=="under5_deaths_fr_15_19_W_per1000" {
        local n 300
        local b -10.704770232
    }
    regress `y' exposed_before15 `macros' i.survey_id i.cohort_id [aw=weighted_women_n] if `y'<. & exposed_before15<. & weighted_women_n>0 & weighted_women_n<., vce(cluster country_id)
    report_s13 B `y' aw weighted_women_n 1 `n' `b'
}
postclose `R'
postclose `C'
use "`out'/s13_wild_results.dta", clear
assert _N==7
assert reps==5000 & feasible==5000 & clusters==16
assert abs(delta_previous)<.0000005
assert abs(delta_fe_b)<.0000001 & abs(delta_fe_se)<.0000001
assert ci_endpoint_max_deviation<.003
format estimate se ci_cluster_low ci_cluster_high ci_wild_low ci_wild_high %12.6f
format p_cluster p_wild %10.6f
sort panel outcome
save "`out'/s13_wild_results.dta", replace
export delimited using "`out'/s13_wild_results_table.txt", delimiter(tab) replace
list panel outcome N estimate ci_wild_low ci_wild_high p_wild ci_parts feasible, noobs abbreviate(40)
use "`out'/s13_wild_intervals.dta", clear
sort panel outcome component
save "`out'/s13_wild_intervals.dta", replace
export delimited using "`out'/s13_wild_intervals_table.txt", delimiter(tab) replace
list, noobs abbreviate(40)
display "S13_WCB_COMPLETE MODELS=7 REPS_PER_MODEL=5000 SEED=20260708 ALL_POINT_ESTIMATES_AND_SAMPLE_SIZES_VERIFIED"
macro drop s13_post_results s13_post_intervals
display "COMBINED_RUN_COMPLETE TABLE_2_MODELS=6 TABLE_S13_MODELS=7 REPS_PER_MODEL=5000"
display "END_DATE=`c(current_date)' END_TIME=`c(current_time)'"
log close combined
