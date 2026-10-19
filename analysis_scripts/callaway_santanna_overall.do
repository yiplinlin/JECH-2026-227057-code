clear all
version 19.5
set more off
set linesize 255
set varabbrev off
set processors 1
args input out
if `"`input'"'=="" local input "data/data.dta"
if `"`out'"'=="" local out "results"
confirm file "`input'"
capture mkdir "`out'"
local output "`out'/callaway_santanna_overall"
capture mkdir "`output'"
local results "`out'"
local dta "`input'"
local requested_reps 5000
local reps 5000
capture log close _all
log using "`output'/callaway_santanna_overall.log", text replace
set maxvar 32767
set processors 4
use data_part primary_sample country_id country_name female_birth_year candidate_reform_year exposed_before15 weight analytic_union analytic_birth18 analytic_birth20 union_before18 birth_before18 birth_before20 using "`input'", clear
keep if primary_sample==1 & data_part==1
keep if !missing(country_id,female_birth_year,candidate_reform_year,exposed_before15,weight)
gen int csa_time=female_birth_year+15
bysort country_id: egen byte csa_any_treated=max(exposed_before15)
gen int csa_cohort=candidate_reform_year if csa_any_treated==1
replace csa_cohort=0 if csa_any_treated==0
gen byte csa_treat=csa_time>=csa_cohort if csa_cohort>0 & csa_time<.
replace csa_treat=0 if csa_cohort==0 & csa_time<.
assert csa_treat==exposed_before15 if exposed_before15<.
tempfile individual_input
save `individual_input'
tempname CS MP
postfile `CS' str8 outcome str8 control long N clusters reps seed double att se p ci_lower ci_upper using "`output'/cs_overall.dta", replace
postfile `MP' str8 outcome str8 control double delta_pp bound_lower_pp bound_upper_pp confidence_lower_pp confidence_upper_pp loss_of_negative_ci_at_delta_pp long reps using "`output'/bounded_variation.dta", replace
foreach control in notyet never {
    forvalues k=1/3 {
        use `individual_input', clear
        local ys union_before18 birth_before18 birth_before20
        local ss analytic_union analytic_birth18 analytic_birth20
        local ids union18 birth18 birth20
        local y : word `k' of `ys'
        local sample : word `k' of `ss'
        local id : word `k' of `ids'
        local seed=20260740+`k'
        display "MODEL_START `id' `control' `c(current_date)' `c(current_time)'"
        hdidregress aipw (`y') (csa_treat) [pweight=weight] if `sample'==1 & `y'<. & weight>0, group(country_id) time(csa_time) usercohort(csa_cohort) controlgroup(`control') basetime(common) vce(cluster country_id) nolog
        
        estat aggregation, overall sci(reps(5000) rseed(`seed'))
        matrix T=r(table)
        post `CS' ("`id'") ("`control'") (e(N)) (e(N_clust)) (5000) (`seed') (T[1,1]) (T[2,1]) (T[4,1]) (T[5,1]) (T[6,1])
        foreach delta in 0 .25 .5 1 2 3 5 {
            post `MP' ("`id'") ("`control'") (`delta') (100*T[1,1]-`delta') (100*T[1,1]+`delta') (100*T[5,1]-`delta') (100*T[6,1]+`delta') (max(0,-100*T[6,1])) (5000)
        }
    }
}
postclose `CS'
postclose `MP'
use "`output'/cs_overall.dta", clear
list, noobs
export delimited using "`output'/cs_overall.csv", replace
use "`output'/bounded_variation.dta", clear
list, noobs
export delimited using "`output'/bounded_variation.csv", replace
log close
