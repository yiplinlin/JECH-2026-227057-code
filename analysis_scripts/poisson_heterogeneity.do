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
local output "`out'/poisson"
capture mkdir "`output'"
local results "`out'"
local dta "`input'"
local requested_reps 5000
local reps 5000
capture log close _all
log using "`output'/poisson_heterogeneity.log", text replace
capture program drop ppml_att
program define ppml_att, rclass
syntax varname, HET(string)
tempvar t g event post expo count cl sfe hcat tmpcat d mu1 hkeep mu0 rd
preserve
quietly keep if data_part==2 & `varlist'<. & weighted_women_n>0 & !missing(country_id,survey_id,female_birth_cohort,candidate_reform_year)
quietly gen int `t'=female_birth_cohort+15
quietly gen int `g'=candidate_reform_year
quietly gen int `event'=`t'-`g'
quietly gen byte `post'=(`g'<. & `event'>=0)
quietly gen double `expo'=weighted_women_n/1000
quietly gen double `count'=`varlist'*`expo'
capture confirm variable __boot_country_id
if _rc {
quietly gen long `cl'=country_id
quietly egen long `sfe'=group(country_id survey_id)
local vceopt "vce(cluster `cl')"
}
else {
quietly gen long `cl'=__boot_country_id
quietly egen long `sfe'=group(__boot_country_id survey_id)
local vceopt "vce(robust)"
}
quietly gen int `hcat'=0
if "`het'"=="event" {
quietly replace `hcat'=`event'+1 if inrange(`event',0,4)
quietly replace `hcat'=6 if `event'>=5 & `event'<.
}
if "`het'"=="cohort" {
quietly egen int `tmpcat'=group(`g') if `post'==1
quietly replace `hcat'=`tmpcat' if `post'==1
}
if "`het'"=="time" {
quietly egen int `tmpcat'=group(`t') if `post'==1
quietly replace `hcat'=`tmpcat' if `post'==1
}
capture quietly ppmlhdfe `count' ib0.`hcat', absorb(`sfe' `t') exposure(`expo') `vceopt' d(`d')
local fitrc=_rc
if !`fitrc' {
capture quietly predict double `mu1', mu
local fitrc=_rc
}
if !`fitrc' {
quietly clonevar `hkeep'=`hcat'
quietly replace `hcat'=0 if `post'==1
capture quietly predict double `mu0', mu
local fitrc=_rc
}
if `fitrc' {
restore
return scalar att=.
return scalar n=0
return scalar rc=`fitrc'
exit
}
quietly replace `hcat'=`hkeep'
quietly gen double `rd'=(`mu1'-`mu0')/`expo'
quietly summarize `rd' [aw=weighted_women_n] if `post'==1 & `rd'<., meanonly
local att=r(mean)
local n=r(N)
restore
return scalar att=`att'
return scalar n=`n'
return scalar rc=0
end
capture program drop ppml_draw
program define ppml_draw, rclass
syntax varname, HET(string)
preserve
capture quietly bsample, cluster(country_id) idcluster(__boot_country_id)
local rc=_rc
if !`rc' {
capture quietly ppml_att `varlist', het(`het')
local rc=_rc
if !`rc' {
local value=r(att)
local rc=r(rc)
}
}
if `rc' local value=.
restore
return scalar att=`value'
return scalar rc=`rc'
end
local target `reps'
local outdir "`output'"
foreach outcome in births infant under5 {
foreach het in event cohort time {
if 1 {
local jobdir "`outdir'/ppml_`outcome'_`het'"
capture mkdir "`jobdir'"
use "`input'", clear
keep if data_part==2 & primary_sample==1
sort country_id survey_id female_birth_cohort
keep if !inlist(country_id,3,4)
local outcomevar "infant_deaths_fr_15_19_W_per1000"
if "`outcome'"=="under5" local outcomevar "under5_deaths_fr_15_19_W_per1000"
if "`outcome'"=="births" local outcomevar "births_15_19_w_per1000"
generate double rate=`outcomevar'
keep if rate<. & weighted_women_n>0 & !missing(country_id,survey_id,female_birth_cohort,candidate_reform_year)
assert rate>=0
local sampleN=_N
local seed=20260745
if "`outcome'"=="births" local seed=20260735
if "`outcome'"=="under5" local seed=20260755
if "`het'"=="cohort" local seed=`seed'+1
if "`het'"=="time" local seed=`seed'+2
display "MODEL_START OUTCOME=`outcome' HET=`het' REPS=`target' SEED=`seed' N=`sampleN'"
quietly ppml_att rate, het(`het')
scalar pointatt=r(att)
scalar pointn=r(n)
assert pointatt<.
display "POINT_ATT=" %16.10f pointatt " POST_N=" pointn
tempname draws summary observed
postfile `draws' long attempt double att int rc byte valid using "`jobdir'/draws.dta", replace
local valid=0
local attempts=0
local failed=0
set seed `seed'
while `valid'<`target' & `attempts'<5*`target' {
capture quietly ppml_draw rate, het(`het')
local rc=_rc
local val=.
if !`rc' {
local val=r(att)
local rc=r(rc)
}
local ++attempts
local ok=(`rc'==0 & `val'<.)
if `ok' local ++valid
else local ++failed
post `draws' (`attempts') (`val') (`rc') (`ok')
if mod(`attempts',250)==0 display "PROGRESS VALID=`valid' FAILED=`failed' ATTEMPTS=`attempts' TARGET=`target'"
}
postclose `draws'
use "`jobdir'/draws.dta", clear
keep if valid==1
save "`jobdir'/valid_draws.dta", replace
matrix `observed'=(pointatt)
matrix colnames `observed'=att
bstat att, stat(`observed') n(`sampleN')
estat bootstrap, normal percentile
ereturn list
quietly summarize att
scalar bootse=r(sd)
scalar bp=2*normal(-abs(pointatt/bootse))
scalar nlo=pointatt-invnormal(.975)*bootse
scalar nhi=pointatt+invnormal(.975)*bootse
quietly _pctile att, p(2.5 97.5)
scalar plo=r(r1)
scalar phi=r(r2)
postfile `summary' str8 outcome str8 heterogeneity double n_groups double n_post double att double se double p_normal double ci_normal_lo double ci_normal_hi double ci_percentile_lo double ci_percentile_hi long seed long valid long failed long attempts using "`jobdir'/summary.dta", replace
post `summary' ("`outcome'") ("`het'") (`sampleN') (pointn) (pointatt) (bootse) (bp) (nlo) (nhi) (plo) (phi) (`seed') (`valid') (`failed') (`attempts')
postclose `summary'
display "MODEL_COMPLETE OUTCOME=`outcome' HET=`het' ATT=" %14.8f pointatt " SE=" %14.8f bootse " P=" %12.8f bp " VALID=`valid' FAILED=`failed'"
if `valid'<`target' {
display as error "TARGET_NOT_REACHED"
exit 459
}
}
}
}

log close
