version 19.5
set more off
set linesize 180
set processors 1
args input folder
if `"`input'"'=="" local input "data/data.dta"
if `"`folder'"'=="" local folder "results"
local reps 5000
capture mkdir "`folder'"
local folder "`folder'/poisson_linear_comparison"
capture mkdir "`folder'"
capture log close _all
log using "`folder'/poisson_comparison.txt", text replace
capture program drop paired_models
program define paired_models, rclass
syntax, [REPORT]
preserve
tempvar cl sfe t event post expo count hcat saved d mu1 mu0 rd sample
capture confirm variable boot_country
if _rc gen long `cl'=country_id
else gen long `cl'=boot_country
egen long `sfe'=group(`cl' survey_id)
gen int `t'=female_birth_cohort+15
gen int `event'=`t'-candidate_reform_year
gen byte `post'=`event'>=0
gen double `expo'=weighted_women_n/1000
gen double `count'=rate*`expo'
tempname answer
matrix `answer'=J(1,8,.)
local k=0
foreach het in constant event cohort time {
local ++k
gen int `hcat'=0
if "`het'"=="constant" replace `hcat'=1 if `post'
if "`het'"=="event" {
replace `hcat'=`event'+1 if inrange(`event',0,4)
replace `hcat'=6 if `event'>=5
}
if "`het'"=="cohort" {
egen int __hc=group(candidate_reform_year) if `post'
replace `hcat'=__hc if `post'
drop __hc
}
if "`het'"=="time" {
egen int __hc=group(`t') if `post'
replace `hcat'=__hc if `post'
drop __hc
}
quietly ppmlhdfe `count' ib0.`hcat', absorb(`sfe' `t') exposure(`expo') d(`d')
gen byte `sample'=e(sample)
quietly count if `sample'
local nfit=r(N)
quietly count if `sample' & `post'
local npost=r(N)
assert `npost'>0
predict double `mu1' if `sample', mu
clonevar `saved'=`hcat'
replace `hcat'=0 if `post'
predict double `mu0' if `sample', mu
replace `hcat'=`saved'
gen double `rd'=(`mu1'-`mu0')/`expo' if `sample'
quietly summarize `rd' [aw=weighted_women_n] if `sample' & `post', meanonly
matrix `answer'[1,`k'+4]=r(mean)
quietly regress rate ib0.`hcat' i.`sfe' i.`t' [aw=weighted_women_n] if `sample', vce(cluster `cl')
quietly levelsof `hcat' if `sample' & `post', local(cats)
quietly summarize weighted_women_n if `sample' & `post', meanonly
local denom=r(sum)
local att=0
foreach c of local cats {
quietly summarize weighted_women_n if `sample' & `hcat'==`c', meanonly
local share=r(sum)/`denom'
local att=`att'+`share'*_b[`c'.`hcat']
}
matrix `answer'[1,`k']=`att'
if "`report'"!="" {
display "PAIRED_MODEL HET=`het' FIT_N=`nfit' POST_N=`npost' LINEAR_ATT=" %16.10f `answer'[1,`k'] " PPML_ATT=" %16.10f `answer'[1,`k'+4]
}
drop `hcat' `saved' `d' `mu1' `mu0' `rd' `sample'
}
restore
forvalues j=1/8 {
return scalar a`j'=`answer'[1,`j']
}
end
use "`input'", clear
keep if primary_sample==1 & data_part==2
gen double rate=births_15_19_w_per1000
keep if !missing(rate,weighted_women_n,country_id,survey_id,female_birth_cohort,candidate_reform_year) & weighted_women_n>0
assert _N==770
sort country_id survey_id female_birth_cohort
tempfile original
save `original'
noisily paired_models, report
matrix observed=J(1,8,.)
forvalues j=1/8 {
matrix observed[1,`j']=r(a`j')
}
assert abs(observed[1,1]-(-50.081071))<.00001
assert abs(observed[1,6]-(-88.446247556))<.00001
assert abs(observed[1,7]-(-66.315718))<.00001
assert abs(observed[1,8]-(-60.939357))<.00001
tempname D S
postfile `D' long attempt byte valid int rc double a1 a2 a3 a4 a5 a6 a7 a8 using "`folder'/poisson_paired_draws.dta", replace
local valid=0
local attempts=0
set seed 20260708
while `valid'<`reps' & `attempts'<5*`reps' {
use `original', clear
bsample, cluster(country_id) idcluster(boot_country)
capture quietly paired_models
local rc=_rc
local ok=(`rc'==0)
forvalues j=1/8 {
local a`j'=.
if !`rc' local a`j'=r(a`j')
if missing(`a`j'') local ok=0
}
local ++attempts
if `ok' local ++valid
post `D' (`attempts') (`ok') (`rc') (`a1') (`a2') (`a3') (`a4') (`a5') (`a6') (`a7') (`a8')
if mod(`attempts',100)==0 display "PAIRED_PROGRESS VALID=`valid' ATTEMPTS=`attempts' TARGET=`reps'"
}
postclose `D'
assert `valid'==`reps'
use "`folder'/poisson_paired_draws.dta", clear
keep if valid
postfile `S' str12 heterogeneity str12 model double estimate se p_normal ci_normal_lo ci_normal_hi ci_pct_lo ci_pct_hi long valid attempts using "`folder'/poisson_comparison_results.dta", replace
local hets constant event cohort time
forvalues j=1/4 {
local het: word `j' of `hets'
foreach model in linear poisson {
local idx=`j'
if "`model'"=="poisson" local idx=`j'+4
local point=observed[1,`idx']
gen double x=a`idx'
quietly summarize x
local se=r(sd)
local p=2*normal(-abs(`point'/`se'))
local lo=`point'-invnormal(.975)*`se'
local hi=`point'+invnormal(.975)*`se'
quietly _pctile x, p(2.5 97.5)
post `S' ("`het'") ("`model'") (`point') (`se') (`p') (`lo') (`hi') (r(r1)) (r(r2)) (`valid') (`attempts')
drop x
}
}
postclose `S'
use "`folder'/poisson_comparison_results.dta", clear
list, noobs clean
display "POISSON_COMPARISON_COMPLETE VALID=`valid' ATTEMPTS=`attempts' SEED=20260708"
log close
