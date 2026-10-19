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
local output "`out'/cutoffs"
capture mkdir "`output'"
local results "`out'"
local dta "`input'"
local requested_reps 5000
local reps 5000
capture log close _all
log using "`output'/cutoffs_and_pretrends.log", text replace
version 19.5
set more off
set linesize 255
set varabbrev off
local cutcsv "`output'/cutoffs.csv"
local precsv "`output'/proximal_pretrends.csv"
local summary "`output'/cutoff_proximal_notes.txt"
confirm file "`dta'"
display "CUTOFFS_AND_PRETRENDS_START"
display "RUN_DATE=`c(current_date)'"
display "REPLICATIONS=`reps' SEED=20260708"
display "INPUT_DTA=`dta'"
which boottest
use "`dta'", clear
keep if primary_sample==1
describe, short
datasignature
capture confirm variable data_part
assert _rc==0
assert inrange(data_part,1,6)
local outcomes "union_before18 birth_before18 birth_before20 births_15_19_w_per1000 infant_deaths_fr_15_19_W_per1000 under5_deaths_fr_15_19_W_per1000"
local samples  "analytic_union analytic_birth18 analytic_birth20 burden burden burden"
local scales   "proportion proportion proportion per_1000_women per_1000_women per_1000_women"
tempfile cutoff_results
tempname CUT
postfile `CUT' int cutoff str45 outcome str22 scale  double N clusters support_n zero_n one_n b se ci_lo ci_hi p_cluster  wcb_p wcb_ci_lo wcb_ci_hi wcb_reps  str244 all_countries str244 support_countries str244 zero_countries str244 one_countries  using "`cutoff_results'", replace
forvalues k=12/17 {
    capture drop exp_cut
    quietly generate byte exp_cut=(candidate_reform_year<=female_birth_cohort+`k')  if inlist(data_part,1,2) & !missing(candidate_reform_year,female_birth_cohort)
    if `k'==15 {
        quietly count if inlist(data_part,1,2) & exp_cut!=exposed_before15  & !missing(exp_cut,exposed_before15)
        display "CUTOFF15_EXPOSURE_MISMATCH_N=" r(N)
        assert r(N)==0
    }
    forvalues m=1/6 {
        local y : word `m' of `outcomes'
        local s : word `m' of `samples'
        local scale : word `m' of `scales'
        if `m'<=3 {
            quietly regress `y' exp_cut i.survey_id i.cohort_id [pweight=weight]  if data_part==1 & `s'==1 & `y'<. & exp_cut<. & weight>0,  vce(cluster country_id)
        }
        else {
            quietly regress `y' exp_cut i.survey_id i.cohort_id [aweight=weighted_women_n]  if data_part==2 & `y'<. & exp_cut<. & weighted_women_n>0,  vce(cluster country_id)
        }
        scalar __N=e(N)
        scalar __b=_b[exp_cut]
        scalar __se=_se[exp_cut]
        scalar __df=e(df_r)
        scalar __crit=invttail(__df,.025)
        scalar __lo=__b-__crit*__se
        scalar __hi=__b+__crit*__se
        scalar __p=2*ttail(__df,abs(__b/__se))
        tempvar esample
        quietly generate byte `esample'=e(sample)
        quietly levelsof country_id if `esample'==1, local(allids)
        local clusters : word count `allids'
        local alliso ""
        local supportiso ""
        local zeroiso ""
        local oneiso ""
        local support_n=0
        local zero_n=0
        local one_n=0
        foreach c of local allids {
            quietly levelsof iso3 if `esample'==1 & country_id==`c', local(thisiso) clean
            local alliso "`alliso' `thisiso'"
            quietly count if `esample'==1 & country_id==`c' & exp_cut==0
            local n0=r(N)
            quietly count if `esample'==1 & country_id==`c' & exp_cut==1
            local n1=r(N)
            if `n0'>0 & `n1'>0 {
                local supportiso "`supportiso' `thisiso'"
                local ++support_n
            }
            if `n0'>0 & `n1'==0 {
                local zeroiso "`zeroiso' `thisiso'"
                local ++zero_n
            }
            if `n0'==0 & `n1'>0 {
                local oneiso "`oneiso' `thisiso'"
                local ++one_n
            }
        }
        local alliso=strtrim("`alliso'")
        local supportiso=strtrim("`supportiso'")
        local zeroiso=strtrim("`zeroiso'")
        local oneiso=strtrim("`oneiso'")
        local bootseed=20260708
        capture noisily boottest exp_cut, cluster(country_id) reps(`reps')  seed(`bootseed') nograph
        if _rc {
            scalar __wp=.
            scalar __wlo=.
            scalar __whi=.
            scalar __wreps=.
            display as error "BOOTTEST_FAILED CUTOFF=`k' OUTCOME=`y' RC=" _rc
        }
        else {
            scalar __wp=r(p)
            matrix __wci=r(CI)
            scalar __wlo=__wci[1,1]
            scalar __whi=__wci[1,2]
            scalar __wreps=r(reps)
        }
        display "CUTOFF_RESULT cutoff=`k' outcome=`y' N=" __N  " clusters=`clusters' support_n=`support_n' b=" %14.9f __b  " se=" %14.9f __se " p_cluster=" %10.7f __p  " wcb_p=" %10.7f __wp
        post `CUT' (`k') ("`y'") ("`scale'") (__N) (`clusters')  (`support_n') (`zero_n') (`one_n') (__b) (__se) (__lo) (__hi)  (__p) (__wp) (__wlo) (__whi) (__wreps)  ("`alliso'") ("`supportiso'") ("`zeroiso'") ("`oneiso'")
        drop `esample'
    }
}
postclose `CUT'
preserve
use "`cutoff_results'", clear
sort cutoff outcome
format %21.15g N clusters support_n zero_n one_n b se ci_lo ci_hi p_cluster  wcb_p wcb_ci_lo wcb_ci_hi wcb_reps
export delimited using "`cutcsv'", replace datafmt
restore
capture drop event15_cat
quietly generate int event15_cat=event15_bin+11 if event15_bin<.
tempfile proximal_results
tempname PRE
postfile `PRE' str10 window str45 outcome str22 scale  double N clusters prepost_n lead_m3_n lead_m2_n  b_m3 se_m3 ci_lo_m3 ci_hi_m3 p_m3 wcb_p_m3 wcb_ci_lo_m3 wcb_ci_hi_m3  b_m2 se_m2 ci_lo_m2 ci_hi_m2 p_m2 wcb_p_m2 wcb_ci_lo_m2 wcb_ci_hi_m2  joint_F joint_p joint_wcb_p joint_wcb_reps  str244 all_countries str244 prepost_countries str244 lead_m3_countries str244 lead_m2_countries  using "`proximal_results'", replace
foreach window in full narrow {
    forvalues m=1/6 {
        local y : word `m' of `outcomes'
        local s : word `m' of `samples'
        local scale : word `m' of `scales'
        local winif "event15_cat<."
        if "`window'"=="narrow" local winif "inrange(event15_bin,-6,8) & event15_cat<."
        if `m'<=3 {
            quietly regress `y' ib10.event15_cat i.survey_id i.cohort_id [pweight=weight]  if data_part==1 & `s'==1 & `winif' & `y'<. & weight>0,  vce(cluster country_id)
        }
        else {
            quietly regress `y' ib10.event15_cat i.survey_id i.cohort_id [aweight=weighted_women_n]  if data_part==2 & `winif' & `y'<. & weighted_women_n>0,  vce(cluster country_id)
        }
        scalar __N=e(N)
        scalar __df=e(df_r)
        scalar __crit=invttail(__df,.025)
        scalar __b3=_b[8.event15_cat]
        scalar __s3=_se[8.event15_cat]
        scalar __l3=__b3-__crit*__s3
        scalar __h3=__b3+__crit*__s3
        scalar __p3=2*ttail(__df,abs(__b3/__s3))
        scalar __b2=_b[9.event15_cat]
        scalar __s2=_se[9.event15_cat]
        scalar __l2=__b2-__crit*__s2
        scalar __h2=__b2+__crit*__s2
        scalar __p2=2*ttail(__df,abs(__b2/__s2))
        tempvar esample
        quietly generate byte `esample'=e(sample)
        quietly levelsof country_id if `esample'==1, local(allids)
        local clusters : word count `allids'
        local alliso ""
        local prepostiso ""
        local lead3iso ""
        local lead2iso ""
        local prepost_n=0
        local lead3_n=0
        local lead2_n=0
        foreach c of local allids {
            quietly levelsof iso3 if `esample'==1 & country_id==`c', local(thisiso) clean
            local alliso "`alliso' `thisiso'"
            quietly count if `esample'==1 & country_id==`c' & event15_bin<0
            local npre=r(N)
            quietly count if `esample'==1 & country_id==`c' & event15_bin>=0
            local npost=r(N)
            quietly count if `esample'==1 & country_id==`c' & event15_bin==-3
            local nm3=r(N)
            quietly count if `esample'==1 & country_id==`c' & event15_bin==-2
            local nm2=r(N)
            if `npre'>0 & `npost'>0 {
                local prepostiso "`prepostiso' `thisiso'"
                local ++prepost_n
            }
            if `nm3'>0 {
                local lead3iso "`lead3iso' `thisiso'"
                local ++lead3_n
            }
            if `nm2'>0 {
                local lead2iso "`lead2iso' `thisiso'"
                local ++lead2_n
            }
        }
        local alliso=strtrim("`alliso'")
        local prepostiso=strtrim("`prepostiso'")
        local lead3iso=strtrim("`lead3iso'")
        local lead2iso=strtrim("`lead2iso'")
        quietly test 8.event15_cat 9.event15_cat
        scalar __jf=r(F)
        scalar __jp=r(p)
        local wcode=cond("`window'"=="full",1,2)
        local seed3=20260708
        local seed2=20260708
        local seedj=20260708
        capture noisily boottest 8.event15_cat, cluster(country_id) reps(`reps')  seed(`seed3') nograph
        if _rc {
            scalar __wp3=.
            scalar __wl3=.
            scalar __wh3=.
        }
        else {
            scalar __wp3=r(p)
            matrix __wci3=r(CI)
            scalar __wl3=__wci3[1,1]
            scalar __wh3=__wci3[1,2]
        }
        capture noisily boottest 9.event15_cat, cluster(country_id) reps(`reps')  seed(`seed2') nograph
        if _rc {
            scalar __wp2=.
            scalar __wl2=.
            scalar __wh2=.
        }
        else {
            scalar __wp2=r(p)
            matrix __wci2=r(CI)
            scalar __wl2=__wci2[1,1]
            scalar __wh2=__wci2[1,2]
        }
        capture noisily boottest 8.event15_cat 9.event15_cat,  cluster(country_id) reps(`reps') seed(`seedj') nograph
        if _rc {
            scalar __jwp=.
            scalar __jwreps=.
        }
        else {
            scalar __jwp=r(p)
            scalar __jwreps=r(reps)
        }
        display "PROXIMAL_RESULT window=`window' outcome=`y' N=" __N  " clusters=`clusters' b_m3=" %14.9f __b3 " p_m3=" %10.7f __p3  " wcb_p_m3=" %10.7f __wp3 " b_m2=" %14.9f __b2  " p_m2=" %10.7f __p2 " wcb_p_m2=" %10.7f __wp2  " joint_p=" %10.7f __jp " joint_wcb_p=" %10.7f __jwp
        post `PRE' ("`window'") ("`y'") ("`scale'") (__N) (`clusters')  (`prepost_n') (`lead3_n') (`lead2_n')  (__b3) (__s3) (__l3) (__h3) (__p3) (__wp3) (__wl3) (__wh3)  (__b2) (__s2) (__l2) (__h2) (__p2) (__wp2) (__wl2) (__wh2)  (__jf) (__jp) (__jwp) (__jwreps)  ("`alliso'") ("`prepostiso'") ("`lead3iso'") ("`lead2iso'")
        drop `esample'
    }
}
postclose `PRE'
preserve
use "`proximal_results'", clear
sort window outcome
format %21.15g N clusters prepost_n lead_m3_n lead_m2_n  b_m3 se_m3 ci_lo_m3 ci_hi_m3 p_m3 wcb_p_m3 wcb_ci_lo_m3 wcb_ci_hi_m3  b_m2 se_m2 ci_lo_m2 ci_hi_m2 p_m2 wcb_p_m2 wcb_ci_lo_m2 wcb_ci_hi_m2  joint_F joint_p joint_wcb_p joint_wcb_reps
export delimited using "`precsv'", replace datafmt
restore
capture drop exp15_check
quietly generate byte exp15_check=(candidate_reform_year<=female_birth_cohort+15)  if inlist(data_part,1,2) & !missing(candidate_reform_year,female_birth_cohort)
quietly count if inlist(data_part,1,2) & exp15_check!=exposed_before15  & !missing(exp15_check,exposed_before15)
assert r(N)==0
local expected_b1=-0.0641641163
local expected_b2=-0.0526707926
local expected_b3=-0.0797277330
local expected_b4=-50.08107109
local expected_b5=-6.978656526
local expected_b6=-10.723971581
local expected_n1=311260
local expected_n2=314592
local expected_n3=314592
local expected_n4=770
local expected_n5=693
local expected_n6=462
local validation_pass=1
forvalues m=1/6 {
    local y : word `m' of `outcomes'
    local s : word `m' of `samples'
    if `m'<=3 {
        quietly regress `y' exp15_check i.survey_id i.cohort_id [pweight=weight]  if data_part==1 & `s'==1 & `y'<. & exp15_check<. & weight>0,  vce(cluster country_id)
    }
    else {
        quietly regress `y' exp15_check i.survey_id i.cohort_id [aweight=weighted_women_n]  if data_part==2 & `y'<. & exp15_check<. & weighted_women_n>0,  vce(cluster country_id)
    }
    scalar __diff=abs(_b[exp15_check]-`expected_b`m'')
    local thispass=(__diff<1e-7 & e(N)==`expected_n`m'')
    if !`thispass' local validation_pass=0
    display "CUTOFF15_VALIDATION outcome=`y' b=" %16.11f _b[exp15_check]  " expected=" %16.11f `expected_b`m'' " absdiff=" %16.11f __diff  " N=" e(N) " expected_N=" `expected_n`m'' " PASS=`thispass'"
}
display "CUTOFF15_ALL_SIX_VALIDATION_PASS=`validation_pass'"
assert `validation_pass'==1
display "CUTOFFS_AND_PRETRENDS_COMPLETE"

log close
