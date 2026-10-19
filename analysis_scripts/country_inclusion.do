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
local output "`out'/countries"
capture mkdir "`output'"
local results "`out'"
local dta "`input'"
local requested_reps 5000
local reps 5000
capture log close _all
log using "`output'/country_inclusion.log", text replace
local csv "`output'/country_inclusion.csv"
tempfile hybrid_results
use "`input'", clear
keep if inlist(data_part,1,2)
assert _N==362401
assert missing(candidate_reform_year) if country_id==3
assert has_candidate_reform==0 if country_id==3
assert exposed_before15==0 if country_id==3
replace analytic_union=inrange(age,20,29) & union_before18<. & weight>0 if country_id==3 & data_part==1
replace analytic_birth18=inrange(age,20,29) & birth_before18<. & weight>0 if country_id==3 & data_part==1
replace analytic_birth20=inrange(age,20,29) & birth_before20<. & weight>0 if country_id==3 & data_part==1
local samples "base16 plus_COG plus_EGY all18"
local outcomes "union_before18 birth_before18 birth_before20 births_15_19_w_per1000 infant_deaths_fr_15_19_W_per1000 under5_deaths_fr_15_19_W_per1000"
local flags "analytic_union analytic_birth18 analytic_birth20 burden burden burden"
local domains "individual individual individual burden burden burden"
local scales "proportion proportion proportion per_1000_women per_1000_women per_1000_women"
local expected_N "311260 314592 314592 770 693 462"
local expected_b "-0.064164116 -0.052670793 -0.079727733 -50.081071091 -6.978656526 -10.723971581"
tempname RES
postfile `RES' str16 construction str10 sample str12 added_countries  str45 outcome str12 domain str22 scale double N double clusters  double support_n double zero_only_n double one_only_n double b  double se double ci_lo double ci_hi double p_cluster double wcb_p  double wcb_ci_lo double wcb_ci_hi double wcb_reps  str244 all_countries str244 support_countries  str244 zero_only_countries str244 one_only_countries  using `hybrid_results', replace
forvalues s=1/4 {
    local sample : word `s' of `samples'
    local keepcond "!inlist(country_id,3,4)"
    local added "none"
    if `s'==2 {
        local keepcond "country_id!=4"
        local added "COG"
    }
    if `s'==3 {
        local keepcond "country_id!=3"
        local added "EGY"
    }
    if `s'==4 {
        local keepcond "country_id<."
        local added "COG_EGY"
    }
    forvalues m=1/6 {
        local y : word `m' of `outcomes'
        local flag : word `m' of `flags'
        local domain : word `m' of `domains'
        local scale : word `m' of `scales'
        if `m'<=3 {
            quietly regress `y' exposed_before15 i.survey_id i.cohort_id  [pweight=weight] if data_part==1 & `flag'==1 &  (`keepcond') & `y'<. & exposed_before15<. & weight>0,  vce(cluster country_id)
        }
        else {
            quietly regress `y' exposed_before15 i.survey_id i.cohort_id  [aweight=weighted_women_n] if data_part==2 &  (`keepcond') & `y'<. & exposed_before15<. &  weighted_women_n>0, vce(cluster country_id)
        }
        scalar __N=e(N)
        scalar __b=_b[exposed_before15]
        scalar __se=_se[exposed_before15]
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
            quietly levelsof iso3 if `esample'==1 & country_id==`c',  local(thisiso) clean
            local alliso "`alliso' `thisiso'"
            quietly count if `esample'==1 & country_id==`c' &  exposed_before15==0
            local n0=r(N)
            quietly count if `esample'==1 & country_id==`c' &  exposed_before15==1
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
        if `s'==1 {
            local gateN : word `m' of `expected_N'
            local gateB : word `m' of `expected_b'
            assert __N==`gateN'
            assert abs(__b-(`gateB'))<0.0000005
            assert `clusters'==16
            display "BASE16_REPRODUCTION_PASS outcome=`y' N=" %12.0f __N  " b=" %14.9f __b
        }
        local bootseed=20260708
        capture noisily boottest exposed_before15, cluster(country_id)  reps(`reps') seed(`bootseed') nograph
        if _rc {
            scalar __wp=.
            scalar __wlo=.
            scalar __whi=.
            scalar __wreps=.
            display as error "BOOTTEST_FAILED sample=`sample' outcome=`y' RC=" _rc
        }
        else {
            scalar __wp=r(p)
            matrix __wci=r(CI)
            scalar __wlo=__wci[1,1]
            scalar __whi=__wci[1,2]
            scalar __wreps=r(reps)
        }
        display "HYBRID_RESULT sample=`sample' added=`added' outcome=`y'"  " N=" %12.0f __N " clusters=`clusters' support_n=`support_n'"  " b=" %14.9f __b " se=" %14.9f __se  " p_cluster=" %10.7f __p " wcb_p=" %10.7f __wp
        post `RES' ("complete") ("`sample'") ("`added'") ("`y'")  ("`domain'") ("`scale'") (__N) (`clusters') (`support_n')  (`zero_n') (`one_n') (__b) (__se) (__lo) (__hi) (__p)  (__wp) (__wlo) (__whi) (__wreps) ("`alliso'")  ("`supportiso'") ("`zeroiso'") ("`oneiso'")
        drop `esample'
    }
}
postclose `RES'
preserve
use `hybrid_results', clear
assert _N==24
assert !missing(N,clusters,b,se,ci_lo,ci_hi,p_cluster,wcb_p,wcb_reps)
assert wcb_reps==`reps'
assert inrange(wcb_p,0,1)
isid sample outcome
display "ALL_24_MODELS_AND_5000_REPLICATIONS_PASS"
sort domain outcome sample
format %21.15g N clusters support_n zero_only_n one_only_n b se ci_lo  ci_hi p_cluster wcb_p wcb_ci_lo wcb_ci_hi wcb_reps
export delimited using "`csv'", replace datafmt
list sample added_countries outcome N clusters support_n b se ci_lo ci_hi  p_cluster wcb_p wcb_reps, noobs abbreviate(32)
restore
display "PERSISTENT_OUTPUT_LOG=`log'"
display "PERSISTENT_OUTPUT_CSV=`csv'"
display "NO_COMBINED_MICRODATA_SAVED=YES"
display "COUNTRY_INCLUSION_COMPLETE"

log close
