version 19.5
clear all
set more off
set linesize 255
set varabbrev off
set processors 1
args input out
if `"`input'"'=="" local input "data/data.dta"
if `"`out'"'=="" local out "results"
capture mkdir "`out'"
local original "`out'/poisson"
local root "`out'/poisson_bounded_variation"
capture mkdir "`root'"
capture log close _all
log using "`root'/analysis.log", text replace
display "POISSON_BOUNDED_DEPARTURE_START"
display "RUN_DATE=`c(current_date)' RUN_TIME=`c(current_time)'"
which ppmlhdfe
confirm file "`input'"
tempname pointout boundout breakout
postfile `pointout' str8 outcome str8 heterogeneity double n_groups double n_fit double n_post double n_countries double theta double rate_treated double rate_untreated double weight_sum double min_weight double max_weight double bootstrap_se double normal_p double normal_lo double normal_hi double percentile_lo double percentile_hi long seed long replications using "`root'/baseline_verification.dta", replace
postfile `boundout' str8 outcome str8 heterogeneity double delta double estimate double bound_lo double bound_hi double normal_lo double normal_hi double percentile_lo double percentile_hi using "`root'/sensitivity_bounds.dta", replace
postfile `breakout' str8 outcome str8 heterogeneity double point_allowance double normal_allowance double percentile_allowance byte normal_zero_at_baseline byte pctile_zero_at_baseline using "`root'/zero_crossing_allowances.dta", replace
foreach outcome in births infant under5 {
    foreach het in event cohort time {
        use data_part primary_sample country_id survey_id female_birth_cohort candidate_reform_year weighted_women_n births_15_19_w_per1000 infant_deaths_fr_15_19_W_per1000 under5_deaths_fr_15_19_W_per1000 using "`input'", clear
        keep if data_part==2 & primary_sample==1
        sort country_id survey_id female_birth_cohort
        keep if !inlist(country_id,3,4)
        local y "births_15_19_w_per1000"
        if "`outcome'"=="infant" local y "infant_deaths_fr_15_19_W_per1000"
        if "`outcome'"=="under5" local y "under5_deaths_fr_15_19_W_per1000"
        keep if `y'<. & weighted_women_n>0 & !missing(country_id,survey_id,female_birth_cohort,candidate_reform_year)
        assert `y'>=0
        local n_groups=_N
        quietly levelsof country_id, local(countries)
        local n_countries: word count `countries'
        assert `n_countries'==16
        assert candidate_reform_year==2003 if country_id==14
        generate int t=female_birth_cohort+15
        generate int g=candidate_reform_year
        generate int event=t-g
        generate byte exposed=(event>=0 & event<.)
        generate double denominator=weighted_women_n/1000
        generate double count=`y'*denominator
        egen long surveyfe=group(country_id survey_id)
        generate int hcat=0
        if "`het'"=="event" {
            replace hcat=event+1 if inrange(event,0,4)
            replace hcat=6 if event>=5 & event<.
        }
        if "`het'"=="cohort" {
            egen int tmpcat=group(g) if exposed==1
            replace hcat=tmpcat if exposed==1
        }
        if "`het'"=="time" {
            egen int tmpcat=group(t) if exposed==1
            replace hcat=tmpcat if exposed==1
        }
        display "POINT_REFIT OUTCOME=`outcome' HETEROGENEITY=`het' N=`n_groups'"
        quietly ppmlhdfe count ib0.hcat, absorb(surveyfe t) exposure(denominator) vce(cluster country_id) d(fixed_effect_sum)
        local n_fit=e(N)
        predict double mu1, mu
        clonevar hkeep=hcat
        replace hcat=0 if exposed==1
        predict double mu0, mu
        replace hcat=hkeep
        generate double rate1=mu1/denominator
        generate double rate0=mu0/denominator
        generate double difference=rate1-rate0
        generate byte aggregate=(exposed==1 & difference<.)
        quietly summarize difference [aw=weighted_women_n] if aggregate, meanonly
        scalar theta=r(mean)
        local n_post=r(N)
        quietly summarize weighted_women_n if aggregate, meanonly
        generate double aggregate_weight=weighted_women_n/r(sum) if aggregate
        quietly summarize aggregate_weight, meanonly
        scalar weight_sum=r(sum)
        scalar min_weight=r(min)
        scalar max_weight=r(max)
        assert abs(weight_sum-1)<1e-12
        assert aggregate_weight>0 if aggregate
        quietly summarize rate1 [aw=weighted_women_n] if aggregate, meanonly
        scalar treated=r(mean)
        quietly summarize rate0 [aw=weighted_women_n] if aggregate, meanonly
        scalar untreated=r(mean)
        assert abs(theta-(treated-untreated))<1e-8
        keep country_id survey_id female_birth_cohort candidate_reform_year weighted_women_n t event exposed rate1 rate0 difference aggregate aggregate_weight
        save "`root'/predictions_`outcome'_`het'.dta", replace
        local saved "`original'/ppml_`outcome'_`het'"
        use "`saved'/summary.dta", clear
        assert _N==1
        assert valid==5000 & failed==0 & attempts==5000
        assert abs(att-theta)<1e-8
        assert n_groups==`n_groups' & n_post==`n_post'
        scalar original_se=se[1]
        scalar original_nlo=ci_normal_lo[1]
        scalar original_nhi=ci_normal_hi[1]
        scalar original_plo=ci_percentile_lo[1]
        scalar original_phi=ci_percentile_hi[1]
        local seed=seed[1]
        use "`saved'/valid_draws.dta", clear
        assert _N==5000
        assert valid==1 & rc==0 & att<.
        quietly summarize att
        scalar bootstrap_se=r(sd)
        scalar normal_p=2*normal(-abs(theta/bootstrap_se))
        scalar nlo=theta-invnormal(.975)*bootstrap_se
        scalar nhi=theta+invnormal(.975)*bootstrap_se
        quietly _pctile att, p(2.5 97.5)
        scalar plo=r(r1)
        scalar phi=r(r2)
        assert abs(bootstrap_se-original_se)<1e-10
        assert abs(nlo-original_nlo)<1e-8 & abs(nhi-original_nhi)<1e-8
        assert abs(plo-original_plo)<1e-8 & abs(phi-original_phi)<1e-8
        post `pointout' ("`outcome'") ("`het'") (`n_groups') (`n_fit') (`n_post') (`n_countries') (theta) (treated) (untreated) (weight_sum) (min_weight) (max_weight) (bootstrap_se) (normal_p) (nlo) (nhi) (plo) (phi) (`seed') (5000)
        foreach delta in 0 .25 .5 1 2 3 5 10 20 30 50 {
            generate double lower_draw=att-`delta'
            generate double upper_draw=att+`delta'
            quietly _pctile lower_draw, p(2.5)
            assert abs(r(r1)-(plo-`delta'))<1e-8
            quietly _pctile upper_draw, p(97.5)
            assert abs(r(r1)-(phi+`delta'))<1e-8
            post `boundout' ("`outcome'") ("`het'") (`delta') (theta) (theta-`delta') (theta+`delta') (nlo-`delta') (nhi+`delta') (plo-`delta') (phi+`delta')
            drop lower_draw upper_draw
        }
        post `breakout' ("`outcome'") ("`het'") (abs(theta)) (max(0,-nhi)) (max(0,-phi)) (nlo<=0 & nhi>=0) (plo<=0 & phi>=0)
        display "MODEL_VERIFIED OUTCOME=`outcome' HETEROGENEITY=`het' ATT=" %14.8f theta " NORMAL_95CI=" %14.8f nlo " to " %14.8f nhi " PERCENTILE_95CI=" %14.8f plo " to " %14.8f phi " REUSED_VALID_DRAWS=5000"
    }
}
postclose `pointout'
postclose `boundout'
postclose `breakout'
use "`root'/baseline_verification.dta", clear
assert _N==9
export delimited using "`root'/baseline_verification.csv", replace
use "`root'/sensitivity_bounds.dta", clear
assert _N==99
bysort outcome heterogeneity (delta): assert normal_lo<=normal_lo[_n-1] & normal_hi>=normal_hi[_n-1] if _n>1
bysort outcome heterogeneity (delta): assert percentile_lo<=percentile_lo[_n-1] & percentile_hi>=percentile_hi[_n-1] if _n>1
export delimited using "`root'/sensitivity_bounds.csv", replace
format estimate bound_lo bound_hi normal_lo normal_hi percentile_lo percentile_hi %9.3f
list outcome heterogeneity delta estimate bound_lo bound_hi normal_lo normal_hi percentile_lo percentile_hi, noobs separator(11) abbreviate(20)
use "`root'/zero_crossing_allowances.dta", clear
assert _N==9
export delimited using "`root'/zero_crossing_allowances.csv", replace
list, noobs abbreviate(32)
display "POISSON_BOUNDED_DEPARTURE_COMPLETE MODELS=9 GRID_ROWS=99 REUSED_DRAWS_PER_MODEL=5000 ALL_CHECKS_PASSED"
log close
