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
local output "`out'/primary_and_event_studies"
capture mkdir "`output'"
local results "`out'"
local dta "`input'"
local requested_reps 5000
local reps 5000
capture log close _all
log using "`output'/primary_and_event_studies.log", text replace
capture program drop showb
program define showb
    args x
    capture scalar __b=_b[`x']
    if _rc {
        display "COEFFICIENT_NOT_AVAILABLE `x'"
    }
    if !_rc {
        scalar __se=_se[`x']
        scalar __t=__b/__se
        capture scalar __p=2*ttail(e(df_r),abs(__t))
        if _rc scalar __p=.
        display "N=" %12.0f e(N) " COEF_`x'=" %12.6f __b " SE=" %12.6f __se " P=" %9.5f __p
        capture scalar drop __b __se __t __p
    }
end
capture program drop runmain
program define runmain
    args y a s label
    capture confirm variable `y'
    if _rc display "MISSING_OUTCOME `label' `y'"
    if !_rc {
        capture confirm variable `a'
        if _rc display "MISSING_EXPOSURE `label' `a'"
        if !_rc {
            capture confirm variable `s'
            if _rc display "MISSING_SAMPLE `label' `s'"
            if !_rc {
                display "MODEL_START `label' OUTCOME `y' EXPOSURE `a' SPEC SURVEY_COHORT_FE"
                capture noisily regress `y' `a' i.survey_id i.cohort_id [pweight=weight] if data_part==1 & `s'==1 & `a'<. & `y'<. & weight>0, vce(cluster country_id)
                if !_rc showb `a'
            }
        }
    }
end
capture program drop runevent
program define runevent
    args y s label
    capture confirm variable `y'
    if _rc display "MISSING_EVENT_OUTCOME `label' `y'"
    if !_rc {
        capture confirm variable `s'
        if _rc display "MISSING_EVENT_SAMPLE `label' `s'"
        if !_rc {
            capture drop event15_cat
            quietly gen event15_cat=event15_bin+11 if event15_bin<.
            display "EVENT_DISTRIBUTION `label' OUTCOME `y' SAMPLE `s'"
            capture noisily tab event15_bin if data_part==1 & `s'==1 & `y'<., missing
            display "EVENT_STUDY_START `label' OUTCOME `y' FULL_WINDOW REF_EVENT_-1"
            capture noisily regress `y' ib10.event15_cat i.survey_id i.cohort_id [pweight=weight] if data_part==1 & `s'==1 & `y'<. & event15_cat<. & weight>0, vce(cluster country_id)
            if !_rc {
                forvalues j=1/21 {
                    if `j'!=10 {
                        local ev=`j'-11
                        capture scalar __b=_b[`j'.event15_cat]
                        if !_rc {
                            scalar __se=_se[`j'.event15_cat]
                            capture scalar __p=2*ttail(e(df_r),abs(__b/__se))
                            if _rc scalar __p=.
                            display "EVENT15=" %6.0f `ev' " COEF=" %12.6f __b " SE=" %12.6f __se " P=" %9.5f __p
                            capture scalar drop __b __se __p
                        }
                    }
                }
                capture test 1.event15_cat 2.event15_cat 3.event15_cat 4.event15_cat 5.event15_cat 6.event15_cat 7.event15_cat 8.event15_cat 9.event15_cat
                if !_rc display "PRETREND_FULL_P=" %9.5f r(p)
            }
            display "EVENT_STUDY_START `label' OUTCOME `y' NARROW_-6_TO_8 REF_EVENT_-1"
            capture noisily regress `y' ib10.event15_cat i.survey_id i.cohort_id [pweight=weight] if data_part==1 & `s'==1 & `y'<. & event15_bin>=-6 & event15_bin<=8 & event15_cat<. & weight>0, vce(cluster country_id)
            if !_rc {
                forvalues j=5/19 {
                    if `j'!=10 {
                        local ev=`j'-11
                        capture scalar __b=_b[`j'.event15_cat]
                        if !_rc {
                            scalar __se=_se[`j'.event15_cat]
                            capture scalar __p=2*ttail(e(df_r),abs(__b/__se))
                            if _rc scalar __p=.
                            display "EVENT15=" %6.0f `ev' " COEF=" %12.6f __b " SE=" %12.6f __se " P=" %9.5f __p
                            capture scalar drop __b __se __p
                        }
                    }
                }
                capture test 5.event15_cat 6.event15_cat 7.event15_cat 8.event15_cat 9.event15_cat
                if !_rc display "PRETREND_NARROW_P=" %9.5f r(p)
            }
        }
    }
end
capture program drop leaveone
program define leaveone
    args y s label
    capture confirm variable `y'
    if !_rc {
        capture confirm variable `s'
        if !_rc {
            levelsof country_id if data_part==1 & `s'==1 & `y'<. & exposed_before15<. & weight>0, local(countries)
            foreach c of local countries {
                display "LEAVE_ONE_COUNTRY_OUT_START `label' OUTCOME `y' DROPPED_COUNTRY_ID=`c'"
                capture noisily regress `y' exposed_before15 i.survey_id i.cohort_id [pweight=weight] if data_part==1 & `s'==1 & `y'<. & exposed_before15<. & weight>0 & country_id!=`c', vce(cluster country_id)
                if !_rc showb exposed_before15
            }
        }
    }
end
capture program drop runburden
program define runburden
    args y s label
    capture confirm variable `y'
    if _rc display "MISSING_BURDEN_OUTCOME `label' `y'"
    if !_rc {
        capture confirm variable `s'
        if _rc display "MISSING_BURDEN_SAMPLE `label' `s'"
        if !_rc {
            count if data_part==2 & `s'==1 & `y'<. & exposed_before15<. & weighted_women_n>0
            if r(N)>0 {
                display "BURDEN_MODEL_START `label' OUTCOME `y' SPEC SURVEY_COHORT_FE"
                capture noisily regress `y' exposed_before15 i.survey_id i.cohort_id [aw=weighted_women_n] if data_part==2 & `s'==1 & `y'<. & exposed_before15<. & weighted_women_n>0, vce(cluster country_id)
                if !_rc showb exposed_before15
            }
            if r(N)==0 display "NO_NONMISSING_BURDEN_OUTCOME `label' `y'"
        }
    }
end
capture program drop runburdenevent
program define runburdenevent
    args y s label
    capture confirm variable `y'
    if _rc display "MISSING_BURDEN_EVENT_OUTCOME `label' `y'"
    if !_rc {
        capture confirm variable `s'
        if _rc display "MISSING_BURDEN_EVENT_SAMPLE `label' `s'"
        if !_rc {
            count if data_part==2 & `s'==1 & `y'<. & exposed_before15<. & event15_bin<. & weighted_women_n>0
            if r(N)>0 {
                capture drop event15_cat
                quietly gen event15_cat=event15_bin+11 if event15_bin<.
                display "BURDEN_EVENT_DISTRIBUTION `label' OUTCOME `y' SAMPLE `s'"
                capture noisily tab event15_bin if data_part==2 & `s'==1 & `y'<. & weighted_women_n>0, missing
                display "BURDEN_EVENT_START `label' OUTCOME `y' FULL_WINDOW REF_EVENT_-1"
                capture noisily regress `y' ib10.event15_cat i.survey_id i.cohort_id [aw=weighted_women_n] if data_part==2 & `s'==1 & `y'<. & event15_cat<. & weighted_women_n>0, vce(cluster country_id)
                if !_rc {
                    forvalues j=1/21 {
                        if `j'!=10 {
                            local ev=`j'-11
                            capture scalar __b=_b[`j'.event15_cat]
                            if !_rc {
                                scalar __se=_se[`j'.event15_cat]
                                capture scalar __p=2*ttail(e(df_r),abs(__b/__se))
                                if _rc scalar __p=.
                                display "EVENT15=" %6.0f `ev' " COEF=" %12.6f __b " SE=" %12.6f __se " P=" %9.5f __p
                                capture scalar drop __b __se __p
                            }
                        }
                    }
                    capture test 1.event15_cat 2.event15_cat 3.event15_cat 4.event15_cat 5.event15_cat 6.event15_cat 7.event15_cat 8.event15_cat 9.event15_cat
                    if !_rc display "BURDEN_PRETREND_FULL_P=" %9.5f r(p)
                }
                display "BURDEN_EVENT_START `label' OUTCOME `y' NARROW_-6_TO_8 REF_EVENT_-1"
                capture noisily regress `y' ib10.event15_cat i.survey_id i.cohort_id [aw=weighted_women_n] if data_part==2 & `s'==1 & `y'<. & event15_bin>=-6 & event15_bin<=8 & event15_cat<. & weighted_women_n>0, vce(cluster country_id)
                if !_rc {
                    forvalues j=5/19 {
                        if `j'!=10 {
                            local ev=`j'-11
                            capture scalar __b=_b[`j'.event15_cat]
                            if !_rc {
                                scalar __se=_se[`j'.event15_cat]
                                capture scalar __p=2*ttail(e(df_r),abs(__b/__se))
                                if _rc scalar __p=.
                                display "EVENT15=" %6.0f `ev' " COEF=" %12.6f __b " SE=" %12.6f __se " P=" %9.5f __p
                                capture scalar drop __b __se __p
                            }
                        }
                    }
                    capture test 5.event15_cat 6.event15_cat 7.event15_cat 8.event15_cat 9.event15_cat
                    if !_rc display "BURDEN_PRETREND_NARROW_P=" %9.5f r(p)
                }
            }
            if r(N)==0 display "NO_NONMISSING_BURDEN_EVENT_OUTCOME `label' `y'"
        }
    }
end
capture program drop runchild
program define runchild
    args y w s label
    capture confirm variable `y'
    if _rc display "MISSING_CHILD_OUTCOME `label' `y'"
    if !_rc {
        capture confirm variable `w'
        if _rc display "MISSING_CHILD_WEIGHT `label' `w'"
        if !_rc {
            capture confirm variable `s'
            if _rc display "MISSING_CHILD_SAMPLE `label' `s'"
            if !_rc {
                count if data_part==3 & `s'==1 & `y'<. & exposed_before15<. & `w'>0
                if r(N)>0 {
                    display "CHILD_MODEL_START `label' OUTCOME `y' SPEC ADJUSTED_WEIGHTED_CELL"
                    capture noisily regress `y' exposed_before15 i.child_sex birth_order_cat multiple_birth maternal_age_at_birth i.survey_id i.cohort_id [aw=`w'] if data_part==3 & `s'==1 & `y'<. & exposed_before15<. & maternal_age_at_birth>=10 & maternal_age_at_birth<=49 & `w'>0, vce(cluster country_id)
                    if !_rc showb exposed_before15
                    display "CHILD_MOTHER_15_24 OUTCOME `y'"
                    capture noisily regress `y' exposed_before15 i.child_sex birth_order_cat multiple_birth maternal_age_at_birth i.survey_id i.cohort_id [aw=`w'] if data_part==3 & `s'==1 & `y'<. & exposed_before15<. & birth_to_mother_15_24==1 & maternal_age_at_birth>=10 & maternal_age_at_birth<=49 & `w'>0, vce(cluster country_id)
                    if !_rc showb exposed_before15
                    display "CHILD_MOTHER_15_19 OUTCOME `y'"
                    capture noisily regress `y' exposed_before15 i.child_sex birth_order_cat multiple_birth maternal_age_at_birth i.survey_id i.cohort_id [aw=`w'] if data_part==3 & `s'==1 & `y'<. & exposed_before15<. & birth_to_mother_15_19==1 & maternal_age_at_birth>=10 & maternal_age_at_birth<=49 & `w'>0, vce(cluster country_id)
                    if !_rc showb exposed_before15
                }
            }
        }
    }
end
use "`input'", clear
keep if primary_sample==1
capture drop overlap_main overlap_analytic_union overlap_analytic_birth18 overlap_analytic_birth20 overlap_sample20_24 burden_overlap child_overlap
generate byte overlap_main=1 if data_part<.
generate byte overlap_analytic_union=(data_part==1 & analytic_union==1)
generate byte overlap_analytic_birth18=(data_part==1 & analytic_birth18==1)
generate byte overlap_analytic_birth20=(data_part==1 & analytic_birth20==1)
generate byte overlap_sample20_24=(data_part==1 & sample_20_24==1)
generate byte burden_overlap=(data_part==2)
generate byte child_overlap=(data_part==3)
tab country_id exposed_before15 if data_part==1, row
mean union_before18 birth_before18 birth_before20 [pw=weight] if data_part==1, over(exposed_before15)
display "SECTION 7 ALTERNATIVE EXPOSURE BEFORE18 IN OVERLAP SAMPLE"
runmain union_before18 exposed_before18 overlap_analytic_union OVERLAP_UNION_EXPOSED18_ALT
runmain birth_before18 exposed_before18 overlap_analytic_birth18 OVERLAP_BIRTH18_EXPOSED18_ALT
runmain birth_before20 exposed_before18 overlap_analytic_birth20 OVERLAP_BIRTH20_EXPOSED18_ALT
display "SECTION 8 AGE 20 TO 24 SENSITIVITY IN OVERLAP SAMPLE"
display "AGE20_24_SENSITIVITY_SPEC=COUNTRY_COHORT_FE_RECOMMENDED_FOR_YOUNG_SUBSAMPLE"
display "MODEL_START OVERLAP_UNION_SAMPLE20_24_PREFERRED OUTCOME union_before18 EXPOSURE exposed_before15 SPEC COUNTRY_COHORT_FE"
capture noisily regress union_before18 exposed_before15 i.country_id i.cohort_id [pweight=weight] if data_part==1 & overlap_sample20_24==1 & union_before18<. & exposed_before15<. & weight>0, vce(cluster country_id)
if !_rc showb exposed_before15
display "MODEL_START OVERLAP_BIRTH18_SAMPLE20_24_PREFERRED OUTCOME birth_before18 EXPOSURE exposed_before15 SPEC COUNTRY_COHORT_FE"
capture noisily regress birth_before18 exposed_before15 i.country_id i.cohort_id [pweight=weight] if data_part==1 & overlap_sample20_24==1 & birth_before18<. & exposed_before15<. & weight>0, vce(cluster country_id)
if !_rc showb exposed_before15
display "MODEL_START OVERLAP_BIRTH20_SAMPLE20_24_PREFERRED OUTCOME birth_before20 EXPOSURE exposed_before15 SPEC COUNTRY_COHORT_FE"
capture noisily regress birth_before20 exposed_before15 i.country_id i.cohort_id [pweight=weight] if data_part==1 & overlap_sample20_24==1 & birth_before20<. & exposed_before15<. & weight>0, vce(cluster country_id)
if !_rc showb exposed_before15
display "SECTION 8B AGE 20 TO 24 SENSITIVITY RURAL_WEALTH_ADJUSTED"
display "MODEL_START OVERLAP_UNION_SAMPLE20_24_RURAL_WEALTH OUTCOME union_before18 EXPOSURE exposed_before15 SPEC SURVEY_COHORT_FE_RURAL_WEALTH"
capture noisily regress union_before18 exposed_before15 i.rural i.wealth_quintile i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_sample20_24==1 & union_before18<. & exposed_before15<. & rural<. & wealth_quintile<. & weight>0, vce(cluster country_id)
if !_rc showb exposed_before15
display "MODEL_START OVERLAP_BIRTH18_SAMPLE20_24_RURAL_WEALTH OUTCOME birth_before18 EXPOSURE exposed_before15 SPEC SURVEY_COHORT_FE_RURAL_WEALTH"
capture noisily regress birth_before18 exposed_before15 i.rural i.wealth_quintile i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_sample20_24==1 & birth_before18<. & exposed_before15<. & rural<. & wealth_quintile<. & weight>0, vce(cluster country_id)
if !_rc showb exposed_before15
display "MODEL_START OVERLAP_BIRTH20_SAMPLE20_24_RURAL_WEALTH OUTCOME birth_before20 EXPOSURE exposed_before15 SPEC SURVEY_COHORT_FE_RURAL_WEALTH"
capture noisily regress birth_before20 exposed_before15 i.rural i.wealth_quintile i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_sample20_24==1 & birth_before20<. & exposed_before15<. & rural<. & wealth_quintile<. & weight>0, vce(cluster country_id)
if !_rc showb exposed_before15
display "SECTION 8C AGE 20 TO 24 NARROW WINDOW RURAL_WEALTH_ADJUSTED"
display "MODEL_START OVERLAP_UNION_SAMPLE20_24_NARROW_RURAL_WEALTH OUTCOME union_before18 EXPOSURE exposed_before15 SPEC SURVEY_COHORT_FE_RURAL_WEALTH_NARROW"
capture noisily regress union_before18 exposed_before15 i.rural i.wealth_quintile i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_sample20_24==1 & event15_bin>=-6 & event15_bin<=8 & union_before18<. & exposed_before15<. & rural<. & wealth_quintile<. & weight>0, vce(cluster country_id)
if !_rc showb exposed_before15
display "MODEL_START OVERLAP_BIRTH18_SAMPLE20_24_NARROW_RURAL_WEALTH OUTCOME birth_before18 EXPOSURE exposed_before15 SPEC SURVEY_COHORT_FE_RURAL_WEALTH_NARROW"
capture noisily regress birth_before18 exposed_before15 i.rural i.wealth_quintile i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_sample20_24==1 & event15_bin>=-6 & event15_bin<=8 & birth_before18<. & exposed_before15<. & rural<. & wealth_quintile<. & weight>0, vce(cluster country_id)
if !_rc showb exposed_before15
display "MODEL_START OVERLAP_BIRTH20_SAMPLE20_24_NARROW_RURAL_WEALTH OUTCOME birth_before20 EXPOSURE exposed_before15 SPEC SURVEY_COHORT_FE_RURAL_WEALTH_NARROW"
capture noisily regress birth_before20 exposed_before15 i.rural i.wealth_quintile i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_sample20_24==1 & event15_bin>=-6 & event15_bin<=8 & birth_before20<. & exposed_before15<. & rural<. & wealth_quintile<. & weight>0, vce(cluster country_id)
if !_rc showb exposed_before15
display "SECTION 10 EVENT STUDY RECOMMENDED OVERLAP SAMPLE"
runevent union_before18 overlap_analytic_union OVERLAP_UNION
runevent birth_before18 overlap_analytic_birth18 OVERLAP_BIRTH18
runevent birth_before20 overlap_analytic_birth20 OVERLAP_BIRTH20
display "SECTION 11 PLACEBO IN OVERLAP SAMPLE"
capture drop placebo_m5
quietly gen placebo_m5=.
quietly replace placebo_m5=1 if data_part==1 & overlap_main==1 & candidate_reform_year-5<=female_birth_year+15 & candidate_reform_year<.
quietly replace placebo_m5=0 if data_part==1 & overlap_main==1 & candidate_reform_year-5>female_birth_year+15 & candidate_reform_year<.
capture noisily regress union_before18 placebo_m5 i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_analytic_union==1 & weight>0, vce(cluster country_id)
if !_rc showb placebo_m5
capture noisily regress birth_before18 placebo_m5 i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_analytic_birth18==1 & weight>0, vce(cluster country_id)
if !_rc showb placebo_m5
capture drop placebo_m10
quietly gen placebo_m10=.
quietly replace placebo_m10=1 if data_part==1 & overlap_main==1 & candidate_reform_year-10<=female_birth_year+15 & candidate_reform_year<.
quietly replace placebo_m10=0 if data_part==1 & overlap_main==1 & candidate_reform_year-10>female_birth_year+15 & candidate_reform_year<.
capture noisily regress union_before18 placebo_m10 i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_analytic_union==1 & weight>0, vce(cluster country_id)
if !_rc showb placebo_m10
capture noisily regress birth_before18 placebo_m10 i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_analytic_birth18==1 & weight>0, vce(cluster country_id)
if !_rc showb placebo_m10
capture drop placebo_p5
quietly gen placebo_p5=.
quietly replace placebo_p5=1 if data_part==1 & overlap_main==1 & candidate_reform_year+5<=female_birth_year+15 & candidate_reform_year<.
quietly replace placebo_p5=0 if data_part==1 & overlap_main==1 & candidate_reform_year+5>female_birth_year+15 & candidate_reform_year<.
capture noisily regress union_before18 placebo_p5 i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_analytic_union==1 & weight>0, vce(cluster country_id)
if !_rc showb placebo_p5
capture noisily regress birth_before18 placebo_p5 i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_analytic_birth18==1 & weight>0, vce(cluster country_id)
if !_rc showb placebo_p5
display "SECTION 12 LEAVE ONE COUNTRY OUT IN OVERLAP SAMPLE"
leaveone union_before18 overlap_analytic_union OVERLAP_UNION
leaveone birth_before18 overlap_analytic_birth18 OVERLAP_BIRTH18
leaveone birth_before20 overlap_analytic_birth20 OVERLAP_BIRTH20
display "SECTION 15 CHILD LEVEL WEIGHTED CELLS IN OVERLAP SAMPLE"
capture noisily summarize child_cell_n child_cell_w_all maternal_age_at_birth neonatal_death_any infant_death_any under5_death_any neonatal_death_fullrisk infant_death_fullrisk under5_death_fullrisk if data_part==3 & overlap_main==1
runchild neonatal_death_any w_neonatal_death_any child_overlap CHILD_NEONATAL_ANY_OVERLAP
runchild infant_death_any w_infant_death_any child_overlap CHILD_INFANT_ANY_OVERLAP
runchild under5_death_any w_under5_death_any child_overlap CHILD_UNDER5_ANY_OVERLAP
runchild neonatal_death_fullrisk w_neonatal_death_fullrisk child_overlap CHILD_NEONATAL_FULLRISK_OVERLAP
runchild infant_death_fullrisk w_infant_death_fullrisk child_overlap CHILD_INFANT_FULLRISK_OVERLAP
runchild under5_death_fullrisk w_under5_death_fullrisk child_overlap CHILD_UNDER5_FULLRISK_OVERLAP
display "SECTION 16 FEMALE COHORT BURDEN V2 OVERLAP SAMPLE"
scalar eff_births1519=.
scalar eff_infant1519=.
scalar eff_under51519=.
scalar eff_births1524=.
scalar eff_infant1524=.
capture noisily summarize women_n weighted_women_n births_15_19_w_per1000 infant_deaths_fr_15_19_W_per1000 under5_deaths_fr_15_19_W_per1000 births_15_24_w_per1000 infant_deaths_fr_15_24_W_per1000 if data_part==2 & overlap_main==1
foreach y in births_15_19_w_per1000 infant_deaths_fr_15_19_W_per1000 under5_deaths_fr_15_19_W_per1000 births_15_24_w_per1000 infant_deaths_fr_15_24_W_per1000 {
    runburden `y' burden_overlap BURDEN_OVERLAP_`y'
    runburdenevent `y' burden_overlap BURDEN_OVERLAP_`y'
}
capture noisily regress births_15_19_w_per1000 exposed_before15 i.survey_id i.cohort_id [aw=weighted_women_n] if data_part==2 & overlap_main==1 & births_15_19_w_per1000<. & exposed_before15<. & weighted_women_n>0, vce(cluster country_id)
if !_rc scalar eff_births1519=-_b[exposed_before15]
capture noisily regress infant_deaths_fr_15_19_W_per1000 exposed_before15 i.survey_id i.cohort_id [aw=weighted_women_n] if data_part==2 & overlap_main==1 & infant_deaths_fr_15_19_W_per1000<. & exposed_before15<. & weighted_women_n>0, vce(cluster country_id)
if !_rc scalar eff_infant1519=-_b[exposed_before15]
capture noisily regress under5_deaths_fr_15_19_W_per1000 exposed_before15 i.survey_id i.cohort_id [aw=weighted_women_n] if data_part==2 & overlap_main==1 & under5_deaths_fr_15_19_W_per1000<. & exposed_before15<. & weighted_women_n>0, vce(cluster country_id)
if !_rc scalar eff_under51519=-_b[exposed_before15]
capture noisily regress births_15_24_w_per1000 exposed_before15 i.survey_id i.cohort_id [aw=weighted_women_n] if data_part==2 & overlap_main==1 & births_15_24_w_per1000<. & exposed_before15<. & weighted_women_n>0, vce(cluster country_id)
if !_rc scalar eff_births1524=-_b[exposed_before15]
capture noisily regress infant_deaths_fr_15_24_W_per1000 exposed_before15 i.survey_id i.cohort_id [aw=weighted_women_n] if data_part==2 & overlap_main==1 & infant_deaths_fr_15_24_W_per1000<. & exposed_before15<. & weighted_women_n>0, vce(cluster country_id)
if !_rc scalar eff_infant1524=-_b[exposed_before15]
display "BURDEN_EFFECT_BIRTHS15_19_PER1000=" %12.6f eff_births1519
display "BURDEN_EFFECT_INFANT15_19_PER1000=" %12.6f eff_infant1519
display "BURDEN_EFFECT_UNDER515_19_PER1000=" %12.6f eff_under51519
display "BURDEN_EFFECT_BIRTHS15_24_PER1000=" %12.6f eff_births1524
display "BURDEN_EFFECT_INFANT15_24_PER1000=" %12.6f eff_infant1524
preserve
keep if data_part==2
tempname DESC
tempfile desc
postfile `DESC' str45 outcome byte exposed double groups women_weight mean using `desc', replace
foreach y in births_15_19_w_per1000 infant_deaths_fr_15_19_W_per1000 under5_deaths_fr_15_19_W_per1000 births_15_24_w_per1000 infant_deaths_fr_15_24_W_per1000 {
forvalues x=0/1 {
quietly summarize `y' [aw=weighted_women_n] if exposed_before15==`x', meanonly
local n=r(N)
local m=r(mean)
quietly summarize weighted_women_n if exposed_before15==`x' & `y'<., meanonly
post `DESC' ("`y'") (`x') (`n') (r(sum)) (`m')
}
}
postclose `DESC'
use `desc', clear
export delimited using "`output'/group_descriptives.csv", replace
restore
preserve
keep if data_part==1
collapse (count) women=weight (sum) exposed=exposed_before15 (sum) total_weight=weight, by(country_id country_name iso3 survey_id survey_year)
export delimited using "`output'/survey_sample_counts.csv", replace
restore
log close
