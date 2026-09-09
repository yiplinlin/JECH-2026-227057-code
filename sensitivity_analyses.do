* The scripts use Stata 19.5 and prepared datasets stored in the data directory.
* The main input is data.dta.
* Country-inclusion analyses additionally require data_women.dta and data_cohort.dta.
* Run the scripts from the repository root.
* Outputs are saved in the results directory.

version 19.5
set more off
set linesize 255
set varabbrev off
capture log close
args task
if "`task'"=="" local task "all"
if !inlist("`task'","all","standard","cutoffs","countries") {
    exit 198
}
local root "`c(pwd)'"
local input "`root'/data/data.dta"
confirm file "`input'"
capture mkdir "`root'/results"
local output "`root'/results/module4"
capture mkdir "`output'"

if inlist("`task'","all","standard") {
    use "`input'", clear
    log using "`output'/standard_sensitivity.log", text replace
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
                    display "MODEL_START `label' OUTCOME `y' EXPOSURE `a' SPEC BASIC"
                    capture noisily regress `y' `a' [pweight=weight] if data_part==1 & `s'==1 & `a'<. & `y'<. & weight>0, vce(cluster country_id)
                    if !_rc showb `a'
                    display "MODEL_START `label' OUTCOME `y' EXPOSURE `a' SPEC COUNTRY_COHORT_FE"
                    capture noisily regress `y' `a' i.country_id i.cohort_id [pweight=weight] if data_part==1 & `s'==1 & `a'<. & `y'<. & weight>0, vce(cluster country_id)
                    if !_rc showb `a'
                    display "MODEL_START `label' OUTCOME `y' EXPOSURE `a' SPEC SURVEY_COHORT_FE"
                    capture noisily regress `y' `a' i.survey_id i.cohort_id [pweight=weight] if data_part==1 & `s'==1 & `a'<. & `y'<. & weight>0, vce(cluster country_id)
                    if !_rc showb `a'
                }
            }
        }
    end
    capture program drop runmainfe
    program define runmainfe
        args y a s label
        capture confirm variable `y'
        if !_rc {
            capture confirm variable `a'
            if !_rc {
                capture confirm variable `s'
                if !_rc {
                    display "MODEL_START `label' OUTCOME `y' EXPOSURE `a' SPEC SURVEY_COHORT_FE"
                    capture noisily regress `y' `a' i.survey_id i.cohort_id [pweight=weight] if data_part==1 & `s'==1 & `a'<. & `y'<. & weight>0, vce(cluster country_id)
                    if !_rc showb `a'
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
                        display "CHILD_MODEL_START `label' OUTCOME `y' SPEC BASIC_WEIGHTED_CELL"
                        capture noisily regress `y' exposed_before15 [aw=`w'] if data_part==3 & `s'==1 & `y'<. & exposed_before15<. & `w'>0, vce(cluster country_id)
                        if !_rc showb exposed_before15
                        display "CHILD_MODEL_START `label' OUTCOME `y' SPEC SURVEY_COHORT_FE_WEIGHTED_CELL"
                        capture noisily regress `y' exposed_before15 i.survey_id i.cohort_id [aw=`w'] if data_part==3 & `s'==1 & `y'<. & exposed_before15<. & `w'>0, vce(cluster country_id)
                        if !_rc showb exposed_before15
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
                    display "BURDEN_MODEL_START `label' OUTCOME `y' SPEC BASIC"
                    capture noisily regress `y' exposed_before15 [aw=weighted_women_n] if data_part==2 & `s'==1 & `y'<. & exposed_before15<. & weighted_women_n>0, vce(cluster country_id)
                    if !_rc showb exposed_before15
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
    capture confirm variable data_part
    if _rc {
        display as error "Required analysis input lacks data_part"
        log close
        exit 498
    }
    capture drop overlap_main
    quietly gen byte overlap_main=1 if data_part<.
    capture drop sparse_exposure_country
    quietly gen byte sparse_exposure_country=0 if data_part<.
    capture drop no_woman_firststage_country
    quietly gen byte no_woman_firststage_country=0 if data_part<.
    capture drop any_birth_under18_h
    quietly gen byte any_birth_under18_h=.
    quietly replace any_birth_under18_h=birth_before18 if data_part==1 & birth_before18<.
    capture drop any_birth_under20_h
    quietly gen byte any_birth_under20_h=.
    quietly replace any_birth_under20_h=birth_before20 if data_part==1 & birth_before20<.
    capture drop low_wealth
    quietly gen byte low_wealth=.
    quietly replace low_wealth=1 if data_part==1 & wealth_quintile<=2 & wealth_quintile<.
    quietly replace low_wealth=0 if data_part==1 & wealth_quintile>2 & wealth_quintile<.
    capture drop high_wealth
    quietly gen byte high_wealth=.
    quietly replace high_wealth=1 if data_part==1 & wealth_quintile>=4 & wealth_quintile<.
    quietly replace high_wealth=0 if data_part==1 & wealth_quintile<4 & wealth_quintile<.
    capture drop country_exposure_share15
    quietly gen double country_exposure_share15=.
    capture drop country_union18_baseline
    quietly gen double country_union18_baseline=.
    levelsof country_id if data_part==1 & country_id<. & weight>0 & exposed_before15<., local(country_levels)
    foreach c of local country_levels {
        quietly summarize exposed_before15 [aw=weight] if data_part==1 & country_id==`c' & exposed_before15<. & weight>0
        quietly replace country_exposure_share15=r(mean) if country_id==`c'
        quietly summarize union_before18 [aw=weight] if data_part==1 & country_id==`c' & union_before18<. & weight>0
        quietly replace country_union18_baseline=r(mean) if country_id==`c'
    }
    capture drop __country_tag
    quietly egen __country_tag=tag(country_id) if data_part==1 & overlap_main==1 & country_union18_baseline<.
    quietly summarize country_union18_baseline if __country_tag==1, detail
    scalar __cm_median=r(p50)
    capture drop high_cm_baseline
    quietly gen byte high_cm_baseline=.
    quietly replace high_cm_baseline=1 if country_union18_baseline>=__cm_median & country_union18_baseline<.
    quietly replace high_cm_baseline=0 if country_union18_baseline<__cm_median
    capture drop __country_tag
    capture drop woman_all
    quietly gen byte woman_all=1 if data_part==1
    quietly replace woman_all=0 if data_part!=1
    capture drop woman_overlap
    quietly gen byte woman_overlap=1 if data_part==1 & overlap_main==1
    quietly replace woman_overlap=0 if data_part!=1 | overlap_main!=1
    capture drop overlap_analytic_union
    quietly gen byte overlap_analytic_union=1 if data_part==1 & overlap_main==1 & analytic_union==1
    quietly replace overlap_analytic_union=0 if data_part!=1 | overlap_main!=1 | analytic_union!=1
    capture drop overlap_analytic_birth18
    quietly gen byte overlap_analytic_birth18=1 if data_part==1 & overlap_main==1 & analytic_birth18==1
    quietly replace overlap_analytic_birth18=0 if data_part!=1 | overlap_main!=1 | analytic_birth18!=1
    capture drop overlap_analytic_birth20
    quietly gen byte overlap_analytic_birth20=1 if data_part==1 & overlap_main==1 & analytic_birth20==1
    quietly replace overlap_analytic_birth20=0 if data_part!=1 | overlap_main!=1 | analytic_birth20!=1
    capture drop overlap_sample20_24
    quietly gen byte overlap_sample20_24=1 if data_part==1 & overlap_main==1 & sample_20_24==1
    quietly replace overlap_sample20_24=0 if data_part!=1 | overlap_main!=1 | sample_20_24!=1
    capture drop burden_overlap
    quietly gen byte burden_overlap=1 if data_part==2 & overlap_main==1
    quietly replace burden_overlap=0 if data_part!=2 | overlap_main!=1
    capture drop burden_all
    quietly gen byte burden_all=1 if data_part==2
    quietly replace burden_all=0 if data_part!=2
    capture drop child_overlap
    quietly gen byte child_overlap=1 if data_part==3 & overlap_main==1
    quietly replace child_overlap=0 if data_part!=3 | overlap_main!=1
    capture drop child_all
    quietly gen byte child_all=1 if data_part==3
    quietly replace child_all=0 if data_part!=3
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
    display "SECTION 14 MACRO ADJUSTED OVERLAP WOMAN MODELS"
    foreach y in union_before18 birth_before18 birth_before20 n_births_15_19 n_births_15_24 n_infant_deaths_fullrisk_15_19 n_under5_deaths_fullrisk_15_19 {
        display "MACRO_ADJUSTED_OVERLAP OUTCOME `y'"
        capture noisily regress `y' exposed_before15 gdp_pc_constant_survey health_exp_gdp_survey urban_pop_pct_survey under5_mortality_rate_survey rule_of_law_survey ucdp_total_deaths_survey dtp3_coverage_survey i.survey_id i.cohort_id [pweight=weight] if data_part==1 & overlap_main==1 & `y'<. & exposed_before15<. & weight>0, vce(cluster country_id)
        if !_rc showb exposed_before15
    }
    display "SECTION 15 CHILD LEVEL WEIGHTED CELLS IN OVERLAP SAMPLE"
    capture noisily summarize child_cell_n child_cell_w_all maternal_age_at_birth neonatal_death_any infant_death_any under5_death_any neonatal_death_fullrisk infant_death_fullrisk under5_death_fullrisk if data_part==3 & overlap_main==1
    runchild neonatal_death_any w_neonatal_death_any child_overlap CHILD_NEONATAL_ANY_OVERLAP
    runchild infant_death_any w_infant_death_any child_overlap CHILD_INFANT_ANY_OVERLAP
    runchild under5_death_any w_under5_death_any child_overlap CHILD_UNDER5_ANY_OVERLAP
    runchild neonatal_death_fullrisk w_neonatal_death_fullrisk child_overlap CHILD_NEONATAL_FULLRISK_OVERLAP
    runchild infant_death_fullrisk w_infant_death_fullrisk child_overlap CHILD_INFANT_FULLRISK_OVERLAP
    runchild under5_death_fullrisk w_under5_death_fullrisk child_overlap CHILD_UNDER5_FULLRISK_OVERLAP
    foreach y in births_15_24_w_per1000 infant_deaths_fr_15_24_w_per1000 under5_deaths_fr_15_24_w_per1000 {
        runburden `y' burden_overlap BURDEN_OVERLAP_`y'
        runburdenevent `y' burden_overlap BURDEN_OVERLAP_`y'
    }
    display "CORE_STANDARD_SENSITIVITY_COMPLETE"
    log close
}

if inlist("`task'","all","cutoffs") {
    log using "`output'/cutoffs_pretrends.log", text replace
    use "`input'", clear
    which boottest
    assert inrange(data_part,1,6)
    local outcomes "union_before18 birth_before18 birth_before20 births_15_19_w_per1000 infant_deaths_fr_15_19_w_per1000 under5_deaths_fr_15_19_w_per1000"
    local samples  "analytic_union analytic_birth18 analytic_birth20 burden burden burden"
    local scales   "proportion proportion proportion per_1000_women per_1000_women per_1000_women"


    local cutcsv "`output'/cutoffs.csv"
    local precsv "`output'/proximal_pretrends.csv"
    tempfile cutoff_results
    tempname CUT
    postfile `CUT' int cutoff str45 outcome str22 scale ///
    double N clusters support_n zero_n one_n b se ci_lo ci_hi p_cluster ///
    wcb_p wcb_ci_lo wcb_ci_hi wcb_reps ///
    str244 all_countries str244 support_countries str244 zero_countries str244 one_countries ///
    using "`cutoff_results'", replace

    forvalues k=12/17 {
        capture drop exp_cut
        quietly generate byte exp_cut=(candidate_reform_year<=female_birth_cohort+`k') ///
        if inlist(data_part,1,2) & !missing(candidate_reform_year,female_birth_cohort)

        if `k'==15 {
            quietly count if inlist(data_part,1,2) & exp_cut!=exposed_before15 ///
            & !missing(exp_cut,exposed_before15)
            display "CUTOFF15_EXPOSURE_MISMATCH_N=" r(N)
            assert r(N)==0
        }

        forvalues m=1/6 {
            local y : word `m' of `outcomes'
            local s : word `m' of `samples'
            local scale : word `m' of `scales'

            if `m'<=3 {
                quietly regress `y' exp_cut i.survey_id i.cohort_id [pweight=weight] ///
                if data_part==1 & `s'==1 & `y'<. & exp_cut<. & weight>0, ///
                vce(cluster country_id)
            }
            else {
                quietly regress `y' exp_cut i.survey_id i.cohort_id [aweight=weighted_women_n] ///
                if data_part==2 & `y'<. & exp_cut<. & weighted_women_n>0, ///
                vce(cluster country_id)
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
            capture noisily boottest exp_cut, cluster(country_id) reps(5000) ///
            seed(`bootseed') nograph
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

            display "CUTOFF_RESULT cutoff=`k' outcome=`y' N=" __N ///
            " clusters=`clusters' support_n=`support_n' b=" %14.9f __b ///
            " se=" %14.9f __se " p_cluster=" %10.7f __p ///
            " wcb_p=" %10.7f __wp

            post `CUT' (`k') ("`y'") ("`scale'") (__N) (`clusters') ///
            (`support_n') (`zero_n') (`one_n') (__b) (__se) (__lo) (__hi) ///
            (__p) (__wp) (__wlo) (__whi) (__wreps) ///
            ("`alliso'") ("`supportiso'") ("`zeroiso'") ("`oneiso'")
            drop `esample'
        }
    }
    postclose `CUT'

    preserve
    use "`cutoff_results'", clear
    sort cutoff outcome
    format %21.15g N clusters support_n zero_n one_n b se ci_lo ci_hi p_cluster ///
    wcb_p wcb_ci_lo wcb_ci_hi wcb_reps
    export delimited using "`cutcsv'", replace datafmt
    restore
    capture drop event15_cat
    quietly generate int event15_cat=event15_bin+11 if event15_bin<.

    tempfile proximal_results
    tempname PRE
    postfile `PRE' str10 window str45 outcome str22 scale ///
    double N clusters prepost_n lead_m3_n lead_m2_n ///
    b_m3 se_m3 ci_lo_m3 ci_hi_m3 p_m3 wcb_p_m3 wcb_ci_lo_m3 wcb_ci_hi_m3 ///
    b_m2 se_m2 ci_lo_m2 ci_hi_m2 p_m2 wcb_p_m2 wcb_ci_lo_m2 wcb_ci_hi_m2 ///
    joint_F joint_p joint_wcb_p joint_wcb_reps ///
    str244 all_countries str244 prepost_countries str244 lead_m3_countries str244 lead_m2_countries ///
    using "`proximal_results'", replace

    foreach window in full narrow {
        forvalues m=1/6 {
            local y : word `m' of `outcomes'
            local s : word `m' of `samples'
            local scale : word `m' of `scales'
            local winif "event15_cat<."
            if "`window'"=="narrow" local winif "inrange(event15_bin,-6,8) & event15_cat<."

            if `m'<=3 {
                quietly regress `y' ib10.event15_cat i.survey_id i.cohort_id [pweight=weight] ///
                if data_part==1 & `s'==1 & `winif' & `y'<. & weight>0, ///
                vce(cluster country_id)
            }
            else {
                quietly regress `y' ib10.event15_cat i.survey_id i.cohort_id [aweight=weighted_women_n] ///
                if data_part==2 & `winif' & `y'<. & weighted_women_n>0, ///
                vce(cluster country_id)
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

            capture noisily boottest 8.event15_cat, cluster(country_id) reps(5000) ///
            seed(`seed3') nograph
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

            capture noisily boottest 9.event15_cat, cluster(country_id) reps(5000) ///
            seed(`seed2') nograph
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

            capture noisily boottest 8.event15_cat 9.event15_cat, ///
            cluster(country_id) reps(5000) seed(`seedj') nograph
            if _rc {
                scalar __jwp=.
                scalar __jwreps=.
            }
            else {
                scalar __jwp=r(p)
                scalar __jwreps=r(reps)
            }

            display "PROXIMAL_RESULT window=`window' outcome=`y' N=" __N ///
            " clusters=`clusters' b_m3=" %14.9f __b3 " p_m3=" %10.7f __p3 ///
            " wcb_p_m3=" %10.7f __wp3 " b_m2=" %14.9f __b2 ///
            " p_m2=" %10.7f __p2 " wcb_p_m2=" %10.7f __wp2 ///
            " joint_p=" %10.7f __jp " joint_wcb_p=" %10.7f __jwp

            post `PRE' ("`window'") ("`y'") ("`scale'") (__N) (`clusters') ///
            (`prepost_n') (`lead3_n') (`lead2_n') ///
            (__b3) (__s3) (__l3) (__h3) (__p3) (__wp3) (__wl3) (__wh3) ///
            (__b2) (__s2) (__l2) (__h2) (__p2) (__wp2) (__wl2) (__wh2) ///
            (__jf) (__jp) (__jwp) (__jwreps) ///
            ("`alliso'") ("`prepostiso'") ("`lead3iso'") ("`lead2iso'")
            drop `esample'
        }
    }
    postclose `PRE'

    preserve
    use "`proximal_results'", clear
    sort window outcome
    format %21.15g N clusters prepost_n lead_m3_n lead_m2_n ///
    b_m3 se_m3 ci_lo_m3 ci_hi_m3 p_m3 wcb_p_m3 wcb_ci_lo_m3 wcb_ci_hi_m3 ///
    b_m2 se_m2 ci_lo_m2 ci_hi_m2 p_m2 wcb_p_m2 wcb_ci_lo_m2 wcb_ci_hi_m2 ///
    joint_F joint_p joint_wcb_p joint_wcb_reps
    export delimited using "`precsv'", replace datafmt
    restore
    display "CORE_CUTOFFS_PRETRENDS_COMPLETE"
    log close
}

if inlist("`task'","all","countries") {
    local wm "`root'/data/data_women.dta"
    local bur "`root'/data/data_cohort.dta"
    confirm file "`wm'"
    confirm file "`bur'"
    local log "`output'/country_inclusion.log"
    local csv "`output'/country_inclusion.csv"
    log using "`log'", text replace
    which boottest
    tempfile women_extra burden_extra hybrid_results

    use "`wm'", clear
    keep if inlist(country_id,3,4)
    keep country_name iso3 country_id survey_id cohort_id survey_year ///
    female_birth_year female_birth_cohort age weight candidate_reform_year ///
    has_candidate_reform exposed_before15 analytic_union analytic_birth18 ///
    analytic_birth20 union_before18 birth_before18 birth_before20

    assert _N==46657
    quietly count if country_id==3
    assert r(N)==6490
    quietly count if country_id==4
    assert r(N)==40167

    assert missing(candidate_reform_year) if country_id==3
    assert has_candidate_reform==0 if country_id==3
    assert exposed_before15==0 if country_id==3
    replace candidate_reform_year=. if country_id==3
    replace has_candidate_reform=0 if country_id==3
    replace exposed_before15=0 if country_id==3
    replace analytic_union=inrange(age,20,29) & union_before18<. & weight>0 ///
    if country_id==3
    replace analytic_birth18=inrange(age,20,29) & birth_before18<. & weight>0 ///
    if country_id==3
    replace analytic_birth20=inrange(age,20,29) & birth_before20<. & weight>0 ///
    if country_id==3

    generate byte data_part=1
    generate str32 source_marker="data_women_COG"
    replace source_marker="data_women_EGY" if country_id==4
    save `women_extra', replace

    use "`bur'", clear
    keep if inlist(country_id,3,4)
    keep country_name iso3 country_id survey_id cohort_id survey_year ///
    female_birth_cohort candidate_reform_year has_candidate_reform ///
    exposed_before15 weighted_women_n births_15_19_w_per1000 ///
    infant_deaths_fr_15_19_w_per1000 under5_deaths_fr_15_19_w_per1000

    assert _N==100
    quietly count if country_id==3
    assert r(N)==20
    quietly count if country_id==4
    assert r(N)==80
    assert missing(candidate_reform_year) if country_id==3
    assert has_candidate_reform==0 if country_id==3
    assert exposed_before15==0 if country_id==3
    replace candidate_reform_year=. if country_id==3
    replace has_candidate_reform=0 if country_id==3
    replace exposed_before15=0 if country_id==3

    generate byte data_part=2
    generate str32 source_marker="data_cohort_COG"
    replace source_marker="data_cohort_EGY" if country_id==4
    save `burden_extra', replace


    use "`input'", clear
    keep if inlist(data_part,1,2)
    keep data_part country_name iso3 country_id survey_id cohort_id survey_year ///
    female_birth_year female_birth_cohort age weight candidate_reform_year ///
    has_candidate_reform exposed_before15 analytic_union analytic_birth18 ///
    analytic_birth20 union_before18 birth_before18 birth_before20 ///
    weighted_women_n births_15_19_w_per1000 ///
    infant_deaths_fr_15_19_w_per1000 under5_deaths_fr_15_19_w_per1000
    assert !inlist(country_id,3,4)
    generate str32 source_marker="data_base16"
    append using `women_extra' `burden_extra'

    quietly count if data_part==1 & source_marker=="data_base16"
    assert r(N)==314874
    quietly count if data_part==2 & source_marker=="data_base16"
    assert r(N)==770
    quietly count if data_part==1 & source_marker=="data_women_COG"
    assert r(N)==6490
    quietly count if data_part==1 & source_marker=="data_women_EGY"
    assert r(N)==40167
    quietly count if data_part==2 & source_marker=="data_cohort_COG"
    assert r(N)==20
    quietly count if data_part==2 & source_marker=="data_cohort_EGY"
    assert r(N)==80

    display "HYBRID_SOURCE_AUDIT"
    tab data_part source_marker, missing


    local samples "base16 plus_COG plus_EGY all18"
    local outcomes "union_before18 birth_before18 birth_before20 births_15_19_w_per1000 infant_deaths_fr_15_19_w_per1000 under5_deaths_fr_15_19_w_per1000"
    local flags "analytic_union analytic_birth18 analytic_birth20 burden burden burden"
    local domains "individual individual individual burden burden burden"
    local scales "proportion proportion proportion per_1000_women per_1000_women per_1000_women"

    local expected_N "311260 314592 314592 770 693 462"
    local expected_b "-0.064164116 -0.052670793 -0.079727733 -50.081071091 -6.978656526 -10.723971581"

    tempname RES
    postfile `RES' str16 construction str10 sample str12 added_countries ///
    str45 outcome str12 domain str22 scale double N double clusters ///
    double support_n double zero_only_n double one_only_n double b ///
    double se double ci_lo double ci_hi double p_cluster double wcb_p ///
    double wcb_ci_lo double wcb_ci_hi double wcb_reps ///
    str244 all_countries str244 support_countries ///
    str244 zero_only_countries str244 one_only_countries ///
    using `hybrid_results', replace

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
                quietly regress `y' exposed_before15 i.survey_id i.cohort_id ///
                [pweight=weight] if data_part==1 & `flag'==1 & ///
                (`keepcond') & `y'<. & exposed_before15<. & weight>0, ///
                vce(cluster country_id)
            }
            else {
                quietly regress `y' exposed_before15 i.survey_id i.cohort_id ///
                [aweight=weighted_women_n] if data_part==2 & ///
                (`keepcond') & `y'<. & exposed_before15<. & ///
                weighted_women_n>0, vce(cluster country_id)
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
                quietly levelsof iso3 if `esample'==1 & country_id==`c', ///
                local(thisiso) clean
                local alliso "`alliso' `thisiso'"
                quietly count if `esample'==1 & country_id==`c' & ///
                exposed_before15==0
                local n0=r(N)
                quietly count if `esample'==1 & country_id==`c' & ///
                exposed_before15==1
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
                display "BASE16_REPRODUCTION_PASS outcome=`y' N=" %12.0f __N ///
                " b=" %14.9f __b
            }

            local bootseed=20260708
            capture noisily boottest exposed_before15, cluster(country_id) ///
            reps(5000) seed(`bootseed') nograph
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

            display "HYBRID_RESULT sample=`sample' added=`added' outcome=`y'" ///
            " N=" %12.0f __N " clusters=`clusters' support_n=`support_n'" ///
            " b=" %14.9f __b " se=" %14.9f __se ///
            " p_cluster=" %10.7f __p " wcb_p=" %10.7f __wp

            post `RES' ("data_expanded") ("`sample'") ("`added'") ("`y'") ///
            ("`domain'") ("`scale'") (__N) (`clusters') (`support_n') ///
            (`zero_n') (`one_n') (__b) (__se) (__lo) (__hi) (__p) ///
            (__wp) (__wlo) (__whi) (__wreps) ("`alliso'") ///
            ("`supportiso'") ("`zeroiso'") ("`oneiso'")
            drop `esample'
        }
    }
    postclose `RES'

    preserve
    use `hybrid_results', clear
    assert _N==24
    assert !missing(N,clusters,b,se,ci_lo,ci_hi,p_cluster,wcb_p,wcb_reps)
    assert wcb_reps==5000
    assert inrange(wcb_p,0,1)
    isid sample outcome
    display "ALL_24_MODELS_AND_5000_REPLICATIONS_PASS"
    sort domain outcome sample
    format %21.15g N clusters support_n zero_only_n one_only_n b se ci_lo ///
    ci_hi p_cluster wcb_p wcb_ci_lo wcb_ci_hi wcb_reps
    export delimited using "`csv'", replace datafmt
    list sample added_countries outcome N clusters support_n b se ci_lo ci_hi ///
    p_cluster wcb_p wcb_reps, noobs abbreviate(32)
    restore

    display "PERSISTENT_OUTPUT_LOG=`log'"
    display "PERSISTENT_OUTPUT_CSV=`csv'"
    display "NO_COMBINED_MICRODATA_SAVED=YES"
    display "CORE_COUNTRY_INCLUSION_COMPLETE"
    log close
}
