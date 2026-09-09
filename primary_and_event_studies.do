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
if !inlist("`task'","all","main","wcb") {
    exit 198
}
local root "`c(pwd)'"
local input "`root'/data/data.dta"
confirm file "`input'"
capture mkdir "`root'/results"
local output "`root'/results/module2"
capture mkdir "`output'"

if inlist("`task'","all","main") {
    use "`input'", clear
    log using "`output'/primary.log", text replace
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
    display "SECTION 6 RECOMMENDED OVERLAP MAIN MODELS"
    runmain union_before18 exposed_before15 overlap_analytic_union OVERLAP_UNION_EXPOSED15
    runmain birth_before18 exposed_before15 overlap_analytic_birth18 OVERLAP_BIRTH18_EXPOSED15
    runmain birth_before20 exposed_before15 overlap_analytic_birth20 OVERLAP_BIRTH20_EXPOSED15
    display "SECTION 10 EVENT STUDY RECOMMENDED OVERLAP SAMPLE"
    runevent union_before18 overlap_analytic_union OVERLAP_UNION
    runevent birth_before18 overlap_analytic_birth18 OVERLAP_BIRTH18
    runevent birth_before20 overlap_analytic_birth20 OVERLAP_BIRTH20
    foreach y in births_15_19_w_per1000 infant_deaths_fr_15_19_w_per1000 under5_deaths_fr_15_19_w_per1000 {
        runburden `y' burden_overlap BURDEN_OVERLAP_`y'
        runburdenevent `y' burden_overlap BURDEN_OVERLAP_`y'
    }
    display "CORE_PRIMARY_COMPLETE"
    log close
}

if inlist("`task'","all","wcb") {
    log using "`output'/primary_wcb.log", text replace
    use "`input'", clear
    which boottest
    assert inrange(data_part,1,6)
    local outcomes "union_before18 birth_before18 birth_before20 births_15_19_w_per1000 infant_deaths_fr_15_19_w_per1000 under5_deaths_fr_15_19_w_per1000"
    local samples  "analytic_union analytic_birth18 analytic_birth20 burden burden burden"
    local scales   "proportion proportion proportion per_1000_women per_1000_women per_1000_women"


    local cutcsv "`output'/primary_wcb.csv"
    tempfile cutoff_results
    tempname CUT
    postfile `CUT' int cutoff str45 outcome str22 scale ///
    double N clusters support_n zero_n one_n b se ci_lo ci_hi p_cluster ///
    wcb_p wcb_ci_lo wcb_ci_hi wcb_reps ///
    str244 all_countries str244 support_countries str244 zero_countries str244 one_countries ///
    using "`cutoff_results'", replace

    forvalues k=15/15 {
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
    display "CORE_PRIMARY_WCB_COMPLETE"
    log close
}
