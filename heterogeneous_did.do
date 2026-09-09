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
if "`task'"=="" local task "cs_notyet"
if !inlist("`task'","cs_notyet","cs_never","ppml") {
    exit 198
}
local root "`c(pwd)'"
local input "`root'/data/data.dta"
confirm file "`input'"
capture mkdir "`root'/results"
local output "`root'/results/`task'"
capture mkdir "`output'"

if inlist("`task'","cs_notyet","cs_never") {
    local control "notyet"
    if "`task'"=="cs_never" local control "never"
    local reps=5000
    confirm file "`input'"
    use "`input'", clear
    capture log close
    log using "`output'/cs_`control'.log", text replace
    version 19.5
    set more off
    foreach v in data_part country_id country_name female_birth_year candidate_reform_year exposed_before15 weight analytic_union analytic_birth18 analytic_birth20 union_before18 birth_before18 birth_before20 {
        capture confirm variable `v'
        if _rc {
            display as error "Required variable not found: `v'"
            log close
            exit 111
        }
    }
    preserve
    keep if data_part==1
    keep if !missing(country_id,female_birth_year,candidate_reform_year,exposed_before15,weight)
    capture drop csa_time csa_any_treated csa_cohort csa_treat csa_event
    gen int csa_time=female_birth_year+15
    bysort country_id: egen byte csa_any_treated=max(exposed_before15)
    gen int csa_cohort=candidate_reform_year if csa_any_treated==1
    replace csa_cohort=0 if csa_any_treated==0
    gen byte csa_treat=csa_time>=csa_cohort if csa_cohort>0 & csa_time<.
    replace csa_treat=0 if csa_cohort==0 & csa_time<.
    gen int csa_event=csa_time-csa_cohort if csa_cohort>0 & csa_time<.
    assert csa_treat==exposed_before15 if exposed_before15<.
    display "FORMAL_CSA_SETUP REPS=`reps' CONTROL=`control' BASETIME=COMMON"
    display "COUNTRY_TREATMENT_STRUCTURE"
    tab country_name csa_treat, row
    tab csa_cohort country_name, missing
    local outcomes "union_before18 birth_before18 birth_before20"
    local analytics "analytic_union analytic_birth18 analytic_birth20"
    local labels "UNION18 BIRTH18 BIRTH20"
    matrix FORMAL_RESULTS=J(3,7,.)
    matrix rownames FORMAL_RESULTS=UNION18 BIRTH18 BIRTH20
    matrix colnames FORMAL_RESULTS=N ATT_PP BOOT_SE_PP P_VALUE SCI95_LL_PP SCI95_UL_PP REPS
    local i=0
    foreach y of local outcomes {
        local ++i
        local a : word `i' of `analytics'
        local lab : word `i' of `labels'
        local seed=20260740+`i'
        quietly count if `a'==1 & `y'<. & weight>0
        local n=r(N)
        display "CSA_MODEL_START OUTCOME=`lab' N=`n' REPS=`reps' SEED=`seed'"
        capture noisily hdidregress aipw (`y') (csa_treat) [pweight=weight] if `a'==1 & `y'<. & weight>0, group(country_id) time(csa_time) usercohort(csa_cohort) controlgroup(`control') basetime(common) vce(cluster country_id) nolog
        local modelrc=_rc
        if `modelrc' {
            display as error "CSA_MODEL_FAILED OUTCOME=`lab' RC=`modelrc'"
        }
        else {
            capture noisily estat ptrends
            local ptrc=_rc
            if `ptrc' display as error "CSA_PTRENDS_FAILED OUTCOME=`lab' RC=`ptrc'"
            capture noisily estat aggregation, overall sci(reps(`reps') rseed(`seed'))
            local aggrc=_rc
            if `aggrc' {
                display as error "CSA_AGGREGATION_FAILED OUTCOME=`lab' RC=`aggrc'"
            }
            else {
                matrix CSA_TABLE=r(table)
                matrix list CSA_TABLE, format(%12.6f)
                return list
                matrix FORMAL_RESULTS[`i',1]=`n'
                matrix FORMAL_RESULTS[`i',2]=100*CSA_TABLE[1,1]
                matrix FORMAL_RESULTS[`i',3]=100*CSA_TABLE[2,1]
                matrix FORMAL_RESULTS[`i',4]=CSA_TABLE[4,1]
                matrix FORMAL_RESULTS[`i',5]=100*CSA_TABLE[5,1]
                matrix FORMAL_RESULTS[`i',6]=100*CSA_TABLE[6,1]
                matrix FORMAL_RESULTS[`i',7]=r(reps)
                display "CSA_FORMAL_RESULT OUTCOME=`lab' N=`n' ATT_PP=" %12.6f FORMAL_RESULTS[`i',2] " BOOT_SE_PP=" %12.6f FORMAL_RESULTS[`i',3] " P=" %9.6f FORMAL_RESULTS[`i',4] " SCI95_LL_PP=" %12.6f FORMAL_RESULTS[`i',5] " SCI95_UL_PP=" %12.6f FORMAL_RESULTS[`i',6] " REPS=" %9.0f FORMAL_RESULTS[`i',7]
            }
        }
        display "CSA_MODEL_END OUTCOME=`lab'"
    }
    display "CSA_FORMAL_SUMMARY_TABLE"
    matrix list FORMAL_RESULTS, format(%12.6f)
    restore
    clear
    svmat double FORMAL_RESULTS, names(col)
    gen str12 outcome="UNION18" in 1
    replace outcome="BIRTH18" in 2
    replace outcome="BIRTH20" in 3
    assert !missing(N,ATT_PP,BOOT_SE_PP,P_VALUE,SCI95_LL_PP,SCI95_UL_PP,REPS)
    assert REPS==`reps'
    export delimited using "`output'/cs_`control'.csv", replace
    display "CS_REPLICATION_COMPLETE control=`control'"
    log close
}

if "`task'"=="ppml" {
    local oldprocessors=c(processors)
    capture set processors 1
    confirm file "`input'"
    use "`input'", clear
    cd "`output'"
    capture log close
    local bootreps=5000
    local maxmult=5
    if "`bootreps'"=="" local bootreps=5000
    if "`maxmult'"=="" local maxmult=5
    capture confirm integer number `bootreps'
    if _rc | `bootreps'<1 {
        display as error "Optional bootstrap replication argument must be a positive integer"
        exit 198
    }
    capture confirm integer number `maxmult'
    if _rc | `maxmult'<1 {
        display as error "Optional maximum-attempt multiplier must be a positive integer"
        exit 198
    }
    if `bootreps'==5000 {
        log using "ppml.log", text replace
    }
    else {
        log using "ppml_smoketest_`bootreps'_result.txt", text replace
    }
    set more off
    capture which ftools
    if _rc {
        exit 499
    }
    capture which reghdfe
    if _rc {
        exit 499
    }
    capture which ppmlhdfe
    if _rc {
        exit 499
    }
    capture program drop burden_ppml_att
    program define burden_ppml_att, rclass
        version 19.5
        syntax varname, HET(string)
        if !inlist("`het'","event","cohort","time") {
            display as error "HET() must be event, cohort, or time"
            exit 198
        }
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
            return scalar ok=0
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
        return scalar ok=1
        return scalar rc=0
    end
    capture program drop burden_ppml_draw
    program define burden_ppml_draw, rclass
        version 19.5
        syntax varname, HET(string)
        preserve
        capture quietly bsample, cluster(country_id) idcluster(__boot_country_id)
        local samplerc=_rc
        if `samplerc' {
            restore
            return scalar att=.
            return scalar ok=0
            return scalar rc=`samplerc'
            exit
        }
        capture quietly burden_ppml_att `varlist', het(`het')
        local drawrc=_rc
        if `drawrc' {
            restore
            return scalar att=.
            return scalar ok=0
            return scalar rc=`drawrc'
            exit
        }
        local drawatt=r(att)
        local drawok=r(ok)
        local fitrc=r(rc)
        restore
        return scalar att=`drawatt'
        return scalar ok=`drawok'
        return scalar rc=`fitrc'
    end
    capture confirm variable __boot_country_id
    if !_rc {
        display as error "Variable __boot_country_id already exists; reload data/data.dta before running"
        log close
        exit 110
    }
    local yi=0
    foreach y in births_15_19_w_per1000 infant_deaths_fr_15_19_w_per1000 under5_deaths_fr_15_19_w_per1000 {
        local ++yi
        preserve
        quietly keep if data_part==2 & `y'<. & weighted_women_n>0 & !missing(country_id,survey_id,female_birth_cohort,candidate_reform_year)
        local sampleN=_N
        quietly levelsof country_id, local(clusterlevels)
        local nclusters: word count `clusterlevels'
        local hi=0
        foreach h in event cohort time {
            local ++hi
            local bseed = 20260724 + 10*`yi' + `hi'
            noisily display "MODEL_START OUTCOME=`y' HETEROGENEITY=`h' REQUESTED_REPS=`bootreps' SEED=`bseed'"
            capture noisily burden_ppml_att `y', het(`h')
            local pointrc=_rc
            if `pointrc' {
                noisily display as error "POINT_MODEL_FAILED OUTCOME=`y' HETEROGENEITY=`h' RC=`pointrc'"
            }
            else if missing(r(att)) {
                noisily display as error "POINT_MODEL_NOT_ESTIMABLE OUTCOME=`y' HETEROGENEITY=`h' FIT_RC=" r(rc)
            }
            else {
                local pointatt = r(att)
                local pointn = r(n)
                noisily return list
                local maxattempts = `maxmult' * `bootreps'
                local valid=0
                local failed=0
                local commandfailed=0
                local notestimable=0
                local attempts=0
                tempname posth observed semat
                tempfile repfile
                postfile `posth' double att using "`repfile'", replace
                set seed `bseed'
                while `valid'<`bootreps' & `attempts'<`maxattempts' {
                    capture quietly burden_ppml_draw `y', het(`h')
                    local drawrc=_rc
                    local drawok=0
                    local drawatt=.
                    if !`drawrc' {
                        local drawok=r(ok)
                        local drawatt=r(att)
                    }
                    local ++attempts
                    if !`drawrc' & `drawok'==1 & `drawatt'<. {
                        post `posth' (`drawatt')
                        local ++valid
                    }
                    else {
                        local ++failed
                        if `drawrc' {
                            local ++commandfailed
                        }
                        else {
                            local ++notestimable
                        }
                    }
                    if mod(`attempts',250)==0 noisily display "BOOTSTRAP_PROGRESS OUTCOME=`y' HETEROGENEITY=`h' VALID=`valid' FAILED=`failed' ATTEMPTED=`attempts' TARGET=`bootreps'"
                }
                postclose `posth'
                local validpct = 100*`valid'/`attempts'
                noisily display "BOOTSTRAP_COUNTS OUTCOME=`y' HETEROGENEITY=`h' VALID=`valid' FAILED=`failed' ATTEMPTED=`attempts' TARGET_VALID=`bootreps' CLUSTERS=`nclusters' VALID_PERCENT=" %6.2f `validpct' " COMMAND_FAILED=`commandfailed' NOT_ESTIMABLE=`notestimable'"
                if `valid'<`bootreps' noisily display as error "BOOTSTRAP_TARGET_NOT_REACHED OUTCOME=`y' HETEROGENEITY=`h' VALID=`valid' TARGET=`bootreps' MAX_ATTEMPTS=`maxattempts'"
                if `validpct'<90 noisily display as error "BOOTSTRAP_LOW_COMPLETION_WARNING OUTCOME=`y' HETEROGENEITY=`h' VALID_PERCENT=" %6.2f `validpct'
                if `valid'<1000 noisily display as error "PERCENTILE_CI_UNSTABLE_WARNING OUTCOME=`y' HETEROGENEITY=`h' VALID_REPS=`valid'"
                if `valid'>=2 {
                    matrix `observed'=(`pointatt')
                    matrix colnames `observed'=att
                    capture noisily bstat att using "`repfile'", stat(`observed') n(`sampleN')
                    local bstatrc=_rc
                    if `bstatrc' {
                        noisily display as error "BSTAT_FAILED OUTCOME=`y' HETEROGENEITY=`h' RC=`bstatrc'"
                    }
                    else {
                        matrix `semat'=e(se)
                        local bootse = `semat'[1,1]
                        local z = `pointatt'/`bootse'
                        local p = 2*normal(-abs(`z'))
                        noisily display "BOOTSTRAP_NORMAL_TEST OUTCOME=`y' HETEROGENEITY=`h' ATT=" %12.6f `pointatt' " SE=" %12.6f `bootse' " Z=" %9.4f `z' " P=" %9.6f `p'
                        capture noisily estat bootstrap, normal percentile
                        local estatrc=_rc
                        if `estatrc' noisily display as error "ESTAT_BOOTSTRAP_FAILED OUTCOME=`y' HETEROGENEITY=`h' RC=`estatrc'"
                    }
                }
                else {
                    noisily display as error "BOOTSTRAP_INSUFFICIENT_VALID_REPS OUTCOME=`y' HETEROGENEITY=`h' VALID=`valid'"
                }
            }
            noisily display "MODEL_END OUTCOME=`y' HETEROGENEITY=`h'"
        }
        restore
    }
    log close

    display "PPML_REPLICATION_COMPLETE"

    cd "`root'"
    capture set processors `oldprocessors'
}
