/*==============================================================================
  SEA-PLM 2019 — School-Level ICC Analysis

  Purpose:
    1. Import the SEA-PLM Regional Students file
    2. Document key variables
    3. Estimate unconditional school-level ICC for each country and subject
    4. Estimate conditional school-level ICC for the Philippines using
       school-level covariates merged from the school questionnaire file

  Method:
    - Two-level null mixed model (students nested in schools)
    - ICC estimated separately for each of the 5 plausible values per subject
    - Final ICC and SE derived using Rubin's (1987) combining rules
    - Student weights (WT2019) scaled to sample size within country to avoid
      inflated denominator degrees of freedom

  Requirements: Stata 15+ (uses estat icc with r(icc2))

  Author:  [your name]
  Date:    [date]
==============================================================================*/

clear all
set more off
set linesize 120


/*------------------------------------------------------------------------------
  GLOBALS
------------------------------------------------------------------------------*/

global data_path  `"/Users/douglasjohnson/Documents/Data/SEA-PLM/2019"'
global out_path   `"/Users/douglasjohnson/Documents/code/sea_plm_analysis/Analysis"'


/*==============================================================================
  SECTION 1 — IMPORT STUDENT DATA
==============================================================================*/

di as txt _n "Loading student data..."
import spss using "$data_path/SEA-PLM_Regional_Students 2.0.SAV", clear

di as txt "Obs: `=_N'   Vars: `=c(k)'"


/*==============================================================================
  SECTION 2 — KEY VARIABLES

  IDENTIFIERS
    CNT        Country code (string): KHM LAO MMR MYS PHL VNM
    StIDStrt   Stratum ID (2-digit)
    StIDSch    School ID within stratum (3-digit)
    ClassID    Class ID (2-digit)
    StIDStd    Student ID (3-digit)

  PLAUSIBLE VALUES  (5 per subject, labelled PV#_X)
    PV1_M – PV5_M   Mathematics
    PV1_R – PV5_R   Reading
    PV1_W – PV5_W   Writing

  PROFICIENCY BANDS  (one per PV, labelled PL_PV#_X)
    PL_PV1_M – PL_PV5_M
    PL_PV1_R – PL_PV5_R
    PL_PV1_W – PL_PV5_W

  WEIGHTS
    WT2019     Final student weight (use for population estimates)
    rwgt1–rwgt95   Replicate weights (Jackknife 2; use for SE via replication)

  STUDENT BACKGROUND
    Gender     1=Female, 2=Male
    S_AGE      Age
    S_LANG     Language spoken at home
    S_HOMLIT   Books in the home
    S_FISCED   Father's highest education
    S_MISCED   Mother's highest education
    SES        Socioeconomic index (nationally standardised)
    PARED      Highest parental education
    HOMERES    Home resources index

  SCHOOL-LEVEL VARIABLES (from school questionnaire, merged in Section 4)
    C_SCHSIZE      Number of students in school
    SCH_SIZE_G5    Grade 5 enrolment
    PCGIRLS_G5     Proportion of Grade 5 girls
    STRATIO_PERM   Permanent teacher-student ratio
    STRATIO_TEMP   Temporary teacher-student ratio
    C_GENDER       Principal gender
    C_ISCED        Principal highest education
    C_AGE          Principal approximate age
    RESOU          Resources in local area index
    HINDER         Issues hindering school capacity
    COMATPA        Community attitudes — parents
    COMATST        Community attitudes — students
    COMATTC        Community attitudes — teachers
    STUDISS        Frequency of student issues
==============================================================================*/
rename cnt CNT

* Confirm key variables exist
foreach v of varlist CNT StIDStrt StIDSch WT2019 ///
    PV1_M PV1_R PV1_W rwgt1 {
    capture confirm variable `v'
    if _rc {
        di as error "Variable `v' not found — check variable names in the SAV file"
        exit 111
    }
}
di as txt "All key variables confirmed."


/*==============================================================================
  SECTION 3 — DATA PREPARATION
==============================================================================*/

*--- Unique school identifier (country + stratum + school) -------------------*
* Concatenate as string to guarantee uniqueness across countries
* (schools in different countries may share numeric IDs)

* Convert numeric IDs to string if needed
capture confirm string variable StIDStrt
if _rc {
    tostring StIDStrt, gen(StIDStrt_s) format("%02.0f")
}
else {
    gen StIDStrt_s = StIDStrt
}

capture confirm string variable StIDSch
if _rc {
    tostring StIDSch, gen(StIDSch_s) format("%03.0f")
}
else {
    gen StIDSch_s = StIDSch
}

* Ensure CNT is string
capture confirm string variable CNT
if _rc {
    decode CNT, gen(CNT_s)
}
else {
    gen CNT_s = CNT
}

gen schoolid_str = CNT_s + StIDStrt_s + StIDSch_s
encode schoolid_str, gen(schoolid)
label variable schoolid "Unique school identifier (CNT+stratum+school)"

* Drop helper variables
drop StIDStrt_s StIDSch_s CNT_s schoolid_str


*--- Scaled weights (sum to sample size within each country) -----------------*
* Using pweights in mixed inflates df when sum(weights) >> N.
* Scaling so mean weight = 1 within country preserves the relative weighting
* while keeping df calibrated to the actual sample.

bysort CNT: egen wt_mean = mean(WT2019)
gen wt_scaled = WT2019 / wt_mean
drop wt_mean
label variable wt_scaled "WT2019 scaled to mean=1 within country"


/*==============================================================================
  SECTION 4 — UNCONDITIONAL SCHOOL-LEVEL ICC BY COUNTRY

  Model: Y_ij = μ + u_j + e_ij
    u_j ~ N(0, σ²_u)   school-level random effect
    e_ij ~ N(0, σ²_e)  student-level residual

  ICC = σ²_u / (σ²_u + σ²_e)

  Combining rules across 5 PVs (Rubin 1987):
    ICC_bar  = mean(ICC_k)
    V_W      = mean(SE²_k)                        [within-imputation variance]
    V_B      = sum((ICC_k - ICC_bar)²) / (5-1)    [between-imputation variance]
    V_total  = V_W + (1 + 1/5) * V_B
    SE_comb  = sqrt(V_total)
==============================================================================*/

di as txt _n(2) "========================================================"
di as txt        "  UNCONDITIONAL ICC BY COUNTRY AND SUBJECT"
di as txt        "========================================================"

* Countries and subjects to loop over
local countries KHM LAO MMR MYS PHL VNM
local subjects  M R W
local subj_lbl_M "Mathematics"
local subj_lbl_R "Reading"
local subj_lbl_W "Writing"

* Results matrix: rows = countries x subjects, cols = ICC, SE, N_students, N_schools
tempname icc_results
local nrow = 18    // 6 countries × 3 subjects
matrix `icc_results' = J(`nrow', 5, .)
matrix colnames `icc_results' = ICC SE N_students N_schools PVs_used
local row = 0

foreach cnt of local countries {
    foreach subj of local subjects {
        local ++row
        di as txt _n "  Country: `cnt'   Subject: `subj_lbl_`subj''"

        * Temporary storage for PV-specific ICC and SE
        matrix pv_icc = J(5,2,.)   // col1=ICC, col2=SE

        forvalues pv = 1/5 {
            quietly mixed PV`pv'_`subj' [pw=wt_scaled] ///
                if CNT == "`cnt'" || schoolid:, mle

            * Extract variance components
            scalar var_u = exp(2 * [lns1_1_1]_cons)
            scalar var_e = exp(2 * [lnsig_e]_cons)
            scalar icc_k = var_u / (var_u + var_e)

            * SE from estat icc (delta method)
            quietly estat icc
            scalar se_k = r(se2)

            matrix pv_icc[`pv', 1] = icc_k
            matrix pv_icc[`pv', 2] = se_k
        }

        * Rubin's combining rules
        scalar icc_bar = 0
        forvalues pv = 1/5 {
            scalar icc_bar = icc_bar + pv_icc[`pv', 1]
        }
        scalar icc_bar = icc_bar / 5

        scalar V_W = 0
        scalar V_B = 0
        forvalues pv = 1/5 {
            scalar V_W = V_W + (pv_icc[`pv', 2])^2
            scalar V_B = V_B + (pv_icc[`pv', 1] - icc_bar)^2
        }
        scalar V_W = V_W / 5
        scalar V_B = V_B / 4
        scalar se_comb = sqrt(V_W + 1.2 * V_B)

        * Sample sizes (from last mixed estimation)
        scalar n_stu    = e(N)
        scalar n_school = e(N_clust)

        di as txt "    ICC = " %6.4f icc_bar ///
                  "   SE = " %6.4f se_comb ///
                  "   95% CI: [" %6.4f (icc_bar - 1.96*se_comb) ///
                  ", " %6.4f (icc_bar + 1.96*se_comb) "]" ///
                  "   N_students=" n_stu "  N_schools=" n_school

        matrix `icc_results'[`row', 1] = icc_bar
        matrix `icc_results'[`row', 2] = se_comb
        matrix `icc_results'[`row', 3] = n_stu
        matrix `icc_results'[`row', 4] = n_school
        matrix `icc_results'[`row', 5] = 5   // all 5 PVs used
    }
}

* Display full results table
di as txt _n(2) "--- Summary: Unconditional School-Level ICC ---"
di as txt "Row order: (KHM_M, KHM_R, KHM_W, LAO_M, ..., VNM_W)"
matrix list `icc_results', format(%8.4f)


/*==============================================================================
  SECTION 5 — MERGE SCHOOL-LEVEL DATA FOR THE PHILIPPINES
==============================================================================*/

di as txt _n(2) "========================================================"
di as txt        "  LOADING AND MERGING SCHOOL DATA (Philippines)"
di as txt        "========================================================"

preserve

    *--- Load school file into tempfile --------------------------------------*
    tempfile school_data

    import spss using "$data_path/SEA-PLM_Regional_Schools 2.0.SAV", clear

    * Standardise country variable name (SPSS import may use lowercase)
    capture confirm variable CNT
    if _rc {
        rename cnt CNT
    }

    * Keep Philippines only
    keep if CNT == "PHL"

    * Keep school identifiers and school-level covariates
    keep CNT StIDStrt StIDSch ///
         C_SCHSIZE SCHSIZE_G5 PCGIRLS_G5 ///
         STRATIO_PERM STRATIO_TEMP ///
         C_GENDER C_ISCED C_AGE ///
         resou hinder comatpa comatst comattc studiss

    * Build same school identifier string
    capture confirm string variable StIDStrt
    if _rc {
        tostring StIDStrt, replace format("%02.0f")
    }
    capture confirm string variable StIDSch
    if _rc {
        tostring StIDSch, replace format("%03.0f")
    }
    capture confirm string variable CNT
    if _rc {
        decode CNT, gen(CNT_temp)
        drop CNT
        rename CNT_temp CNT
    }

    gen schoolid_str = CNT + StIDStrt + StIDSch

    drop StIDStrt StIDSch CNT

    * Confirm no duplicates on school key
    duplicates report schoolid_str

    save `school_data', replace

restore

*--- Merge back onto student data (Philippines only) -------------------------*

* Rebuild the schoolid_str in the student data for merge
capture confirm string variable StIDStrt
if _rc {
    tostring StIDStrt, gen(StIDStrt_s) format("%02.0f")
}
else {
    gen StIDStrt_s = StIDStrt
}

capture confirm string variable StIDSch
if _rc {
    tostring StIDSch, gen(StIDSch_s) format("%03.0f")
}
else {
    gen StIDSch_s = StIDSch
}

capture confirm string variable CNT
if _rc {
    decode CNT, gen(CNT_s)
}
else {
    gen CNT_s = CNT
}

gen schoolid_str = CNT_s + StIDStrt_s + StIDSch_s
drop StIDStrt_s StIDSch_s CNT_s

merge m:1 schoolid_str using `school_data', keep(master match) gen(_merge_school)

di as txt "School merge results for Philippines:"
tab _merge_school if CNT == "PHL"

drop schoolid_str _merge_school

* Standardise continuous school-level covariates (mean 0, SD 1) among PHL students
* to aid interpretation of conditional ICC
foreach v of varlist C_SCHSIZE SCHSIZE_G5 PCGIRLS_G5 ///
                     STRATIO_PERM STRATIO_TEMP C_AGE ///
                     resou hinder comatpa comatst comattc studiss {
    quietly summarize `v' if CNT == "PHL"
    if r(sd) > 0 {
        gen z_`v' = (`v' - r(mean)) / r(sd) if CNT == "PHL"
        label variable z_`v' "`v' standardised (PHL)"
    }
}


/*==============================================================================
  SECTION 6 — CONDITIONAL ICC FOR THE PHILIPPINES

  Two model families per subject:

  Model A — School resources & size
    Covariates: z_C_SCHSIZE  z_STRATIO_PERM  z_STRATIO_TEMP  z_PCGIRLS_G5

  Model B — Principal characteristics & community context
    Covariates: C_GENDER  z_C_AGE  C_ISCED
                z_RESOU  z_HINDER
                z_COMATPA  z_COMATST  z_COMATTC  z_STUDISS

  Conditional ICC = σ²_u|X / (σ²_u|X + σ²_e|X)

  Variance explained at school level (VES):
    VES = 1 - (σ²_u|X / σ²_u_null)
==============================================================================*/

di as txt _n(2) "========================================================"
di as txt        "  CONDITIONAL ICC — PHILIPPINES"
di as txt        "========================================================"

* Unconditional school-level variances for Philippines (needed for VES)
* (These were estimated in Section 4; we re-estimate here for easy reference)

local cov_A "z_C_SCHSIZE z_STRATIO_PERM z_STRATIO_TEMP z_PCGIRLS_G5"
local cov_B "i.C_GENDER z_C_AGE i.C_ISCED z_resou z_hinder z_comatpa z_comatst z_comattc z_studiss"
local model_names "A_Resources B_Principal_Community"

foreach subj of local subjects {
    di as txt _n "  --- Subject: `subj_lbl_`subj'' ---"

    *--- Step 1: null model variance for Philippines (for VES denominator) ---
    matrix pv_icc_null = J(5,1,.)
    scalar var_u_null_avg = 0

    forvalues pv = 1/5 {
        quietly mixed PV`pv'_`subj' [pw=wt_scaled] ///
            if CNT == "PHL" || schoolid:, mle
        scalar var_u_null_`pv' = exp(2 * [lns1_1_1]_cons)
        scalar var_e_null_`pv' = exp(2 * [lnsig_e]_cons)
        scalar icc_null_`pv'   = var_u_null_`pv' / (var_u_null_`pv' + var_e_null_`pv')
        scalar var_u_null_avg  = var_u_null_avg + var_u_null_`pv'
    }
    scalar var_u_null_avg = var_u_null_avg / 5

    * Combine null ICC across PVs
    scalar icc_null_bar = 0
    scalar V_W_null = 0
    scalar V_B_null = 0
    forvalues pv = 1/5 {
        scalar icc_null_bar = icc_null_bar + icc_null_`pv'
    }
    scalar icc_null_bar = icc_null_bar / 5

    * SE for null model
    forvalues pv = 1/5 {
        quietly mixed PV`pv'_`subj' [pw=wt_scaled] ///
            if CNT == "PHL" || schoolid:, mle
        quietly estat icc
        scalar V_W_null = V_W_null + (r(se2))^2
        scalar V_B_null = V_B_null + (icc_null_`pv' - icc_null_bar)^2
    }
    scalar V_W_null = V_W_null / 5
    scalar V_B_null = V_B_null / 4
    scalar se_null  = sqrt(V_W_null + 1.2 * V_B_null)

    di as txt _n "    Null model (unconditional):"
    di as txt "      ICC = " %6.4f icc_null_bar "  SE = " %6.4f se_null

    *--- Step 2: conditional models A and B ---------------------------------*
    local m = 0
    foreach model of local model_names {
        local ++m

        if "`model'" == "A_Resources" {
            local covars "`cov_A'"
            local desc "Model A: School resources & size"
        }
        else {
            local covars "`cov_B'"
            local desc "Model B: Principal characteristics & community"
        }

        di as txt _n "    `desc':"

        matrix pv_icc_cond = J(5, 2, .)
        scalar var_u_cond_avg = 0

        forvalues pv = 1/5 {
            capture mixed PV`pv'_`subj' `covars' [pw=wt_scaled] ///
                if CNT == "PHL" || schoolid:, mle

            if _rc {
                di as txt "      PV`pv': model did not converge — skipping"
                matrix pv_icc_cond[`pv', 1] = .
                matrix pv_icc_cond[`pv', 2] = .
                continue
            }

            scalar var_u_c = exp(2 * [lns1_1_1]_cons)
            scalar var_e_c = exp(2 * [lnsig_e]_cons)
            scalar icc_c_k = var_u_c / (var_u_c + var_e_c)

            quietly estat icc
            scalar se_c_k  = r(se2)

            matrix pv_icc_cond[`pv', 1] = icc_c_k
            matrix pv_icc_cond[`pv', 2] = se_c_k
            scalar var_u_cond_avg = var_u_cond_avg + var_u_c
        }
        scalar var_u_cond_avg = var_u_cond_avg / 5

        * Combine across PVs (skip missing)
        scalar icc_cond_bar = 0
        scalar npv_ok = 0
        forvalues pv = 1/5 {
            if pv_icc_cond[`pv', 1] != . {
                scalar icc_cond_bar = icc_cond_bar + pv_icc_cond[`pv', 1]
                scalar npv_ok = npv_ok + 1
            }
        }
        scalar icc_cond_bar = icc_cond_bar / npv_ok

        scalar V_W_c = 0
        scalar V_B_c = 0
        forvalues pv = 1/5 {
            if pv_icc_cond[`pv', 1] != . {
                scalar V_W_c = V_W_c + (pv_icc_cond[`pv', 2])^2
                scalar V_B_c = V_B_c + (pv_icc_cond[`pv', 1] - icc_cond_bar)^2
            }
        }
        scalar V_W_c = V_W_c / npv_ok
        scalar V_B_c = V_B_c / (npv_ok - 1)
        scalar se_cond = sqrt(V_W_c + (1 + 1/npv_ok) * V_B_c)

        * Variance explained at school level (VES)
        scalar ves = 1 - (var_u_cond_avg / var_u_null_avg)

        di as txt "      ICC = " %6.4f icc_cond_bar "  SE = " %6.4f se_cond ///
                  "   95% CI: [" %6.4f (icc_cond_bar - 1.96*se_cond) ///
                  ", " %6.4f (icc_cond_bar + 1.96*se_cond) "]"
        di as txt "      Variance explained at school level (VES) = " %6.4f ves
    }
}


/*==============================================================================
  END OF DO FILE
==============================================================================*/

di as txt _n(2) "Analysis complete."
