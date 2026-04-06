/*==============================================================================
  SEA-PLM 2019 — School-Level ICC Analysis

  Purpose:
    1. Import the SEA-PLM Regional Students file
    2. Estimate unconditional school-level ICC for each country and subject
    3. Estimate conditional school-level ICC for the Philippines using
       school-level covariates merged from the school questionnaire file

  Method:
    Unconditional ICC: loneway (one-way ANOVA estimator, fast)
    Conditional ICC:   mixed followed by estat icc (r(icc2), r(se2))
    Multiple imputation combining: Rubin (1987) rules across 5 plausible values

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

capture log close
log using "$out_path/icc_analysis.log", replace text


/*==============================================================================
  SECTION 1 — IMPORT STUDENT DATA
==============================================================================*/

di as txt _n "Loading student data..."
import spss using "$data_path/SEA-PLM_Regional_Students 2.0.SAV", clear
rename cnt CNT
di as txt "Obs: `=_N'   Vars: `=c(k)'"


/*==============================================================================
  SECTION 2 — DATA PREPARATION
==============================================================================*/

*--- Unique school identifier (country + stratum + school) -------------------*
capture confirm string variable StIDStrt
if _rc tostring StIDStrt, gen(StIDStrt_s) format("%02.0f")
else   gen StIDStrt_s = StIDStrt

capture confirm string variable StIDSch
if _rc tostring StIDSch, gen(StIDSch_s) format("%03.0f")
else   gen StIDSch_s = StIDSch

capture confirm string variable CNT
if _rc decode CNT, gen(CNT_s)
else   gen CNT_s = CNT

gen schoolid_str = CNT_s + StIDStrt_s + StIDSch_s
encode schoolid_str, gen(schoolid)
label variable schoolid "Unique school identifier (CNT+stratum+school)"
drop StIDStrt_s StIDSch_s CNT_s schoolid_str

*--- Scaled weights (mean = 1 within country) --------------------------------*
bysort CNT: egen wt_mean = mean(WT2019)
gen wt_scaled = WT2019 / wt_mean
drop wt_mean
label variable wt_scaled "WT2019 scaled to mean=1 within country"


/*==============================================================================
  SECTION 3 — UNCONDITIONAL ICC BY COUNTRY AND SUBJECT  (via loneway)

  loneway uses the one-way ANOVA estimator: much faster than mixed.
  Returns r(rho) = ICC, r(se_rho) = SE, r(N) = students, r(ngroup) = schools.

  Combining rules across 5 PVs (Rubin 1987):
    ICC_bar  = mean(ICC_k)
    V_W      = mean(SE²_k)
    V_B      = sum((ICC_k − ICC_bar)²) / 4
    V_total  = V_W + (1 + 1/5) * V_B
    SE_comb  = sqrt(V_total)
==============================================================================*/

di as txt _n(2) "========================================================"
di as txt        "  UNCONDITIONAL ICC BY COUNTRY AND SUBJECT"
di as txt        "========================================================"

local countries KHM LAO MMR MYS PHL VNM
local subjects  M R W
local subj_lbl_M "Mathematics"
local subj_lbl_R "Reading"
local subj_lbl_W "Writing"

tempname icc_results
matrix `icc_results' = J(18, 4, .)   // 6 countries × 3 subjects
matrix colnames `icc_results' = ICC SE N_students N_schools
local row = 0

foreach cnt of local countries {
    foreach subj of local subjects {
        local ++row
        di as txt _n "  Country: `cnt'   Subject: `subj_lbl_`subj''"

        * Collect PV-specific ICC (loneway returns r(rho), r(N), r(N_g); no SE)
        forvalues pv = 1/5 {
            quietly loneway PV`pv'_`subj' schoolid if CNT == "`cnt'"
            scalar icc_`pv' = r(rho)
        }

        * Store sample sizes from last PV (stable across PVs)
        scalar n_stu    = r(N)
        scalar n_school = r(N_g)

        * Rubin's combining rules (between-PV variance only; loneway has no SE)
        scalar icc_bar = (icc_1 + icc_2 + icc_3 + icc_4 + icc_5) / 5
        scalar V_B = ((icc_1 - icc_bar)^2 + (icc_2 - icc_bar)^2 + ///
                      (icc_3 - icc_bar)^2 + (icc_4 - icc_bar)^2 + ///
                      (icc_5 - icc_bar)^2) / 4
        scalar se_comb = sqrt(1.2 * V_B)

        di as txt "    ICC = " %6.4f icc_bar ///
                  "   SE = " %6.4f se_comb ///
                  "   95% CI: [" %6.4f (icc_bar - 1.96*se_comb) ///
                  ", " %6.4f (icc_bar + 1.96*se_comb) "]" ///
                  "   N_students=" n_stu "  N_schools=" n_school

        matrix `icc_results'[`row', 1] = icc_bar
        matrix `icc_results'[`row', 2] = se_comb
        matrix `icc_results'[`row', 3] = n_stu
        matrix `icc_results'[`row', 4] = n_school
    }
}

di as txt _n "--- Summary: Unconditional School-Level ICC ---"
di as txt "Row order: KHM_M KHM_R KHM_W  LAO_M ... VNM_W"
matrix list `icc_results', format(%8.4f)


/*==============================================================================
  SECTION 4 — MERGE SCHOOL-LEVEL DATA FOR THE PHILIPPINES
==============================================================================*/

di as txt _n(2) "========================================================"
di as txt        "  LOADING AND MERGING SCHOOL DATA (Philippines)"
di as txt        "========================================================"

preserve

    tempfile school_data
    import spss using "$data_path/SEA-PLM_Regional_Schools 2.0.SAV", clear

    capture confirm variable CNT
    if _rc rename cnt CNT

    keep if CNT == "PHL"

    keep CNT StIDStrt StIDSch ///
         C_SCHSIZE SCHSIZE_G5 PCGIRLS_G5 ///
         STRATIO_PERM STRATIO_TEMP ///
         C_GENDER C_ISCED C_AGE ///
         resou hinder comatpa comatst comattc studiss

    capture confirm string variable StIDStrt
    if _rc tostring StIDStrt, replace format("%02.0f")
    capture confirm string variable StIDSch
    if _rc tostring StIDSch, replace format("%03.0f")
    capture confirm string variable CNT
    if _rc {
        decode CNT, gen(CNT_temp)
        drop CNT
        rename CNT_temp CNT
    }

    gen schoolid_str = CNT + StIDStrt + StIDSch
    drop CNT StIDStrt StIDSch
    duplicates report schoolid_str
    save `school_data', replace

restore

* Rebuild schoolid_str in student data for merge
capture confirm string variable StIDStrt
if _rc tostring StIDStrt, gen(StIDStrt_s) format("%02.0f")
else   gen StIDStrt_s = StIDStrt

capture confirm string variable StIDSch
if _rc tostring StIDSch, gen(StIDSch_s) format("%03.0f")
else   gen StIDSch_s = StIDSch

capture confirm string variable CNT
if _rc decode CNT, gen(CNT_s)
else   gen CNT_s = CNT

gen schoolid_str = CNT_s + StIDStrt_s + StIDSch_s
drop StIDStrt_s StIDSch_s CNT_s

merge m:1 schoolid_str using `school_data', keep(master match) gen(_merge_school)
di as txt "School merge results (Philippines):"
tab _merge_school if CNT == "PHL"
drop schoolid_str _merge_school

* Standardise continuous school covariates (PHL only)
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
  SECTION 5 — CONDITIONAL ICC FOR THE PHILIPPINES  (via mixed + estat icc)

  After each mixed call, estat icc returns:
    r(icc2)  = school-level ICC
    r(se2)   = SE of school-level ICC (delta method)

  Variance explained at school level (VES):
    VES = 1 − (mean σ²_u|X  /  mean σ²_u_null)
  School-level variance σ²_u is still drawn from the mixed coefficient
  [lns1_1_1]_cons since estat icc does not return it directly.

  Model A: School resources & size
  Model B: Principal characteristics & community context
==============================================================================*/

di as txt _n(2) "========================================================"
di as txt        "  CONDITIONAL ICC — PHILIPPINES"
di as txt        "========================================================"

local cov_A "z_C_SCHSIZE z_STRATIO_PERM z_STRATIO_TEMP z_PCGIRLS_G5"
local cov_B "i.C_GENDER z_C_AGE i.C_ISCED z_resou z_hinder z_comatpa z_comatst z_comattc z_studiss"

* Re-declare subjects (defensive; ensures the macro is set for this section)
local subjects M R W

foreach subj of local subjects {
    di as txt _n "  ===== Subject: `subj_lbl_`subj'' ====="

    *--- Null model (unconditional, PHL only) ---------------------------------*
    scalar var_u_null_sum = 0
    forvalues pv = 1/5 {
        quietly mixed PV`pv'_`subj' [pw=wt_scaled] if CNT == "PHL" || schoolid:, mle
        quietly estat icc
        scalar icc_`pv' = r(icc2)
        scalar se_`pv'  = r(se2)
        scalar var_u_null_sum = var_u_null_sum + exp(2 * [lns1_1_1]_cons)
    }
    scalar var_u_null_avg = var_u_null_sum / 5

    scalar icc_bar = (icc_1 + icc_2 + icc_3 + icc_4 + icc_5) / 5
    scalar V_W = ((se_1)^2 + (se_2)^2 + (se_3)^2 + (se_4)^2 + (se_5)^2) / 5
    scalar V_B = ((icc_1 - icc_bar)^2 + (icc_2 - icc_bar)^2 + ///
                  (icc_3 - icc_bar)^2 + (icc_4 - icc_bar)^2 + ///
                  (icc_5 - icc_bar)^2) / 4
    scalar se_comb = sqrt(V_W + 1.2 * V_B)

    di as txt _n "    Null model (unconditional):"
    di as txt "      ICC = " %6.4f icc_bar "  SE = " %6.4f se_comb

    *--- Conditional models ---------------------------------------------------*
    foreach spec in A B {
        local covars "`cov_A'"
        local desc   "Model A: School resources & size"
        if "`spec'" == "B" {
            local covars "`cov_B'"
            local desc   "Model B: Principal characteristics & community"
        }

        di as txt _n "    `desc':"

        scalar var_u_cond_sum = 0
        local npv_ok = 0

        forvalues pv = 1/5 {
            scalar icc_`pv' = .
            scalar se_`pv'  = .
            capture mixed PV`pv'_`subj' `covars' [pw=wt_scaled] ///
                if CNT == "PHL" || schoolid:, mle
            if _rc di as txt "      PV`pv': model did not converge -- skipping"
            if !_rc {
                quietly estat icc
                scalar icc_`pv' = r(icc2)
                scalar se_`pv'  = r(se2)
                scalar var_u_cond_sum = var_u_cond_sum + exp(2 * [lns1_1_1]_cons)
                local ++npv_ok
            }
        }

        * Rubin's rules (converged PVs only; cond() avoids nested if-block)
        scalar icc_bar = 0
        forvalues pv = 1/5 {
            scalar icc_bar = icc_bar + cond(icc_`pv' == ., 0, icc_`pv')
        }
        scalar icc_bar = icc_bar / `npv_ok'

        scalar V_W = 0
        scalar V_B = 0
        forvalues pv = 1/5 {
            scalar V_W = V_W + cond(icc_`pv' == ., 0, (se_`pv')^2)
            scalar V_B = V_B + cond(icc_`pv' == ., 0, (icc_`pv' - icc_bar)^2)
        }
        scalar V_W = V_W / `npv_ok'
        scalar V_B = V_B / (`npv_ok' - 1)
        scalar se_comb = sqrt(V_W + (1 + 1/`npv_ok') * V_B)

        scalar var_u_cond_avg = var_u_cond_sum / `npv_ok'
        scalar ves = 1 - (var_u_cond_avg / var_u_null_avg)

        di as txt "      ICC = " %6.4f icc_bar "  SE = " %6.4f se_comb ///
                  "   95% CI: [" %6.4f (icc_bar - 1.96*se_comb) ///
                  ", " %6.4f (icc_bar + 1.96*se_comb) "]"
        di as txt "      Variance explained at school level (VES) = " %6.4f ves
    }
}


/*==============================================================================
  END
==============================================================================*/

di as txt _n(2) "Analysis complete."
capture log close
