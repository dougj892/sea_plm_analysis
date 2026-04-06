

* this file entirely written by me, Doug
import spss using "/Users/douglasjohnson/Documents/Data/SEA-PLM/2019/SEA-PLM_Regional_Students 2.0.SAV", clear

keep cnt StIDStrt StIDSch PV*_M
keep if cnt == "PHL"

* generate a unique school ID
gen schoolid_str = StIDStrt + StIDSch
encode schoolid_str, gen(schoolid)

levelsof StIDStrt, local(strata)

local rho_sum = 0
local rho_count = 0

foreach region of local strata {
    quietly loneway PV2_M schoolid if StIDStrt =="`region'"
    local rho `r(rho)'
    display "`region' ICC: `rho'"
    local rho_sum = `rho_sum' + `rho'
    local rho_count = `rho_count' + 1
}

local rho_mean = `rho_sum' / `rho_count'
display "Average ICC across regions: `rho_mean'"


* estimate the share of variance explained by school level variables
preserve
    tempfile school_data
    import spss using "/Users/douglasjohnson/Documents/Data/SEA-PLM/2019/SEA-PLM_Regional_Schools 2.0.SAV", clear
    keep if cnt == "PHL"
    keep StIDStrt StIDSch ///
         C_SCHSIZE SCHSIZE_G5 PCGIRLS_G5 ///
         STRATIO_PERM STRATIO_TEMP ///
         C_GENDER C_ISCED C_AGE ///
         resou hinder comatpa comatst comattc studiss
    gen schoolid_str = StIDStrt + StIDSch
    duplicates report schoolid_str
    save `school_data', replace

restore
merge m:1 schoolid_str using `school_data', keep(master match) gen(_merge_school)

regress PV1_M C_SCHSIZE SCHSIZE_G5 PCGIRLS_G5 ///
         STRATIO_PERM STRATIO_TEMP ///
         C_GENDER C_ISCED C_AGE ///
         resou hinder comatpa comatst comattc studiss

