# SEA-PLM ICC Analysis

This repository contains Stata code analyzing school-level intraclass correlations (ICCs) using data from the **SEA-PLM 2019** (Southeast Asia Primary Learning Metrics) assessment. SEA-PLM is a regional assessment of Grade 5 students in six countries — Cambodia, Laos, Myanmar, Malaysia, the Philippines, and Vietnam — covering Mathematics, Reading, and Writing. The data include student-level files with five plausible values per subject and sampling weights, as well as a school-level questionnaire file with principal and school characteristics.

## Analysis Files

**[Analysis/icc_analysis.do](Analysis/icc_analysis.do)** — More complicated analysis written by Claude Code. Estimates unconditional school-level ICCs for all six countries and three subjects using `loneway`, combining estimates across five plausible values via Rubin's (1987) rules. Then, for the Philippines, estimates conditional ICCs using two sets of school-level covariates (school resources/size; principal characteristics and community context) via multilevel models (`mixed` + `estat icc`), and reports the share of school-level variance explained by each covariate set.

**Results:**
According to this do file, the ICC for the Philippines for math is .43.


**[Analysis/icc_analysis_doug.do](Analysis/icc_analysis_doug.do)** — I created this simpler do file because I wanted to estimate ICC at the regional level and also didn't trust the conditional ICC estimates. Unlike the other do file:

- I just use data from the Philippines
- I just use the first plausible value for math (rather than estimating the ICCs for each of the 5 plausible values and combining using the multiple imputation formula)
- I use the simpler `loneway' command
- I estimate conditional ICC by getting the r squared from a regression of math scores on school level variables and then manually adjusting

I am fairly confident in all of these simplifications. 

**Results**
The ICC of the first plausible value for math for region 5 is .22 but the ICC for other regions are higher. The average across regions is .36

If we go with an ICC of .25 for region 

