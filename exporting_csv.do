***************************************************
*Project: Assignment 2 -Micrometrics
*Purpose: Setting the panel & attrition
*Last modified: 14th October 2020 by Ishwara Hegde
***************************************************


***************************************************
*00. Preamble
***************************************************
clear all
set maxvar 30000
set more off

*Set user and the directory
local users  "Ishwara" // "Jonathan" "Aleksei"

	if "`users'" == "Federico" {
	global dir "Enter your directory here"
}

	if "`users'" == "Ishwara" {
		global dir "C:/Users/user/OneDrive - City University of Hong Kong/University/Cemfi/year_2/Micrometrics/Data-Code_AEJApp_MicrocreditMorocco/Output"
	}

***************************************************
*1. Preparing data for R
***************************************************
*Loading endline survey
use "$dir/endline_minienquete_outcomes.dta", clear

*Keeping only the necessary variables from the endline 
keep output_total profit_total id*
save "$dir/endline_vars.dta",replace

*Loading baseline 
use "$dir/baseline_minienquete_outcomes.dta", clear
keep id* m1 nadults_resid a7_11 d2_6 i1 a3_1 paire
merge 1:1 ident using "$dir/endline_vars.dta"

*keeping only the matched data
drop if _merge!=3
export delimited using "$dir/clean_data.csv", replace
save "$dir/clean_data.dta",replace

