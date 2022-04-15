log using "Z:\Azathioprine\Dickson\persistence\ACKR1\AIM\Revision\Final_analysis.log", replace

/*
*******RACE, GENOTYPE, AND AZATHIOPRINE DISCONTINUATION: A COHORT STUDY*********
This study examines whether genotype predicts discontinuation of azathioprine attributed 
to hematopoietic toxicity; secondary outcomes include weight-adjusted last dose, 
last white blood cell count (WBC), last neutrophil count, and delta white blood cell count; 
the validation cohort assesses weight-adjusted dose by genotype


Primary Exposure: 
	ackr = rs2814778-CC genotype
Primary Outcome: 
	dose_depend = discontinuation of azathioprine attributed to leukopenia, neutropenia, 
	pancytopenia, anemia, and/or thrombocytopenia
Secondary Outcomes: 
	last_wbc = final WBC (closest to last dose in -30 to +3 days); 
	last_neutro = final neutrophil count (closest to last dose in -30 to +3 days); 
	delta_wbc = change in WBC from initial dose to last dose; 
	last_ratio = final weight-adjusted dose during follow=up; 
	avgdi = 6-month average dose intensity relative to target dose [validation cohort]

***************Table of Contents (follows order of appearance in manuscript)**************
Primary Cohort
	I. Information reported or referenced in manuscript, not in Table or Figure
	II. Table 1
	III. Supplement Table 3
	IV. Primary Analysis
	V. Supplement Table 4
	VI. Supplement Table 5
	VII. Figure 2
	VIII. Supplement Table 6
	IX. Supplement Figure 1
	X. Figure 3
	XI. Table 2
Responses to Reviewers
Validation Cohort
	XIII. Figure 4
	XIV. Supplement Figure 2


*******************************************************************************************
*******************************PRIMARY COHORT: VUMC COHORT*********************************
*******************************************************************************************/
use "Z:\Azathioprine\Dickson\persistence\ACKR1\final_ACKR1_dataset.dta"
describe
label list

*************I. Information in manuscript without reference to Table or Figure*************
*primary exposure
tab ackr

*overall cohort follow-up, age, race, and sex 
sum time, d
sum age, d
tab race
tab sex
tab ancestry race, col row

*assess follow-up (in days) by genotpe
bysort ackr: sum time, d
ranksum time, by(ackr)

*transform follow-up to months and assess by genotype
gen months=days/30.4167
sum months, d
bysort ackr: sum months, d
ranksum months, by (ackr)

*************II. Table 1 - BASELINE CHARACTERISTICS BY ackr (GENOTYPE)*******************
***RACE
tab ackr race, col row
*Standardized Mean Difference (SMD) -by manual calculation
quietly ci proportions race if ackr == 1, wilson
gen prop0 = r(proportion)
quietly ci proportions race if ackr == 0, wilson 
gen prop1 = r(proportion)
gen smd=(prop1-prop0)/sqrt( (prop1*(1-prop1) + prop0*(1-prop0)) /2)
display smd
drop prop0 prop1 smd
*Quality Control (QC) - stddiff package to check SMD, also check Fisher's and Chi-square
*ssc install stddiff
stddiff race, by(ackr)
tab ackr race, exact
tab ackr race, chi

***SEX
tab ackr sex, col row
*SMD - by manual calculation
quietly ci proportions sex if ackr == 1, wilson
gen prop0 = r(proportion)
quietly ci proportions sex if ackr == 0, wilson 
gen prop1 = r(proportion)
gen smd=(prop1-prop0)/sqrt( (prop1*(1-prop1) + prop0*(1-prop0)) /2)
display smd
drop prop0 prop1 smd
*QC - stddiff package to check SMD, also check Fisher's and Chi-square
stddiff sex, by(ackr)
tab ackr sex, exact
tab ackr sex, chi

***AGE
bysort ackr: sum age_at_baseline, d
*SMD -by manual calculation
ttest age, by(ackr)
gen smd = (  r(mu_1) - r(mu_2))/sqrt(( r(sd_1)^2 + r(sd_2)^2 )/2)
display smd
drop smd
*QC - stddiff package and Wilcoxon ranksum
stddiff age, by(ackr)
ranksum age_at, by (ackr)

***INDICATION (primary is categorical variable; sle, ibd, octd are the dummy binary variables)
*categorical 
tab ackr primary, col row
*QC - Fisher's exact and Chi-square
tab ackr primary, exact
tab ackr primary, chi

*SLE only
tab ackr sle, col row
*SMD - by manual calculation
quietly ci proportions sle if ackr == 1, wilson
gen prop0 = r(proportion)
quietly ci proportions sle if ackr == 0, wilson 
gen prop1 = r(proportion)
gen smd=(prop1-prop0)/sqrt( (prop1*(1-prop1) + prop0*(1-prop0)) /2)
display smd
drop prop0 prop1 smd
*QC - stddiff package
stddiff sle, by(ackr)

*OCTD only
tab ackr octd, col row
*SMD - by manual calculation
quietly ci proportions octd if ackr == 1, wilson
gen prop0 = r(proportion)
quietly ci proportions octd if ackr == 0, wilson 
gen prop1 = r(proportion)
gen smd=(prop1-prop0)/sqrt( (prop1*(1-prop1) + prop0*(1-prop0)) /2)
display smd
drop prop0 prop1 smd
*QC - stddiff package
stddiff octd, by(ackr)

*IBD only
tab ackr ibd, col row
*SMD - by manual calculation
quietly ci proportions ibd if ackr == 1, wilson
gen prop0 = r(proportion)
quietly ci proportions ibd if ackr == 0, wilson 
gen prop1 = r(proportion)
gen smd=(prop1-prop0)/sqrt( (prop1*(1-prop1) + prop0*(1-prop0)) /2)
display smd
drop prop0 prop1 smd
*QC - stddiff package
stddiff ibd, by(ackr)

***INITIAL DOSE
bysort ackr: sum initial_dose, d
*SMD -by manual calculation
ttest initial_dose, by(ackr)
gen smd = (  r(mu_1) - r(mu_2))/sqrt(( r(sd_1)^2 + r(sd_2)^2 )/2)
display smd
drop smd
*QC - stddiff package and Wilcoxon ranksum
stddiff initial_dose, by(ackr)
ranksum initial_dose, by(ackr)

***CALENDAR YEAR OF INITIAL DOSE
bysort ackr: sum cal_yr, d
*SMD -by manual calculation
ttest cal_yr, by(ackr)
gen smd = (  r(mu_1) - r(mu_2))/sqrt(( r(sd_1)^2 + r(sd_2)^2 )/2)
display smd
drop smd
*QC - stddiff package and Wilcoxon ranksum
stddiff cal_yr, by(ackr)
ranksum cal_yr, by(ackr)

***TESTED FOR TPMT BEFORE INITIATION
tab ackr tpmt_before, col row
*SMD - by manual calculation
quietly ci proportions tpmt_before if ackr == 1, wilson
gen prop0 = r(proportion)
quietly ci proportions tpmt_before if ackr == 0, wilson 
gen prop1 = r(proportion)
gen smd=(prop1-prop0)/sqrt( (prop1*(1-prop1) + prop0*(1-prop0)) /2)
display smd
drop prop0 prop1 smd
*QC - stddiff package, Fishers exact and Chi-square
tab ackr tpmt_before, exact
tab ackr tpmt_before, chi
stddiff tpmt_before, by(ackr)

***TPMT/NUDT15 METABOLIZER PHENOTYPE (overall is the categorical variable; normal, 
*intermediate, poor, and indeterminate are binary dummy variables)
tab ackr overall, col row
*QC - Fishers exact and Chi-square
tab ackr overall, exact
tab ackr overall, chi

*Normal
tab ackr normal, col row
*SMD - by manual calculation
quietly ci proportions normal if ackr == 1, wilson
gen prop0 = r(proportion)
quietly ci proportions normal if ackr == 0, wilson 
gen prop1 = r(proportion)
gen smd=(prop1-prop0)/sqrt( (prop1*(1-prop1) + prop0*(1-prop0)) /2)
display smd
drop prop0 prop1 smd
*QC - stddiff package
stddiff normal, by(ackr)

*Intermediate
tab ackr intermediate, col row
*SMD -by manual calculation
quietly ci proportions intermediate if ackr == 1, wilson
gen prop0 = r(proportion)
quietly ci proportions intermediate if ackr == 0, wilson 
gen prop1 = r(proportion)
gen smd=(prop1-prop0)/sqrt( (prop1*(1-prop1) + prop0*(1-prop0)) /2)
display smd
drop prop0 prop1 smd
*QC -stddiff package
stddiff intermediate, by(ackr)

*Poor
tab ackr poor, col row
*SMD -by manual calculation
quietly ci proportions poor if ackr == 1, wilson
gen prop0 = r(proportion)
quietly ci proportions poor if ackr == 0, wilson 
gen prop1 = r(proportion)
gen smd=(prop1-prop0)/sqrt( (prop1*(1-prop1) + prop0*(1-prop0)) /2)
display smd
drop prop0 prop1 smd
*QC - stddiff package
stddiff poor, by(ackr)

*Indeterminate
tab ackr indeterminate, col row
*SMD -by manual calculation
quietly ci proportions indeterminate if ackr == 1, wilson
gen prop0 = r(proportion)
quietly ci proportions indeterminate if ackr == 0, wilson 
gen prop1 = r(proportion)
gen smd=(prop1-prop0)/sqrt( (prop1*(1-prop1) + prop0*(1-prop0)) /2)
display smd
drop prop0 prop1 smd
*QC - stddiff package
stddiff indeterminate, by(ackr)

***BASELINE WBC
bysort ackr: sum baseline_wbc, d
*SMD - by manual calculation
ttest baseline_wbc, by(ackr)
gen smd = (  r(mu_1) - r(mu_2))/sqrt(( r(sd_1)^2 + r(sd_2)^2 )/2)
display smd
drop smd
*QC - stddiff package and Wilcoxon ranksum
stddiff baseline_wbc, by(ackr)
ranksum baseline_wbc, by(ackr)



************III. Supplement Table 3: Genotype by EHR-reported race and genetic ancestry****
tab ancestry race
tab ancestry ackr if race==0
tab ancestry ackr if race==1



***************IV. Primary Analysis (also reported in Supplement Table 4)****************
***Establish COMPETING RISK 
*prepare competing risk variable (discontinued for reasons other than hematopoietic toxicities, 
*including other side effects (e.g., fever) and non-side effects (e.g., costs))
gen other=.
replace other=1 if non_hem_excl==1|non_side_excl==1
label var other "Discontinued for reasons other than hematopoietic toxicity"
tab other dose_depend
tab other non_hem_excl
tab other non_side_excl

***PRIMARY ANALYSIS - outcome is discontinuation for hematopoetic toxicity (dose_depend==1), 
*competing risk is other discontinuation (other==1), censoring is all other end of followup
stset year, failure(dose_depend)
*Competing risk - unadjusted
stcrreg ackr, compete(other==1)

***QC for competing risk on unadjusted model
stcurve, cif at1(ackr=0) at2(ackr=1)
graph save "Graph" "Z:\Azathioprine\Dickson\persistence\ACKR1\competing risk cif.gph", replace
stcurve, cumhaz at1(ackr=0) at2(ackr=1)
graph save "Graph" "Z:\Azathioprine\Dickson\persistence\ACKR1\competing risk cumulative subhazard.gph", replace

*Competing risk - Model 1 (adjusted for age and sex)
stcrreg ackr sex age_at_baseline, compete(other==1)

*Competing risk - Model 2 (adjusted for age, sex, and grouped TPMT/NUDT15 metabolizer status)
*pheno is collapsed overall TPMT metabolizer status
tab pheno overall
stcrreg ackr sex age_at_baseline pheno, compete(other==1)

*Competing risk - Model 3 (first without, then with propensity score)(adjusted for age, sex, 
*grouped TPMT/NUDT15 metabolizer status, tertile of calendar year at initiation, and SLE indication)
*i.tert is tertile of calendar year at initiation
bysort tert: sum cal_yr
stcrreg ackr sex age_at_baseline pheno i.tert sle, compete(other==1)
*generate propensity score (PS) to be used as adjustor: ps_adj
logistic ackr sex age_at_baseline pheno i.tert sle
predict ps_adj
stcox ackr ps_adj
stcrreg ackr ps_adj, compete (other==1)


******QC checks of PS
*Check original PS for summary stats by genotype to confirm good overlap and no values <0.01 or >0.99
bysort ackr: sum ps_adj, d
*Histogram of original PS by genotype
twoway (histogram ps_adj if ackr==0, percent color(red%30))        ///
       (histogram ps_adj if ackr==1, percent color(blue%30))       ///
     , legend(order(1 "TT/TC" 2 "CC" ) col(1) ring(0) position(1)  ///
           subtitle(Patient Genotype))                             ///
       xlabels(0 (.05) .35)                                        ///
       xtitle("Propensity Score by Genotype")
*Confirm expected problem of extreme weighting
gen iptw_check =  ackr/ps_adj + (1-ackr)/(1-ps_adj)
sum iptw_check, d
bysort ackr: sum iptw_check, d
drop ps_adj iptw_check

***use linear predictor instead of untransformed propensity score as adjustor [reported in manuscript](Model 3)
logistic ackr sex age_at_baseline pheno i.tert sle
predict logodds, xb
stcox ackr logodds
stcrreg ackr logodds, compete (other==1)

***Additional QC for Primary analysis [incidence rates reported in manuscript]
*IRR
ir dose_depend ackr year
*Poisson (unadjusted & Model 3)
poisson dose_depend ackr1 , exposure(year) irr
poisson dose_depend ackr1 sex age_at sle i.tert pheno, exposure(year) irr
logistic ackr sex age_at_baseline pheno i.tert sle
predict ps_adj
poisson dose_depend ackr ps_adj, exposure(year) irr
drop ps_adj
*Poisson additionally adjusted for race (unadjusted & Model 3)
poisson dose_depend ackr1 race, exposure(year) irr
poisson dose_depend ackr1 race sex age_at sle i.tert pheno, exposure(year) irr
logistic ackr sex age_at_baseline pheno i.tert sle
predict ps_adj
poisson dose_depend ackr race ps_adj, exposure(year) irr
drop ps_adj

**********V. Supplement Table 4: Sensitivity Analyses (primary results also reported in this table)************
*Sensitivity/Additional Analyses 
* 1) analysis for those with available baseline WBC (unadjusted, Model 1, Model 2, Model 3)
stset year, failure (dose_depend)
stcrreg ackr if baseline_wbc!=., compete(other==1)
stcrreg ackr baseline_wbc, compete(other==1)
stcrreg ackr sex age_at_baseline baseline_wbc, compete(other==1)
stcrreg ackr sex age_at_baseline pheno baseline_wbc, compete(other==1)
stcrreg ackr sex age_at_baseline pheno i.tert sle baseline_wbc, compete(other==1)
logistic ackr sex age_at_baseline pheno i.tert sle baseline_wbc
predict ps_adj
stcrreg ackr ps_adj, compete (other==1)
drop ps_adj

* 2) restricting to 18 yo or older (unadjusted, Model 1, Model 2, Model 3)
stset year, failure (dose_depend)
stcrreg ackr if age_at>=18, compete(other==1)
stcrreg ackr sex age_at_baseline if age_at>=18, compete(other==1)
stcrreg ackr sex age_at_baseline pheno if age_at>=18, compete(other==1)
stcrreg ackr sex age_at_baseline pheno i.tert sle if age_at>=18, compete(other==1)
logistic ackr sex age_at_baseline pheno i.tert sle if age_at>=18
predict ps_adj
stcrreg ackr ps_adj if age_at>=18, compete (other==1)
drop ps_adj

* 3) stratified by TPMT status (not enough in indeterminate and poor for results)
stset year, failure (dose_depend)
bysort overall: stcrreg ackr, compete(other==1)
*NOTE: Model 1 has no convergence for poor or indeterminate metabolizers 
bysort overall: stcrreg ackr sex age_at, compete(other==1)
*NOTE: Model 2 has too much collinearity for meaningful difference of additional adjustment by overall metabolizer status
*bysort overall: stcrreg ackr sex age_at pheno, compete(other==1)  [too much collinearity] 
*NOTE: Model 3 has no convergence for poor or indeterminate metabolizers
*normal 
stcrreg ackr sex age_at_baseline pheno i.tert sle baseline_wbc if overall==0, compete(other==1)
logistic ackr sex age_at_baseline pheno i.tert sle baseline_wbc if overall==0
predict ps_adj
stcrreg ackr ps_adj if overall==0, compete (other==1)
drop ps_adj
*intermediate
stcrreg ackr sex age_at_baseline pheno i.tert sle baseline_wbc if overall==1, compete(other==1)
logistic ackr sex age_at_baseline pheno i.tert sle baseline_wbc if overall==1
predict ps_adj
stcrreg ackr ps_adj if overall==1, compete (other==1)
drop ps_adj
*poor (have to use capture because single event causes errors)
capture stcrreg ackr sex age_at_baseline pheno i.tert sle baseline_wbc if overall==2, compete(other==1)
capture logistic ackr sex age_at_baseline pheno i.tert sle baseline_wbc if overall==2
capture predict ps_adj
capture stcrreg ackr ps_adj if overall==2, compete (other==1)
capture drop ps_adj
*indeterminate (have to use capture because single event causes errors)
capture stcrreg ackr sex age_at_baseline pheno i.tert sle baseline_wbc if overall==3, compete(other==1)
capture logistic ackr sex age_at_baseline pheno i.tert sle baseline_wbc if overall==3
capture predict ps_adj
capture stcrreg ackr ps_adj if overall==3, compete (other==1)
capture drop ps_adj

* 4) without competing risk
stset year, failure(dose_depend)
stcox ackr
stcox ackr sex age_at_baseline
stcox ackr sex age_at_baseline pheno
stcox ackr sex age_at_baseline pheno i.tert sle
logistic ackr sex age_at_baseline pheno i.tert sle
predict ps_adj
stcox ackr ps_adj
drop ps_adj

* 5) alternative outcomes  - 5A & 5B

* 5A) more restrictive outcome (only neutropenia, leukopenia, and pancytopenia); 
*anemia and thrombocytopenia included as competing risks
gen restricted=0
replace restricted=1 if reason_leuko==1|reason_neutro==1|reason_pancy==1
gen compete=0
replace compete=1 if other==1 | (dose_depend==1 & reason_leuko==0 & reason_neutro==0 & reason_pancy==0)
stset year, failure(restricted)
stcrreg ackr, compete(compete==1)
stcrreg ackr sex age_at, compete(compete==1)
stcrreg ackr sex age_at pheno, compete(compete==1)
stcrreg ackr sex age_at pheno i.tert sle, compete(compete==1)
logistic ackr sex age_at pheno i.tert sle
predict ps_adj
stcrreg ackr ps_adj, compete(compete==1)
drop ps_adj

* 5B) other (other==1) discontinuation as primary outcome (dose_depend==1 as competing risk)
stset year, failure(other)
*Competing risk - unadjusted
stcrreg ackr, compete(dose_depend==1)
*Competing risk - Model 1 
stcrreg ackr sex age_at_baseline, compete(dose_depend==1)
*Competing risk - Model 2
stcrreg ackr sex age_at_baseline pheno, compete(dose_depend==1)
*Competing risk - Model 3 (first without, then with propensity score)
stcrreg ackr sex age_at_baseline pheno i.tert sle, compete(dose_depend==1)
logistic ackr sex age_at_baseline pheno i.tert sle
*generate propensity score to be used as adjustor
predict ps_adj
stcox ackr ps_adj
stcrreg ackr ps_adj, compete (dose_depend==1)
drop ps_adj


**************VI. Supplement Table 5 - Propensity score characteristics weighted and unweighted 
*(uses programs "propwt" and "pbalchk" from the "propensity" suite at  
*http://personalpages.manchester.ac.uk/staff/mark.lunt/; using command 
*"net from http://personalpages.manchester.ac.uk/staff/mark.lunt")*********************
*recreate primary analysis Model 3 propensity score
logistic ackr sex age_at_baseline pheno i.tert sle 
predict ps_adj
propwt ackr ps_adj, ipt smr
pbalchk ackr sex age_at_baseline pheno tert sle
pbalchk ackr sex age_at_baseline pheno tert sle , wt(ipt_wt)
drop ps_adj ipt smr


***************VII. Figure 2: Discontinuation by genotype in first 24 months*********
gen dose_depend_all=dose_depend
replace dose_depend_all=0 if months>24
gen month_graph=months
replace month_graph=24 if months>24 
stset month_graph, failure(dose_depend_all)
stcox ackr
mylabels 0(2.0)14.0, myscale(@/100) local(myla) 
sts graph, cumhaz by(ackr) risktable(0(6)24, size (small) order(1 "TT/TC" 2 "CC")) ///
     plot1opts(lpattern("-")) plot2opts(color(forest_green))  title("") ///
	 ytitle("Cumulative probability of treatment" "discontinuation due to hematopoietic toxicity (%)") ///
	 ylabel(`myla', angle(horizontal) ) ymtick(##2) xtitle("Time since start of azathioprine therapy (months)") ///
	 xlabel(0(6)24) xmtick(##5) legend(order(2 "rs2814778-CC" 1 "rs2814778-TT/rs2814778-TC") ///
	 cols(1) region(fcolor(none) lcolor(none)) position(11) ring(0)) title(, color(none) ring(0)) name(graph, replace)
graph export "Z:\Azathioprine\Dickson\persistence\ACKR1\ackr1_myelotoxicity_genotypeonly.pdf", replace
graph save "graph" "Z:\Azathioprine\Dickson\persistence\ackr1\genotypeonly.gph", replace




**************VIII. Supplement Table 6: Race and ACKR1 in competing risk survival analysis - 
*why ACKR is the important category*******************
* 1) Race as exposure, rather than ackr (unadjusted, Model 1, Model 2, and Model 3)
stset year, failure (dose_depend)
stcrreg race, compete (other==1)
stcrreg race sex age_at, compete (other==1)
stcrreg race sex age_at_baseline pheno, compete(other==1)
stcrreg race sex age_at_baseline pheno i.tert sle, compete(other==1)
logistic race sex age_at_baseline pheno i.tert sle
predict ps_adj
stcrreg race ps_adj, compete (other==1)
drop ps_adj

* 2) race adjusted by genotype (unadjusted, Model 1, Model 2, Model 3)
stcrreg race ackr1, compete(other==1)
stcrreg race ackr1 sex age_at, compete (other==1)
stcrreg race ackr sex age_at_baseline pheno, compete(other==1)
stcrreg race ackr sex age_at_baseline pheno i.tert sle, compete(other==1)
logistic race sex age_at_baseline pheno i.tert sle
predict ps_adj
stcrreg race ackr ps_adj, compete (other==1)
drop ps_adj

* 3) genotyped adjusted by race (unadjusted, Model 1, Model 2, Model 3)
stcrreg ackr1 race, compete(other==1)
stcrreg ackr1 race sex age_at, compete (other==1)
stcrreg ackr race sex age_at_baseline pheno, compete(other==1)
stcrreg ackr race sex age_at_baseline pheno i.tert sle, compete(other==1)
logistic ackr sex age_at_baseline pheno i.tert sle
predict ps_adj
stcrreg ackr race ps_adj, compete (other==1)
drop ps_adj

* 4) by genetic ancestry (>60% African genetic ancestry versus >60% European genetic ancestry  
*(8 people in neither category excluded)(unadjusted, Model 1, Model 2, Model 3)
stset year, failure(dose_depend)
stcrreg AA_60 if ancestry!=2, compete (other==1)
stcrreg AA_60 sex age_at if ancestry!=2, compete(other==1)
stcrreg AA_60 sex age_at pheno if ancestry!=2, compete(other==1)
stcrreg AA_60 sex age_at pheno i.tert sle if ancestry!=2, compete(other==1)
logistic AA_60 sex age_at pheno i.tert sle if ancestry!=2
predict ps_adj
stcrreg AA_60 ps_adj if ancestry!=2, compete(other==1)
drop ps_adj

* 5) stratified by genotype and race (not reporting white CC because too few for 
*significance; n=5)(unadjusted, Model 1, Model 2, Model 3)
*generate stratification
gen ehr_ackr = .
replace ehr_ackr=0 if race==0 & ackr==0
replace ehr_ackr=1 if race==0 & ackr==1
replace ehr_ackr=2 if race==1 & ackr==0
replace ehr_ackr=3 if race==1 & ackr==1
label var ehr_ackr "Stratified by EHR-reported race and ACKR1"
tab ehr_ackr dose_depend
*remove White CC from analysis
replace ehr_ack=. if race==0 & ackr==1

*Black-TT/TC and Black-CC relative to White-TT/TC
stset year, failure(dose_depend)
stcrreg i.ehr_ackr, compete(other==1)
stcrreg i.ehr_ackr sex age_at, compete(other==1)
stcrreg i.ehr_ackr sex age_at pheno, compete(other==1)
stcrreg i.ehr_ackr sex age_at pheno i.tert sle, compete(other==1)
logistic ehr_ackr sex age_at pheno i.tert sle if ehr_ackr==0|ehr_ackr==2
predict ps_adj
stcrreg i.ehr_ackr ps_adj if ehr_ackr==0|ehr_ackr==2, compete(other==1)
drop ps_adj
logistic ehr_ackr sex age_at pheno i.tert sle if ehr_ackr==0|ehr_ackr==3
predict ps_adj
stcrreg i.ehr_ackr ps_adj if ehr_ackr==0|ehr_ackr==3, compete(other==1)
drop ps_adj

* 6) other SNPs associated with African ancestries (sickle cell)(unadjusted, Model 1, Model 2, Model 3)
stset year, failure (dose_depend)
tab sickle
stcrreg sickle, compete(other==1)
stcrreg sickle sex age_at, compete(other==1)
stcrreg sickle sex age_at pheno, compete(other==1)
stcrreg sickle sex age_at pheno i.tert sle, compete(other==1)
logistic sickle sex age_at pheno i.tert sle
predict ps_adj
stcrreg sickle ps_adj, compete(other==1)
drop ps_adj

* 7) by genotype in Black race only (7A - with and 7B - without competing risk)
* 7A) with competing risk (unadjusted, Model 1, Model 2, Model 3)
stset year, failure(dose_depend)
stcrreg ackr if race==1, compete (other==1)
stcrreg ackr sex age_at if race==1, compete(other==1)
stcrreg ackr sex age_at pheno if race==1, compete(other==1)
stcrreg ackr sex age_at pheno i.tert sle if race==1, compete(other==1)
logistic ackr sex age_at pheno i.tert sle if race==1
predict ps_adj
stcrreg ackr ps_adj if race==1, compete(other==1)
drop ps_adj

* 7B) without competing risk (unadjusted, Model 1, Model 2, Model 3)
stcox ackr if race==1
stcox ackr sex age_at_baseline if race==1
stcox ackr sex age_at_baseline pheno if race==1
stcox ackr sex age_at_baseline pheno i.tert sle if race==1
logistic ackr sex age_at_baseline pheno i.tert sle if race==1
predict ps_adj
stcox ackr ps_adj if race==1
drop ps_adj


*************************IX. Supplement Figure 1 *******************************
stset year, failure (dose_depend)
stcox race
stcox race sex age_at_baseline
stcox race sex age_at_baseline pheno
logistic race sex age_at_baseline pheno i.tert sle
predict ps_adj
stcox race ps_adj
drop ps_adj

stset year, failure (other)
stcox race
stcox race sex age_at_baseline
stcox race sex age_at_baseline pheno
logistic race sex age_at_baseline pheno i.tert sle
predict ps_adj
stcox race ps_adj
drop ps_adj

*************X. Figure 3: Discontinuation by genotype and race in first 24 months****************
*EHR race and ACKR1 - no white CC (n=5)

stset year, failure (dose_depend)
stcox i.ehr_ack
stcrreg i.ehr_ack, compete(other==1)
stset month_graph, failure (dose_depend_all)
stcox i.ehr_ack
mylabels 0(2.0)14.0, myscale(@/100) local(myla) 
sts graph, cumhaz by(ehr_ackr) risktable(0(6)24, size (small) ///
     order(1 "White-TT/TC" 2 "Black-TT/TC" 3 "Black-CC")) plot1opts(lpattern("-")) plot2opts(lpattern("-..")) title("") ///
	 ytitle("Cumulative probability of treatment" "discontinuation due to hematopoietic toxicity (%)") ///
	 ylabel(`myla', angle(horizontal) ) ymtick(##2) xtitle("Time since start of azathioprine therapy (months)") ///
	 xlabel(0(6)24) xmtick(##5) legend(order(3 "Black-CC" 2 "Black-TT/TC" 1 "White-TT/TC") ///
	 cols(1) region(fcolor(none) lcolor(none)) position(11) ring(0)) title(, color(none) ring(0)) name(graph, replace)
graph export "Z:\Azathioprine\Dickson\persistence\ACKR1\ackr1_myelotoxicity_race_genotype.pdf", replace
graph save "graph" "Z:\Azathioprine\Dickson\persistence\ackr1\race_genotype.gph", replace

*************************SECONDARY OUTCOMES***********************************
***********XI. Table 2: Last WBC, Last Neutrophil Count, Change in WBC (initial to last), 
*and Weight-adjusted Last Dose*******************
*create change in WBC variable
gen delta_wbc=last_wbc-baseline_wbc if last_wbc!=. & baseline_wbc!=.
bysort ackr1: sum last_wbc last_neutro delta_wbc last_ratio, d
ranksum last_wbc, by (ackr1)
ranksum last_neutro, by (ackr1)
ranksum delta_wbc, by (ackr1)
ranksum last_ratio, by (ackr)
bysort ackr1: sum last_ratio if overall==0, d
ranksum last_ratio if overall==0, by (ackr)




********************************************************************************
************RESPONSES TO REVIEWERS (not reported in manuscript)*****************
********************************************************************************
*reviewer 1, comment 3
bysort primary: sum age_at, d
kwallis age, by (primary)

bysort primary: sum age_at if ackr==0, d
kwallis age if ackr==0, by (primary)

bysort primary: sum age_at if ackr==1, d
kwallis age if ackr==1, by (primary)

bysort primary: sum age_at if race==0, d
kwallis age if race==0, by (primary)

bysort primary: sum age_at if race==1, d
kwallis age if race==1, by (primary)

*reviewer 2, comment 3
bysort ackr: tab overall dose_depend, col row

*reviewer 2, comment 5
stset year, failure(dose_depend)
bysort primary: stcrreg ackr, compete(other==1)
bysort primary: stcrreg ackr sex age_at, compete(other==1)
bysort primary: stcrreg ackr sex age_at pheno, compete(other==1)
bysort primary: stcrreg ackr sex age_at pheno i.tert, compete(other==1)
logistic ackr sex age_at pheno i.tert if primary==1
predict ps_adj
stcrreg ackr ps_adj if primary==1, compete(other==1)
drop ps_adj
logistic ackr sex age_at pheno i.tert if primary==2
predict ps_adj
stcrreg ackr ps_adj if primary==2, compete(other==1)
drop ps_adj
logistic ackr sex age_at pheno i.tert if primary==3
predict ps_adj
stcrreg ackr ps_adj if primary==3, compete(other==1)
drop ps_adj

clear


********************************************************************************
******************************VALIDATION COHORT*********************************
********************************************************************************
*import st. jude data
import delimited "Z:\Azathioprine\ACKR1\aall03n1_black_akcr1_data_0930.txt"

*prep for STATA use
tab ethnicity
tab snp
tab avgdi
replace avgdi="" if avgdi=="NA"
destring avg, replace
tab rtgdi
replace rtg="" if rtg=="NA"
destring rtg, replace
tab mmpndi
replace mmpndi="" if mmpn=="NA"
destring mmpn, replace
tab afr
replace afr="" if afr=="NA"
destring afr, replace

*drop patients without dose intensity (outcome)
list pat_id if avgdi==.
drop if avgdi==.

*generate primary exposure variable: ackr
gen ackr=0
replace ackr=1 if rs281!="CC"
label def ackr1 0 "CC" 1 "TT or TC"
label values ackr ackr1
tab ackr
tab ackr ethnicity
tab ackr snp
tab ethnicity snp if tpmt==0

**REPORTED IN PAPER: all patients with predominantly African ancestry, normal TPMT
*average dose intensity
tab ackr if snp=="African" & tpmt==0
bysort ackr: sum avgdi if snp=="African" & tpmt==0, d
ranksum avg if snp=="African" & tpmt==0, by(ackr)

*metabolites
bysort ackr: sum rtg if snp=="African" & tpmt==0, d
ranksum rtg if snp=="African" & tpmt==0, by(ackr)

bysort ackr: sum mmpn if snp=="African" & tpmt==0, d
ranksum mmpn if snp=="African" & tpmt==0, by(ackr)

*******************XIII. Figure 3**********************************************
by ackr: egen med=median(avgdi)
by ackr: egen lqt=pctile(avgdi), p(25)
by ackr: egen uqt=pctile(avgdi), p(75)
by ackr: egen iqr=iqr(avgdi)
by ackr: egen mean=mean(avgdi)
by ackr: egen ls = min(max(avgdi, lqt-1.5*iqr))
by ackr: egen us = max(min(avgdi, uqt+1.5*iqr))
gen outliers = avgdi if(avgdi<=lqt-1.5*iqr | avgdi>=uqt+1.5*iqr)

twoway rbar lqt med ack, pstyle(p1) barw(.65) || rbar med uqt ack,pstyle(p1) barw(.65) || ///
     rspike lqt ls ack,pstyle(p1) || rspike uqt us ack, pstyle(p1) ||   ///
	 rcap ls ls ack, msize(*6) || rcap us us ack, msize(*6) pstyle(p1) ||  ///
	 scatter outliers ack, pstyle(p1)   legend(off) xsize(3) xtitle("") ///
	 xlabel( 0 "CC" 1 "TT/TC") ytitle(Average Dose Intensity (relative to 75 mg/m2))
graph save "Graph" "Z:\Azathioprine\Dickson\persistence\ACKR1\st jude.gph", replace
graph export "Z:\Azathioprine\Dickson\persistence\ACKR1\st_jude_dose.pdf", replace
*QC Figure 3
graph box avgdi if snp=="African" & tpmt==0, over(ackr) ytitle(Average Dose Intensity (relative to 75 mg/m2))


*******************XIV. Supplement Figure 2**************************************
*2A
graph box rtg if snp=="African" & tpmt==0, over(ackr) ytitle(Dose-adjusted log2 rTG (pmol/ 8 × 10^8 erythrocytes))
graph save "Graph" "Z:\Azathioprine\Dickson\persistence\ACKR1\st jude rtg.gph", replace
graph export "Z:\Azathioprine\Dickson\persistence\ACKR1\st_jude_rtg.pdf", replace
*2B
graph box mmpn if snp=="African" & tpmt==0, over(ackr) ytitle(Dose-adjusted log2 MMPN (pmol/ 8 × 10^8 erythrocytes))
graph save "Graph" "Z:\Azathioprine\Dickson\persistence\ACKR1\st jude mmpn.gph", replace
graph export "Z:\Azathioprine\Dickson\persistence\ACKR1\st_jude_mmpn.pdf", replace


clear

log close
exit
