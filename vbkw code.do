**************************************************************************
**VBKW                                        							**
**Version: November 2019       		      								**
**Code by: Jessica Lum 													**
**Details of methods available in Garrido MM, Lum J, Pizer SD.          **
**Vector-based kernel weighting: A simple estimator for improving       **
**precision and bias of average treatment effects in multiple           **
**treatment settings. Stat Med. 2021; 40(5): 1204-1223.                 **
**																		**
**Contact jessica.lum2@va.gov or melissa.garrido@va.gov with questions  **
**************************************************************************




*Some info: 

*In this code, there are three treatment groups:
* medgrp == 1 for NSAID + Coxib
* medgrp == 3 for Opioid
* medgrp == 5 for Aceta

*Outcome variable of interest: edvisitspost

*The names of the propensity scores variables in this file: p_nsaid, p_opioid, p_aceta

*This file allows you to estimate ATT 1 vs. 3 | T = 1, ATT 1 vs. 3 | T = 3, and ATE 1 vs. 3
********************************************************************************
*The easiest way to use this file for your purposes is to keep the naming convention in this file
*as is, and just rename things later. 

*To do this: 
*1. If your treatment variable is named treat, with values 1, 2, and 3, change them to match the labels in this file so treatment levels are 1, 3, and 5, respectively. 
*2. Rename your treatment variable as medgrp 
*3. Rename your outcome variable as edvisitspost


********************************************************************************
*Please replace this with the dataset you're using:
use "data.dta", clear

*estimate propensity scores: Please replace x1 x2 x3... with the list of confounders you'd like balanced. 
mlogit medgrp x1 x2 x3... 
 
*After this, there should be nothing else for you to change in the dataset you're working with. You can run the rest of the code as is.
*Line 337 in this code is "sum vbkw_ATE13 vbkw_ATT13_1 vbkw_ATT13_3" and will give you the estimated effects, but not standard errors.
*Line 349 shows you how to get the standard errors from estimation of ATE 1 vs. 3 from GLM. 
 
*predicted probabilities:  
predict p_nsaid p_opioid p_aceta
 
********************************************************************************

g byte _support = 1 

sum p_nsaid if medgrp == 1 
scalar min1 = r(min) 
scalar max1 = r(max) 
sum p_nsaid if medgrp == 3 
scalar min3 = r(min) 
scalar max3 = r(max)
sum p_nsaid if medgrp == 5 
scalar min5 = r(min) 
scalar max5 = r(max)

scalar maxofmin = max(min1, min3, min5)
scalar minofmax = min(max1, max3, max5)

replace _support = 0 if p_nsaid < maxofmin | p_nsaid > minofmax


sum p_opioid if medgrp == 1 
scalar min1 = r(min) 
scalar max1 = r(max) 
sum p_opioid if medgrp == 3 
scalar min3 = r(min) 
scalar max3 = r(max)
sum p_opioid if medgrp == 5 
scalar min5 = r(min) 
scalar max5 = r(max)

scalar maxofmin = max(min1, min3, min5)
scalar minofmax = min(max1, max3, max5)

replace _support = 0 if p_opioid < maxofmin | p_opioid > minofmax



sum p_aceta if medgrp == 1 
scalar min1 = r(min) 
scalar max1 = r(max) 
sum p_aceta if medgrp == 3 
scalar min3 = r(min) 
scalar max3 = r(max)
sum p_aceta if medgrp == 5 
scalar min5 = r(min) 
scalar max5 = r(max)

scalar maxofmin = max(min1, min3, min5)
scalar minofmax = min(max1, max3, max5)

replace _support = 0 if p_aceta < maxofmin | p_aceta > minofmax



tab _support


 
*Drop units outside the common support 
keep if _support == 1 

*generate an id variable to keep track of sort order.  
g id = _n 

********************************************************************************


mata: mata clear
set matastrict on 
mata: 


void match(string viewvars, string strategy, real scalar N_ref, real scalar N_comp, real scalar bwidth, real scalar bwidth2, real scalar bwidth3)
{

		real scalar Nobs, pscore_ref, pscore_ref2, pscore_ref3, i, j
		real colvector pscore_tref, pscore_tref2, pscore_tref3, pscore_tcomp, pscore_tcomp2, pscore_tcomp3, dif, dif2, dif3, weight, _y
		real matrix X
		
		
		Nobs = pscore_ref = pscore_ref2 = pscore_ref3 = i = j = .
		pscore_tref = pscore_tref2 = pscore_tref3 = pscore_tcomp = pscore_tcomp2 = pscore_tcomp3 = dif = dif2 = dif3 = weight = _y = .

		if (strategy == "kw") 
		{
		st_view(X = ., ., tokens(viewvars), 0)
		}
		else if (strategy == "vbkw" & tokens(viewvars)[1] == "p_nsaid") {
		st_view(X = ., ., tokens(viewvars + " p_opioid" + " p_aceta"), 0)
		}
		else if (strategy == "vbkw" & tokens(viewvars)[1] == "p_aceta") {
		st_view(X = ., ., tokens(viewvars + " p_nsaid" + " p_opioid"), 0)
		}
		else if (strategy == "vbkw" & tokens(viewvars)[1] == "p_opioid") {
		st_view(X = ., ., tokens(viewvars + " p_nsaid" + " p_aceta"), 0)
		}
		
		Nobs = rows(X)
		pscore_tref = X[(1..N_ref), 1]
		pscore_tcomp = X[((N_ref + 1)..Nobs), 1]

		
		if (strategy == "vbkw") {
		pscore_tref2 = X[(1..N_ref), 7]
		pscore_tref3 = X[(1..N_ref), 8]
		pscore_tcomp2 = X[((N_ref + 1)..Nobs), 7]
		pscore_tcomp3 = X[((N_ref + 1)..Nobs), 8]
		}

		
		
		for (i = 1; i <= N_ref; i = i + 1) {
				pscore_ref = pscore_tref[i, 1]
				dif = abs(pscore_tcomp :- pscore_ref)
				if (strategy == "vbkw") {
				pscore_ref2 = pscore_tref2[i, 1]
				pscore_ref3 = pscore_tref3[i, 1]
				dif2 = abs(pscore_tcomp2 :- pscore_ref2) 
				dif3 = abs(pscore_tcomp3 :- pscore_ref3) 
				}
				weight = J(N_comp, 1, .)
				weight = (3/4) :* (1 :- (dif :/ bwidth) :^2)
				
				if (strategy == "kw"){
					dif = dif :<= bwidth
					weight = weight :* dif
					X[((N_ref + 1)..Nobs), 6] = X[((N_ref + 1)..Nobs), 6] :+ dif
					X[i, 6] = colsum(dif)
				}
				else {

				dif = dif :<= bwidth
				dif2 = dif2 :<= bwidth2
				dif3 = dif3 :<= bwidth3
				dif = dif :* dif2 :* dif3
				
				weight = weight :* dif
				X[((N_ref + 1)..Nobs), 6] = X[((N_ref + 1)..Nobs), 6] :+ dif
				X[i, 6] = colsum(dif)				
				}				
				
		if (mean(weight) == 0) X[i, 2] = 0
		if (mean(weight) == 0) X[i, 4] = 0
		if (mean(weight) != 0) weight = weight / sum(weight)
		X[((N_ref + 1)..Nobs), 4] = X[((N_ref + 1)..Nobs), 4] :+ weight
		_y = X[((N_ref + 1)..Nobs), 5] :* weight
		_editvalue(_y, 0, .)
		if (mean(weight) == 0) X[i, 3] = .
		if (mean(weight) != 0) X[i, 3] = sum(_y)
		}
		
}


end 

 
********************************************************************************
*generate vars needed for pairwise comparisons of t1 and t3: nsaid vs opioid 

cap drop medgrp1
gen medgrp1= 1 if medgrp==1
replace medgrp1 = 0 if medgrp==3

capture drop _AT*support* _AT*weight* _AT*trueY*

gen byte _ATTsupport = .
replace _ATTsupport = (medgrp1 <=1)
gen double _ATTtrueY = 0 if _ATTsupport 
 
gen byte _ATUsupport = .
replace _ATUsupport = (medgrp1 <=1)
gen double _ATUtrueY = 0 if _ATUsupport 

gen double _ATUweight = 1 - medgrp1 if _ATUsupport == 1 
gen double _ATTweight = medgrp1 if _ATTsupport == 1
********************************************************************************
*VBKW ATT 1 v 3 | T = 3 estimation:
cap drop n_used
g double n_used = 0 if medgrp1 !=. 

sort medgrp1 id 
local varlist p_opioid _ATUsupport _ATUtrueY _ATUweight edvisitspost n_used
local method vbkw 
count if medgrp == 3 
local nref = r(N)
count if medgrp == 1 
local ncomp = r(N)
sum p_opioid
local bandwidth = .2*r(sd)
sum p_nsaid 
local bandwidth2 = .2*r(sd)
sum p_aceta
local bandwidth3 = .2*r(sd)
mata: match("`varlist'", "`method'", `nref', `ncomp', `bandwidth', `bandwidth2', `bandwidth3')
replace _ATUsupport = 0 if _ATUweight == 0 | _ATUweight == . 


*for reference subjects, n_used gives us the number of comparison subjects used to calculate
*the counterfactual outcomes for each reference subject. 
sum n_used if _ATUsupport == 1 & medgrp == 3, detail
g double vbkw_ATT133_compused_mu = r(mean)
g double vbkw_ATT133_compused_p50 = r(p50)
g double vbkw_ATT133_compused_min = r(min)
g double vbkw_ATT133_compused_max = r(max)


*for comparison subjects, n_used gives the number of reference subjects each comparison subject
*has been matched to. 
sum n_used if _ATUsupport == 1 & medgrp == 1, detail
g double vbkw_ATT133_refmatched_mu = r(mean)
g double vbkw_ATT133_refmatched_p50 = r(p50)
g double vbkw_ATT133_refmatched_min = r(min)
g double vbkw_ATT133_refmatched_max = r(max)


*of the medgrp = 3 reference subjects on overall common support, how many of them were able to find matches? 
count if _ATUsupport == 1 & medgrp == 3
scalar on = r(N)
count if medgrp == 3
scalar total = r(N)

g double vbkw_133_matched_pct = (on/total)*100 
*ATT 1 v 3 | T = 3: 
sum _ATUtrueY if medgrp == 3 & _ATUsupport == 1 
scalar m1u = r(mean)
local N2 = r(N)
sum edvisitspost if medgrp == 3 & _ATUsupport == 1 
scalar m2u = r(mean)
scalar k_ATU = m1u - m2u 

g double vbkw_ATT13_3 = k_ATU



********************************************************************************
*VBKW ATT 1 v 3 | T = 1 estimation:
cap drop n_used
g double n_used = 0 if medgrp1 !=. 


gsort -medgrp1 id 
local varlist p_nsaid _ATTsupport _ATTtrueY _ATTweight edvisitspost n_used
local method vbkw 
count if medgrp == 1 
local nref = r(N)
count if medgrp == 3 
local ncomp = r(N)
sum p_nsaid
local bandwidth = .2*r(sd)
sum p_opioid
local bandwidth2 = .2*r(sd)
sum p_aceta
local bandwidth3 = .2*r(sd)
mata: match("`varlist'", "`method'", `nref', `ncomp', `bandwidth', `bandwidth2', `bandwidth3')
replace _ATTsupport = 0 if _ATTweight == 0 | _ATTweight == . 

*for reference subjects, n_used gives us the number of comparison subjects used to calculate
*the counterfactual outcomes for each reference subject. 
sum n_used if _ATTsupport == 1 & medgrp == 1, detail
g double vbkw_ATT131_compused_mu = r(mean)
g double vbkw_ATT131_compused_p50 = r(p50)
g double vbkw_ATT131_compused_min = r(min)
g double vbkw_ATT131_compused_max = r(max)


*for comparison subjects, n_used gives the number of reference subjects each comparison subject
*has been matched to. 
sum n_used if _ATTsupport == 1 & medgrp == 3, detail
g double vbkw_ATT131_refmatched_mu = r(mean)
g double vbkw_ATT131_refmatched_p50 = r(p50)
g double vbkw_ATT131_refmatched_min = r(min)
g double vbkw_ATT131_refmatched_max = r(max)


*of the medgrp = 1 reference subjects on overall common support, how many of them were able to find matches? 
count if _ATTsupport == 1 & medgrp == 1
scalar on = r(N)
count if medgrp == 1
scalar total = r(N)

g double vbkw_131_matched_pct = (on/total)*100 

*ATT 1 v 3 | T = 1: 
sum edvisitspost if medgrp == 1 & _ATTsupport == 1 
scalar m1t = r(mean)
local N1 = r(N)
sum _ATTtrueY if medgrp == 1 & _ATTsupport == 1 
scalar m2t = r(mean)
scalar k_ATT = m1t - m2t

g double vbkw_ATT13_1 = k_ATT 

********************************************************************************
*VBKW ATE 1 vs. 3: 

g double vbkw_ATE13 = (k_ATT*`N1'/(`N1'+`N2')) +  (k_ATU*`N2'/(`N1'+`N2')) 
gen double _ATEweight = _ATTweight + _ATUweight if medgrp1 <= 1


sum vbkw_ATE13 vbkw_ATT13_1 vbkw_ATT13_3
 
********************************************************************************

*This just renames all these variables to vbkw13_ATTsupport, vbkw13_ATTtrueY, etc. in case
*you'd like to run other pairwise comparisons. 
foreach i in _ATTsupport _ATTtrueY _ATUsupport _ATUtrueY _ATUweight _ATTweight _ATEweight{
	rename `i' vbkw13`i'
}

save vbkw_1vs3
********************************************************************************

 *The coefficient will give you ATE 1 vs. 3, which should match vbkw_ATE13 from line 337 in this code. 
glm edvisitspost ib3.medgrp [pw=vbkw13_ATEweight] if _support == 1 & inlist(medgrp, 1, 3), family(gamma) link(log) 

/*
. glm edvisitspost ib3.medgrp [pw=_ATEweight] if _support == 1 & inlist(medgrp, 1, 3), family(gamma) link(log) 

Iteration 0:   log pseudolikelihood = -68702.403  
Iteration 1:   log pseudolikelihood = -16600.125  
Iteration 2:   log pseudolikelihood = -15875.989  
Iteration 3:   log pseudolikelihood =   -15874.3  
Iteration 4:   log pseudolikelihood =   -15874.3  

Generalized linear models                         No. of obs      =     35,221
Optimization     : ML                             Residual df     =     35,219
                                                  Scale parameter =   15.67045
Deviance         =  70656.04259                   (1/df) Deviance =   2.006191
Pearson          =  551897.4089                   (1/df) Pearson  =   15.67045

Variance function: V(u) = u^2                     [Gamma]
Link function    : g(u) = ln(u)                   [Log]

                                                  AIC             =   .9015247
Log pseudolikelihood = -15874.30044               BIC             =  -298065.7

--------------------------------------------------------------------------------
               |               Robust
  edvisitspost |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
---------------+----------------------------------------------------------------
        medgrp |
NSAID + Coxib  |  -.0785505   .0377019    -2.08   0.037     -.152445   -.0046561
         _cons |   -.735372   .0334509   -21.98   0.000    -.8009345   -.6698094
--------------------------------------------------------------------------------
*/

********************************************************************************
*If you went from treat values of 1, 2, and 3 to medgrp values of 1, 3, and 5: 
*treat = 1 -> medgrp = 1
*treat = 2 -> medgrp = 3 
*treat = 3 -> medgrp = 5
*the estimated vbkw_ATT13_1 (see line 337) would actually represent ATT treat = 1 vs. treat = 2 amongst treat = 1 subjects
*the estimated vbkw_ATT13_3 (see line 337) would actually represent ATT treat = 1 vs. treat = 2 amongst treat = 2 subjects
*the estimated vbkw_ATE13 (see line 337) would actually represent ATE treat = 1 vs. treat = 2

*If you'd like to estimate comparisons representing treat = 1 vs. treat = 3, in line 18, the instruction would be changed to: 
	*1. If your treatment variable is named treat, with values 1, 2, and 3, change them to match the labels in this file so treatment levels are 1, 5, and 3, respectively. 
		*This makes it so that: 
			*treat = 1 -> medgrp = 1
			*treat = 2 -> medgrp = 5 
			*treat = 3 -> medgrp = 3
		*the estimated vbkw_ATT13_1 would actually represent ATT treat = 1 vs. treat = 3 amongst treat = 1 subjects
		*the estimated vbkw_ATT13_3 would actually represent ATT treat = 1 vs. treat = 3 amongst treat = 3 subjects
		*the estimated vbkw_ATE13 would actually represent ATE treat = 1 vs. treat = 3















