

**********************************************************************************************************************************************************
**#HAZARD RATIOS AND ABSOLUTE VALUES FIGURE 1 -- ABSOLUTE VALUES FIGURE 2*********************************************************************************
**********************************************************************************************************************************************************

*******************
*FIGURE 1**********
*******************

*#** --> WOMEN | LINEAR NO INTERACTIONS: HR (AGE TIME SCALE) 95%CI FOR SELECTED CENTILES COMPARISONS <-- ***
clear all
use "database1_volint_R2", replace

keep if male == 0
tab currempl, gen(dcurr)
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)
tab smok,     gen(dsmok)
gen cage = age + tdied

foreach j in 10 25 50 75 90 {
	qui sum vol, d
	local v = `r(p`j')'
	gen vol`j' = `v'
	qui sum tin, d
	local t = `r(p`j')'
	gen tin`j' = `t'
}

qui stset cage, failure(died==1) enter(age)
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"     
stpm2 `adj1' vol tin, scale(hazard) df(4)
estimate store wadj1

preserve
clear
tempfile HR_w
save `HR_w', emptyok replace
restore

local rv10 = vol10[1]
local rt10 = tin10[1]

foreach j in 10 25 50 75 90 {
	
	foreach w in 10 25 50 75 90 {
		
		local vol`j' = vol`j'[1]
		local tin`w' = tin`w'[1]
		
		qui estimate restore wadj1

		qui nlcom (lnhr: ( (`vol`j'')*(_b[xb0:vol]) + (`tin`w'')*(_b[xb0:tin]) ) - ( (`rv10')*(_b[xb0:vol]) + (`rt10')*(_b[xb0:tin]) ) ), post
		
		preserve
		qui parmest, fast
		qui gen vol  = `vol`j''
		qui gen c_vol    = `j'
		qui gen tin  = `tin`w''
		qui gen c_tin    = `w'
		qui append using `HR_w'
		qui save `HR_w', replace
		restore
	}
}
use `HR_w', replace
gen modelling 	= "HR"
gen sex 	= "Women"
gen model 	= "Adj1"
order modelling parm model sex vol c_vol tin c_tin 
drop stderr z p
sort  vol tin
tempfile R1_HR_W
save `R1_HR_W', replace


*#** --> MEN | LINEAR WITH INTERACTIONS: HR (AGE TIME SCALE) 95%CI FOR SELECTED CENTILES COMPARISONS <-- ***
clear all
use "database1_volint_R2", replace

keep if male == 1
tab currempl, gen(dcurr)
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)
tab smok,     gen(dsmok)
gen cage = age + tdied
gen voltin = vol*tin

foreach j in 10 25 50 75 90 {
	qui sum vol, d
	local v = `r(p`j')'
	gen vol`j' = `v'
	qui sum tin, d
	local t = `r(p`j')'
	gen tin`j' = `t'
}
foreach j in 10 25 50 75 90 {
	foreach p in 10 25 50 75 90 {
		gen vol`j'tin`p' = vol`j'*tin`p'
	}
}

qui stset cage, failure(died==1) enter(age)
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"     
stpm2 `adj1' vol tin voltin, scale(hazard) df(4)
estimate store madj1

preserve
clear
tempfile HR_m
save `HR_m', emptyok replace
restore

local rv10  = vol10[1]
local rt10  = tin10[1]
local irt10 = vol10tin10[1]

foreach j in 10 25 50 75 90 {
	
	foreach w in 10 25 50 75 90 {
		
		local vol`j' = vol`j'[1]
		local tin`w' = tin`w'[1]
		local vol`j'tin`w' = vol`j'tin`w'[1]
		
		qui estimate restore madj1

		qui nlcom (lnhr: ( (`vol`j'')*(_b[xb0:vol]) + (`tin`w'')*(_b[xb0:tin]) + (`vol`j'tin`w'')*(_b[xb0:voltin]) ) - ( (`rv10')*(_b[xb0:vol]) + (`rt10')*(_b[xb0:tin]) + (`irt10')*(_b[xb0:voltin]) ) ), post
		
		preserve
		qui parmest, fast
		qui gen vol    = `vol`j''
		qui gen c_vol  = `j'
		qui gen tin    = `tin`w''
		qui gen c_tin  = `w'
		qui gen voltin = `vol`j'tin`w''
		qui gen inter  = "v`j'_t`w'"
		qui append using `HR_m'
		qui save `HR_m', replace
		restore

	}
}
use `HR_m', replace
gen modelling 	= "HR"
gen sex 	= "Men"
gen model 	= "Adj1"
order modelling parm model sex vol c_vol tin c_tin voltin inter
drop stderr z p
sort  vol tin
tempfile R1_HR_M
save `R1_HR_M', replace


*#** --> COMBINED WOMEN AND MEN | HR (AGE TIME SCALE) 95%CI FOR SELECTED CENTILES COMPARISONS <-- ***
use `R1_HR_W', replace
append using `R1_HR_M'
gen hr = exp(estimate)
gen lb95 = exp(min95)
gen ub95 = exp(max95)
drop modelling parm estimate min95 max95
replace model = "Adjusted" if model == "Adj1"
tempfile HR_wm
save `HR_wm', replace


*#** --> COMBINED WITH FIGURE 1 DATA | HR (AGE TIME SCALE) 95%CI FOR SELECTED CENTILES COMPARISONS <-- ***
use "database_yll_graphs.dta", clear
keep if agec == 60 & model == "Adjusted"
keep model sex c_vol c_tin auc_act* yld*
foreach var of varlist auc_act* {
    replace `var' = 60 + `var'
}
merge 1:1 model sex c_vol c_tin using `HR_wm', update
drop _merge
gen plotid = c_vol
replace c_vol = . if c_tin != 50
foreach var of varlist hr lb95 ub95 {
	tostring `var', replace format(%5.2f) force
}
gen HR = hr + " (" + lb95 + ", " + ub95 + ")"
replace HR = "Ref" if HR == "1.00 (1.00, 1.00)"
drop hr lb95 ub95
label variable HR "Hazard ratio (95% CI)"
tostring vol, replace format(%5.1f) force
tostring tin, replace format(%5.2f) force
foreach var of varlist c_vol c_tin {
	tostring `var', replace
}
replace c_vol = c_vol + " (" + vol + ")" if c_vol != "."
replace c_tin = c_tin + " (" + tin + ")"
label variable c_tin "Intensity centile (IG)"
label variable c_vol "Volume centile (mg)"
replace c_vol = "" if c_vol == "."

foreach sex in Men Women {
	forestplot auc_act auc_act_lci auc_act_uci if sex == "`sex'", effect("Life expectancy (years)") lcols(c_vol c_tin HR)   	///
	nonull nonames noov nosu nowt dp(1) classic boxscale(60) astext(40) textsize(105) xlabel(88(1)98, labsize(7pt) nogrid) 		///
	spacing(2) yline(5.5(5)20.5, lwidth(vthin) lpattern(vshortdash)) xtitle("Life expectancy at 60 years", size(9pt))			///
	leftjustify ciopts(lwidth(vthin)) plotid(plotid)         																	///
	box1opts(mcolor(red))    ci1opts(lcolor(red))    box2opts(mcolor(blue))   ci2opts(lcolor(blue))                  			///
	box3opts(mcolor(black))  ci3opts(lcolor(black))  box4opts(mcolor(orange)) ci4opts(lcolor(orange))               			///
	box5opts(mcolor(purple)) ci5opts(lcolor(purple)) 																			///
	title("`sex'", size(small)) name("Fig1_HR_`sex'", replace) xsize(6) ysize(4) scale(0.8) nodraw
}
graph combine Fig1_HR_Women Fig1_HR_Men, ycommon cols(2) xsize(8.5) ysize(6) nocopies scale(0.9) name("R1_Fig1", replace) nodraw
graph save "R1_Fig1" "R1_Fig1.gph", replace


*******************
*FIGURE 2**********
*******************
foreach sex in Men Women {
	qui forestplot yld yld_lb yld_ub if sex == "`sex'", effect("Difference (years)") lcols(c_vol c_tin)                        		///
    nonull nonames noov nosu nowt dp(1) classic boxscale(60) astext(40) textsize(105) xlabel(0(1)6, labsize(7pt) nogrid) 			///
	spacing(2) yline(5.5(5)20.5, lwidth(vthin) lpattern(vshortdash)) xtitle("Life expectancy difference at 60 years", size(9pt))	///
	leftjustify ciopts(lwidth(vthin)) plotid(plotid)         																		///
	box1opts(mcolor(red))    ci1opts(lcolor(red))    box2opts(mcolor(blue))   ci2opts(lcolor(blue))                  				///
	box3opts(mcolor(black))  ci3opts(lcolor(black))  box4opts(mcolor(orange)) ci4opts(lcolor(orange))               				///
	box5opts(mcolor(purple)) ci5opts(lcolor(purple)) 																				///
	xline(0, lcolor(black) lpattern(solid) lwidth(vthin)) 																			///
	title("`sex'", size(vsmall)) name("Fig2_`sex'", replace) xsize(6) ysize(4) scale(0.8) nodraw
}		
graph combine Fig2_Women Fig2_Men, ycommon cols(2) xsize(8) ysize(5.5) nocopies scale(1.15) name("R1_Fig2", replace) nodraw
graph save "R1_Fig2" "R1_Fig2.gph", replace


********************************************************************************************************************************************************
**#MULTIPLE IMPUTATIONS*********************************************************************************************************************************
********************************************************************************************************************************************************

**************
***CENTILES***
**************
clear all
use "database1_volint_R2_missing", replace
mdesc spring winter age male tws ndrugs nnca alcf pmeat currempl smok gsleep rmeat lsd fvscore
gen cage = age + tdied
stset cage, failure(died==1) enter(age)
sts gen H = na, by(male)

sencode currempl, gsort(currempl) replace
sencode fvscore, gsort(fvscore) replace
gen sln     = 1 if gsleep == "<7"
replace sln = 2 if gsleep == ">=7 to <=8"
replace sln = 3 if gsleep == ">8"
sencode gsleep, gsort(sln) replace
drop sln

local adj0 = "spring winter tws ndrugs i.currempl nnca i.lsd rmeat pmeat i.fvscore i.alcf i.gsleep i.smok"   
stpm2 `adj0' vol tin, scale(hazard) df(4)

mi set flong
mi register imputed tws alcf pmeat currempl smok gsleep rmeat lsd fvscore   
mi register regular spring winter age male ndrugs nnca vol tin
mi register passive _rcs* _d_rcs* _s0_rcs*

set seed 21722
mi impute chained (reg) tws pmeat rmeat (mlogit, augment) currempl (ologit) alcf smok gsleep fvscore (logit) lsd = spring winter age male ndrugs nnca vol tin H _d, add(10)
save "database1_volint_R2_imputed", replace


**#*** --> WOMEN | LINEAR NO INTERACTIONS <-- ***
clear all
use "database1_volint_R2_imputed", replace
keep if male == 0

foreach j in 10 25 50 75 90 {
	qui sum vol, d
	local v = `r(p`j')'
	gen vol`j' = `v'
	qui sum tin, d
	local t = `r(p`j')'
	gen tin`j' = `t'
}

tab currempl, gen(dcurr)
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)
tab smok,     gen(dsmok)

local maxage = 100								                        
local steps  = 0.1 								                       
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"                                       

mata: mata clear
mata:
	function udiff(at) {
		return(at[4]/at[3] - at[2]/at[1])		
	}
end

tempfile aucdiff_w
cap postclose stats
postfile stats agec model men maxage steps vol c_vol tin c_tin auc_act auc_act_lci auc_act_uci auc_con auc_con_lci auc_con_uci auc_sdiff auc_sdiff_lci auc_sdiff_uci imputed using `aucdiff_w'

foreach tcon in 60 {									

	foreach acov in 1 {
		
		forvalues mi = 1/10 {
		
			preserve
			
			mi extract `mi', clear
			
			local minage = `tcon'						
			
			qui sum male, meanonly
			local q = r(mean)

			qui stset cage, failure(died==1) enter(age)
			qui stpm2 `adj`acov'' vol tin, scale(hazard) df(4) 						

			qui range tage `minage' `maxage' `= 1 + (`maxage'-`minage')/`steps''
			qui gen tcon = `tcon' if tage !=.
			
			foreach j in 10 25 50 75 90 {
				
				foreach w in 10 25 50 75 90 {
					
						local m`j' = vol`j'[1]
						local u`w' = tin`w'[1]
						
						cap drop auc* sdiff* act* con*	
				
						qui standsurv, 					 						///		
							at1(vol = vol10 tin = tin10, attimevar(tcon)) 		///
							at2(vol = vol10 tin = tin10, attimevar(tage)) 		///								  
							contrast(ratio) atvar(resc*) ci contrastvars(con)
							qui drop resc*
								
						qui standsurv, 				 	 						///		
							at1(vol = vol`j' tin = tin`w', attimevar(tcon)) 	///
							at2(vol = vol`j' tin = tin`w', attimevar(tage)) 	///								  
							contrast(ratio) atvar(resa*) ci contrastvars(act)
							qui drop resa*
								
						qui standsurv, 					 						///		
							at1(vol = vol10  tin = tin10,  attimevar(tcon)) 	///
							at2(vol = vol10  tin = tin10,  attimevar(tage))    	///	
							at3(vol = vol`j' tin = tin`w', attimevar(tcon))   	///	
							at4(vol = vol`j' tin = tin`w', attimevar(tage))	    ///						  
							atvar(standres*) ci userfunction(udiff) userfunctionvar(sdiff)
							qui drop standres*

						foreach var of varlist con* act* sdiff* {
							qui integ `var' tage, gen(auc_`var') 
						}
						
						foreach nm in auc_act auc_act_lci auc_act_uci auc_con auc_con_lci auc_con_uci auc_sdiff auc_sdiff_lci auc_sdiff_uci {
							qui sum `nm' if tage==`maxage', meanonly
							local v`nm' = r(mean)
						}
						
						post stats (`tcon') (`acov') (`q') (`maxage') (`steps') (`m`j'') (`j') (`u`w'') (`w') 						///
								    (`vauc_act') (`vauc_act_lci') (`vauc_act_uci') (`vauc_con') (`vauc_con_lci') (`vauc_con_uci') 	///
								    (`vauc_sdiff') (`vauc_sdiff_lci') (`vauc_sdiff_uci') (`mi')
							
						di " --- AgeCond = `tcon' | Model = `acov' | Men = `q' | CVol = `j' | CTin = `w' | Imputed = `mi' --- $S_TIME  $S_DATE"
					}
				}	
			restore 
		}
	}
}
postclose stats
use `aucdiff_w', replace
compress
gen modelling = "Imputations_R2"
save "R1_estimatesYLL_DIF_WOMEN_R2_imputations", replace

	
**#*** --> MEN | LINEAR WITH INTERACTIONS: YLL (AGE TIME SCALE) 95%CI FOR SELECTED CENTILES COMPARISONS <-- ***
clear all
use "database1_volint_R2_imputed", replace

keep if male == 1
tab currempl, gen(dcurr)
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)
tab smok,     gen(dsmok)
gen voltin = vol*tin

foreach j in 10 25 50 75 90 {
	qui sum vol, d
	local v = `r(p`j')'
	gen vol`j' = `v'
	qui sum tin, d
	local t = `r(p`j')'
	gen tin`j' = `t'
}
foreach j in 10 25 50 75 90 {
	foreach p in 10 25 50 75 90 {
		gen vol`j'tin`p' = vol`j'*tin`p'
	}
}

local maxage = 100								                            
local steps  = 0.1 								                           
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3" 					                            
 
mata: mata clear
mata:
	function udiff(at) {
		return(at[4]/at[3] - at[2]/at[1])		
	}
end

tempfile aucdiff_m
cap postclose stats
postfile stats agec model men maxage steps vol c_vol tin c_tin auc_act auc_act_lci auc_act_uci auc_con auc_con_lci auc_con_uci auc_sdiff auc_sdiff_lci auc_sdiff_uci imputed using `aucdiff_m'

foreach tcon in 60 {									

	foreach acov in 1 {
		
		forvalues mi = 1/10 {
		
			preserve
			
			mi extract `mi', clear
			
			local minage = `tcon'						
			
			qui sum male, meanonly
			local q = r(mean)

			qui stset cage, failure(died==1) enter(age)
			qui stpm2 `adj`acov'' vol tin voltin, scale(hazard) df(4)

			qui range tage `minage' `maxage' `= 1 + (`maxage'-`minage')/`steps''
			qui gen tcon = `tcon' if tage !=.
			
			foreach j in 10 25 50 75 90 {
				
				foreach w in 10 25 50 75 90 {
					
						local m`j' = vol`j'[1]
						local u`w' = tin`w'[1]
			
						cap drop auc* sdiff* act* con*	
						
						qui standsurv, 					 											///		
							at1(vol = vol10 tin = tin10 voltin = vol10tin10, attimevar(tcon)) 		///
							at2(vol = vol10 tin = tin10 voltin = vol10tin10, attimevar(tage)) 		///								  
							contrast(ratio) atvar(resc*) ci contrastvars(con)
							qui drop resc*
								
						qui standsurv, 				 	 											///		
							at1(vol = vol`j' tin = tin`w' voltin = vol`j'tin`w', attimevar(tcon)) 	///
							at2(vol = vol`j' tin = tin`w' voltin = vol`j'tin`w', attimevar(tage)) 	///								  
							contrast(ratio) atvar(resa*) ci contrastvars(act)
							qui drop resa*
								
						qui standsurv, 					 											///		
							at1(vol = vol10  tin = tin10  voltin = vol10tin10,   attimevar(tcon)) 	///
							at2(vol = vol10  tin = tin10  voltin = vol10tin10,   attimevar(tage)) 	///	
							at3(vol = vol`j' tin = tin`w' voltin = vol`j'tin`w', attimevar(tcon)) 	///
							at4(vol = vol`j' tin = tin`w' voltin = vol`j'tin`w', attimevar(tage)) 	///						  
							atvar(standres*) ci userfunction(udiff) userfunctionvar(sdiff)
							qui drop standres*

						foreach var of varlist con* act* sdiff* {
							qui integ `var' tage, gen(auc_`var') 
						}
						
						foreach nm in auc_act auc_act_lci auc_act_uci auc_con auc_con_lci auc_con_uci auc_sdiff auc_sdiff_lci auc_sdiff_uci {
							qui sum `nm' if tage==`maxage', meanonly
							local v`nm' = r(mean)
						}
						
						post stats (`tcon') (`acov') (`q') (`maxage') (`steps') (`m`j'') (`j') (`u`w'') (`w') 						///
								    (`vauc_act') (`vauc_act_lci') (`vauc_act_uci') (`vauc_con') (`vauc_con_lci') (`vauc_con_uci') 	///
								    (`vauc_sdiff') (`vauc_sdiff_lci') (`vauc_sdiff_uci') (`mi')
							
						di " --- AgeCond = `tcon' | Model = `acov' | Men = `q' | CVol = `j' | CTin = `w' | Imputed = `mi' --- $S_TIME  $S_DATE"
				}
			}	
		restore
		}
	}
}
postclose stats
use `aucdiff_m', replace
compress
gen modelling = "Imputations_R2"
save "R1_estimatesYLL_DIF_MEN_R2_imputations", replace


*****************************
***MODELLING BRISK WALKING***
*****************************
clear all
use "database2_volint_R2_missing", replace
mdesc spring winter age male tws ndrugs nnca alcf pmeat currempl smok gsleep rmeat lsd fvscore
gen cage = age + tdied
stset cage, failure(died==1) enter(age)
sts gen H = na, by(male)

sencode currempl, gsort(currempl) replace
sencode fvscore, gsort(fvscore) replace
gen sln     = 1 if gsleep == "<7"
replace sln = 2 if gsleep == ">=7 to <=8"
replace sln = 3 if gsleep == ">8"
sencode gsleep, gsort(sln) replace
drop sln

local adj0 = "spring winter tws ndrugs i.currempl nnca i.lsd rmeat pmeat i.fvscore i.alcf i.gsleep i.smok"   
stpm2 `adj0' vol tin, scale(hazard) df(4)

mi set flong
mi register imputed alcf pmeat gsleep rmeat lsd fvscore   
mi register regular spring winter age male tws ndrugs nnca currempl smok vol tin   
mi register passive _rcs* _d_rcs* _s0_rcs*

set seed 21722
mi impute chained (reg) pmeat rmeat (ologit) alcf gsleep fvscore (logit) lsd = spring winter age male tws ndrugs nnca currempl smok vol tin H _d, add(10)
save "database2_volint_R2_imputed", replace


**#*** --> WOMEN | LINEAR NO INTERACTIONS: YLL (AGE TIME SCALE) 95%CI FOR MODELLING OVERALL EFFECT OF BRISK WALKING <-- ***
clear all
use "database2_volint_R2_imputed", replace
keep if male == 0
tab currempl, gen(dcurr)
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)
tab smok,     gen(dsmok)

local maxage = 100								                            
local steps  = 0.1 								                            
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"

mata: mata clear
mata:
	function udiff(at) {
		return(at[4]/at[3] - at[2]/at[1])		
	}
end

tempfile aucbw_w
cap postclose stats
postfile stats Npeople ndeaths Cpeople CVpeople CTpeople agec model men maxage steps eff auc_sdiff auc_sdiff_lci auc_sdiff_uci imputed using `aucbw_w'

foreach z in 10 25 50 {

	foreach tcon in 60 {												

		foreach acov in 1 {
			
				foreach f in 10 30 {
					
					forvalues mi = 1/10 {
						
						preserve
						
						mi extract `mi', clear
						
						local minage = `tcon'								
						
						qui sum male, meanonly
						local q = r(mean)
						
						qui _pctile vol, p(`z')								
						local Vsub`z' = r(r1)
			    
						qui _pctile tin, p(`z')
						local Tsub`z' = r(r1)

						qui stset cage, failure(died==1) enter(age)
						qui stpm2 `adj`acov'' vol tin, scale(hazard) df(4)

						qui keep if vol <= `Vsub`z''							
						qui keep if tin <= `Tsub`z''
				
						qui sum died
						qui local tot = r(N)
						qui local ndeaths = r(sum)
						
						qui range tage `minage' `maxage' `= 1 + (`maxage'-`minage')/`steps''
						qui gen tcon = `tcon' if tage !=.
						
						qui standsurv, 					 						///
							at1(vol = vol0 	  tin = tin0,    attimevar(tcon))	///	
							at2(vol = vol0 	  tin = tin0,    attimevar(tage))	///
							at3(vol = wvol`f' tin = wtin`f', attimevar(tcon)) 	///	
							at4(vol = wvol`f' tin = wtin`f', attimevar(tage)) 	///									  
							atvar(standres*) ci userfunction(udiff) userfunctionvar(sdiff)
						
						qui drop standres*

						foreach var of varlist sdiff* {
							qui integ `var' tage, gen(auc_`var')
						}
						
						foreach nm in auc_sdiff auc_sdiff_lci auc_sdiff_uci {
							qui sum `nm' if tage==`maxage', meanonly
							local v`nm' = r(mean)
						}

						post stats (`tot') (`ndeaths') (`z') (`Vsub`z'') (`Tsub`z'') (`tcon') (`acov') (`q') (`maxage') (`steps') (`f') (`vauc_sdiff') (`vauc_sdiff_lci') (`vauc_sdiff_uci') (`mi')
								
						di " --- Npeople = `tot' | ndeaths = `ndeaths' | CPeople = `z' | AgeCond = `tcon' | Model = `acov' | Men = `q' | Effect = Brisk_Walking-`f' | Imputed = `mi' --- $S_TIME  $S_DATE"
						
						restore
				}
			}
		}	
	}
}
postclose stats
use `aucbw_w', replace
compress
gen modelling = "Imputations_Brisk_Walking_R2"
save "R1_estimatesYLL_BW_WOMEN_R2_imputations", replace


**#*** --> MEN | LINEAR WITH INTERACTIONS: YLL (AGE TIME SCALE) 95%CI FOR MODELLING OVERALL EFFECT OF BRISK WALKING <-- ***
clear all
use "database2_volint_R2_imputed", replace

keep if male == 1
tab currempl, gen(dcurr)
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)
tab smok,     gen(dsmok)

gen voltin = vol*tin
gen wvol0tin0 = vol0*tin0
foreach j in 10 30 {
	gen wvol`j'tin`j' = wvol`j'*wtin`j'
}

local maxage = 100								                            
local steps  = 0.1 								                            
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"					                            

mata: mata clear
mata:
	function udiff(at) {
		return(at[4]/at[3] - at[2]/at[1])		
	}
end

tempfile aucbw_m
cap postclose stats
postfile stats Npeople ndeaths Cpeople CVpeople CTpeople agec model men maxage steps eff auc_sdiff auc_sdiff_lci auc_sdiff_uci imputed using `aucbw_m'

foreach z in 10 25 50 {

	foreach tcon in 60 {												

		foreach acov in 1 {
			
				foreach f in 10 30 {
					
					forvalues mi = 1/10 {

						preserve
						
						mi extract `mi', clear
						
						local minage = `tcon'								
						
						qui sum male, meanonly
						local q = r(mean)
						
						qui _pctile vol, p(`z')								
						local Vsub`z' = r(r1)
			    
						qui _pctile tin, p(`z')
						local Tsub`z' = r(r1)

						qui stset cage, failure(died==1) enter(age)
						qui stpm2 `adj`acov'' vol tin voltin, scale(hazard) df(4)

						qui keep if vol <= `Vsub`z''						
						qui keep if tin <= `Tsub`z''
			    
						qui sum died
						qui local tot = r(N)
						qui local ndeaths = r(sum)

						qui range tage `minage' `maxage' `= 1 + (`maxage'-`minage')/`steps''
						qui gen tcon = `tcon' if tage !=.
						
						qui standsurv, 	at1(vol = vol0    tin = tin0    voltin = wvol0tin0,     attimevar(tcon))	/// 	
								at2(vol = vol0    tin = tin0    voltin = wvol0tin0,     attimevar(tage))			///	
								at3(vol = wvol`f' tin = wtin`f' voltin = wvol`f'tin`f', attimevar(tcon)) 			///
								at4(vol = wvol`f' tin = wtin`f' voltin = wvol`f'tin`f', attimevar(tage))			///						  
								atvar(standres*) ci userfunction(udiff) userfunctionvar(sdiff)
						
						qui drop standres*

						foreach var of varlist sdiff* {
							qui integ `var' tage, gen(auc_`var')
						}
						
						foreach nm in auc_sdiff auc_sdiff_lci auc_sdiff_uci {
							qui sum `nm' if tage==`maxage', meanonly
							local v`nm' = r(mean)
						}

						post stats (`tot') (`ndeaths') (`z') (`Vsub`z'') (`Tsub`z'') (`tcon') (`acov') (`q') (`maxage') (`steps') (`f') (`vauc_sdiff') (`vauc_sdiff_lci') (`vauc_sdiff_uci') (`mi')
								
						di " --- Npeople = `tot' | ndeaths = `ndeaths' | CPeople = `z' | AgeCond = `tcon' | Model = `acov' | Men = `q' | Effect = Brisk_Walking-`f' | Imputed = `mi' --- $S_TIME  $S_DATE"
						
						restore
				}
			}
		}	
	}
}
postclose stats
use `aucbw_m', replace
compress
gen modelling = "Imputations_Brisk_Walking_R2"
save "R1_estimatesYLL_BW_MEN_R2_imputations", replace


**********************************
***SUMMARY MULTIPLE IMPUTATIONS***
**********************************

**************
***Centiles***
use "R1_estimatesYLL_DIF_WOMEN_R2_imputations", clear
append using "R1_estimatesYLL_DIF_MEN_R2_imputations"
keep men c_vol c_tin auc_act-imputed
drop auc_con auc_con_lci auc_con_uci
gen le1 = auc_act- auc_act_lci
gen ue1 = auc_act_uci- auc_act
gen re1 = auc_act/ auc_act_lci
gen re2 = auc_act_uci/ auc_act
drop le1-re2
tostring c_vol, replace
tostring c_tin, replace
gen est = "Vol(" + c_vol + ") -- " + "Int(" + c_tin + ")"
sencode est, replace
order men est imputed, first
drop c_tin
sort men est imputed

*Median absolute difference [BMC Medical Research Methodology 2009, 9:57 doi:10.1186/1471-2288-9-57]*
drop auc_act_lci auc_act_uci auc_sdiff_lci auc_sdiff_uci
by men est, sort : egen float med_act = median(auc_act)
by men est, sort : egen float med_sdiff = median(auc_sdiff)
gen dmed_act = abs(auc_act - med_act)
gen dmed_sdiff = abs(auc_sdiff - med_sdiff)
drop auc_act-med_sdiff
tostring men, replace
replace men = "Men" if men == "1"
replace men = "Women" if men == "0"
sencode men, replace gsort(-men)
sum dmed_act dmed_sdiff

stripplot dmed_act, over(est) by(men, note("") xrescale) yla(1/25, valuelabel labsize(small)) yla(, ticks) ysc(reverse) 						///
		    xlab(#5, format(%7.3f) labsize(small)) ytitle("") xtitle("Absolute deviance from median, life expectancy (years)")					///
		    ms(circle) separate(c_vol) name("dmed_act", replace) nodraw

stripplot dmed_sdiff, over(est) by(men, note("") xrescale) yla(1/25, valuelabel labsize(small)) yla(, ticks) ysc(reverse) 						///
		    xlab(#5, format(%7.3f) labsize(small)) ytitle("") xtitle("Absolute deviance from median, difference in life expectancy (years)")	///
		    ms(circle) separate(c_vol) name("dmed_sdiff", replace) nodraw

graph combine dmed_act dmed_sdiff, cols(1) scale(0.8) ysize(6) xsize(5) name("MI_centiles", replace)

****************************
***Modelling Walking Pace***
use "R1_estimatesYLL_BW_WOMEN_R2_imputations", clear
append using "R1_estimatesYLL_BW_MEN_R2_imputations"
keep men Cpeople eff auc_* imputed
tostring Cpeople, replace
tostring eff, replace
gen est = "People(" + Cpeople + ") -- " + "BW(" + eff + ")"
sencode est, replace
order men est imputed, first
drop eff
sort men est imputed
drop auc_sdiff_lci auc_sdiff_uci
by men est, sort : egen float med_sdiff = median(auc_sdiff)
gen dmed_sdiff = abs(auc_sdiff - med_sdiff)
drop auc_sdiff med_sdiff
tostring men, replace
replace men = "Men" if men == "1"
replace men = "Women" if men == "0"
sencode men, replace gsort(-men)
sum dmed_sdiff

stripplot dmed_sdiff, over(est) by(men, note("") xrescale) yla(1/6, valuelabel labsize(small)) yla(, ticks) ysc(reverse) 		///
		    xlab(#5, format(%7.3f) labsize(small)) ytitle("") ms(circle) separate(Cpeople)										///
		    xtitle("Absolute deviance from median, difference in life expectancy (years) associated with brisk walking")		///
		    name("dmed_sdiff_bw", replace)

		    

**********************************************************************************************************************************************************
**#SPLINES INTERNAL KNOTS NUMBER 4************************************************************************************************************************
**********************************************************************************************************************************************************

*#** --> BIC FOR SPLINES INTERNAL KNOTS NUMBER 4<-- ***
clear all
use "database1_volint_R2", replace

gen cage = age + tdied
sencode currempl, replace
sencode fvscore, replace
sencode gsleep, replace

forvalues k = 0/1 {
	
	preserve
	
	keep if male == `k'
	
	qui rcsgen vol, gen(vs) df(4) orthog				
	qui rcsgen tin, gen(ts) df(4) orthog

	forvalues i = 1/3 {
		forvalues j = 1/3 {
			qui gen vs`i'_ts`j' = vs`i' * ts`j'
		}
	}   

	qui stset cage, failure(died==1) enter(age)
	
	qui stpm2 spring winter tws ndrugs ib(1).currempl nnca ib(0).lsd rmeat pmeat ib(2).fvscore ib(1).alcf ib(3).gsleep ib(1).smok vol tin, scale(hazard) df(4) allbaselevels
	estimate store male`k'_l_no
	qui stpm2 spring winter tws ndrugs ib(1).currempl nnca ib(0).lsd rmeat pmeat ib(2).fvscore ib(1).alcf ib(3).gsleep ib(1).smok c.vol##c.tin, scale(hazard) df(4)
	estimate store male`k'_l_yes
	qui stpm2 spring winter tws ndrugs ib(1).currempl nnca ib(0).lsd rmeat pmeat ib(2).fvscore ib(1).alcf ib(3).gsleep ib(1).smok vs1 vs2 vs3 ts1 ts2 ts3, scale(hazard) df(4)
	estimate store male`k'_nl_no
	qui stpm2 spring winter tws ndrugs ib(1).currempl nnca ib(0).lsd rmeat pmeat ib(2).fvscore ib(1).alcf ib(3).gsleep ib(1).smok vs1 vs2 vs3 ts1 ts2 ts3 vs1_ts1-vs3_ts3, scale(hazard) df(4)
	estimate store male`k'_nl_yes
	
	restore
}

estimates stats male0_* male1_*, bicdetail



**********************************************************************************************************************************************************
**#DIFFERENT Tmax*****************************************************************************************************************************************
**********************************************************************************************************************************************************

*******************************
**CENTILES*********************
*******************************

*#** --> Max age at baseline and death/censoring <-- **
clear all
use "database1_volint_R2", replace
gen cage = age + tdied
tabstat age cage, statistics(max) by(male) columns(statistics)


*#** --> WOMEN | LINEAR NO INTERACTIONS: YLL (AGE TIME SCALE) 95%CI FOR SELECTED CENTILES COMPARISONS <-- **
clear all
use "database1_volint_R2", replace

keep if male == 0
tab currempl, gen(dcurr)
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)
tab smok,     gen(dsmok)
gen cage = age + tdied

foreach j in 10 25 50 75 90 {
	qui sum vol, d
	local v = `r(p`j')'
	gen vol`j' = `v'
	qui sum tin, d
	local t = `r(p`j')'
	gen tin`j' = `t'
}

local steps = 0.1
local adj1  = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"                                       

mata: mata clear
mata:
	function udiff(at) {
		return(at[4]/at[3] - at[2]/at[1])		
	}
end

tempfile aucdiff_w
cap postclose stats
postfile stats agec model men maxage steps vol c_vol tin c_tin auc_act auc_act_lci auc_act_uci auc_con auc_con_lci auc_con_uci auc_sdiff auc_sdiff_lci auc_sdiff_uci using `aucdiff_w'

foreach maxage in 105 {
	
	foreach tcon in 60 {								

		foreach acov in 1 {
			
				preserve
				
				local minage = `tcon'						
				
				qui sum male, meanonly
				local q = r(mean)

				qui stset cage, failure(died==1) enter(age)
				qui stpm2 `adj`acov'' vol tin, scale(hazard) df(4) 						

				qui range tage `minage' `maxage' `= 1 + (`maxage'-`minage')/`steps''
				qui gen tcon = `tcon' if tage !=.
				
				foreach j in 10 25 50 75 90 {	
					
					foreach w in 10 25 50 75 90 {
						
							local m`j' = vol`j'[1]
							local u`w' = tin`w'[1]
							
							cap drop auc* sdiff* act* con*	
					
							qui standsurv, 					 						///				
								at1(vol = vol10 tin = tin10, attimevar(tcon)) 		///
								at2(vol = vol10 tin = tin10, attimevar(tage)) 		///								  
								contrast(ratio) atvar(resc*) ci contrastvars(con)
								qui drop resc*
									
							qui standsurv, 				 	 						///				
								at1(vol = vol`j' tin = tin`w', attimevar(tcon)) 	///
								at2(vol = vol`j' tin = tin`w', attimevar(tage)) 	///								  
								contrast(ratio) atvar(resa*) ci contrastvars(act)
								qui drop resa*
									
							qui standsurv, 					 						///				
								at1(vol = vol10  tin = tin10,  attimevar(tcon)) 	///
								at2(vol = vol10  tin = tin10,  attimevar(tage))    	///	
								at3(vol = vol`j' tin = tin`w', attimevar(tcon))   	///	
								at4(vol = vol`j' tin = tin`w', attimevar(tage))	    ///						  
								atvar(standres*) ci userfunction(udiff) userfunctionvar(sdiff)
								qui drop standres*

							foreach var of varlist con* act* sdiff* {
								qui integ `var' tage, gen(auc_`var') 
							}
							
							foreach nm in auc_act auc_act_lci auc_act_uci auc_con auc_con_lci auc_con_uci auc_sdiff auc_sdiff_lci auc_sdiff_uci {
								qui sum `nm' if tage==`maxage', meanonly
								local v`nm' = r(mean)
							}
							
							post stats (`tcon') (`acov') (`q') (`maxage') (`steps') (`m`j'') (`j') (`u`w'') (`w') 					///
									   (`vauc_act') (`vauc_act_lci') (`vauc_act_uci') (`vauc_con') (`vauc_con_lci') (`vauc_con_uci') 	///
									   (`vauc_sdiff') (`vauc_sdiff_lci') (`vauc_sdiff_uci')
								
							di " --- Maxage = `maxage' | AgeCond = `tcon' | Model = `acov' | Men = `q' | CVol = `j' | CTin = `w' --- $S_TIME  $S_DATE"
					}
				}	
			restore 
		}
	}
}

postclose stats
use `aucdiff_w', replace
compress
gen modelling = "Tmax"
save "R1_Results/R1_WTmax", replace


*#** --> MEN | LINEAR WITH INTERACTIONS: YLL (AGE TIME SCALE) 95%CI FOR SELECTED CENTILES COMPARISONS <-- **
clear all
use "database1_volint_R2", replace

keep if male == 1
tab currempl, gen(dcurr)
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)
tab smok,     gen(dsmok)
gen cage = age + tdied
gen voltin = vol*tin

foreach j in 10 25 50 75 90 {
	qui sum vol, d
	local v = `r(p`j')'
	gen vol`j' = `v'
	qui sum tin, d
	local t = `r(p`j')'
	gen tin`j' = `t'
}
foreach j in 10 25 50 75 90 {
	foreach p in 10 25 50 75 90 {
		gen vol`j'tin`p' = vol`j'*tin`p'
	}
}

local steps  = 0.1 								                           
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"                                       
   
mata: mata clear
mata:
	function udiff(at) {
		return(at[4]/at[3] - at[2]/at[1])		
	}
end

tempfile aucdiff_m
cap postclose stats
postfile stats agec model men maxage steps vol c_vol tin c_tin auc_act auc_act_lci auc_act_uci auc_con auc_con_lci auc_con_uci auc_sdiff auc_sdiff_lci auc_sdiff_uci using `aucdiff_m'

foreach maxage in 105 {

	foreach tcon in 60 {								

		foreach acov in 1 {
			
				preserve
				
				local minage = `tcon'						
				
				qui sum male, meanonly
				local q = r(mean)

				qui stset cage, failure(died==1) enter(age)
				qui stpm2 `adj`acov'' vol tin voltin, scale(hazard) df(4)

				qui range tage `minage' `maxage' `= 1 + (`maxage'-`minage')/`steps''
				qui gen tcon = `tcon' if tage !=.
				
				foreach j in 10 25 50 75 90 {
					
					foreach w in 10 25 50 75 90 {
						
							local m`j' = vol`j'[1]
							local u`w' = tin`w'[1]
				
							cap drop auc* sdiff* act* con*	
							
							qui standsurv, 					 											///		
								at1(vol = vol10 tin = tin10 voltin = vol10tin10, attimevar(tcon)) 		///
								at2(vol = vol10 tin = tin10 voltin = vol10tin10, attimevar(tage)) 		///								  
								contrast(ratio) atvar(resc*) ci contrastvars(con)
								qui drop resc*
									
							qui standsurv, 				 	 											///		
								at1(vol = vol`j' tin = tin`w' voltin = vol`j'tin`w', attimevar(tcon)) 	///
								at2(vol = vol`j' tin = tin`w' voltin = vol`j'tin`w', attimevar(tage)) 	///								  
								contrast(ratio) atvar(resa*) ci contrastvars(act)
								qui drop resa*
									
							qui standsurv, 					 											///		
								at1(vol = vol10  tin = tin10  voltin = vol10tin10,   attimevar(tcon)) 	///
								at2(vol = vol10  tin = tin10  voltin = vol10tin10,   attimevar(tage)) 	///	
								at3(vol = vol`j' tin = tin`w' voltin = vol`j'tin`w', attimevar(tcon)) 	///
								at4(vol = vol`j' tin = tin`w' voltin = vol`j'tin`w', attimevar(tage)) 	///						  
								atvar(standres*) ci userfunction(udiff) userfunctionvar(sdiff)
								qui drop standres*

							foreach var of varlist con* act* sdiff* {
								qui integ `var' tage, gen(auc_`var') 
							}
							
							foreach nm in auc_act auc_act_lci auc_act_uci auc_con auc_con_lci auc_con_uci auc_sdiff auc_sdiff_lci auc_sdiff_uci {
								qui sum `nm' if tage==`maxage', meanonly
								local v`nm' = r(mean)
							}
							
							post stats (`tcon') (`acov') (`q') (`maxage') (`steps') (`m`j'') (`j') (`u`w'') (`w') 						///
									    (`vauc_act') (`vauc_act_lci') (`vauc_act_uci') (`vauc_con') (`vauc_con_lci') (`vauc_con_uci') 	///
									    (`vauc_sdiff') (`vauc_sdiff_lci') (`vauc_sdiff_uci')
								
							di " --- Maxage = `maxage' | AgeCond = `tcon' | Model = `acov' | Men = `q' | CVol = `j' | CTin = `w' --- $S_TIME  $S_DATE"
					}
				}	
			restore
		}
	}
}

postclose stats
use `aucdiff_m', replace
compress
gen modelling = "Tmax"
save "R1_Results/R1_MTmax", replace


*******************************
**MODELLING BRISK WALKING******
*******************************

**#*** --> WOMEN | LINEAR NO INTERACTIONS: YLL (AGE TIME SCALE) 95%CI FOR MODELLING OVERALL EFFECT OF BRISK WALKING <-- ***
clear all
use "database2_volint_R2", replace

keep if male == 0
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)

local maxage = 105								                            
local steps  = 0.1 								                            
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"                                       

mata: mata clear
mata:
	function udiff(at) {
		return(at[4]/at[3] - at[2]/at[1])		
	}
end

tempfile aucbw_w
cap postclose stats
postfile stats Npeople ndeaths Cpeople CVpeople CTpeople agec model men maxage steps eff auc_sdiff auc_sdiff_lci auc_sdiff_uci using `aucbw_w'

foreach z in 10 25 50 {

	foreach tcon in 60 {											

		foreach acov in 1 {
			
				foreach f in 10 30 {

					preserve
					
					local minage = `tcon'								
					
					qui sum male, meanonly
					local q = r(mean)
					
					qui _pctile vol, p(`z')								
					local Vsub`z' = r(r1)
                    
					qui _pctile tin, p(`z')
					local Tsub`z' = r(r1)

					qui stset cage, failure(died==1) enter(age)
					qui stpm2 `adj`acov'' vol tin, scale(hazard) df(4)

					qui keep if vol <= `Vsub`z''							
					qui keep if tin <= `Tsub`z''
			
					qui sum died
					qui local tot = r(N)
					qui local ndeaths = r(sum)
					
					qui range tage `minage' `maxage' `= 1 + (`maxage'-`minage')/`steps''
					qui gen tcon = `tcon' if tage !=.
					
					qui standsurv, 					 									///
						at1(vol = vol0 	  tin = tin0,    attimevar(tcon))				///	
						at2(vol = vol0 	  tin = tin0,    attimevar(tage))				///
						at3(vol = wvol`f' tin = wtin`f', attimevar(tcon)) 				///	
						at4(vol = wvol`f' tin = wtin`f', attimevar(tage)) 				///									  
						atvar(standres*) ci userfunction(udiff) userfunctionvar(sdiff)
					
					qui drop standres*

					foreach var of varlist sdiff* {
						qui integ `var' tage, gen(auc_`var')
					}
					
					foreach nm in auc_sdiff auc_sdiff_lci auc_sdiff_uci {
						qui sum `nm' if tage==`maxage', meanonly
						local v`nm' = r(mean)
					}

					post stats (`tot') (`ndeaths') (`z') (`Vsub`z'') (`Tsub`z'') (`tcon') (`acov') (`q') (`maxage') (`steps') (`f') (`vauc_sdiff') (`vauc_sdiff_lci') (`vauc_sdiff_uci')
							
					di " --- Npeople = `tot' | ndeaths = `ndeaths' | CPeople = `z' | AgeCond = `tcon' | Model = `acov' | Men = `q' | Effect = Brisk_Walking-`f' --- $S_TIME  $S_DATE"
					
					restore
			}
		}
	}	
}
postclose stats
use `aucbw_w', replace
compress
gen modelling = "BW_Tmax"
save "R1_WTmax_BW", replace


**#*** --> MEN | LINEAR WITH INTERACTIONS: YLL (AGE TIME SCALE) 95%CI FOR MODELLING OVERALL EFFECT OF BRISK WALKING <-- ***
clear all
use "database2_volint_R2", replace

keep if male == 1
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)

gen voltin = vol*tin
gen wvol0tin0 = vol0*tin0
foreach j in 10 30 {
	gen wvol`j'tin`j' = wvol`j'*wtin`j'
}

local maxage = 105								                            
local steps  = 0.1 								                            
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"                                       

mata: mata clear
mata:
	function udiff(at) {
		return(at[4]/at[3] - at[2]/at[1])		
	}
end

tempfile aucbw_m
cap postclose stats
postfile stats Npeople ndeaths Cpeople CVpeople CTpeople agec model men maxage steps eff auc_sdiff auc_sdiff_lci auc_sdiff_uci using `aucbw_m'

foreach z in 10 25 50 {

	foreach tcon in 60 {												

		foreach acov in 1 {
			
				foreach f in 10 30 {

					preserve
					
					local minage = `tcon'								
					
					qui sum male, meanonly
					local q = r(mean)
					
					qui _pctile vol, p(`z')								
					local Vsub`z' = r(r1)
                    
					qui _pctile tin, p(`z')
					local Tsub`z' = r(r1)

					qui stset cage, failure(died==1) enter(age)
					qui stpm2 `adj`acov'' vol tin voltin, scale(hazard) df(4)

					qui keep if vol <= `Vsub`z''							
					qui keep if tin <= `Tsub`z''
                    
					qui sum died
					qui local tot = r(N)
					qui local ndeaths = r(sum)

					qui range tage `minage' `maxage' `= 1 + (`maxage'-`minage')/`steps''
					qui gen tcon = `tcon' if tage !=.
					
					qui standsurv, 	at1(vol = vol0    tin = tin0    voltin = wvol0tin0,     attimevar(tcon))	///
							at2(vol = vol0    tin = tin0    voltin = wvol0tin0,     attimevar(tage))			///	
							at3(vol = wvol`f' tin = wtin`f' voltin = wvol`f'tin`f', attimevar(tcon)) 			///	
							at4(vol = wvol`f' tin = wtin`f' voltin = wvol`f'tin`f', attimevar(tage))			///						  
							atvar(standres*) ci userfunction(udiff) userfunctionvar(sdiff)
					
					qui drop standres*

					foreach var of varlist sdiff* {
						qui integ `var' tage, gen(auc_`var')
					}
					
					foreach nm in auc_sdiff auc_sdiff_lci auc_sdiff_uci {
						qui sum `nm' if tage==`maxage', meanonly
						local v`nm' = r(mean)
					}

					post stats (`tot') (`ndeaths') (`z') (`Vsub`z'') (`Tsub`z'') (`tcon') (`acov') (`q') (`maxage') (`steps') (`f') (`vauc_sdiff') (`vauc_sdiff_lci') (`vauc_sdiff_uci')
							
					di " --- Npeople = `tot' | ndeaths = `ndeaths' | CPeople = `z' | AgeCond = `tcon' | Model = `acov' | Men = `q' | Effect = Brisk_Walking-`f' --- $S_TIME  $S_DATE"
					
					restore
			}
		}
	}	
}
postclose stats
use `aucbw_m', replace
compress
gen modelling = "BW_Tmax"
save "R1_MTmax_BW", replace


*******************************
**GRAPHS***********************
*******************************

***CENTILES****
use "R1_WTmax", replace
append using "R1_MTmax"
drop steps vol tin 
renames auc_sdiff auc_sdiff_lci	auc_sdiff_uci \ yld yld_lb yld_ub
tostring men, replace
replace men = "Women" if men == "0"
replace men = "Men"   if men == "1"
rename men sex
tostring model, replace
replace model = "Adjusted" if model == "1"
order agec model sex c_vol c_tin
tempfile dtmax
save `dtmax', replace

use "database_yll_graphs.dta", clear
keep if agec == 60 & model == "Adjusted"
gen maxage = 100, after(sex)
tempfile dmain
save `dmain', replace

use `dtmax', clear
append using `dmain'
sort agec model sex c_vol c_tin maxage

keep model sex c_vol c_tin auc_act* yld* maxage
foreach var of varlist auc_act* {
    replace `var' = 60 + `var'
}
by sex c_vol c_tin, sort : egen float seq = seq()
replace c_vol = . if seq == 2
replace c_vol = . if c_tin != 50
replace c_tin = . if seq == 2
label variable c_tin "Intensity centile (IG)"
label variable c_vol "Volume centile (mg)"

foreach sex in Men Women {
	forestplot auc_act auc_act_lci auc_act_uci if sex == "`sex'", effect("Life expectancy (years)") lcols(c_vol c_tin)   	///
	nonull nonames noov nosu nowt dp(1) classic boxscale(30) astext(30) textsize(130) xlabel(88(1)101, labsize(7pt) nogrid) ///
	spacing(2) yline(10.5(10)40.5, lwidth(vthin) lpattern(vshortdash)) xtitle("Life expectancy at 60 years", size(8pt))		///
	leftjustify ciopts(lwidth(vthin)) plotid(seq)         																	///
	box1opts(mcolor(forest_green)) ci1opts(lcolor(forest_green)) box2opts(mcolor(gold)) ci2opts(lcolor(gold))              	///
	title("`sex'", size(small)) name("Fig_Tmax_`sex'", replace) xsize(8) ysize(4) scale(0.9) nodraw
}
graph combine Fig_Tmax_Women Fig_Tmax_Men, ycommon cols(2) xsize(8.5) ysize(6) nocopies scale(1.2) name("Fig_Tmax", replace) nodraw
graph save "Fig_Tmax" "Fig_TmaxC.gph", replace


***DIFFERENCES CENTILES****
foreach sex in Men Women {
	qui forestplot yld yld_lb yld_ub if sex == "`sex'", effect("Difference (years)") lcols(c_vol c_tin)                        		///
    nonull nonames noov nosu nowt dp(1) classic boxscale(60) astext(40) textsize(160) xlabel(0(1)8, labsize(7pt) nogrid) 			///
	spacing(2) yline(10.5(10)40.5, lwidth(vthin) lpattern(vshortdash)) xtitle("Life expectancy difference at 60 years", size(9pt))	///
	leftjustify ciopts(lwidth(vthin)) plotid(seq)         																			///
	box1opts(mcolor(forest_green)) ci1opts(lcolor(forest_green)) box2opts(mcolor(gold)) ci2opts(lcolor(gold))              			///
	xline(0, lcolor(black) lpattern(solid) lwidth(vthin)) 																			///
	title("`sex'", size(vsmall)) name("Fig_TmaxD_`sex'", replace) xsize(6) ysize(4) scale(0.8) nodraw
}		
graph combine Fig_TmaxD_Women Fig_TmaxD_Men, ycommon cols(2) xsize(8.5) ysize(6) nocopies scale(1.2) name("Fig_TmaxD", replace) nodraw
graph save "Fig_TmaxD" "Fig_TmaxCD.gph", replace


***MODELLING BRISK WALKING****
use "R1_WTmax_BW", replace
append using "R1_MTmax_BW"
drop Npeople ndeaths CVpeople CTpeople steps
renames Cpeople agec men eff auc_sdiff auc_sdiff_lci auc_sdiff_uci \ subset cond_age sex mins yld yld_lb yld_ub
tostring sex, replace
replace sex = "Women" if sex == "0"
replace sex = "Men"   if sex == "1"
tostring model, replace
replace model = "Adjusted" if model == "1"
order subset cond_age model sex mins maxage yld yld_lb yld_ub modelling
tempfile bwtmax
save `bwtmax', replace

use "database_BWyll_graphs.dta", clear
keep if cond_age == 60 & model == "Adjusted"
gen maxage = 100, after(mins)
tempfile bwmain
save `bwmain', replace

use `bwtmax', clear
append using `bwmain'
sort sex subset mins maxage

by sex subset: egen float seq = seq()
replace mins = . if seq == 2 | seq == 4
replace subset = . if seq != 1
label variable subset "Population (bottom centile)"
label variable mins "Activity (min)"

foreach sex in Men Women {
	qui forestplot yld yld_lb yld_ub if sex == "`sex'", effect("Difference (years)") lcols(subset mins)         							///
    nonull nonames noov nosu nowt dp(1) classic boxscale(110) astext(40) textsize(120) xlabel(0(0.5)4, labsize(10pt) nogrid format(%3.1f)) 	///
	spacing(2) yline(2.5(2)10.5, lwidth(vthin) lpattern(vshortdash)) xtitle(" Life expectancy difference at 60 years", size(12pt))			///
	leftjustify ciopts(lwidth(vthin)) plotid(seq)         																					///
	box1opts(mcolor(forest_green)) ci1opts(lcolor(forest_green)) box2opts(mcolor(gold)) ci2opts(lcolor(gold))              					///
	box3opts(mcolor(forest_green)) ci3opts(lcolor(forest_green)) box4opts(mcolor(gold)) ci4opts(lcolor(gold))              					///
	xline(0, lcolor(black) lpattern(solid) lwidth(vthin)) 																					///
	title("`sex'", size(med)) name("Fig_TmaxWP_`sex'", replace) xsize(6) ysize(4) scale(0.8) nodraw
}		
graph combine Fig_TmaxWP_Women Fig_TmaxWP_Men, ycommon cols(2) xsize(8.5) ysize(6) nocopies scale(0.9) name("Fig_TmaxWP", replace) nodraw
graph save "Fig_TmaxWP" "Fig_TmaxWP.gph", replace



**********************************************************************************************************************************************************
**#CAUSE-SPECIFIC DEATHS**********************************************************************************************************************************
**********************************************************************************************************************************************************
use "db_causedeath_volint_R1", replace
keep id icd10* dod
mdesc
drop icd10_2
merge 1:1 id using "database1_volint_R2", update
keep if _merge == 3    
drop _merge
tab died
groups icd10_1, missing
list id centre died icd10_1 dod if icd10_1 != "" & died == 0
replace icd10_1 = "" if died == 0
icd10 clean icd10_1, generate(check)
gen checkn = 1 if check != ""
tab checkn died, m
gen term     = "Cancer" if strpos(check, "C")>0
replace term = "CVD"    if strpos(check, "I")>0
replace term = "Other"  if term == "" & check != ""
tab term died, m
tab icd10_1 if died == 1 & term == "", m
tab term male if died == 1, m
tab term male if died == 1, col
drop if term == "" & died == 1
tab term male if died == 1, m col

gen cancer = 1 if strpos(check, "C")>0
gen cvd    = 1 if strpos(check, "I")>0
gen other  = 1 if cancer != 1 & cvd != 1 & died == 1
foreach c in cancer cvd other {
	replace `c' = 0 if `c' == .
	di "--------------------------------------------------------------------"
	tab died `c', m
	tab `c' male, col
	di "--------------------------------------------------------------------"
}
tempfile csd
save `csd', replace

**************
***WOMEN******
keep if male == 0
tab currempl, gen(dcurr)
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)
tab smok,     gen(dsmok)
gen cage = age + tdied

foreach c in cancer cvd other {
	stset cage, f(`c' == 1) enter(age)
	stpm2 spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3 vol tin, scale(hazard) df(4)
	estimate store `c'
}	
			
