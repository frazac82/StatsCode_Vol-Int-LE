
********************************************************************************************************************************************************
*CONFOUNDING ADJUSTMENT AND SENSITIVITY*****************************************************************************************************************
********************************************************************************************************************************************************

**#*** --> WOMEN | LINEAR NO INTERACTIONS: YLL (AGE TIME SCALE) 95%CI FOR SELECTED CENTILES COMPARISONS <-- ***
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

local maxage = 100								                        
local steps  = 0.1
local adj0   = "spring winter"    
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"                                       
local adj2   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3 bmi"                                       
local adj3   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3 bmi bdia bckd bia"   

mata: mata clear
mata:
	function udiff(at) {
		return(at[4]/at[3] - at[2]/at[1])		
	}
end

tempfile aucdiff_w
cap postclose stats
postfile stats agec model men maxage steps vol c_vol tin c_tin auc_act auc_act_lci auc_act_uci auc_con auc_con_lci auc_con_uci auc_sdiff auc_sdiff_lci auc_sdiff_uci using `aucdiff_w'

foreach tcon in 50 60 70 {							

	forvalues acov = 0/3 {
		
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
						
						post stats (`tcon') (`acov') (`q') (`maxage') (`steps') (`m`j'') (`j') (`u`w'') (`w') 						///
								   (`vauc_act') (`vauc_act_lci') (`vauc_act_uci') (`vauc_con') (`vauc_con_lci') (`vauc_con_uci') 	///
								   (`vauc_sdiff') (`vauc_sdiff_lci') (`vauc_sdiff_uci')
							
						di " --- AgeCond = `tcon' | Model = `acov' | Men = `q' | CVol = `j' | CTin = `w' --- $S_TIME  $S_DATE"
				}
			}	
		restore 
	}
}

postclose stats
use `aucdiff_w', replace
compress
gen modelling = "Centiles_R2"
save "R0_estimatesYLL_DIF_WOMEN", replace


**#*** --> MEN | LINEAR WITH INTERACTIONS: YLL (AGE TIME SCALE) 95%CI FOR SELECTED CENTILES COMPARISONS <-- ***
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

local maxage = 100								                            
local steps  = 0.1 								                           
local adj0   = "spring winter"    
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"                                       
local adj2   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3 bmi"                                       
local adj3   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3 bmi bdia bckd bia"  
   
mata: mata clear
mata:
	function udiff(at) {
		return(at[4]/at[3] - at[2]/at[1])		
	}
end

tempfile aucdiff_m
cap postclose stats
postfile stats agec model men maxage steps vol c_vol tin c_tin auc_act auc_act_lci auc_act_uci auc_con auc_con_lci auc_con_uci auc_sdiff auc_sdiff_lci auc_sdiff_uci using `aucdiff_m'

foreach tcon in 50 60 70 {								

	forvalues acov = 0/3 {
		
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
							
						di " --- AgeCond = `tcon' | Model = `acov' | Men = `q' | CVol = `j' | CTin = `w' --- $S_TIME  $S_DATE"
				}
			}	
		restore
	}
}

postclose stats
use `aucdiff_m', replace
compress
gen modelling = "Centiles_R2"
save "R0_estimatesYLL_DIF_MEN", replace


**#*** --> WOMEN | LINEAR NO INTERACTIONS: YLL (AGE TIME SCALE) 95%CI FOR MODELLING OVERALL EFFECT OF BRISK WALKING <-- ***
clear all
use "database2_volint_R2", replace

keep if male == 0
tab lsd,      gen(disab)
tab fvscore,  gen(dfvs)
tab alcf,     gen(deh)
tab gsleep,   gen(dsle)

local maxage = 100								                            
local steps  = 0.1 								                            
local adj0   = "spring winter"    
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"                                       
local adj2   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3 bmi"                                       
local adj3   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3 bmi bdia bckd bia"   

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

	foreach tcon in 50 60 70 {											

		forvalues acov = 0/3 {
			
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
					
					qui standsurv, 					 								///
						at1(vol = vol0 	  tin = tin0,    attimevar(tcon))			///	
						at2(vol = vol0 	  tin = tin0,    attimevar(tage))			///
						at3(vol = wvol`f' tin = wtin`f', attimevar(tcon)) 			///	
						at4(vol = wvol`f' tin = wtin`f', attimevar(tage)) 			///									  
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
gen modelling = "Brisk_Walking_R2"
save "R0_estimatesYLL_BW_WOMEN", replace


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

local maxage = 100								                            
local steps  = 0.1 								                            
local adj0   = "spring winter"    
local adj1   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3"                                       
local adj2   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3 bmi"                                       
local adj3   = "spring winter tws ndrugs dcurr2-dcurr8 nnca disab2 rmeat pmeat dfvs2 deh2-deh6 dsle2 dsle3 dsmok2 dsmok3 bmi bdia bckd bia" 

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

	foreach tcon in 50 60 70 {											

		forvalues acov = 0/3 {
			
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
gen modelling = "Brisk_Walking_R2"
save "R0_estimatesYLL_BW_MEN", replace


*********************************************************
*********BASELINE TABLES WITH CONFOUNDERS****************
*********************************************************
use "database1_volint_R2", replace
label variable bdia "Diabetes"
label variable bckd "Chronic kidney disease"
label variable bia  "Inflammatory arthritis"

*#*[TABLE 1] DESCRIPTIVE***
/*baseline tables by sex*/																								
baselinetable   							/*
*/	age(cts tab("p50 (p25, p75)")) 			/*
*/	tws(cts tab("p50 (p25, p75)")) 			/*
*/	ndrugs(cts tab("p50 (p25, p75)")) 		/*
*/	currempl(cat countformat(%15.0fc))		/*
*/  nnca(cts tab("p50 (p25, p75)"))     	/*
*/  lsd(cat countformat(%15.0fc) value(1))  /*
*/  bdia(cat countformat(%15.0fc) value(1)) /*
*/  bckd(cat countformat(%15.0fc) value(1)) /*
*/  bia(cat countformat(%15.0fc) value(1))  /*
*/  rmeat(cts tab("p50 (p25, p75)"))    	/*
*/  pmeat(cts tab("p50 (p25, p75)"))    	/*
*/  fvscore(cat countformat(%15.0fc))   	/*
*/  alcf(cat countformat(%15.0fc))      	/*
*/  gsleep(cat countformat(%15.0fc))    	/*   
*/	smok(cat countformat(%15.0fc))			/*
*/	bmi(cts tab("p50 (p25, p75)"))			/*
*/	v_days(cts tab("p50 (p25, p75)"))		/*
*/	vol(cts tab("p50 (p25, p75)")) 			/*
*/	tin(cts tab("p50 (p25, p75)") medianformat(%5.3f)) 									/*
*/	, by(male, total) countformat(%15.0fc) notable meanformat(%5.2f) sdformat(%5.1f) 	/*
*/	exportexcel("R0_Table1", replace)

*#*[TABLE RATES - COMBINED IN TABLE 1] MORTALITY RATES***
/*outcome rates*/																									
stset tdied, f(died==1) id(id)
sum tdied, d
strate, per(1000) output("table_rates_all", replace)
strate male, per(1000) output("table_rates_sex", replace)
use "table_rates_all", replace
append using "table_rates_sex"
sdecode male, replace
replace male = "All" if male == ""
renames male _D _Y _Rate _Lower _Upper \ Sex Events PYRs Rate1000 Lb95 Ub95
gen out = "All_cause"
foreach var of varlist Rate1000-Ub95 {
	tostring `var', replace format(%5.1f) force
}
gen rate = Rate1000 + " (" + Lb95 + ", " + Ub95 + ")"
drop out
replace PYRs = round(PYRs*1000)
format %15.0f PYRs
drop Rate1000-Ub95
tostring Events, replace
tostring PYRs, replace
order Sex, first
sxpose, clear firstnames
gen var = "Deaths" in 1
replace var = "PYRs" in 2
replace var = "Rate1000" in 3
order var Female Male All
export delimited using "R0_Table1_rates.csv", datafmt replace
erase "table_rates_all.dta"
erase "table_rates_sex.dta"

*#*[TABLE STRATIFIED] DESCRIPTIVE***
use "database1_volint_R2", replace
label variable bdia "Diabetes"
label variable bckd "Chronic kidney disease"
label variable bia  "Inflammatory arthritis"

xtilew vol3 = vol, within(male) nq(3)
xtilew tin3 = tin, within(male) nq(3)

tostring vol3, replace
tostring tin3, replace
gen sex = "Men" if male == 1
replace sex = "Women" if male == 0

gen combv = sex + "_" + vol3
gen combt = sex + "_" + tin3

*#*Baseline tables by sex_vol3*
baselinetable   							/*
*/	age(cts tab("p50 (p25, p75)")) 			/*
*/	tws(cts tab("p50 (p25, p75)")) 			/*
*/	ndrugs(cts tab("p50 (p25, p75)")) 		/*
*/	currempl(cat countformat(%15.0fc))		/*
*/  nnca(cts tab("p50 (p25, p75)"))     	/*
*/  lsd(cat countformat(%15.0fc))       	/*
*/  bdia(cat countformat(%15.0fc) value(1)) /*
*/  bckd(cat countformat(%15.0fc) value(1)) /*
*/  bia(cat countformat(%15.0fc) value(1))  /*
*/  rmeat(cts tab("p50 (p25, p75)"))    	/*
*/  pmeat(cts tab("p50 (p25, p75)"))    	/*
*/  fvscore(cat countformat(%15.0fc))   	/*
*/  alcf(cat countformat(%15.0fc))      	/*
*/  gsleep(cat countformat(%15.0fc))    	/*   
*/	smok(cat countformat(%15.0fc))			/*
*/	bmi(cts tab("p50 (p25, p75)"))			/*
*/	v_days(cts tab("p50 (p25, p75)"))		/*
*/	vol(cts tab("p50 (p25, p75)")) 			/*
*/	tin(cts tab("p50 (p25, p75)") medianformat(%5.3f)) 									/*
*/	, by(combv, total) countformat(%15.0fc) notable meanformat(%5.2f) sdformat(%5.1f) 	/*
*/	exportexcel("R0_TableV3", replace)

*#*Baseline tables by sex_tin3*
baselinetable   							/*
*/	age(cts tab("p50 (p25, p75)")) 			/*
*/	tws(cts tab("p50 (p25, p75)")) 			/*
*/	ndrugs(cts tab("p50 (p25, p75)")) 		/*
*/	currempl(cat countformat(%15.0fc))		/*
*/  nnca(cts tab("p50 (p25, p75)"))     	/*
*/  lsd(cat countformat(%15.0fc))       	/*
*/  bdia(cat countformat(%15.0fc) value(1)) /*
*/  bckd(cat countformat(%15.0fc) value(1)) /*
*/  bia(cat countformat(%15.0fc) value(1))  /*
*/  rmeat(cts tab("p50 (p25, p75)"))    	/*
*/  pmeat(cts tab("p50 (p25, p75)"))    	/*
*/  fvscore(cat countformat(%15.0fc))   	/*
*/  alcf(cat countformat(%15.0fc))      	/*
*/  gsleep(cat countformat(%15.0fc))    	/*   
*/	smok(cat countformat(%15.0fc))			/*
*/	bmi(cts tab("p50 (p25, p75)"))			/*
*/	v_days(cts tab("p50 (p25, p75)"))		/*
*/	vol(cts tab("p50 (p25, p75)")) 			/*
*/	tin(cts tab("p50 (p25, p75)") medianformat(%5.3f)) 									/*
*/	, by(combt, total) countformat(%15.0fc) notable meanformat(%5.2f) sdformat(%5.1f) 	/*
*/	exportexcel("R0_TableI3", replace)

*#*[TABLE RATES - COMBINED IN STRATIFIFED TABLES] MORTALITY RATES***
sencode combv, replace
sencode combt, replace
stset tdied, f(died==1) id(id)
strate combv, per(1000) output("table_rates_vol3", replace)
strate combt, per(1000) output("table_rates_tin3", replace)
use "table_rates_vol3", replace
append using "table_rates_tin3"
gen par = "vol" if combv !=.
replace par = "tin" if combt !=.
sdecode combv, replace
sdecode combt, replace
gen group = combv if combv != ""
replace group = combt if combt != ""
order par group, first
drop comb*
sort par group
split group, p("_")
drop group
renames _D _Y _Rate _Lower _Upper group1 group2 \ Events PYRs Rate1000 Lb95 Ub95 sex tertile
order par sex tertile, first
foreach var of varlist Rate1000-Ub95 {
	tostring `var', replace format(%5.1f) force
}
gen rate = Rate1000 + " (" + Lb95 + ", " + Ub95 + ")"
replace PYRs = round(PYRs*1000)
format %15.0f PYRs
drop Rate1000-Ub95
tostring Events, replace
tostring PYRs, replace
order par sex tertile, first
sxpose, clear firstnames
gen var = "Sex" in 1
replace var = "Tertile" in 2
replace var = "Deaths" in 3
replace var = "PYRs" in 4
replace var = "Rate1000" in 5
order var, first
set obs 6
replace var = "Par" if var == ""
foreach var of varlist tin-_var6 {
    replace `var' = "tin" if `var' == ""
}
foreach var of varlist vol-_var12 {
    replace `var' = "vol" if `var' == ""
}
export delimited using "R0_TableV3I3_rates.csv", datafmt replace
erase "table_rates_vol3.dta"
erase "table_rates_tin3.dta"