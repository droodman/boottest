cd "D:\OneDrive\Documents\Work\Econometrics\Wild cluster"
set more off
set rmsg off
set trace off
set linesize 200
cap set processors 4

cap program drop myprobit
program myprobit // custom likelihood evaluator
  args lnf theta
  quietly replace `lnf' = lnnormal((2*$ML_y1-1) * `theta')
end

cap log close
qui log using "D:\OneDrive\Documents\Macros\boottest\unit tests.log", replace

version 13

set seed 0193284710

use collapsed, clear

qui regress hasinsurance selfemployed post post_self, cluster(year)
boottest post_self=.04, nogr // wild bootstrap, Rademacher weights, null imposed, 999 replications
boottest post_self=.04, weight(webb) noci // wild bootstrap, Webb weights, null imposed, 999 replications, no graph or CI
scoretest post_self=.04, nogr                   // Rao score/Lagrange multipler test of same

boottest post_self post, reps(999) weight(webb) nogr // wild bootstrap test of joint null, Webb weights, null imposed, 9,999 replications
boottest (post_self) (post), reps(999) weight(webb) nogr // same
boottest {post_self=.04} {post}, nogr // separate tests, no correction for multiple hypotheses
boottest {(post) (post_self=.04)} {(post) (post_self=.08)}, madj(sidak) nogr  // separate tests, Sidak correction for  multiple hypotheses

use nlsw88

qui regress wage tenure ttl_exp collgrad, cluster(industry)
boottest tenure, svmat nogr        // wild bootstrap test of joint null, Rademacher weights, null imposed, saving  simulated distribution

constraint 1 ttl_exp = .2
qui cnsreg wage tenure ttl_exp collgrad, constr(1) cluster(industry)
boottest tenure, nogr // wild bootstrap test of tenure=0, conditional on ttl_exp=2, Rademacher weights, null imposed, 999 replications

regress wage tenure ttl_exp collgrad south#union, cluster(industry)
margins south
boottest, margins nogr // bootstrap CI of average predicted wage for south = 0 and 1
margins, dydx(south)
boottest, margins graphopt(xtitle(Average effect of south)) nogr // bootstrap CI of average impact in sample of changing south from 0 to 1

qui ivregress 2sls wage ttl_exp collgrad (tenure = union), cluster(industry)
boottest tenure, ptype(equaltail) seed(987654321) nogr  // Wald test, wild restricted efficient bootstrap, Rademacher weights, null imposed, 999 reps
boottest tenure, ptype(equaltail) seed(987654321) stat(c) nogr // same but bootstrap-c
boottest tenure, ptype(equaltail) seed(987654321) stat(c) gridmin(-2) gridmax(2) nogr // same but limit graphing range
boottest, ar nogr // same bootstrap, but Anderson-Rubin test (much faster)
scoretest tenure, nogr  // Rao/LM test of same
waldtest tenure, nogr  // Wald test of same

qui ivregress liml wage (tenure = collgrad ttl_exp), cluster(industry)
boottest tenure, noci // WRE bootstrap, Rademacher weights, 999 replications
qui cmp (wage = tenure) (tenure = collgrad ttl_exp), ind(1 1) qui nolr cluster(industry)
boottest tenure // reasonable match on test statistic and p value

qui ivreg2 wage collgrad smsa race age (tenure = union married), cluster(industry) fuller(1)
boottest tenure, nograph // Wald test, WRE bootstrap, Rademacher weights, 999 replications
boottest, nograph ar // same, but Anderson-Rubin (faster, but CI misleading if instruments invalid)

qui ivregress liml wage (collgrad tenure = ttl_exp union), cluster(industry)
boottest, ar nogr // Anderson-Rubin test, with contour plot of p value surface
boottest collgrad tenure, gridpoints(10 10) nogr // WRE boostrap also with contour plot

qui regress wage tenure ttl_exp collgrad, robust  // no clustering
boottest tenure, nogr

qui ivregress liml wage (collgrad tenure = ttl_exp union), robust  // no clustering
boottest, ar nogr
boottest collgrad tenure, gridpoints(10 10) nogr

qui regress wage ttl_exp collgrad tenure, cluster(industry)
waldtest collgrad tenure, cluster(industry age) nogr // multi-way-clustered tests after estimation command not offering such
boottest tenure, cluster(industry age) bootcluster(industry) gridmin(-.2) gridmax(.2) nogr

qui areg wage ttl_exp collgrad tenure [aw=hours] if occupation<., cluster(age) absorb(industry)
boottest tenure, cluster(age occupation) bootcluster(occupation) seed(999) nograph  // override estimate's clustering
qui reg wage ttl_exp collgrad tenure i.industry [aw=hours] if occupation<., cluster(age)
boottest tenure, cluster(age occupation) bootcluster(occupation) seed(999) nograph  // should match previous result

qui probit c_city tenure wage ttl_exp collgrad, cluster(industry)
boottest tenure, nogr                          // score bootstrap, Rademacher weights, null imposed, 999 replications
boottest tenure, cluster(industry age) bootcluster(industry) small nogr  // multi-way-clustered, finite-sample-corrected test with score bootstrap

qui gsem (c_city <- tenure wage ttl_exp collgrad), vce(cluster industry) probit  // same probit estimate as previous
boottest tenure                                                             // requires Stata 14.0 or later
boottest tenure, cluster(industry age) bootcluster(industry) small          // requires Stata 14.0 or later

sysuse auto, clear
ml model lf myprobit (foreign = mpg weight)  // define model
qui ml max // estimate
boottest mpg, cmdline(ml model lf myprobit (foreign = mpg weight)) 

use pixel-level-baseline-final, clear
global pix lnkm pixpetro pixdia pixwaterd pixcapdist pixmal pixsead pixsuit pixelev pixbdist
global geo lnwaterkm lnkm2split mean_elev mean_suit malariasuit petroleum diamondd
global poly capdistance1 seadist1 borderdist1
encode pixwbcode, gen(ccode)  // make numerical country identifier

qui areg lnl0708s centr_tribe lnpd0 $pix $geo $poly, absorb(ccode)
waldtest centr_tribe, nogr cluster(ccode pixcluster)
set seed 2309487  // for exact reproducibility of results
boottest centr_tribe, nogr reps(999) clust(ccode pixcluster) bootcluster(ccode)
boottest centr_tribe, nogr reps(999) clust(ccode pixcluster) bootcluster(pixcluster)
boottest centr_tribe, nogr reps(999) clust(ccode pixcluster) bootcluster(ccode pixcluster)

infile coll merit male black asian year state chst using regm.raw, clear
set seed 3541641
qui regress coll merit male black asian i.year i.state, cluster(state)	
generate individual = _n  // unique ID for each observation
boottest merit, nogr reps(999) gridpoints(10)  // defaults to bootcluster(state)
boottest merit, nogr reps(999) gridpoints(10) nonull
boottest merit, nogr reps(999) gridpoints(10) bootcluster(state year)
boottest merit, nogr reps(999) gridpoints(10) nonull bootcluster(state year)
boottest merit, nogr reps(999) gridpoints(10) bootcluster(individual)
boottest merit, nogr reps(999) gridpoints(10) nonull bootcluster(individual)

qui regress coll merit male black asian i.year i.state if !inlist(state,34,57,59,61,64,71,72,85,88), cluster(state)
foreach bootcluster in state "state year" individual {
  foreach nonull in "" nonull {
    boottest merit, reps(999) gridpoints(10) bootcl(`bootcluster') `nonull' nogr
  }
}

use Levitt, clear
set seed 8723419
foreach crimevar in Violent Property {
  qui ivregress 2sls D.l`crimevar'pop (DL.lpris_totpop = ibnL.stage#i(1/3)L.substage)  ///
    D.(lincomepop unemp lpolicepop metrop black a*pop) i.year i.state, robust
  boottest DL.lpris_totpop, cluster(state year) bootcluster(year) ptype(equaltail) reps(999) gridmin(-2) gridmax(2) nogr
}

qui ivregress 2sls D.lPropertypop (DL.lpris_totpop = ibnL.stage#i(1/3)L.substage) D.(lincomepop unemp lpolicepop metrop black a*pop) i.year i.state, robust small
boottest DL.lpris_totpop=-1/3, cluster(state year) bootcluster(year) ptype(equaltail) reps(999) seed(1231) nogr gridmin(-2) gridmax(2)
qui ivreghdfe D.lPropertypop (DL.lpris_totpop = ibnL.stage#i(1/3)L.substage) D.(lincomepop unemp lpolicepop metrop black a*pop) i.year, absorb(state) small
boottest DL.lpris_totpop=-1/3, cluster(state year) bootcluster(year) ptype(equaltail) reps(999) seed(1231) nogr gridmin(-2) gridmax(2)
qui xtivreg   D.lPropertypop (DL.lpris_totpop = ibnL.stage#i(1/3)L.substage) D.(lincomepop unemp lpolicepop metrop black a*pop) i.year, fe small
boottest DL.lpris_totpop=-1/3, cluster(state year) bootcluster(year) ptype(equaltail) reps(999) seed(1231) nogr gridmin(-2) gridmax(2)

qui ivreg2 D.lPropertypop (DL.lpris_totpop = ibnL.stage#i(1/3)L.substage) i.year i.state, robust
boottest DL.lpris_totpop=-1/3, cluster(state year) bootcluster(year) ptype(equaltail) reps(199) noci seed(1231)
qui ivreg2 D.lPropertypop (DL.lpris_totpop = ibnL.stage#i(1/3)L.substage) i.year i.state, robust partial(i.year i.state)
boottest DL.lpris_totpop=-1/3, cluster(state year) bootcluster(year) ptype(equaltail) reps(199) noci seed(1231)

qui log close