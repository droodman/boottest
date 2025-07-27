cap cd "D:\OneDrive\Documents\Work\Econometrics\Wild cluster"
cap cd "/mnt/d/OneDrive/Documents/Work/Econometrics/Wild cluster"
cap cd "/Users/davidroodman/Library/CloudStorage/OneDrive-Personal/Documents/Macros/boottest"

set more off
set rmsg off
set trace off
set linesize 200
cap set processors 4

cap program drop myprobit
program myprobit  // custom likelihood evaluator
  args lnf theta
  quietly replace `lnf' = lnnormal((2*$ML_y1-1) * `theta')
end

cap log close
qui log using "unit tests.log", replace

version 13

set seed 0193284710

foreach julia in "" julia {
  use test/collapsed, clear

  qui regress hasinsurance selfemployed post post_self, cluster(year)
  boottest post_self=.04, `julia' nogr
  boottest post_self=.04, `julia' weight(webb) noci
  boottest post_self=.04, `julia' weight(webb) jk nogr
  boottest post_self=.04, `julia' weight(webb) jk nogr nonull
  scoretest post_self=.04, `julia' nogr

  boottest post_self post, `julia' reps(999) weight(webb) nogr  // wild bootstrap test of joint null, Webb weights, null imposed, 9,999 replications
  boottest (post_self) (post), `julia' reps(999) weight(webb) nogr  // same
  boottest {post_self=.04} {post}, `julia' nogr  // separate tests, no correction for multiple hypotheses
  boottest {(post) (post_self=.04)} {(post) (post_self=.08)}, `julia' madj(sidak) nogr  // separate tests, Sidak correction for  multiple hypotheses

  use test/nlsw88

  qui regress wage tenure ttl_exp collgrad, cluster(industry)
  boottest tenure, `julia' svmat nogr  // wild bootstrap test of joint null, Rademacher weights, null imposed, saving  simulated distribution

  constraint 1 ttl_exp = .2
  qui cnsreg wage tenure ttl_exp collgrad, constr(1) cluster(industry)
  boottest tenure, `julia' nogr  // wild bootstrap test of tenure=0, conditional on ttl_exp=2, Rademacher weights, null imposed, 999 replications

  regress wage tenure ttl_exp collgrad south#union, cluster(industry)
  margins south
  boottest, `julia' margins nogr  // bootstrap CI of average predicted wage for south = 0 and 1
  margins, dydx(south)
  boottest, `julia' margins graphopt(xtitle(Average effect of south)) nogr  // bootstrap CI of average impact in sample of changing south from 0 to 1

  local julia julia
  qui ivregress 2sls wage ttl_exp collgrad (tenure = union), cluster(industry)
  boottest tenure, `julia' ptype(equaltail) seed(987654321) nogr  // Wald test, wild restricted efficient bootstrap, Rademacher weights, null imposed, 999 reps
  boottest tenure, `julia' ptype(equaltail) stat(c) nogr  // same but bootstrap-c
  boottest tenure, `julia' ptype(equaltail) stat(c) gridmin(-2) gridmax(2) nogr // same but limit graphing range
  boottest, `julia' ar nogr  // same bootstrap, but Anderson-Rubin test (much faster)
  scoretest tenure, `julia' nogr  // Rao/LM test of same
  waldtest tenure, `julia' nogr  // Wald test of same

  qui ivregress liml wage (tenure = collgrad ttl_exp), cluster(industry)
  boottest tenure, `julia' noci  // WRE bootstrap, Rademacher weights, 999 replications
  boottest tenure, `julia' noci jk // WRE bootstrap, Rademacher weights, 999 replications
  qui cmp (wage = tenure) (tenure = collgrad ttl_exp), ind(1 1) qui nolr cluster(industry)
  boottest tenure, `julia'  // reasonable match on test statistic and p value

  qui ivreg2 wage collgrad smsa race age (tenure = union married), cluster(industry) fuller(1)
  boottest tenure, `julia' nograph  // Wald test, WRE bootstrap, Rademacher weights, 999 replications
  boottest tenure, `julia' nograph jk // Wald test, WRE bootstrap, Rademacher weights, 999 replications
  boottest, `julia' nograph ar  // same, but Anderson-Rubin (faster, but CI misleading if instruments invalid)

  qui ivregress liml wage (collgrad tenure = ttl_exp union), cluster(industry)
  boottest, `julia' ar nogr  // Anderson-Rubin test, with contour plot of p value surface
  boottest, `julia' ar nogr jk  // Anderson-Rubin test, with contour plot of p value surface
  boottest collgrad tenure, `julia' gridpoints(10 10) nogr  // WRE boostrap also with contour plot

  qui regress wage tenure ttl_exp collgrad, robust  // no clustering
  boottest tenure, `julia' nogr
  boottest tenure, `julia' nogr jk

  qui ivregress liml wage (collgrad tenure = ttl_exp union), robust  // no clustering
  boottest, `julia' ar nogr
  boottest collgrad tenure, `julia' gridpoints(10 10) nogr

  qui regress wage ttl_exp collgrad tenure, cluster(industry)
  waldtest collgrad tenure, cluster(industry age) nogr // multi-way-clustered tests after estimation command not offering such
  boottest tenure, `julia' cluster(industry age) bootcluster(industry) gridmin(-.2) gridmax(.2) nogr

  qui areg wage ttl_exp collgrad tenure [aw=hours] if occupation<., cluster(age) absorb(industry)
  boottest tenure, `julia' cluster(age occupation) bootcluster(occupation) seed(999) nograph  // override estimate's clustering
  boottest tenure, `julia' cluster(age occupation) bootcluster(occupation) seed(999) nograph jk // override estimate's clustering
  qui reg wage ttl_exp collgrad tenure i.industry [aw=hours] if occupation<., cluster(age)
  boottest tenure, `julia' cluster(age occupation) bootcluster(occupation) seed(999) nograph  // should match previous result
  boottest tenure, `julia' cluster(age occupation) bootcluster(occupation) seed(999) nograph jk // should match previous result

  qui probit c_city tenure wage ttl_exp collgrad, cluster(industry)
  boottest tenure, `julia' nogr                          // score bootstrap, Rademacher weights, null imposed, 999 replications
  boottest tenure, `julia' cluster(industry age) bootcluster(industry) small nogr  // multi-way-clustered, finite-sample-corrected test with score bootstrap

  qui gsem (c_city <- tenure wage ttl_exp collgrad), vce(cluster industry) probit  // same probit estimate as previous
  boottest tenure, `julia'                                                              // requires Stata 14.0 or later
  boottest tenure, `julia' cluster(industry age) bootcluster(industry) small          // requires Stata 14.0 or later

  sysuse auto, clear
  ml model lf myprobit (foreign = mpg weight)  // define model
  qui ml max // estimate
  boottest mpg, `julia' cmdline(ml model lf myprobit (foreign = mpg weight)) 

  probit foreign i.mpg
  scoretest 14.mpg

  use test/collapsed, clear

  qui regress hasinsurance selfemployed post post_self, cluster(year)
  boottest post_self=.04, `julia' weight(webb) nogr
  boottest post_self=.04, `julia' weight(webb) reps(9999999) noci
  boottest post_self=.04, `julia' weight(normal) reps(9999) noci
  boottest post_self=.04, `julia' weight(gamma) reps(9999) noci svv
  boottest post_self=.04, `julia' weight(mammen) reps(9999) noci
  boottest post_self=.04, `julia' weight(mammen) reps(9999) boottype(score) nogr

  qui regress hasinsurance selfemployed post post_self, robust
  boottest post_self=.04, `julia' weight(webb) nogr

  qui regress hasinsurance selfemployed post post_self, cluster(year)
  boottest (post_self=.05) (post=-.02), `julia' reps(9999) weight(webb) nogr
  boottest (post_self=.05) (post=-.02) (selfemployed=-.15), `julia' reps(9999) weight(webb) nogr

  qui regress hasinsurance selfemployed post post_self
  boottest post_self=.04, `julia' weight(webb) nogr
  boottest (post_self=.05) (post=-.02), `julia' reps(9999) weight(webb) nogr
  scoretest (post_self=.05), `julia' nogr
  scoretest (post_self=.05) (post=-.02), `julia' nogr
  boottest (post_self=.08), `julia' boottype(score) reps(9999) nogr
  boottest (post_self=.05) (post=-.02), `julia' boottype(score) reps(9999) nogr

  use test/nlsw88, clear
  constraint 1 ttl_exp = .2
  qui cnsreg wage tenure ttl_exp collgrad, constr(1) cluster(industry)
  boottest tenure, `julia' nogr

  keep if e(sample)
  gen id = _n - cond(_n>1000, 1000, 0)
  qui cnsreg wage tenure ttl_exp collgrad, constr(1) cluster(id)  // granular but not pure robust
  boottest tenure, `julia' reps(9999) nogr

  qui areg wage tenure ttl_exp collgrad, cluster(id) a(industry)
  boottest tenure, `julia' reps(9999) nogr

  use test/nlsw88, clear
  qui ivregress liml wage ttl_exp collgrad (tenure = union), cluster(industry)
  boottest tenure, `julia' ptype(equaltail) reps(9999) nogr
  boottest tenure, `julia' nonull reps(99999) matsize(.1) nogr
  boottest tenure, `julia' ptype(upper) svmat(t) reps(9999) nogr
  boottest tenure, `julia' ptype(lower) svmat(numer) reps(9999) nogr

  qui ivregress liml wage ttl_exp collgrad (tenure = union), cluster(industry)
  boottest tenure, `julia' ptype(equaltail) reps(9999) nogr

  qui ivregress liml wage ttl_exp collgrad (tenure = union) if industry<., robust
  boottest tenure, `julia' ptype(equaltail) reps(99) noci

  qui ivregress liml wage ttl_exp collgrad (tenure = union) if industry<., robust
  boottest tenure, ptype(equaltail) reps(99) noci
  boottest collgrad tenure, `julia' ptype(equaltail) reps(99) noci

  qui ivregress 2sls wage ttl_exp collgrad (tenure = union) if industry<.
  boottest tenure, `julia' ptype(equaltail) reps(99) noci
  boottest tenure collgrad, `julia' ptype(equaltail) reps(99) noci

  qui ivregress 2sls wage ttl_exp collgrad (tenure = union), cluster(industry)
  boottest tenure, `julia' ptype(equaltail) weight(webb) stat(c) gridmin(-5) gridmax(5) gridpoints(100) nogr
  boottest tenure, `julia' ptype(equaltail) weight(webb) stat(c) gridmin(-5) gridmax(5) gridpoints(100) matsize(.01) nogr

  qui ivregress 2sls wage ttl_exp collgrad (tenure = union) if industry<., robust
  boottest tenure, `julia' ptype(equaltail) matsize(.005) noci weight(webb)

  preserve
  keep if e(sample)
  gen id = _n - cond(_n>1000, 1000, 0)
  boottest tenure, `julia' cluster(id) ptype(equaltail) matsize(.005) noci weight(webb)
  restore

  qui ivregress 2sls wage ttl_exp collgrad (tenure = union), cluster(industry)
  boottest, `julia' ar nogr
  boottest, `julia' ar nonull nogr
  scoretest tenure, `julia' nogr
  waldtest tenure, `julia' ptype(upper) nogr

  qui ivregress liml wage (tenure = collgrad ttl_exp), cluster(industry)
  boottest tenure, `julia' nogr

  qui ivreg2 wage collgrad smsa race age (tenure = union married), cluster(industry) fuller(1)
  boottest tenure, `julia' nograph weight(webb) reps(9999)
  qui gen individual = _n
  boottest tenure, `julia' noci bootcluster(individual) weight(webb)
  boottest tenure, `julia' nograph bootcluster(collgrad) cluster(collgrad industry) weight(webb) reps(9999)

  qui areg wage ttl_exp collgrad tenure [aw=hours] if occupation<. & grade<. & union<., cluster(age) absorb(industry)
  boottest tenure, `julia' nograph cluster(age occupation) bootcluster(occupation)

  qui areg wage ttl_exp collgrad tenure if occupation<. & grade<. & union<. & hours<., robust absorb(industry)
  boottest tenure, `julia' nograph

  qui areg wage ttl_exp collgrad tenure [aw=hours] if occupation<. & grade<. & union<., robust absorb(industry)
  boottest tenure, `julia' nograph 

// ivreghdfe crashing
//   qui ivreghdfe wage ttl_exp collgrad tenure (occupation = union married) [aw=hours] if grade<., liml cluster(industry) absorb(industry)
//   boottest tenure, `julia' nograph
//   boottest occupation, `julia' nograph
//
//   qui ivreghdfe wage ttl_exp collgrad tenure (occupation = union married) [aw=hours] if grade<., liml cluster(industry) absorb(age)
//   boottest tenure, `julia' nograph
//   boottest collgrad tenure, `julia' nograph
//   boottest occupation, `julia' gridmin(-1) gridmax(1) nograph
//
//   constraint 1 [wage]collgrad
//   qui ivreghdfe wage ttl_exp /*collgrad*/ tenure (occupation = union married) [aw=hours], liml cluster(industry) absorb(age)  // approximate contrained LIML with reghdfejl
//   boottest tenure, `julia' nograph

  use test/abdata, clear
  qui areg n w k, absorb(ind)
  boottest k, `julia' cluster(id year) nograph
  qui areg n w k [aw=ys], absorb(ind)
  boottest k, `julia' cluster(id year) nograph

  use test/pixel-level-baseline-final, clear
  global pix lnkm pixpetro pixdia pixwaterd pixcapdist pixmal pixsead pixsuit pixelev pixbdist
  global geo lnwaterkm lnkm2split mean_elev mean_suit malariasuit petroleum diamondd
  global poly capdistance1 seadist1 borderdist1
  qui encode pixwbcode, gen(ccode)  // make numerical country identifier
  qui areg lnl0708s centr_tribe lnpd0 \$pix \$geo \$poly, absorb(ccode)
  boottest centr_tribe, `julia' nogr reps(999) clust(ccode pixcluster) bootcluster(ccode)
  boottest centr_tribe, `julia' nogr reps(999) clust(ccode pixcluster) bootcluster(pixcluster)
  boottest centr_tribe, `julia' nogr reps(999) clust(ccode pixcluster) bootcluster(ccode pixcluster)

  infile coll merit male black asian year state chst using test/regm.raw, clear
  qui regress coll merit male black asian i.year i.state if !inlist(state,34,57,59,61,64,71,72,85,88), cluster(state)	
  generate individual = _n  // unique ID for each observation
  boottest merit, `julia' nogr reps(999) gridpoints(10)  // defaults to bootcluster(state)
  boottest merit, `julia' nogr reps(999) gridpoints(10) nonull
  boottest merit, `julia' nogr reps(999) gridpoints(10) bootcluster(state year)
  boottest merit, `julia' nogr reps(999) gridpoints(10) nonull bootcluster(state year)
  boottest merit, `julia' nogr reps(999) gridpoints(10) bootcluster(individual)
  boottest merit, `julia' nogr reps(999) gridpoints(10) nonull bootcluster(individual)
  boottest merit, `julia' nogr reps(999) gridpoints(10) nonull bootcluster(individual) matsize(.1)
}
qui log close
