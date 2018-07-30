*! boottest 2.1.6 29 July 2018
*! Copyright (C) 2015-18 David Roodman

* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.

*! Version history at bottom

cap program drop boottest
program define boottest
	version 11
	
	cap noi _boottest `0'
	local rc = _rc
	constraint drop `anythingh0'
	cap mata mata drop _boottestp
	cap mata mata drop _boottestC
	exit `rc'
end

cap program drop _boottest
program define _boottest, rclass sortpreserve
	version 11

	mata st_local("StataVersion", boottestStataVersion()); st_local("CodeVersion", boottestVersion())
	if `StataVersion' != c(stata_version) | "`CodeVersion'" < "02.01.00" {
		cap findfile "lboottest.mlib"
		while !_rc {
			erase "`r(fn)'"
			cap findfile "lboottest.mlib"
		}
		qui findfile "boottest.mata"
		run "`r(fn)'"
	}

	local   cmd = cond(substr("`e(cmd)'", 1, 6)=="ivreg2", "ivreg2", "`e(cmd)'")
	local ivcmd = cond("`cmd'"=="reghdfe", "`e(subcmd)'", cond("`cmd'"=="xtivreg2", "ivreg2", "`cmd'"))
	
	if "`e(prefix)'" == "svy" {
		di as err "Doesn't work after {cmd:svy}."
		exit 198
	}
	if inlist("`cmd'", "mvreg", "sureg") {
		di as err "Doesn't work after {cmd:`e(cmd)'}."
		exit 198
	}
	if inlist("`cmd'", "xtreg", "xtivreg") & "`e(model)'"!="fe" {
		di as err "Doesn't work after {`cmd', `e(model)'}."
		exit 198
	}
	if "`cmd'"=="xtivreg2" & "`e(xtmodel)'"!="fe" {
		di as err "Doesn't work after {`cmd', `e(xtmodel)'}."
		exit 198
	}
	if "`cmd'"=="reghdfe" & `:word count `e(absvars)''>1 {
		di as err "Doesn't work after {cmd:reghdfe) with more than one set of fixed effects."
		exit 198
	}
	if inlist("`cmd'", "sem", "gsem") & c(stata_version) < 14 {
		di as err "Requires Stata version 14.0 or later to work with {cmd:`e(cmd)'}."
		exit 198
	}
	if ("`ivcmd'"=="ivreg2") & e(partial_ct)>0 & e(partial_ct)<. {
		di as err "Doesn't work after {cmd:ivreg2} with the {cmd:partial()} option."
		exit 198
	}

	tokenize `"`0'"', parse(",")
	if `"`1'"' != "," {
		local h0s `h0s' `1'
		macro shift
	}
	local 0 `*'
	syntax, [h0(numlist integer >0) Reps(integer 999) seed(string) BOOTtype(string) CLuster(string) Robust BOOTCLuster(string) noNULl QUIetly WEIGHTtype(string) Ptype(string) NOCI Level(real `c(level)') SMall SVMat ///
						noGRaph gridmin(string) gridmax(string) gridpoints(string) graphname(string asis) graphopt(string asis) ar MADJust(string) CMDline(string) MATSIZEgb(string) *]

	if `"`matsizegb'"' != "" & `"`matsizegb'"' != "." {
		cap confirm numeric `matsizegb'
		if _rc {
			di as err "{cmdab:maxmat:size()} option must be a positive number (of gigabytes)."
			exit 198
		}
		else if `matsizegb' <= 0 {
			di as err "{cmdab:maxmat:size()} option must be a positive number (of gigabytes)."
			exit 198
		}
	}
	else local matsizegb .

	if `reps'==0 local svmat
		else {
			local _svmat = cond("`svmat'"=="", "", "t")
			local 0, `options'
			syntax, [SVMat(string) *]
			if !inlist(`"`svmat'"', "numer", "t", "") {
				di as err "option " as inp "svmat(" `svmat' ")" as err " not allowed."
				error 198
			}
			if "`svmat'" == "" local svmat `_svmat'
		}

	if `"`options'"' != "" {
		if `"`options'"'=="ci" di as err "Option {cmd:ci} is obsolete, because it is now the default."
			else di as err `"Option `options' not allowed."'
		exit 198
	}

	if `reps' < 0 {
		di as err "{cmdab:r:eps()} option must be non-negative."
		exit 198
	}

	local 0, `madjust'
	syntax, [Bonferroni Sidak]
	local madjust `bonferroni'`sidak'

	if `"`graphname'"' != "" {
		local 0 `graphname'
		syntax [anything(name=graphname)], [replace]
	}
	else local graphname Graph

	if `"`gridpoints'"'=="" local gridpoints .
	else {
		cap confirm integer number `gridpoints'
		local rc = _rc
		if !_rc {
			local rc = `gridpoints' <= 0
		}
		if `rc' {
			di as err "{cmd:gridpoints()} not a positive integer."
			exit 198
		}
	}
	if `"`gridmin'"'=="" local gridmin .
	else {
		cap confirm number `gridmin'
		if _rc {
			di as err "{cmd:gridmin()} not numeric."
			exit 198
		}
	}
	if `"`gridmax'"'=="" local gridmax .
	else {
		cap confirm number `gridmax'
		if _rc {
			di as err "{cmd:gridmax()} not numeric."
			exit 198
		}
	}
	
	if `"`robust'`cluster'"' != "" {
		local hasrobust 1
		local clustvars `cluster'
		if `"`clustvars'"'!="" confirm var `clustvars'
		local override 1
		di as txt _n "Overriding estimator's cluster/robust settings with " as res `"`=cond("`clustvars'"=="", "robust", "cluster(`clustvars')")'"'
	}
	else {
		local hasrobust = "`e(vcetype)'" == "Robust"
		local clustvars `e(clustvar)'
	}

	local hasclust = `"`clustvars'"' != ""

	local wtype `e(wtype)'
	local null = "`null'" == ""
	if "`svmat'"!="" tempname dist
	local ar = "`ar'" != ""
	if `ar' & `"`h0s'`h0'"' == "" local h0s `e(instd)'

	if `:word count `weighttype'' > 1 {
		di as err "The {cmd:weight:type} option must be {cmdab:rad:emacher}, {cmdab:mam:men}, {cmdab:nor:mal}, or {cmdab:web:b}."
		exit 198
	}
	if "`weighttype'"'=="" local weighttype rademacher
	else {
		local 0, `weighttype'
		syntax, [RADemacher MAMmen NORmal WEBb GAMma]
		local weighttype `rademacher'`mammen'`normal'`webb'`gamma'
	}

	if `:word count `ptype'' > 1 {
		di as err "The {cmd:wp:type} option must be {cmdab:sym:metric}, {cmdab:eq:qualtail}, {cmd:lower}, or {cmd:upper}."
		exit 198
	}
	if "`ptype'"'=="" local ptype symmetric
	else {
		local 0, `ptype'
		syntax, [SYMmetric EQualtail LOWer UPper]
		local ptype `symmetric'`equaltail'`lower'`upper'
	}

	local ML = e(converged) != .
	local IV = "`e(instd)'`e(endogvars)'" != ""
	local LIML = ("`cmd'"=="ivreg2" & "`e(model)'"=="liml") | ("`cmd'"=="ivregress" & "`e(estimator)'"=="liml")| ("`cmd'"=="reghdfe" & strpos("`e(title)'", "LIML"))
	local WRE = `"`boottype'"'!="score" & `IV' & `reps'
	local small = e(df_r) != . | "`small'" != ""

	local fuller `e(fuller)' // "" if missing
	local K = e(kclass) // "." if missing

	local tz = cond(`small', "t", "z")
	local symmetric Prob>|`tz'|
	local equaltail 2 * min(Prob>|`tz'|, Prob<-|`tz'|)
	local lower     Prob<`tz'
	local upper     Prob>`tz'

	if `ar' & !`IV' {
		di as err "Anderson-Rubin test is only for IV models."
		exit 198
	}

	if "`boottype'"'=="" local boottype = cond(`ML', "score", "wild")
	else {
		local 0, `boottype'
		syntax, [Wild SCore]
		local boottype `score'`wild'
		if "`boottype'" == "wild" & `ML' {
			di as err "{cmd:boottype(wild)} not accepted after GMM or Maximum Likelihood-based estimation."
			exit 198
		}
		if "`boottype'"=="score" & "`fuller'`K'" != "." {
			di as err "{cmd:boottype(score)} not accepted after Fuller LIML or k-class estimation."
			exit 198
		}
	}
	local scoreBS = "`boottype'"=="score"
	
	local FEname = cond(inlist("`cmd'","xtreg","xtivreg","xtivreg2"), "`e(ivar)'", "`e(absvar)'`e(absvars)'")

	if `"`seed'"'!="" set seed `seed'

	tempname p padj se stat df df_r hold C C0 CC0 b V keepC keepW repsname repsFeasname
	mat `b' = e(b)
	local k = colsof(`b') // column count before possible removal of _cons term in FE models

	if "`e(wtype)'" != "" {
		tempname wtname
		qui gen double `wtname' `e(wexp)'
	}

	if `"`h0s'"' != "" {
		if "`h0'" != "" {
			di as err "Specify hypotheses before comma or with {cmd:h0()} option, but not both."
			exit 198
		}
		local multiple = strpos(`"`h0s'"', "{")
		if !`multiple' local h0s {`h0s'}
		local N_h0s 0 // number of nulls
		while `"`h0s'"' != "" {
			gettoken h0 h0s: h0s, parse("{}")
			if !inlist(`"`h0'"', "{", "}") {
				local ++N_h0s
				while `"`h0'"' != "" {
					gettoken cns h0: h0, parse("()") match(m)
					constraint free
					constraint `r(free)' `cns'
					local h0_`N_h0s' `h0_`N_h0s'' `r(free)'
					local constraints  `constraints' `r(free)'
					c_local anythingh0 `constraints'
				}
			}
		}
		local anythingh0 1
	}
	else {
		local N_h0s 1 // number of nulls
		if "`h0'" == "" {
			di as txt _n "({cmd:h0(1)} assumed)"
			local h0 1
		}
		foreach c of numlist `h0' {
			if `"`: constraint `c''"' == "" di as res "Constraint `c' not found and will be skipped."
		}
		local h0_1 `h0'
	}

	if `N_h0s'==1 local madjust

	makecns
	if "`e(Cns)'" != "" {
		local hascns 1
		mat `C' = e(Cns)
	}
	cap _estimates drop `hold'
	if !`ML' {
		if ("`ivcmd'"=="ivreg2" & !inlist("`e(model)'", "ols", "iv", "gmm2s", "liml")) | ("`ivcmd'"=="ivregress" & !inlist("`e(estimator)'", "liml", "2sls", "gmm")) {
			di as err "Only for use with OLS, 2SLS, single-equation GMM, and ML-based estimators."
			exit 198
		}

		local Ynames `e(depvar)'

		local colnames: colnames e(b)
		local _cons _cons
		local colnames: list colnames - _cons
		if `IV' {
			local Xnames_endog `e(instd)'
			local Xnames_exog: list colnames - Xnames_endog
			local cons = inlist("`cmd'`e(cons)'`e(constant)'", "ivreg21", "ivregress")
			local ZExclnames `e(insts)'
			local ZExclnames `:list ZExclnames - Xnames_exog'
		}
		else {
			local Xnames_exog: colnames e(b)
			local _cons _cons
			local cons = "`:list Xnames_exog & _cons'"!=""
			local Xnames_exog: list Xnames_exog - _cons
		}
		if !`cons' local _cons

		local GMM = ("`ivcmd'"=="ivreg2" & "`e(model)'"=="gmm2s") | ("`cmd'"=="ivregress" & "`e(estimator)'"=="gmm")
		if `GMM' {
			if `reps' di as txt _n "Bootstrapping purely with final GMM moment weighting matrix."
			tempname W
			mat `W' = e(W)
			local colnamesW: colnames `W'
			local colnamesW: list colnamesW - _cons
			mata _boottestp = J(`cons',1,`k') \ order(tokens("`colnamesW'")', 1)[invorder(order(tokens("`Xnames_exog' `ZExclnames'")', 1))]
			mata st_matrix("`W'", st_matrix("`W'")[_boottestp,_boottestp]) // ensure weight matrix in order cons-other included exog-excluded exog
		}

		mata _boottestp = J(`cons',1,`k') \ order(tokens("`colnames'")', 1)[invorder(order(tokens("`Xnames_exog' `Xnames_endog'")', 1))] \ `k'+1

		if "`FEname'"!="" & `cons' {
			mata _boottestp = _boottestp[|2\.|]
			local cons 0
			local _cons
		}
		
		if 0`hascns' mata st_matrix("`C'" , st_matrix("`C'" )[,_boottestp])

		if `cons' {
			mat `keepC' = 1
			mat `keepW' = 1
		}
		local colnamesC `_cons' `Xnames_exog' `Xnames_endog'
		local colnamesW `_cons' `Xnames_exog' `ZExclnames'
		foreach varlist in Ynames Xnames_exog Xnames_endog ZExclnames {
			fvrevar ``varlist'' if e(sample)
			local revarlist `r(varlist)'
			local _revarlist
			forvalues i=1/`:word count `revarlist'' {
				local var: word `i' of ``varlist''
				_ms_parse_parts `var'
				if !r(omit) {
					local _revarlist `_revarlist' `:word `i' of `revarlist''
					if inlist("`varlist'", "Xnames_exog", "Xnames_endog") {
						local pos: list posof "`var'" in colnamesC
						if `pos' mat `keepC' = nullmat(`keepC'), `pos'
					}
					if inlist("`varlist'", "Xnames_exog", "ZExclnames") & `GMM' {
						local pos: list posof "`var'" in colnamesW
						if `pos' mat `keepW' = nullmat(`keepW'), `pos'
					}
				}
			}
			local `varlist' `_revarlist'
		}

		if 0`hascns' mata _boottestC = st_matrix("`C'")[,(st_matrix("`keepC'"),cols(st_matrix("`C'")))]; _boottestC = select(_boottestC,rowsum(_boottestC:!=0)); st_matrix("`C'" , _boottestC)

		if `GMM' mata st_matrix("`W'", st_matrix("`W'" )[st_matrix("`keepW'"), st_matrix("`keepW'")])
		if `cons' local Xnames_exog `hold' `Xnames_exog' // add constant term
	}

	if `hasclust' {
		if 0`override' {
			cap assert !missing(`:subinstr local clustvars " " ",", all') if e(sample), `=cond(c(stata_version)>=12.1, "fast", "")'
			if _rc {
				di as err "A clustering variable has a missing value for at least one observation."
				exit 9
			}
		}

		if `"`bootcluster'"' == "" {
			local bootcluster `clustvars'
			if `reps' & `:word count `clustvars''>1 di as txt "({cmdab:bootcl:uster(`clustvars')} assumed)"
		}
		local    clustvars `:list clustvars & bootcluster' `:list clustvars - bootcluster'
		local allclustvars `:list bootcluster - clustvars' `clustvars' // put bootstrapping clusters first, error clusters last, overlap in middle

		sort `clustvars' `:list bootcluster - clustvars', stable

		foreach clustvar in `allclustvars' {
			cap confirm numeric var `clustvar'
			if _rc {
				tempvar IDname
				qui egen long `IDname' = group(`clustvar') if e(sample)
				local _allclustvars `_allclustvars' `IDname'
			}
			else local _allclustvars `_allclustvars' `clustvar'
		}
		local allclustvars `_allclustvars'
	}

	forvalues h=1/`N_h0s' { // loop over multiple independent constraints
		_estimates hold `hold', restore
		ereturn post `b'
		mat `b' = e(b)
		// process hypothesis constraints into e(Cns)
		qui makecns `h0_`h''
		if "`h0_`h''" != "`r(clist)'" {
			local clist `r(clist)'
			local clist: list h0_`h' - clist
			foreach c in `clist' {
				di as txt `"(note: constraint `=cond(0`anythingh0', `"{res}`:constraint `c''{txt}"', "number `c'")' caused error {search r(111)})"'
			}
			exit 111
		}
		
		cap mat `C0' = e(Cns)
		if _rc exit 111
		_estimates unhold `hold'
		if r(k_autoCns) mat `C0' = `C0'[r(k_autoCns)+1...,1...]
		
		if "`FEname'" != "" {
			mata st_local("rc", strofreal(any(0:!=select(st_matrix("`C0'"),st_matrixcolstripe("`C0'")[,2]':=="_cons"))))

			if `rc' {
				di as err "In fixed-effect models, null hypotheses may not involve constant term."
				exit 111
			}
		}
		
		local plotmat
		local peakmat
		local cimat
		if "`noci'"=="" & !`ML' & rowsof(`C0')==1 {
			if !(`reps'==0 & `scoreBS' & !`null' & "`graph'"!="") tempname plotmat peakmat // don't request plot if doing simple Wald with no graph
			if `level'<100 tempname cimat
		}

		if rowsof(`C0')>1 & "`ptype'"!="symmetric" di as txt "Note: {cmd:ptype(`ptype')} ignored for multi-constraint null hypotheses."

		if `ML' {
			local K .

			if `null' {
				if "`e(cmdline)'"=="" & `"`cmdline'"'=="" {
					di as err "Original estimation command line needed. Include it in a {cmd:cmdline()} option."
					exit 198
				}
				if "`cmd'"=="tobit" & c(stata_version)<15 {
					di as err "Cannot impose null after {help tobit} in Stata versions below 15. Try re-fitting the model with {stata ssc describe cmp:cmp}, {help intreg}, or {help gsem}."
					exit 198
				}

				mat `CC0' = `C0' \ nullmat(`C') // add null to model constraints

				* get rid of any remaining o. and b. constraints; r(k_autcns) doesn't capture all
				mata _boottestC = st_matrixcolstripe("`CC0'")[,2]
				mata _boottestC = !(strmatch(_boottestC, "*b*.*") :| strmatch(_boottestC, "*o*.*"))'; _boottestC[cols(_boottestC)] = 0 // don't look at r column
				mata st_matrix("`CC0'", select(st_matrix("`CC0'"), rowsum(select(st_matrix("`CC0'"), _boottestC) :!= 0)))
				cap mat `CC0'[1,1] = `CC0'[1,1]
				if _rc {
					di as error _n "Null constraint applies only to omitted variables or base levels of factor variables."
					continue
				}

				`quietly' di as res _n "Re-running regression with null imposed." _n
				if `"`cmdline'"'=="" local 0 `e(cmdline)'
				                else local 0 `cmdline'
				syntax [anything] [aw pw fw iw] [if] [in], [CONSTraints(passthru) from(passthru) INIt(passthru) ITERate(passthru) CLuster(passthru) Robust vce(string) *]
				if "`e(cmd)'"=="ml" local max max // tack on max option in case ml run in interactive model and cmdline() is missing it

				cap _estimates drop `hold'
				_estimates hold `hold', restore
				local coleq: coleq `b'
				if "`coleq''"=="_" local coleq
				local colnames: colnames `b'

				forvalues r=1/`=rowsof(`CC0')' {
					local terms 0
					local _constraint
					forvalues c=1/`=colsof(`CC0')-1' {
						if `CC0'[`r',`c'] {
							if `terms++' local _constraint `_constraint' +
							local _constraint `_constraint' `=cond(`CC0'[`r',`c']!=1,"`=`CC0'[`r',`c']' * ","")' `=cond("`coleq'"=="","","[`:word `c' of `coleq'']")'`:word `c' of `colnames''
						}
					}
					local _constraint `_constraint' = `=`CC0'[`r',`=colsof(`CC0')']'
					constraint free
					local _constraints `_constraints' `r(free)'
					constraint `r(free)' `_constraint'
				}

				cap `=cond("`quietly'"=="", "noisily", "")' `anything' if `hold' `=cond("`weight'"=="","",`"[`weight'`exp']"')', constraints(`_constraints') `from' `init' `iterate' `max' `options'
				local rc = _rc
				constraint drop `_constraints'
				if e(converged)==0 {
					di as err "Could not impose null."
					exit 430
				}
				if !`rc' {
					tempname b0
					mat `b0' = e(b)
					cap `anything' if e(sample), `=cond(inlist("`cmd'","cmp","ml"),"init(`b0')","from(`b0',skip)")' iterate(0) `options' `max'
					local rc = _rc
				}
				if `rc' {
					di as err "Error imposing null. Perhaps {cmd:`cmd'} does not accept the {cmd:constraints()}, {cmd:from()}, and {cmd:iterate()} options, as needed."
					exit `rc'
				}
				if "`quietly'"=="" {
					_ms_omit_info e(b)
					tempname t
					mata st_numscalar("`t'", !all(diagonal(st_matrix("e(V)")) :| st_matrix("r(omit)")'))
					if `t' {
						di _n as res _n "{res}Warning: Negative Hessian under null hypothesis not positive definite. Results may be unreliable." _n
						di "This may indicate that the constrained optimum is far from the unconstrained one, i.e., that the hypothesis is incompatible with the data." _n
					}
				}
			}
			
			mat `V' = e(V`=cond("`e(V_modelbased)'"!="", "_modelbased", "")')

			tempvar sc
			cap predict double `sc'* if e(sample), score
			if _rc {
				di as err "Could not generate scores from `e(cmd)' with {cmd:predict}."
				exit _rc
			}

			unab scnames: `sc'*
			if `:word count `scnames'' == e(k_eq) { // generate score for each parameter from those for each equation
				local colnames: colnames `b'
				tokenize `:coleq `b''
				local lasteq
				local eq 0
				local _sc
				qui forvalues i=1/`k' {
					if "``i''" != "`lasteq'" | "``i''"=="/" {
						local ++eq
						local lasteq ``i''
					}
					local var: word `i' of `colnames'
					if "`var'"=="_cons" | strmatch("`var'","*._cons") | "``i''"=="/" local _sc `_sc' `:word `eq' of `scnames'' // constant term or free parameter
					else {
						tempvar t
						gen double `t' = `:word `eq' of `scnames'' * `var' if e(sample)
						local _sc `_sc' `t'
					}
				}
				local scnames `_sc'
			}

			if !`null' {
				cap _estimates drop `hold'
				_estimates hold `hold', restore // rename estimation as `hold' instead of generating new sample marker
			}
		}
		else { // not ML
			cap _estimates drop `hold'
			_est hold `hold', restore

			mata st_matrix("`C0'", st_matrix("`C0'")[,_boottestp]) // ensure constraint matrix in order cons-other exog-endog
			mata st_matrix("`C0'", st_matrix("`C0'" )[,(st_matrix("`keepC'"),`=colsof(`C0')')]) // drop columns for 0-variance vars

			if `ar' {
				if !`IV' | `GMM' {
					di as err "Anderson-Rubin test only available for non-GMM instrumental variables models."
					exit 198
				}
				local kEnd: word count `Xnames_endog'
				local kEx : word count `Xnames_exog'
				mata st_numscalar("`p'", rows(st_matrix("`C0'"))!=`kEnd' | (`kEx'? any(st_matrix("`C0'")[|.,.\.,`kEx'|]) : 0))
				if `p' {
					di as err "For Anderson-Rubin test, null hypothesis must constrain all the instrumented variables and no others."
					exit 198
				}			
			}
		}

		return local seed = cond("`seed'"!="", "`seed'", "`c(seed)'")

		mata boottest_stata("`stat'", "`df'", "`df_r'", "`p'", "`padj'", "`cimat'", "`plotmat'", "`peakmat'", `level', `ML', `LIML', 0`fuller', `K', `ar', `null', `scoreBS', "`weighttype'", "`ptype'", ///
												"`madjust'", `N_h0s', "`Xnames_exog'", "`Xnames_endog'", 0`cons', ///
												"`Ynames'", "`b'", "`V'", "`W'", "`ZExclnames'", "`hold'", "`scnames'", `hasrobust', "`allclustvars'", `:word count `bootcluster'', `:word count `clustvars'', ///
												"`FEname'", "`wtname'", "`wtype'", "`C'", "`C0'", `reps', "`repsname'", "`repsFeasname'", `small', "`svmat'", "`dist'", ///
												`gridmin', `gridmax', `gridpoints', `matsizegb')
		_estimates unhold `hold'
		
		local reps = `repsname' // in case reduced to 2^G
		if !`ML' & `reps' & `:word count `clustvars''>1 {
			return scalar repsFeas = `repsFeasname'
			if `repsFeasname' < `reps' di _n "Warning: " `reps' - `repsFeasname' " replications returned an infeasible test statistic and were deleted from the bootstrap distribution."
		}
		
		if `N_h0s'>1 local _h _`h'

		if `df'==1 {
			return scalar `=cond(`small', "t", "z")'`_h' = `stat'
		}
		
		di
		if `reps' di as txt strproper("`boottype'") " bootstrap, null " cond(0`null', "", "not ") "imposed, " as txt `reps' as txt " replications, " _c
		di as txt cond(`ar', "Anderson-Rubin ", "") cond(!`reps' & `null' & "`boottype'"=="score", "Rao score (Lagrange multiplier)", "Wald") " test" _c
		if "`cluster'"!="" di ", clustering by " as inp "`cluster'" _c
		if ("`bootcluster'"!="" | `:word count `clustvars'' > 1) & `reps' di as txt ", bootstrap clustering by " as inp "`bootcluster'" _c
		if `reps'	di as txt ", " strproper("`weighttype'") " weights" _c
		di as txt ":"
		
		foreach c in `h0_`h'' {
			di "  " _c
			if !0`anythingh0' di as txt %4.0g `c' ": " _c
			di as res "`:constraint `c''"
		}

		if `df'==1 {
			local line1 = cond(`small', "t(`=`df_r'')", "z")
			di _n as txt _col(`=45-strlen("`line1'")') "`line1' = " as res %10.4f `stat' _n _col(`=45-strlen("``ptype''")') as txt "``ptype'' = " as res %10.4f `p'
			local `stat' = `stat' * `stat'
		}
		else {
			if `small' {
				local line1 F(`=`df'', `=`df_r'')
				local line2 Prob > F
			}
			else {
				local line1 chi2(`=`df'')
				local line2 Prob > chi2
			}
			di _n as txt _col(`=21-strlen("`line1'")') "`line1' = " as res %10.4f `stat' _n as txt _col(`=21-strlen("`line2'")') "`line2' = " as res %10.4f `p'
		}

		if "`madjust'" != "" {
			di _col(`=strlen("bonferroni")-strlen("`madjust'")+1') as txt strproper("`madjust'") "-adjusted prob = " as res %10.4f `padj'
			return scalar padj`_h' = `padj'
		}

		if "`plotmat'" != "" {
			tempvar X Y _plotmat
			mat `_plotmat' = `plotmat'
			return matrix plot`_h' = `plotmat'
			mat `plotmat' = `_plotmat'

			if "`graph'"=="" {
				cap confirm matrix `peakmat'
				if !_rc {
					mata st_matrix("`plotmat'", st_matrix("`plotmat'") \ st_matrix("`peakmat'"))
					return matrix peak`_h' = `peakmat'
				}
				cap confirm matrix `cimat'
				if _rc mat `_plotmat' = `plotmat'
				else {
					mata st_matrix("`_plotmat'", st_matrix("`plotmat'") \ ((st_matrix("`cimat'")[,1] \ st_matrix("`cimat'")[,2]), J(2*`=rowsof(`cimat')', 1, 1-`level'/100)))
				}
				mat colnames `_plotmat' = `X' `Y'
				qui svmat `_plotmat', names(col)
				label var `Y' " " // in case user turns on legend
				format `X' %5.0g
				if `"`graphname'"'=="Graph" cap graph drop Graph`_h'
				local 0, `graphopt'
				syntax, [LPattern(passthru) LWidth(passthru) LColor(passthru) LStyle(passthru) *]
				line `Y' `X', sort(`X') `lpattern' `lwidth' `lcolor' `lstyle' || scatter `Y' `X' if _n>rowsof(`plotmat'), mlabel(`X') mlabpos(6) xtitle("") ///
					ytitle(`"`=strproper("`madjust'")+ cond("`madjust'"!="","-adjusted r", "R")'ejection p value"') ///
					yscale(range(0 .)) name(`graphname'`_h', `replace') ///
					`=cond(`level'<100, "yline(`=1-`level'/100') ylabel(#6 `=1-`level'/100')", "")' legend(off) `options'
			}
		}

		if "`cimat'" != "" {
			di _n as res `level' "%" as txt " confidence set for null hypothesis expression: " _c
			forvalues i=1/`=rowsof(`cimat')' {
				if `i'>1 di " U " _c
				di  as txt "[" as res strofreal(`cimat'[`i',1], "%10.4g") as txt ", " as res strofreal(`cimat'[`i',2], "%10.4g") as txt "]" _c
			}
			di

			if inlist("`ptype'", "symmetric", "equaltail") {
				tempname t
				mata st_numscalar("`t'", anyof(st_matrix("`cimat'"), .))
				if `t' di "(A confidence interval could not be bounded. Try widening the search range with the {cmd:gridmin()} and {cmd:gridmax()} options.)"
			}
			mat colnames `cimat' = lo hi
			return matrix CI`_h' = `cimat'
		}
		
		return scalar `=cond(`small', "F", "chi2")'`_h' = cond(`df'==1, `stat'*`stat', `stat')
		if `small' return scalar df_r`_h' = `df_r'
		return scalar df`_h' = `df'
		return scalar p`_h' = `p'
		if `hasrobust' return local robust robust
		if "`svmat'"!="" return matrix dist`_h' = `dist'
	}
	return scalar level = `level'
	return local weighttype `weighttype'
	return local boottype `boottype'
	return local clustvars `clustvars'
	return scalar null = `null'
	return scalar reps = `repsname'
end

* Version history
* 2.1.6 Changed to faster computation for WCR with many clusters, but not quite WB. In multiway only applies to intersection of all clusters.
* 2.1.4 Fixed CI scaling issue introduced in 2.0.5 that affects scoretest, waldtest
* 2.1.3 Fixed crash on testing of multiple independent hypotheses on ML models
*       Added more return values and Roodman et al. cite to help file. Blocked warning about \alpha(B+1) being an integer for Rademacher with <=12 groups.
* 2.1.2 Fixed error in removing half-counting of ties in 2.0.6
*       Stopped crashes after mlogit (and probably other mlogit and mprobit commands)
* 2.1.1 Fixed failure to detect FE variable after xtreg, fe cluster()
* 2.1.0 Added matsizegb feature.
*       Fixed 2.0.6 failure to subtract 1 from Mammen, Webb weights in WRE non-AR
*       Fixed failure in subcluster bootstrap to sort data by error clusterings before bootstrap clustering
*       Avoided creating diagonal/sparse crosstab matrix
*       Terminate search for CI bounds when bracketing p values are within 1/reps (2/reps for equaltail), with 1 last linear interpolation, rather than find precise step-up point, which not so meaningful
*       Fixed minor new bugs when using Mata interface
* 2.0.6 Stopped (half-)counting ties. Changed default reps from 1000 to 999. Fixed swapped labeling of equal-tail and symmetric p-values(!).
* 2.0.5 Fixed subcluster bootstrap bugs: need to sort data even when c=1; don't falsely flag pure-robust case in WB subcluster
*       Fixed possible failure to find graph bounds when many replications infeasible and bounds not manually set
*       Fixed crash in FE estimation when FE cluster = error cluster
*       Accelerated generation of some wild weights by exploiting fact that results are invariant to rescaling of the weights
*       Finally optimized CI construction for nonull case.
*       Made default plot bounds more symmetric.
* 2.0.4 Made "unrestricted WRE" (WUE?) work.
* 2.0.3 Added automatic reporting of any infeasible replication statistics in multi-way clustering. Made r(reps) return value reflect possible reduction to 2^G.
* 2.0.2 Dropped citations from output but added reporting of weight type. Added warning if alpha*(B+1) not integer. Sped up Webb weight generation.
*       If seed() not specified, then return c(seed) in r(seed). For waldtest, nograph just compute CI analytically.
*       Prevented WRE crash when no clustering.
* 2.0.1 Reworked info matrix construction and organization to fix 2.0.0 bug in "subcluster" bootstrapping
* 2.0.0 Implemented code optimized for pure robust case. Allowed bootstrapping clusters to be chosen arbitrarily, independent of error clusterings.
* 1.9.7 Fixed crash on score bootstrap without observation weights. Improved run time when clusters are many by avoiding computation of Q'Q.
*       Fixed failure to recenter score test (not score bootstrap); bug introduced circa 1.9.0. Fixed failure to square t/z to make r(F)/r(chi2).
*       Fixed 1.9.6 bug causing normal weights to be replace by Mammen weights.
* 1.9.6 Added Gamma(4, .5) - 2 wild weight distribution option per Liu (1988)
* 1.9.5 Fixed score test bugs from 1.9.0, and bugs after ML estimation in Stata 15 because of new free parameter matrix label system
* 1.9.4 Cleaned up display of results for symmetric, equal-tail, etc.
* 1.9.3 Tweaked to no longer explode wild weights in multiway clustering; not needed since 1.9.0
* 1.9.2 Fixed crash with FE and omitted dummies for other vars. Fixed 1.9.0 crash in old Stata versions.
* 1.9.1 Fixed crash with FE and df>1. Stopped waldtest and scoretest ignoring small
* 1.9.0 Added fixed effect support
* 1.8.3 Added svmat(numer) option
* 1.8.2 Fixed bug after ML: was using V from unconstrained instead of constrained fit
* 1.8.1 Fixed bugs in handling robust non-cluster
* 1.8.0 Reworked multiway clustering to first collapse data to one obs per all-cluster-var intersections.
*       Reworked test stat computation for df>1 to mostly iterate over constraints rather than replications. Speeds AR test too.
* 1.7.1 Changed residual dof for multi-way clustered, small-sample-corrected models to smallest number of groups across grouping variables
* 1.7.0 Made bootcluster() accept more than one variable. Fixed error causing it to always bootstrap on combination of all vars in multi-way clustered models.
* 1.6.2 Fixed ado bug in 1.6.1
* 1.6.1 Fixed AR test crash. Dropped nowarning in favor of capture because commands such as poisson don't accept it. Changed left and right to lower and upper. Fixed bugs.
*       Suppressed non-concavity warning when imposing null after ML that was incorrectly triggered by omitted factor variables.
* 1.6.0 Added left and right p value types. Added cmdline option. Added nowarning option to ml, iter(0) call to suppress non-convergence warning.
* 1.5.7 renamed _selectindex() to boottest_selectindex() to reduce conflicts with old cmp versions
* 1.5.6 Fixed crash on waldtest after ML
* 1.5.5 Fixed bug in determining confidence intervals when some test results is missing.
* 1.5.4 Fixed two bugs causing crashes after GMM estimation
* 1.5.3 Simplify _selectindex(). Switch from invt() to invttail() since invt() added in Stata 13. Work around Mata garbage-collecting bug.
* 1.5.2 Switch from "[]" to "{}" in multiple hypothesis syntax. Prevent crashes if hypothesis applies to constrained/dropped parameter, or in non-robust ML.
* 1.5.1 Make wrapper accept level() and graphname(). Handle graphopt() options relevant for format of line.
* 1.5.0 Added Anderson-Rubin test, new null hypothesis syntax, multiple hypothesis testing and corrections thereof
* 1.4.0 Added WRE for LIML, sped up WRE, k-class and Fuller support after ivreg2
* 1.3.1 Bug fixes in WRE with observation weights and no clusters; and in clustered non-ML estimation with time series operators
*       Added level(100) option to disable plotting of CI's
* 1.3.0 Implemented Davidson & MacKinnon WRE. Bug fixes affecting non-robust fweight, weighted OLS/2SLS/GMM when null imposed
* 1.2.5 Added stable option to sort by cluster, for stability of results
* 1.2.4 More thorough clean-up of fv handling.
* 1.2.3 Fixed bug in search for upper bound of CI
* 1.2.2 Fixed bug in handling of fv base vars in non-ML estimation
* 1.2.1 LIML bug fixes.
* 1.2.0 Added svmat option. LIML support. Fixed bug in handling of weights. Replaced ci option with noci.
* 1.1.6 Fixed minor 1.1.4 and 1.1.5 bugs.
* 1.1.5 After ML estimation, don't restrict to sample until after scores generated, in case of ts ops
* 1.1.4 Erase old mlib before compiling new one. Don't recenter Z'e inside sandwich if not bootstrapping.
* 1.1.3 Fixed bug in multiway cluster ci's. Made compatible with estimation commands whose constraints() option only accept numlists. Drop out-of-sample obs *after* predicting scores.
* 1.1.2 Force use of score bootstrap after 2SLS/GMM. Thanks to Federico Belotti.
* 1.1.1 Added support for single-equation linear GMM with ivreg2 and ivregress.
* 1.1.0 Fixed 1.0.1 bug in observation weight handling. Added multiway clustering, robust, cluster(), and bootcluster() options. Added waldtest wrapper.
* 1.0.1 Added check for empty constraints. Fixed crash on use of weights after clustered estimation in Stata >13.
