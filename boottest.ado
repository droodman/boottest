*! boottest 1.5.7 19 June 2017
*! Copyright (C) 2015-17 David Roodman

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
	if `StataVersion' != c(stata_version) | "`CodeVersion'" < "01.05.00" {
		cap findfile "lboottest.mlib"
		while !_rc {
			erase "`r(fn)'"
			cap findfile "lboottest.mlib"
		}
		qui findfile "boottest.mata"
		run "`r(fn)'"
	}

	local cmd = cond(substr("`e(cmd)'", 1, 6)=="ivreg2", "ivreg2", "`e(cmd)'")

	if "`e(prefix)'" == "svy" {
		di as err "Doesn't work after {cmd:svy}."
		exit 198
	}
	if inlist(`"`cmd'"', "xtreg", "areg", "mvreg", "sureg") {
		di as err "Doesn't work after {cmd:`e(cmd)'}."
		exit 198
	}
	if inlist(`"`cmd'"', "sem", "gsem") & c(stata_version) < 14 {
		di as err "Requires Stata version 14.0 or later to work with {cmd:`e(cmd)'}."
		exit 198
	}
	if "`cmd'"=="ivreg2" & e(partial_ct)>0 & e(partial_ct)<. {
		di as err "Doesn't work after {cmd:ivreg2} with the {cmd:partial()} option."
		exit 198
	}

	tokenize `"`0'"', parse(",")
	if `"`1'"' != "," {
		local h0s `h0s' `1'
		macro shift
	}
	local 0 `*'
	syntax, [h0(numlist integer >0) Reps(integer 1000) seed(string) BOOTtype(string) CLuster(string) Robust BOOTCLuster(string) noNULl QUIetly WEIGHTtype(string) Ptype(string) NOCI Level(real `c(level)') SMall SVMat ///
						noGRaph gridmin(string) gridmax(string) gridpoints(string) graphname(string asis) graphopt(string asis) ar MADJust(string) *]

	if `reps'==0 local svmat

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
		confirm var `clustvars'
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

	if `:word count `clustvars'' <= 1 & `"`bootcluster'"' != "" {
		di as err "{cmdab:bootcl:uster()} option only accepted for multi-way clustered estimation."
		exit 198
	}

	if `:word count `weighttype'' > 1 {
		di as err "The {cmd:weight:type} option must be {cmdab:rad:emacher}, {cmdab:mam:men}, {cmdab:nor:mal}, or {cmdab:web:b}."
		exit 198
	}
	if "`weighttype'"'=="" local weighttype rademacher
	else {
		local 0, `weighttype'
		syntax, [RADemacher MAMmen NORmal WEBb]
		local weighttype `rademacher'`mammen'`normal'`webb'
	}

	if `:word count `ptype'' > 1 {
		di as err "The {cmd:wp:type} option must be {cmdab:sym:metric} or {cmdab:eq:qualtail}."
		exit 198
	}
	if "`ptype'"'=="" local ptype symmetric
	else {
		local 0, `ptype'
		syntax, [SYMmetric EQualtail]
		local ptype `symmetric'`equaltail'
	}

	local ML = e(converged) != .
	local IV = "`e(instd)'" != ""
	local LIML = ("`cmd'"=="ivreg2" & "`e(model)'"=="liml") | ("`cmd'"=="ivregress" & "`e(estimator)'"=="liml")
	local WRE = `"`boottype'"'!="score" & `IV' & `reps'
	local small = e(df_r) != . | "`small'" != ""

	local fuller `e(fuller)' // "" if missing
	local K = e(kclass) // "." if missing

	if "`boottype'"'=="" local boottype = cond(`ML', "score", "wild")
	else {
		local 0, `boottype'
		syntax, [SCore Wild]
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

	if `"`seed'"'!="" set seed `seed'

	tempname p padj se stat df df_r hold C C0 CC0 b V keepC keepW
	mat `b' = e(b)
	mat `V' = e(V)
	local k = colsof(`b')

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
	local k_autoCns = r(k_autoCns)
	if "`e(Cns)'" != "" mat `C' = e(Cns)
	cap _estimates drop `hold'

	if !`ML' {
		if ("`cmd'"=="ivreg2" & !inlist("`e(model)'", "ols", "iv", "gmm2s", "liml")) |  ("`cmd'"=="ivregress" & !inlist("`e(estimator)'", "liml", "2sls", "gmm")) {
			di as err "Only for use with OLS, 2SLS, single-equation GMM, and ML-based estimators."
			exit 198
		}

		local Ynames `e(depvar)'

		if `IV' {
			local Xnames_exog = cond("`cmd'"=="ivreg2", "`e(inexog)'", "`e(exogr)'")
			local Xnames_endog `e(instd)'
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
		
		local _cons = cond(`cons', "_cons", "")

		local GMM = ("`cmd'"=="ivreg2" & "`e(model)'"=="gmm2s") | ("`cmd'"=="ivregress" & "`e(estimator)'"=="gmm")
		if `GMM' {
			if `reps' di as txt _n "Bootstrapping purely with final GMM weight matrix."
			tempname W
			mat `W' = e(W)
			mata _boottestp = order(tokens("`:colnames `W''")', 1)[invorder(order(tokens("`_cons' `Xnames_exog' `ZExclnames'")', 1))]
			mata st_matrix("`W'", st_matrix("`W'")[_boottestp,_boottestp]) // ensure weight matrix in order cons-other included exog-excluded exog
		}

		mata _boottestp = order(tokens("`:colnames e(b)'")', 1)[invorder(order(tokens("`_cons' `Xnames_exog' `Xnames_endog'")', 1))] \ `=`k'+1'
		cap mata st_matrix("`C'" , st_matrix("`C'" )[.,_boottestp])
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
			forvalues i = 1/`:word count `revarlist'' {
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

		cap mata _boottestC = st_matrix("`C'" )[,(st_matrix("`keepC'"),`=`k'+1')]; _boottestC = select(_boottestC,rowsum(_boottestC:!=0)); st_matrix("`C'" , _boottestC)
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

		if `:word count `clustvars'' > 1 {
			if `"`bootcluster'"' == "" {
				if `reps' {
					di as err "{cmdab:bootcl:uster()} option required after multi-way clustered estimation."
					exit 198
				}
				else local bootcluster: word 1 of `clustvars'
			}
			cap confirm var `bootcluster'
			if _rc {
				di as err `"`bootcluster' not found."'
				exit 198
			}
			if `"`:list clustvars & bootcluster'"' == "" {
				di as txt "`bootcluster'" as err " not among the estimate's clustering variables."
				exit 198
			}
			
			local clustvars `bootcluster' `:list clustvars - bootcluster'
		}

		sort `clustvars', stable

		foreach clustvar in `clustvars' {
			cap confirm numeric var `clustvar'
			if _rc {
				tempname IDname
				qui egen long `IDname' = group(`clustvar') if e(sample)
				local _clustvars `_clustvars' `IDname'
			}
			else local _clustvars `_clustvars' `clustvar'
		}
		local clustvars `_clustvars'
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
		if `k_autoCns' mat `C0' = `C0'[`k_autoCns'+1...,1...]
		if "`noci'"=="" & !`ML' & rowsof(`C0')==1 {
			tempname plotmat
			if `level'<100 tempname cimat
		}

		if rowsof(`C0')>1 & "`ptype'"=="equaltail" di as txt "Note: {cmd:ptype(equaltail)} ignored for multi-constraint null hypotheses."

		if `ML' {
			local K .

			if `null' {
				if "`e(cmdline)'"=="" {
					di as err "Can only impose null after Maximum Likelihood-based estimation commands that provide the return value e(cmdline)."
					exit 198
				}
				if "`cmd'"=="tobit" {
					di as err "Cannot impose null after {help tobit}. Try re-fitting the model with {stata ssc describe cmp:cmp}, {help intreg}, or {help gsem}."
					exit 198
				}

				mat `CC0' = `C0' \ nullmat(`C') // add null to model constraints

				`quietly' di as res _n "Re-running regression with null imposed." _n
				local 0 `e(cmdline)'
				syntax [anything] [aw pw fw iw] [if] [in], [CONSTraints(passthru) from(passthru) INIt(passthru) ITERate(passthru) CLuster(passthru) Robust vce(string) *]

				cap _estimates drop `hold'
				_estimates hold `hold', restore
				local coleq: coleq `b'
				if "`coleq''"=="_" local coleq
				local colnames: colnames `b'
				forvalues r=1/`=rowsof(`CC0')' {
					local terms 0
					local _constraint
					forvalues c=1/`k' {
						if `CC0'[`r',`c'] {
							if `terms++' local _constraint `_constraint' +
							local _constraint `_constraint' `=cond(`CC0'[`r',`c']!=1,"`=`CC0'[`r',`c']' * ","")' `=cond("`coleq'"=="","","[`:word `c' of `coleq'']")'`:word `c' of `colnames''
						}
					}
					local _constraint `_constraint' = `=`CC0'[`r',`=`k'+1']'
					constraint free
					local _constraints `_constraints' `r(free)'
					constraint `r(free)' `_constraint'
				}
				
				cap `=cond("`quietly'"=="", "noisily", "")' `anything' if `hold' `=cond("`weight'"=="","",`"[`weight'`exp']"')', constraints(`_constraints') `from' `init' `iterate' `options'
				local rc = _rc
				constraint drop `_constraints'
				if e(converged)==0 {
					di as err "Could not impose null."
					exit 430
				}
				if !`rc' {
					tempname b0
					mat `b0' = e(b)
					cap noi `anything' if e(sample), `=cond(inlist("`cmd'", "cmp","ml"),"init","from")'(`b0') iterate(0) `options'
					local rc = _rc
				}
				if `rc' {
					di as err "Error imposing null. Perhaps {cmd:`e(cmd)'} does not accept the {cmd:constraints()}, {cmd:from()}, and {cmd:iterate()} options, as needed."
					exit `rc'
				}
				mat `V' = e(V)
				if diag0cnt(`V') {
					`quietly' di _n as res "Warning: Negative Hessian under null hypothesis not positive definite. Results may be unreliable."
					`quietly' di           "This may indicate that the constrained optimum is far from the unconstrained one, i.e., that the hypothesis is incompatible with the data." _n
				}
			}
			else mat `V' = e(V`=cond("`e(V_modelbased)'"!="", "_modelbased", "")')

			tempvar sc
			cap noi predict double `sc'* if e(sample), score
			if _rc {
				di as err "Could not generate scores from `e(cmd)' with {cmd:predict}."
				exit _rc
			}

			unab scnames: `sc'*
			if `:word count `scnames'' == e(k_eq) { // generate score for each parameter from those for each equation
				tempname beq
				local coleq: coleq `b'
				tokenize `:list uniq coleq'
				qui forvalues eq = 1/`:word count `scnames'' {
					mat `beq' = `b'[1, "``eq'':"]
					foreach var in `:colnames `beq'' {
						if "`var'" == "_cons" local _sc `_sc' `sc'`eq'
						else {
							tempvar t
							gen double `t' = `sc'`eq' * `var' if e(sample)
							local _sc `_sc' `t'
						}
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

			mata st_matrix("`C0'", st_matrix("`C0'")[.,_boottestp]) // ensure constraint matrix in order cons-other exog-endog

			if `ar' {
				if !`IV' | `GMM' {
					di as err "Anderson-Rubin test only available for non-GMM instrumental variables models."
					exit 198
				}
				local kEnd: word count `Xnames_endog'
				mata st_numscalar("`p'", rows(st_matrix("`C0'"))!=`kEnd' | (`k'==`kEnd'? 0 : any(st_matrix("`C0'")[|.,.\.,`k'-`kEnd'|])))
				if `p' {
					di as err "For Anderson-Rubin test, null hypothesis must constrain all the instrumented variables and no others."
					exit 198
				}			
			}
			mata _boottestC = st_matrix("`C0'")[,(st_matrix("`keepC'"),`=`k'+1')]; _boottestC = select(_boottestC,rowsum(_boottestC:!=0)); st_matrix("`C0'", _boottestC)
			cap confirm matrix `C0'
			if _rc {
				di as err "Hypothesis refers to constrained parameter."
				exit 111
			}
		}

		mata boottest_stata("`stat'", "`df'", "`df_r'", "`p'", "`padj'", "`cimat'", "`plotmat'", `level', `ML', `LIML', 0`fuller', `K', `ar', `null', `scoreBS', "`weighttype'", "`ptype'", ///
												"`madjust'", `N_h0s', "`Xnames_exog'", "`Xnames_endog'", 0`cons', ///
												"`Ynames'", "`b'", "`V'", "`W'", "`ZExclnames'", "`hold'", "`scnames'", `hasrobust', "`clustvars'",  "`wtname'", "`wtype'", "`C'", "`C0'", `reps', `small', "`dist'", ///
												`gridmin', `gridmax', `gridpoints')

		_estimates unhold `hold'
		
		if `N_h0s'>1 local _h _`h'

		if `df'==1 & !`ar' {
			return scalar `=cond(`small', "t", "z")'`_h' = `stat'
			scalar `stat' = `stat' * `stat'
		}
		
		di _n
		if `reps' di as txt strproper("`boottype'") " bootstrap (" cond(`scoreBS',"Kline and Santos 2012",cond(`IV',"Davidson & MacKinnon 2010","Wu 1986")) "), null " cond(0`null', "", "not ") "imposed, " as txt `reps' as txt " replications, " _c
		di as txt cond(`ar', "Anderson-Rubin ", "") cond(!`reps' & `null' & "`boottype'"=="score", "Rao score (Lagrange multiplier)", "Wald") " test:" _n
		
		foreach c in `h0_`h'' {
			di "  " _c
			if !0`anythingh0' di as txt %4.0g `c' ": " _c
			di as res "`:constraint `c''"
		}

		if `small' di _n as txt    "     " cond(`reps', "quasi-","      ")    "F(" %3.0f `df' "," %6.0f `df_r' ") = " as res %10.4f `stat' _n as txt "                Prob > F = " as res %10.4f `p'
			else     di _n as txt "        " cond(`reps', "quasi-","      ") "chi2(" %3.0f `df'                  ") = " as res %10.4f `stat' _n as txt "            Prob > chi2 = "  as res %10.4f `p'
		
		if "`madjust'" != "" {
			di _col(`=strlen("bonferroni")-strlen("`madjust'")+1') as txt strproper("`madjust'") "-adjusted prob = " as res %10.4f `padj'
			return scalar padj`_h' = `padj'
		}

		if "`graph'"=="" & "`plotmat'"!="" {
			tempvar X Y _plotmat

			if "`cimat'"!="" mat `_plotmat' = `plotmat' \ ((`cimat'[1...,1] \ `cimat'[1...,2]), J(2*rowsof(`cimat'), 1, 1-`level'/100))
				else           mat `_plotmat' = `plotmat'
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

			return matrix plot`_h' = `plotmat'
		}

		if "`cimat'" != "" {
			di _n as res `level' "%" as txt " confidence set for null hypothesis expression: " _c
			forvalues i=1/`=rowsof(`cimat')' {
				if `i'>1 di " U " _c
				di  as txt "[" as res strofreal(`cimat'[`i',1], "%10.4g") as txt ", " as res strofreal(`cimat'[`i',2], "%10.4g") as txt "]" _c
			}
			di

			tempname t
			mata st_numscalar("`t'", anyof(st_matrix("`cimat'"), .))
			if `t' di "(A confidence interval could not be bounded. Try widening the search range with the {cmd:gridmin()} and {cmd:gridmax()} options.)"

			return matrix CI`_h' = `cimat'
		}
		
		return scalar `=cond(`small', "F", "chi2")'`_h' = `stat'
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
	return scalar reps = `reps'
	return local seed `seed'
end

* Version history
* 1.5.7 renamed _selectindex() to boottest_selectindex() to reduce conflicts with old cmp versions
* 1.5.6 Fixed crash on waldtest after ML
* 1.5.5 Fixed bug in determining confidence intervals when some test results is missing.
* 1.5.4 Fixed two bugs causing crashes after GMM estimation
* 1.5.3 Simplify _selectindex(). Switch from invt() to invttail() since invt() added in Stata 13. Work around Mata garbage-cleaning bug.
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
