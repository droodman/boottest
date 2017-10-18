*! boottest 1.9.0 18 October 2017
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
	if `StataVersion' != c(stata_version) | "`CodeVersion'" < "01.09.00" {
		cap findfile "lboottest.mlib"
		while !_rc {
			erase "`r(fn)'"
			cap findfile "lboottest.mlib"
		}
		qui findfile "boottest.mata"
		run "`r(fn)'"
	}

	local cmd = cond(substr("`e(cmd)'", 1, 6)=="ivreg2", "ivreg2", "`e(cmd)'")
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
	syntax, [h0(numlist integer >0) Reps(integer 1000) seed(string) BOOTtype(string) CLuster(string) Robust BOOTCLuster(string) noNULl QUIetly WEIGHTtype(string) Ptype(string) NOCI Level(real `c(level)') SMall SVMat ///
						noGRaph gridmin(string) gridmax(string) gridpoints(string) graphname(string asis) graphopt(string asis) ar MADJust(string) CMDline(string) *]

	if `reps'==0 local svmat
		else {
			local _svmat = cond("`svmat'"=="", "", "t")
			local 0, `options'
			syntax, [SVMat(string) *]
			if !inlist(`"`svmat'"', "numer", "t", "") {
				di as err "option " as inp "svmat(" `svmat' ")" as err " not allowed"
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
	local LIML = ("`ivcmd'"=="ivreg2" & "`e(model)'"=="liml") | ("`cmd'"=="ivregress" & "`e(estimator)'"=="liml")
	local WRE = `"`boottype'"'!="score" & `IV' & `reps'
	local small = e(df_r) != . | "`small'" != ""

	local fuller `e(fuller)' // "" if missing
	local K = e(kclass) // "." if missing

	if `ar' & !`IV' {
		di as err "Anderson-Rubin test is only for IV models."
		exit 198
	}

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
	
	local FEname = cond(inlist("`cmd'","xtivreg","xtivreg2"), "`e(ivar)'", "`e(absvar)'`e(absvars)'")

	if `"`seed'"'!="" set seed `seed'

	tempname p padj se stat df df_r hold C C0 CC0 b V keepC keepW
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
	local k_autoCns = r(k_autoCns)
	if "`e(Cns)'" != "" mat `C' = e(Cns)
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
		
		cap mata st_matrix("`C'" , st_matrix("`C'" )[,_boottestp])
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
				local bootcluster `clustvars'
				if `reps' di as txt "({cmdab:bootcl:uster(`clustvars')} assumed)"
			}
			else {
				confirm var `bootcluster'
				if `"`:list bootcluster - clustvars'"' != "" {
					di as txt "{cmdab:bootcl:uster()} option includes variables not among the estimate's clustering variables."
					exit 198
				}
			}
			
			local clustvars `bootcluster' `:list clustvars - bootcluster' // sort, by bootstrapping clusters first
		}
		else local bootcluster `clustvars'

		sort `clustvars', stable

		foreach clustvar in `clustvars' {
			cap confirm numeric var `clustvar'
			if _rc {
				tempvar IDname
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
			tempname plotmat peakmat
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
				if "`cmd'"=="tobit" {
					di as err "Cannot impose null after {help tobit}. Try re-fitting the model with {stata ssc describe cmp:cmp}, {help intreg}, or {help gsem}."
					exit 198
				}

				mat `CC0' = `C0' \ nullmat(`C') // add null to model constraints

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
					cap `anything' if e(sample), `=cond(inlist("`cmd'", "cmp","ml"),"init","from")'(`b0') iterate(0) `options' `max'
					local rc = _rc
				}
				if `rc' {
					di as err "Error imposing null. Perhaps {cmd:`e(cmd)'} does not accept the {cmd:constraints()}, {cmd:from()}, and {cmd:iterate()} options, as needed."
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

		mata boottest_stata("`stat'", "`df'", "`df_r'", "`p'", "`padj'", "`cimat'", "`plotmat'", "`peakmat'", `level', `ML', `LIML', 0`fuller', `K', `ar', `null', `scoreBS', "`weighttype'", "`ptype'", ///
												"`madjust'", `N_h0s', "`Xnames_exog'", "`Xnames_endog'", 0`cons', ///
												"`Ynames'", "`b'", "`V'", "`W'", "`ZExclnames'", "`hold'", "`scnames'", `hasrobust', "`clustvars'", `:word count `bootcluster'', ///
												"`FEname'", "`wtname'", "`wtype'", "`C'", "`C0'", `reps', `small', "`svmat'", "`dist'", ///
												`gridmin', `gridmax', `gridpoints')

		_estimates unhold `hold'
		
		if `N_h0s'>1 local _h _`h'

		if `df'==1 {
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

		if "`plotmat'"!="" {
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
* 1.9.0 Added FE support
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
