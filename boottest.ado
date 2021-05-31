*! boottest 3.2.2 31 May 2021
*! Copyright (C) 2015-21 David Roodman

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
  cap version 13.1
  if _rc {
     di as err "This version of {cmd:boottest} requires Stata version 13 or later. An older version compatible with Stata `c(stata_version)'"
     di as err "is at https://github.com/droodman/boottest/releases/tag/v2.6.0."
     exit _rc
  }
	cap noi _boottest `0'
	local rc = _rc
	constraint drop `anythingh0'
	cap mata mata drop _boottestp
	cap mata mata drop _boottestC
	cap mata mata drop _boottestkeepC
	cap mata mata drop _boottestm
	cap mata mata drop _boottestt
	cap drop `contourvars'
	exit `rc'
end

cap program drop _boottest
program define _boottest, rclass sortpreserve
	version 11

	local   cmd = cond(substr("`e(cmd)'", 1, 6)=="ivreg2", "ivreg2", "`e(cmd)'")
	local ivcmd = cond(inlist("`cmd'","reghdfe","ivreghdfe"), cond("`e(model)'"=="iv", "ivreg2", ""), cond("`cmd'"=="xtivreg2", "ivreg2", "`cmd'"))

	if "`e(cmd)'" == "" {
		di as err "No estimates detected."
		error 198
	}
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
	if inlist("`cmd'","reghdfe","ivreghdfe") {
		if `:word count `e(extended_absvars)''>1 {
			di as err "Doesn't work after {cmd:`cmd'} with more than one set of absorbed fixed effects or with absorbed interaction terms."
			exit 198
		}
		cap _fv_check_depvar `e(extended_absvars)'
		if _rc {
			di as err "Doesn't work after {cmd:`cmd'} with absorbed interaction terms."
			exit 198
		}
	}
	if inlist("`cmd'", "sem", "gsem") & c(stata_version) < 14 {
		di as err "Requires Stata version 14.0 or later to work with {cmd:`e(cmd)'}."
		exit 198
	}
	if "`e(cmd)'"=="regress" & "`e(model)'" == "iv" {
		di as err "Doesn't support the undocumented IV syntax of {cmd:regress}."
		exit 198
	}
  
  if ("`ivcmd'"=="ivreg2" & "`e(model)'"=="gmm2s") | ("`cmd'"=="ivregress" & "`e(estimator)'"=="gmm") {
  	di as err "GMM no longer supported."
    exit 198
  }

	tokenize `"`0'"', parse(",")
	if `"`1'"' != "," {
		local h0s `h0s' `1'
		macro shift
	}
	local 0 `*'
	syntax, [h0(numlist integer >0) Reps(integer 999) seed(string) BOOTtype(string) CLuster(string) Robust BOOTCLuster(string) noNULl QUIetly WEIGHTtype(string) Ptype(string) STATistic(string) NOCI Level(real `c(level)') SMall SVMat ///
						noGRaph gridmin(string) gridmax(string) gridpoints(string) graphname(string asis) graphopt(string asis) ar MADJust(string) CMDline(string) MATSIZEgb(integer 1000000) PTOLerance(real 1e-6) svv MARGins *]

  local margins = "`margins'" != ""
  if `margins' & `"`h0s'`h0'"' != "" {
    di as err "Include the {cmd:margins} option or state null hyptheses, but don't do both."
    exit 198
  }

	if `matsizegb'==1000000 local matsizegb .
  
  if "`svv'" != "" tempname svv

	if `reps'==0 local svmat
		else {
			local _svmat = cond("`svmat'"=="", "", "t")
			local 0, `options'
			syntax, [SVMat(string) *]
			if !inlist(`"`svmat'"', "numer", "t", "") {
				di as err "option " as res "svmat(" `svmat' ")" as err " not allowed."
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

	if `"`gridpoints'"'!="" {
		foreach g of local gridpoints {
			cap confirm integer number `g'
			local rc = _rc
			if !_rc {
				local rc = `g' <= 0
			}
			if `rc' {
				di as err "{cmd:gridpoints()} entry not a positive integer."
				exit 198
			}
		}
	}

	foreach macro in gridmin gridmax {
		if `"``macro''"'!="" {
			foreach g of local `macro' {
				if `"`g'"' != "." {
					cap confirm number `g'
					if _rc {
						di as err "{cmd:`macro'()} entry not numeric."
						exit 198
					}
				}
			}
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
		local clustvars `e(clustvar)'
		local hasrobust = "`e(vcetype)'" == "Robust" | "`clustvars'" != ""
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
 
	if `"`statistic'"'=="" local statistic t
	else if !inlist(`"`statistic'"', "t", "c") {
		di as err "The {cmd:stat:istic} option must be {cmd:t} or {cmd:c}."
		exit 198
	}
	else if "`statistic'"=="c" & `reps'==0 {
		di as err "{cmdab:stat:istic(c)} not allowed with non-bootstrap tests."
		exit 198
	}

	local ML = e(converged) != .
	local IV = "`e(instd)'`e(endogvars)'" != ""
	local LIML = ("`cmd'"=="ivreg2" & "`e(model)'"=="liml") | ("`cmd'"=="ivregress" & "`e(estimator)'"=="liml")| (inlist("`cmd'","reghdfe","ivreghdfe") & strpos("`e(title)'", "LIML"))
	local WRE = `"`boottype'"'!="score" & `IV' & `reps'
	local small = e(df_r) != . | "`small'" != "" | e(cmd)=="cgmreg"
  local partial = 0`e(partial_ct)' & "`e(cmd)'"=="ivreg2"

  local DID = inlist("`cmd'", "didregress", "xtdidregress")
  local treatment `e(treatment)'
  local DDD = `:word count `e(groupvars)'' == 2

	local fuller `e(fuller)'  // "" if missing
	local K = e(kclass)  // "." if missing

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
	
  if `DDD' {
    tempname FEname
    qui egen long `FEname' = group(`e(groupvars)') if e(sample)
    sum `FEname', meanonly
    local NFE = r(max)
  }
  else {
    local  NFE    = cond(inlist("`cmd'","xtreg","xtivreg","xtivreg2"),   e(N_g)   , cond(`DID', 0`e(N_clust)', cond("`cmd'"=="areg", 1+e(df_a), /*reghdfe*/ max(0`e(K1)',0`e(df_a_initial)'))))
    local _FEname = cond(inlist("`cmd'","xtreg","xtivreg","xtivreg2"), "`e(ivar)'", cond(`DID', "`e(groupvars)'", "`e(absvar)'`e(extended_absvars)'"))
    if `"`_FEname'"' != "" {
      cap confirm numeric var `_FEname'
      if _rc {
        tempvar FEname
        qui egen long `FEname' = group(`_FEname') if e(sample)
      }
      else local FEname `_FEname'
    }
  }

  local FEdfadj = !inlist("`cmd'","xtreg","xtivreg","xtivreg2","xtdidregress")  // these commands don't count time FE in DOF adjustment
  if "`cmd'" == "xtreg" {  // exception: xtreg, fe. dfadj. https://stata.com/statalist/archive/2013-01/msg00540.html
    local 0 `e(cmdline)'
    syntax [anything], [dfadj *]
    local FEdfadj = "`dfadj'" != ""
  }

	if `"`seed'"'!="" set seed `seed'

	tempname p padj se teststat df df_r hold C1 C R1 R r1 r1r R1R r b V b0 V0 keepC repsname repsFeasname t NBootClustname marginsH0
	mat `b' = e(b)
	if "`e(wtype)'" != "" {
		tokenize `e(wexp)'
		cap confirm var `2'
		if _rc {
			tempname wtname
			qui gen double `wtname' `e(wexp)'
		}
		else local wtname `2'
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
	else if `margins' {
  	if strlen("`r(Jacobian)'") {
    	mata _boottestp = selectindex(rowsum(st_matrix("r(Jacobian)"):!=0))  // skip all-zero rows
      mata st_local("marginsnames", invtokens(st_matrixrowstripe("r(Jacobian)")[_boottestp,2]'))
      mata st_matrix("`marginsH0'", st_matrix("r(Jacobian)")[_boottestp,])
      mata st_local("N_h0s", strofreal(length(_boottestp)))
      if `N_h0s'==0 {
        di as err "No valid {cmd:margins} results not found."
        error 303
      }
      scalar `df' = 1  // always when working with margins results, df = 1 and constant term = 0
      mat `r' = 0
    }
    else {
    	di as err "{cmd:margins} results not found."
      error 303
    }
  }
  else {
		local N_h0s 1 // number of nulls
		if "`h0'" == "" {
			di as txt _n "({cmd:h0(1)} assumed)"
			local h0 1
		}
		foreach c of numlist `h0' {
			`quietly' if `"`: constraint `c''"' == "" di as res "Constraint `c' not found and will be skipped."
		}
		local h0_1 `h0'
	}

	if `N_h0s'==1 local madjust

	makecns
	if "`e(Cns)'" != "" {
		local hascns 1
		mat `C1' = e(Cns)
    mat `R1' = `C1'[1...,1..colsof(`C1')-1]
    mat `r1' = `C1'[1...,   colsof(`C1')  ]
	}

	cap _estimates drop `hold'
	if !`ML' {
		if ("`ivcmd'"=="ivreg2" & !inlist("`e(model)'", "ols", "iv", "gmm2s", "liml")) | ("`ivcmd'"=="ivregress" & !inlist("`e(estimator)'", "liml", "2sls", "gmm")) {
			di as err "Only for use with OLS, 2SLS, single-equation GMM, and ML-based estimators."
			exit 198
		}

		local Ynames `e(depvar)'

		local colnames: colnames e(b)
    if `partial' local colnames `colnames' `e(partial1)'
    local k = `:word count `colnames'' + (`partial' & 0`e(partialcons)')
		local _cons _cons
		local colnames: list colnames - _cons
		if `IV' {
			local Xnames_endog `e(instd)'
			local Xnames_exog: list colnames - Xnames_endog
			local cons = inlist("`cmd'`e(cons)'`e(constant)'", "ivreg21", "ivregress") | e(partialcons)
			local ZExclnames `e(insts)'
			local ZExclnames: list ZExclnames - Xnames_exog
		}
		else {
			local Xnames_exog: colnames e(b)
      if `DID' local Xnames_exog: subinstr local Xnames_exog "r1vs0." ""
			local cons = "`:list Xnames_exog & _cons'" != ""
      local Xnames_exog: list Xnames_exog - _cons
		}

		local cons = `cons' & !0`NFE'  // no constant in FE model
		if !`cons' local _cons

		mata _boottestp = J(`cons', 1, `k') \ order(tokens("`colnames'")', 1)[invorder(order(tokens("`Xnames_exog' `Xnames_endog'")', 1))]  // for putting vars in cons-exog-endog order
		if `cons' mat `keepC' = 1
		local colnames `_cons' `Xnames_exog' `Xnames_endog'

		foreach varlist in Xnames_exog Ynames Xnames_endog ZExclnames {
			local _revarlist
			forvalues i=1/`:word count ``varlist''' {
				local var: word `i' of ``varlist''
				_ms_parse_parts `var'
				if !r(omit) {
        	fvrevar `:word `i' of ``varlist''' if e(sample)
					local _revarlist `_revarlist' `r(varlist)'
					if inlist("`varlist'", "Xnames_exog", "Xnames_endog") {
						local pos: list posof "`var'" in colnames
						if `pos' mat `keepC' = nullmat(`keepC'), `pos'
					}
				}
			}
			local `varlist' `_revarlist'
		}
		mata _boottestkeepC = st_matrix("`keepC'"); _boottestp = cols(_boottestkeepC)? _boottestp[_boottestkeepC] : J(0,1,0)

		if 0`hascns' {
      if `partial' mat `R1' = `R1', J(rowsof(`R1'), `:word count `e(partial1)'', 0)  // add blanks for partialled-out 
      mata st_matrix("`R1'", st_matrix("`R1'")[,_boottestp])  // put cols in standardized order & drop those for omitted regressors
      mata _boottestkeepC = rowsum(st_matrix("`R1'"):!=0)
      mata st_matrix("`R1'", select(st_matrix("`R1'"), _boottestkeepC))  // drop rows corresponding purely to omitted variables
      mata st_matrix("`r1'", select(st_matrix("`r1'"), _boottestkeepC))
    }    

		if `cons' local Xnames_exog `hold' `Xnames_exog'  // add constant term
	}

	local Nclustvars: word count `clustvars'

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
			if `reps' & `Nclustvars'>1 di as txt "({cmdab:bootcl:uster(`clustvars')} assumed)"
		}
		local    clustvars `:list clustvars & bootcluster' `:list clustvars - bootcluster'
		local allclustvars `:list bootcluster - clustvars' `clustvars' // put bootstrapping clusters first, error clusters last, overlap in middle

		sort `clustvars' `:list bootcluster - clustvars', stable  // bootstrapping-only clusters left out of this sort. Effect is that in subcluster bootstrap, bootstrapping clusters are intersections of all clusterings

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
    if `margins' mat `R' = `marginsH0'[`h', 1...]
    else {
    	_estimates hold `hold', restore
      ereturn post `b'
      mat `b' = e(b)

    	// process hypothesis constraints into e(Cns)
      qui makecns `h0_`h''

      if "`h0_`h''" != "`r(clist)'" {
        local clist `r(clist)'
        local clist: list h0_`h' - clist
        foreach c in `clist' {
          di as txt `"note: constraint `=cond(0`anythingh0', `"{res}`:constraint `c''{txt}"', "number `c'")' caused error {search r(111)}"'
          if `DID' di `"If your test relates to the treatment effect, you may need to reference "{cmd:r1vs0.`treatment'}" instead of "{cmd:`treatment'}"."' _n
        }
        exit 111
      }
      
      cap mat `C' = e(Cns)
      if _rc exit 111
      _estimates unhold `hold'
      if r(k_autoCns) mat `C' = `C'[r(k_autoCns)+1...,1...]
      scalar `df' = rowsof(`C')
      mat `R' = `C'[1...,1..colsof(`C')-1]
      mat `r' = `C'[1...,   colsof(`C')  ]
    }

		if `NFE' {
			mata st_local("rc", strofreal(any(0:!=select(st_matrix("`R'"),st_matrixcolstripe("`R'")[,2]':=="_cons"))))

			if `rc' {
				di as err "In fixed-effect models, null hypotheses may not involve constant term."
				exit 111
			}
		}
		
		local plotmat
		local peakmat
		local cimat
		if "`noci'"=="" & !`ML' & (`df'<=1 | `df'==2 & "`graph'"=="") {
			if `reps' | !`scoreBS' | `null' | "`graph'"=="" tempname plotmat peakmat // don't request plot if doing simple Wald with no graph
			if `level'<100 & `df'==1 tempname cimat
		}

		if `df'>1 & "`ptype'"!="symmetric" di as txt "Note: {cmd:ptype(`ptype')} ignored for multi-constraint null hypotheses."

		if `df'<=2 { // a bit duplicative of code in the if-`ML' block just below...
			if  "`graph'"=="" {  // construct axis label(s)
				if `margins' local constraintLHS1 "Margins of `:word `h' of `marginsnames''"
        else {
        	local coleq: coleq `b'
          if "`:word 1 of `coleq''"=="_" local coleq
          local colnames: colnames `b'
          forvalues i=1/`=`df'' {
            local terms 0
            local constraintLHS`i'
            forvalues j=1/`=colsof(`R')' {
              if `R'[`i',`j'] {
                if `terms++' local constraintLHS`i' `constraintLHS`i''+
                local constraintLHS`i' `constraintLHS`i'' `=cond(`R'[`i',`j']!=1,"`=`R'[`i',`j']'*","")'`=cond("`coleq'"=="","","[`:word `j' of `coleq'']")'`:word `j' of `colnames''
              }
            }
					}
				}
			}
			
			foreach option in gridmin gridmax gridpoints {
				if "``option''" == "" local `option' = `df' * ". "
				else if `df' != `:word count ``option''' {
					di as err "{cmd:`option'()} option needs " cond(`df'==1, "one entry", "two entries")
					exit 199
				}
			}
		}

		if `ML' {
			local K .
      local k = colsof(`b')

			if `null' {
				if "`e(cmdline)'"=="" & `"`cmdline'"'=="" {
					di as err "Original estimation command line needed. Include it in a {cmd:cmdline()} option."
					exit 198
				}
				if "`cmd'"=="tobit" & c(stata_version)<15 {
					di as err "Cannot impose null after {help tobit} in Stata versions below 15. Try re-fitting the model with {stata ssc describe cmp:cmp}, {help intreg}, or {help gsem}."
					exit 198
				}

				mat `R1R' = `R' \ nullmat(`R1') // add null to model constraints
				mat `r1r' = `r' \ nullmat(`r1')

				* get rid of any remaining o. and b. constraints; r(k_autcns) doesn't capture all
				mata _boottestC = st_matrixcolstripe("`R1R'")[,2]
				mata _boottestC = rowsum(select(st_matrix("`R1R'"), !(strmatch(_boottestC, "*b*.*") :| strmatch(_boottestC, "*o*.*"))') :!= 0)
				mata st_matrix("`R1R'", select(st_matrix("`R1R'"), _boottestC))
				mata st_matrix("`r1r'", select(st_matrix("`r1r'"), _boottestC))
				cap mat `R1R'[1,1] = `R1R'[1,1]
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
				if "`:word 1 of `coleq''"=="_" local coleq
				local colnames: colnames `b'
				local _constraints
				forvalues i=1/`=rowsof(`R1R')' {
					local terms 0
					local _constraint
					forvalues j=1/`=colsof(`R1R')' {
						if `R1R'[`i',`j'] {
							if `terms++' local _constraint `_constraint' +
							local _constraint `_constraint' `=cond(`R1R'[`i',`j']!=1,"`=`R1R'[`i',`j']' * ","")' `=cond("`coleq'"=="","","[`:word `j' of `coleq'']")'`:word `j' of `colnames''
						}
					}
					local _constraint `_constraint' = `=`r1r'[`i',1]'
					constraint free
					local _constraints `_constraints' `r(free)'
					constraint `r(free)' `_constraint'
				}

				cap `=cond("`quietly'"=="", "noisily", "")' `=cond("`cmd'"=="tobit", "version 15:", "")' `anything' if `hold' `=cond("`weight'"=="","",`"[`weight'`exp']"')', constraints(`_constraints') `from' `init' `iterate' `max' `options'

				local rc = _rc
				constraint drop `_constraints'
				if e(converged)==0 {
					di as err "Could not impose null."
					exit 430
				}
				if !`rc' {
					mat `b0' = e(b)
					cap `=cond("`cmd'"=="tobit", "version 15:", "")' `anything' if e(sample), `=cond(inlist("`cmd'","cmp","ml"),"init(`b0')","from(`b0',skip)")' iterate(0) `options' `max'
					local rc = _rc
				}
				if `rc' {
					di as err "Error imposing null. Perhaps {cmd:`cmd'} does not accept the {cmd:constraints()}, {cmd:from()}, and {cmd:iterate()} options, as needed."
					exit `rc'
				}
				if "`quietly'"=="" {
					mata st_numscalar("`t'", any(!diagonal(st_matrix("e(V)")) :& st_matrix("e(b)")'))
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
      if `partial' mat `R'  = `R', J(`df', `:word count `e(partial1)'' + e(partialcons), 0)  // add blanks for partialled-out
      mata st_matrix("`R'" , st_matrix("`R'" )[,_boottestp])  // put cols in standardized order & drop those for 0-variance vars

 			cap _estimates drop `hold'
			_est hold `hold', restore

			if `ar' {
				local kEnd: word count `Xnames_endog'
				local kEx : word count `Xnames_exog'
				mata st_numscalar("`p'", rows(st_matrix("`R'"))!=`kEnd' | (`kEx'? any(st_matrix("`R'")[|.,.\.,`kEx'|]) : 0))
				if `p' {
					di as err "For Anderson-Rubin test, null hypothesis must constrain all the instrumented variables and no others."
					exit 198
				}			
			}
		}

		return local seed = cond("`seed'"!="", "`seed'", "`c(seed)'")

		mata boottest_stata("`teststat'", "`df'", "`df_r'", "`p'", "`padj'", "`cimat'", "`plotmat'", "`peakmat'", `level', `ptolerance', ///
                        `ML', `LIML', 0`fuller', `K', `ar', `null', `scoreBS', "`weighttype'", "`ptype'", "`statistic'", ///
												"`madjust'", `N_h0s', "`Xnames_exog'", "`Xnames_endog'", 0`cons', ///
												"`Ynames'", "`b'", "`V'", "`ZExclnames'", "`hold'", "`scnames'", `hasrobust', "`allclustvars'", `:word count `bootcluster'', `Nclustvars', ///
												"`FEname'", 0`NFE', `FEdfadj', "`wtname'", "`wtype'", "`R1'", "`r1'", "`R'", "`r'", `reps', "`repsname'", "`repsFeasname'", `small', "`svmat'", "`dist'", ///
												"`gridmin'", "`gridmax'", "`gridpoints'", `matsizegb', "`quietly'"!="", "`b0'", "`V0'", "`svv'", "`NBootClustname'")
		_estimates unhold `hold'

    `quietly' if 2^`NBootClustname' < `reps' & inlist("`weighttype'", "rademacher", "mammen") {
      di _n "Warning: with " `NBootClustname' " boostrap clusters, the number of replications, `reps', exceeds the universe of " strproper("`weighttype'") " draws, 2^"`NBootClustname' " = " 2^`NBootClustname' ". " _c
      if "`weighttype'"=="rademacher" di "Sampling each once." _c
      di _n "Consider Webb weights instead, using {cmd:weight(webb)}."
    }
    local reps = `repsname'  // in case reduced to 2^G
    `quietly' if ((`NBootClustname'>12 | "`weighttype'" != "rademacher") & floor(`level'/100 * (`reps'+1)) != `level'/100 * (`reps'+1)) {
      di _n "Note: The bootstrap usually performs best when the confidence level (here, `level'%)" _n "      times the number of replications plus 1 (`reps'+1=" `reps'+1 ") is an integer."
    }

		if !`ML' & `reps' & `Nclustvars'>1 & `teststat'<. {
			return scalar repsFeas = `repsFeasname'
			if `repsFeasname' < `reps' di _n "Warning: " `reps' - `repsFeasname' " replications returned an infeasible test statistic and were deleted from the bootstrap distribution."
		}
		
		if `N_h0s'>1 local _h _`h'

		if `df'==1 {
			local tz = cond(`small', "t", "z")
			return scalar `tz'`_h' = `teststat'
		}

		di
		if `reps' di as txt strproper("`boottype'") " bootstrap-`statistic', null " cond(0`null', "", "not ") "imposed, " as txt `reps' as txt " replications, " _c
		di as txt cond(`ar', "Anderson-Rubin ", "") cond(!`reps' & `null' & "`boottype'"=="score", "Rao score (Lagrange multiplier)", "Wald") " test" _c
		if "`cluster'"!="" di ", clustering by " as res "`cluster'" _c
		if ("`bootcluster'"!="" | `Nclustvars' > 1) & `reps' di as txt ", bootstrap clustering by " as res "`bootcluster'" _c
		if `reps'	di as txt ", " strproper("`weighttype'") " weights" _c
		di as txt ":"
		
		if `margins' di `"  `:word `h' of "`marginsnames'""'
    else {
    	foreach c in `h0_`h'' {
        di "  " _c
        if !0`anythingh0' di as txt %4.0g `c' ": " _c
        di as res "`:constraint `c''"
      }
    }

		if `df'==1 {
			local line1 = cond(`small', "t(`=`df_r'')", "z")
			di _n as txt _col(`=33-strlen("`line1'")') "`line1' = " as res %10.4f `teststat' _n _col(`=33-strlen("``ptype''")') as txt "``ptype'' = " as res %10.4f `p'
			local `teststat' = `teststat' * `teststat'
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
			di _n as txt _col(`=27-strlen("`line1'")') "`line1' = " as res %10.4f `teststat' _n as txt _col(`=27-strlen("`line2'")') "`line2' = " as res %10.4f `p'
		}

		if "`madjust'" != "" {
			di _col(`=strlen("bonferroni")-strlen("`madjust'")+3') as txt strproper("`madjust'") "-adjusted prob = " as res %10.4f `padj'
			return scalar padj`_h' = `padj'
		}

		if `Nclustvars' > 1 {
			cap mat `t' = syminv(`V0')
			cap noi if diag0cnt(`t') di _n "Test statistic could not be computed. The multiway-clustered variance estimate " cond(`df'==1, "", "matrix ") "is not positive" cond(`df'==1, "", "-definite") "."
		}
		cap return matrix b`_h' = `b0'
		cap return matrix V`_h' = `V0'

		cap confirm mat `plotmat'
		if _rc == 0 {
			tempvar X1 Y _plotmat
			if `df'==2 tempvar X2
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
				mat colnames `_plotmat' = `X1' `X2' `Y'
				local _N = _N
        qui svmat `_plotmat', names(col)
				label var `Y' " " // in case user turns on legend
				if `"`graphname'"'=="Graph" cap graph drop Graph`_h'

				if `df'==1 {
        	cap mata st_local("nonmiss", strofreal(nonmissing(st_matrix("`cimat'"))))
          if 0`nonmiss' > 1 {
            mata _boottestm = (-min(-st_matrix("`cimat'")) - min(st_matrix("`cimat'"))) / 2 / (`nonmiss' + 1)  // margin
            mata _boottestt = min(select(st_matrix("`plotmat'")[,1], st_matrix("`plotmat'")[,1] :> max(st_matrix("`cimat'")) + _boottestm)); st_local("plotmax", rows(_boottestt)? strofreal(_boottestt * (1 + sign(_boottestt)*1e-6)) : .)
            mata _boottestt = max(select(st_matrix("`plotmat'")[,1], st_matrix("`plotmat'")[,1] :< min(st_matrix("`cimat'")) - _boottestm)); st_local("plotmin", rows(_boottestt)? strofreal(_boottestt * (1 - sign(_boottestt)*1e-6)) : .)
          }
          else {
          	local plotmin .
            local plotmax .
          }
          
					format `X1' `X2' %5.0g

					local 0, `graphopt'
					syntax, [LPattern(passthru) LWidth(passthru) LColor(passthru) LStyle(passthru) *]

					line `Y' `X1', sort(`X1') `lpattern' `lwidth' `lcolor' `lstyle' || scatter `Y' `X1' if _n>rowsof(`plotmat'), mlabel(`X1') mlabpos(6) xtitle("`constraintLHS1'") ///
            || if (`gridmin'<. | -`X1'<=-`plotmin') & (`gridmax'<. | `X1'<=`plotmax'), ///
						ytitle(`"`=strproper("`madjust'")+ cond("`madjust'"!="","-adjusted ", "")' p value"') ///
						yscale(range(0 .)) name(`graphname'`_h', `replace') ///
						`=cond(`level'<100, "yline(`=1-`level'/100') ylabel(#6 `=1-`level'/100')", "")' legend(off) `options'
				}
				else {
					foreach var in Y X1 X2 {
						if "``var''" != "" {
							cap drop _boottest``var''
              ren ``var'' _boottest``var'' // work-around for weird bug in twoway contour, rejecting temp vars
							local `var' _boottest``var''
						}
					}
					c_local contourvars `Y' `X1' `X2' // pass these to calling program for dropping in case twoway contour crashes
					
					twoway contour `Y' `X2' `X1' if _n<=rowsof(`plotmat'), xtitle("`constraintLHS1'") ytitle("`constraintLHS2'") ///
						ztitle(`"`=strproper("`madjust'")+ cond("`madjust'"!="","-adjusted ", "")' p value"', orient(rvert)) ///
						name(`graphname'`_h', `replace') crule(linear) scolor(gs5) ecolor(white) ccut(0(.05)1) plotregion(margin(zero)) /// // defaults copied from weakiv
						`graphopt'
				}
        qui keep in 1/`_N'
			}
		}

		cap confirm mat `cimat'
		if _rc==0 {
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
		
		if "`statistic'" == "c" & `df'==1 {
			di "Note: denominator for `tz' statistic computed from the bootstrap replications of the numerator."
		}
		
		return scalar `=cond(`small', "F", "chi2")'`_h' = cond(`df'==1, `teststat'*`teststat', `teststat')
		if `small' return scalar df_r`_h' = `df_r'
		return scalar df`_h' = `df'
		return scalar p`_h' = `p'
		if `hasrobust' return local robust robust
		if "`svmat'"!="" return matrix dist`_h' = `dist'
    if "`svv'" != "" return matrix v`_h' = `svv'
	}  // loop over independent H0s

  cap mat_put_rr `C'  // can fail in boottest, margins
	return scalar level = `level'
	return scalar ptol = `ptolerance'
	return local statistic `statistic'
	return local weighttype `weighttype'
	return local boottype `boottype'
	return local clustvars `clustvars'
	return scalar null = `null'
	return scalar reps = `repsname'
end

* Version history
* 3.2.2 Add didregress, xtdidregress support. After xtXXX estimation, emulate those commands by in counting FE in dof adjustment, unless "xtreg, dfadj"
* 3.2.1 Prevent it from expanding data set when number of points in graph exceed # of rows in data set
* 3.2.0 Added margins option
* 3.1.4 Fixed crashes with matsizegb() and with "granular" (many-clustered) FE estimates with FE & cluster groups coherent
* 3.1.3 Fixed crash on stat(c) after OLS
* 3.1.2 Incorporated small-sample factor in r(dist)
* 3.1.1 Minor bug fixes and speed-ups
* 3.1.0 Complete overhaul of WRE for ~200X speed gain. Dropped GMM support. Added support for ivreg2's partial().
* 3.0.2 Dropped "KK" calculation (last expression in eq 60 in paper) because inefficient when interpolating. Refined plotting to minimize interpolation anchor resets. Refined criterion to use "granular"-optimized code (many small clusters).
* 3.0.1 Recompiled in Stata 13
* 3.0.0 Exploit linearity/quadratic form in denominators too. ~10X speed-up over 2018 version for inverting tests after OLS.
* 2.8.1 Fix 2.8.0 bugs. Even more fully exploit linearity in test stat numerators, for ~2X speed gain in inverting test after OLS.
* 2.8.0 More fully exploit linearity of numerators and K matrices as in appendix A.2. Recenter classical score test even for non-OLS, reversing 2.4.1 change. Added ptolerance() option.
* 2.7.4 Fixed errors affecting results after GMM
* 2.7.3 Prevented crash on WRE bootstrap-t with svmat(numer)
* 2.7.2 Added error check for absorbed interaction terms in (iv)reghdfe
* 2.7.1 Fixed crash after IV without clustering. Added check for regress with undocumented IV syntax. Fixed crash on 2-D test with "nonull" but not "nograph".
* 2.7.0 Require Stata 13 or later, so no need to work around possible lack of panelsum()
* 2.6.0 Added svv option.
* 2.5.7 If lusolve() fails on non-invertible matrix, resort to invsym() for generalized inverse.
* 2.5.6 Fixed 2.4.0 bug: look in e(df_a_initial) rather than e(df_a). Matters if clustering on the absorbed var.
* 2.5.5 Fixed wrong results when absorbed variable is type string
* 2.5.4 Added support for ivreghdfe
* 2.5.3 Fixed crash in score test (including waldtest) after "robust" estimation without observation weights
* 2.5.2 More graceful handling of degenerate cases: multiway t stat = .; test hypothesis refers to dropped/constrained variable
* 2.5.1 Fixed 2.5.0 bug after "robust" estimation
* 2.5.0 Added bootstrap-c. Fixed bug affecting results from WCU when using Mammen/Rademacher/Webb and Rademacher universe not exhausted; doubles test stat in published Conley & Tabor reanalysis.
* 2.4.3 minor bug fixes and edits
* 2.4.2 Fixed 2.4.1 bug. Added r(b) and r(V) return values.
* 2.4.1 Optimized classical tests; removed bug in score test after FE est (wrongly droped term 2 in (63) in non-classical use of score test)
* 2.4.0 After reghdfe look for FE count in e(df_a) as well as e(K1)
* 2.3.9 Prevented crash in pure "robust" non-WRE
* 2.3.8 Prevented crash when it can't recompile boottest.mata; instead issues an explanatory warning
* 2.3.7 Fixed crash after tobit estimation in Stata 15
* 2.3.6 Fixed crash in score test/bootstrap with multiple independent hypotheses
* 2.3.5 Fixed stupid 2.3.4 crash
* 2.3.4 Dropped "Rejection" from axis labels. Added check for right number of entries in gridmin(), gridmax(), gridpoints().
* 2.3.3 Eliminated false warning that neg Hessian not pos def when a parameter is constrained to 0 in model
* 2.3.2 Fixed 2.2.0 crash when errors are non-robust
* 2.3.1 Fixed 2.2.0 bug--left behind temporary variables
* 2.3.0 Removed optimization hacks from WRE code because they created matrices with 1 row per obs and 1 col per replication
* 2.2.2 Allowed quietly option in ado interface to suppress dots. Made sorts in Mata code stable.
*       For LIML, reverted to finding eigenvalue of I-TT/TPZT instead of TT/TPZT; seems to avoid instances of eigensystem() returning all missing
* 2.2.1 Fixed failure to detect # of FE after areg in Stata version < 15
* 2.2.0 Added contour plotting for 2-D tests.
* 2.1.9 Work-around for Stata crash when number of fixed effects is very large: require # of FE as input, and don't represent them as linked list.
* 2.1.8 Fixed 2.1.6 crash with FE. Fixed parsing of matsizegb() option.
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
