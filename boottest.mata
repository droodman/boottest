*! boottest 1.6.1 10 July 2017
*! Copyright (C) 2015-17 David Roodman

* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.

mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

string scalar boottestStataVersion() return("`c(stata_version)'")
string scalar      boottestVersion() return("01.05.00")

struct smatrix {
	real matrix M
}

struct boottest_clust {
	real scalar N, multiplier
	real rowvector cols
	real colvector ClustShare, order
	real matrix info
}

class AnalyticalModel { // class for analyitcal OLS, 2SLS, LIML, GMM estimation--everything but iterative ML
	real scalar LIML, YY, ee, eec, Fuller, AR, K, k
	real matrix ZX, ZXEnd, ZZ, H_2SLS, invH, A, VR0, ZVR0, e2, numer, Splus, pi, XXEnd, dbetads, Ze2
	real colvector sAll, e, beta, beta0, Ze, ZY
	real rowvector YXEnd
	pointer(real colvector) scalar pY
	pointer(real matrix) scalar pXEnd, pXX, pXY, pZExclY, pXExY, pZExclXEnd, pV, pW, pinvZZ, pXExXEx, pZXEx, pXExXEnd, pH, pR0, pS
	pointer (class boottestModel scalar) scalar parent
	pointer (class AnalyticalModel scalar) scalar DGP
	
	void new(), InitExog(), InitEndog(), InitTestDenoms(), SetDGP(), SetS(), InitEstimate(), Estimate(), setParent(), SetLIMLFullerK(), SetAR()
}

class boottestModel {
	real scalar scoreBS, reps, small, wildtype, null, dirty, initialized, Neq, ML, Nobs, _Nobs, k, kEx, el, sumwt, Nclust, robust, weights, REst, multiplier, quietly, sqrt, cons, LIML, Fuller, K, IV, WRE, WREnonAR, ptype, twotailed, gridstart, gridstop, gridpoints, df, df_r, AR, d, cuepoint, willplot, NumH0s, p
	pointer (real matrix) scalar pZExcl, pR, pR0, pID, pXEnd, _pXEnd, pXEx
	pointer (real colvector) scalar pr, pr0, pY, pSc, pwt, pW, pV
	real matrix numer, u, U, S, SAR, SAll, LAll_invRAllLAll, plot, CI
	string scalar wttype, madjtype
	real colvector Dist, DistCDR, s, sAR, plotX, plotY, sAll, beta
	real rowvector peak
	struct boottest_clust colvector clust
	class AnalyticalModel scalar M_DGP
	pointer (class AnalyticalModel scalar) scalar pM_Repl, pM

	void new(), set_dirty(), set_sqrt(), boottest(), make_DistCDR(), plot()
	real scalar r0_to_p(), search(), get_p(), get_padj(), get_stat(), get_df(), get_df_r()
	real matrix combs(), count_binary()
	real colvector get_dist()
}

void AnalyticalModel::new()
	AR = 0

void AnalyticalModel::setParent(class boottestModel scalar B)
	parent = &B

void AnalyticalModel::SetLIMLFullerK(real scalar _LIML, real scalar _Fuller, real scalar _K) {
	LIML = _LIML
	Fuller = _Fuller
	K = _K
}

void AnalyticalModel::SetAR(real scalar _AR) {
	if (AR = _AR) LIML = Fuller = K = 0
}

// stuff that can be done before r0 set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
void AnalyticalModel::InitExog() {
	real matrix ZExclXEx

	pXExXEx = &cross(*parent->pXEx, *parent->pwt, *parent->pXEx)
	if (cols(*parent->pZExcl)) { // GMM, 2SLS, LIML
		ZExclXEx = cross(*parent->pZExcl, *parent->pwt, *parent->pXEx)
		pZXEx = &(*pXExXEx \ ZExclXEx)
		if (parent->IV)
			pinvZZ = &invsym((ZZ = *pZXEx, (ZExclXEx' \ cross(*parent->pZExcl, *parent->pwt, *parent->pZExcl))))
	} else
		pXX = pXExXEx
}

void AnalyticalModel::SetDGP(class AnalyticalModel scalar _DGP) {
	DGP = &_DGP
}
void AnalyticalModel::SetS(real matrix S) {
	pS = &S // DGP==NULL means this is the DGP
	if (LIML) Splus = blockdiag(S, 1)  // add an entry to S for the dep var
}
// stuff that can be done before S & r0 set, but depend on endogenous variables, which are bootstrapped in WRE
void AnalyticalModel::InitEndog(pointer (real colvector) scalar _pY, pointer (real matrix) scalar _pXEnd, | ///
		pointer (real colvector) scalar _pZExclY, pointer (real rowvector) scalar _pXExY, real scalar _YY, pointer (real matrix) scalar _pZExclXEnd, pointer (real matrix) scalar _pXExXEnd) {

	pY = _pY; pXEnd = _pXEnd

	pXExY = _pXExY==NULL? &cross(*parent->pXEx, *parent->pwt, *parent->pY) : _pXExY
	if (K | AR)
		pZExclY = _pZExclY==NULL? &cross(*parent->pZExcl, *parent->pwt, *pY) : _pZExclY
	if (K) {
		pXExXEnd   = _pXExXEnd  ==NULL? &cross(*parent->pXEx  , *parent->pwt, *pXEnd) : _pXExXEnd
		pZExclXEnd = _pZExclXEnd==NULL? &cross(*parent->pZExcl, *parent->pwt, *pXEnd) : _pZExclXEnd
		ZXEnd = *pXExXEnd \ *pZExclXEnd
		ZX = *pZXEx, ZXEnd
		XXEnd = *pXExXEnd \ cross(*pXEnd, *parent->pwt, *pXEnd)
		pXX = &(*pXExXEx,  *pXExXEnd \ XXEnd')
		YXEnd = cross(*pY, *parent->pwt, *pXEnd)
		ZY = *pXExY \ *pZExclY
		pXY = &(*pXExY \ YXEnd')
		if (LIML | !(parent->robust | parent->scoreBS))
			YY = _YY==.? cross(*pY, *parent->pwt, *pY) : _YY

		if (parent->IV) // if GMM weight matrix not provided, prepare 2SLS one
			A = (I(parent->kEx) \ J(parent->el-parent->kEx, parent->kEx, 0)), *pinvZZ * ZXEnd // 2SLS is (A' ZX)^-1 * (A'ZY). Also apparently used in k-class and LIML robust VCV by Stata convention
		else
			A = *parent->pW * ZX
		H_2SLS = A ' ZX // Hessian
	} else { // OLS / AR
		pXY = pXExY
		if (AR) {
			pXY = &(*pXY \ *pZExclY)
			pXX = &ZZ
		}
	}
	k = cols(*pXX)
}

// stuff that can be done before r0 set but depends on S and endogenous variables
void AnalyticalModel::InitEstimate() {
	real rowvector val
	real matrix _ZY, TT, TPZT, vec
	pointer (real matrix) scalar pbetadenom
	pragma unset vec; pragma unset val

	if (LIML)
		if (parent->el == k) // exactly identified LIML = 2SLS
			K = 1
		else {
			_ZY = ZXEnd, ZY
			TT = *pXX, (*pXExY \ YXEnd') \ *pXExY', YXEnd, YY  // where T = *pXEx, *pXEnd, *pY
			(TPZT = TT)[|parent->kEx+1,parent->kEx+1\.,.|] = _ZY ' (*pinvZZ) * _ZY

			A = *pinvZZ * ZX
			H_2SLS = A ' ZX // Hessian

			if (rows(*pS)) { // includes H0 if pSAll really points to SAll rather than S
				TT = Splus ' TT * Splus
				TPZT = Splus ' TPZT * Splus
			}
			eigensystemselecti( I(rows(TT)) - invsym(TT) * TPZT, 1\1, vec, val)
			K = 1/Re(val) - Fuller / (parent->_Nobs - parent->el)   // sometimes a tiny imaginary component sneaks in
		}

	pH = K? (K==1? &H_2SLS : &((1-K)* *pXX + K*H_2SLS)) : pXX

	if (rows(*pS)) {
		pbetadenom = &(*pS * invsym(*pS ' (*pH) * *pS) * *pS')
		invH = J(0,0,0)
	} else
		pbetadenom = &(invH = invsym(*pH))

	if (K)
		if (K==1) { // 2SLS
			beta0 = *pbetadenom * A ' ZY
			dbetads = I(rows(beta0)) - *pbetadenom * A ' ZX
		} else { // k-class, LIML
			beta0 = *pbetadenom * (K * A ' ZY + (1-K) * *pXY)
			dbetads = I(rows(beta0)) -  *pbetadenom * (K * A ' ZX + (1-K) * *pXX)
		}
	else { // OLS / AR
		beta0 = *pbetadenom * *pXY
		dbetads = I(rows(beta0)) - *pbetadenom * *pXX
	}
}

// stuff that doesn't depend on r0, for test stat denominators in replication regressions
void AnalyticalModel::InitTestDenoms(real matrix S) {
	real matrix AVR0
	
	if (rows(S))
		pV = &(S * invsym(S ' (*pH) * S) * S')
	else
		pV = &(rows(invH)? invH : invsym(*pH))
	VR0 = *pV * *parent->pR0'

	if (parent->scoreBS | (parent->robust & !(parent->WREnonAR & parent->Nclust==1))) {
		if (K) {
			AVR0 = A * VR0
			ZVR0 = *parent->pZExcl * AVR0[|parent->kEx+1,.\.,.|]; if (parent->kEx) ZVR0 = ZVR0 + *parent->pXEx * AVR0[|.,.\parent->kEx,.|]
		} else if (AR) {
			ZVR0 = *parent->pZExcl * VR0[|cols(*parent->pXEx)+1,.\.,.|]
			if (cols(*parent->pXEx))
				ZVR0 = ZVR0 + *parent->pXEx * VR0[|.,.\cols(*parent->pXEx),.|]
		} else
			ZVR0 = *parent->pXEx * VR0
	}
}

// stuff that depends on r0 and endogenous variables: compute beta and residuals
void AnalyticalModel::Estimate(real colvector s) {
	real matrix invZMeZ; real colvector negZeinvee

	if (parent->null | !rows(beta)) { // don't need to recompute if we're not actually imposing the null
		beta = rows(s)? beta0 + dbetads * s : beta0

		if (AR) {
			e = *pY - *parent->pZExcl * beta[|cols(*parent->pXEx)+1\.|]
			if (cols(*parent->pXEx))
				e =  e - *parent->pXEx * beta[|.\cols(*parent->pXEx)|]
		} else if (parent->IV)
			if (parent->kEx == k)
				e = *pY - *pXEnd * beta[|parent->kEx+1\.|]
			else
				e = *pY - *pXEnd * beta[|parent->kEx+1\.|] - *parent->pXEx * beta[|.\parent->kEx|]
		else
				e = *pY                                    - *parent->pXEx * beta[|.\parent->kEx|]

		if (!(parent->robust | parent->scoreBS) | (DGP==NULL & LIML)) // useful in non-robust, residual-based bootstrap, and in computing e2 in LIML (just below)
			ee = YY - 2 * *pXY ' beta + beta ' (*pXX) * beta
		if (!(parent->robust | parent->scoreBS))
			eec = parent->cons? ee : ee - (parent->weights? cross(e, *parent->pwt) : sum(e))^2 / parent->_Nobs // sum of squares after centering, N * Var

		if (DGP==NULL & LIML) {
			Ze = ZY - ZX * beta
			negZeinvee = Ze / -ee
			invZMeZ = invsym(ZZ + negZeinvee * Ze')
			pi = (invZMeZ * negZeinvee) * (YXEnd - beta ' XXEnd) + invZMeZ * ZXEnd  // coefficients in reduced-form equations; Davidson & MacKinnon (2010), eq 15

			e2 = *pXEnd - *parent->pZExcl * pi[|parent->kEx+1,.\.,.|]; if (parent->kEx) e2 = e2 - *parent->pXEx * pi[|.,.\parent->kEx,.|]
			if (parent->AR)
				Ze2 = ZXEnd - ZZ * pi
		}
	}
}

void boottestModel::new() {
	AR = LIML = Fuller = WRE = small = scoreBS = wildtype = Neq = ML = initialized = quietly = sqrt = cons = IV = ptype = robust = willplot = 0
	twotailed = null = dirty = 1
	cuepoint = .
	pXEnd = pXEx = pZExcl = pY = pSc = pID = pR = pR0 = pwt = &J(0,0,0)
	pr = pr0 = &J(0,1,0)
}

void boottestModel::set_dirty(real scalar _dirty) {
	dirty = _dirty
	if (_dirty)
		initialized = 0
}
void boottestModel::set_sqrt(real scalar _sqrt) {
	if (_sqrt < sqrt)
		if (!dirty) Dist = Dist :* Dist
	else
		set_dirty(1)
	sqrt = _sqrt
}
void boottest_set_dirty   (class boottestModel scalar M                      ) {
	 M.set_dirty(1)
}
void boottest_set_ptype(class boottestModel scalar M, string scalar ptype)     {
	real scalar p
	p = cross( (strtrim(strlower(ptype)) :== ("symmetric"\"equaltail"\"lower"\"upper")), 1::4 ) - 1
	if (p<0) 
		_error(198, `"p-value type must be "symmetric", "equaltail", "lower", or "upper.""')
	M.ptype = p
	M.twotailed = p<=1
}
void boottest_set_XEnd    (class boottestModel scalar M, real matrix X       ) {
	M.pXEnd  = &X; M.set_dirty(1)
}
void boottest_set_XEx    (class boottestModel scalar M, real matrix X        ) {
	M.pXEx  = &X; M.set_dirty(1)
}
void boottest_set_Y       (class boottestModel scalar M, real matrix Y       ) {
	M.pY  = &Y; M.set_dirty(1)
}
void boottest_set_ZExcl   (class boottestModel scalar M, real matrix Z       ) {
	M.pZExcl  = &Z; M.set_dirty(1)
}
void boottest_set_wt       (class boottestModel scalar M, real matrix wt     ) {
	M.pwt  = &wt; M.set_dirty(1)
}
void boottest_set_sc      (class boottestModel scalar M, real matrix Sc) {
	M.pSc  = &Sc
	M.set_dirty(1)
}
void boottest_set_ML      (class boottestModel scalar M, real scalar ML) {
	M.ML  = ML; M.set_dirty(1)
	if (ML) boottest_set_scoreBS(M, 1)
}
void boottest_set_LIML    (class boottestModel scalar M, real scalar LIML) {
	M.LIML = LIML; M.set_dirty(1)
}
void boottest_set_AR    (class boottestModel scalar M, real scalar AR) {
	M.AR = AR; M.set_dirty(1)
}
void boottest_set_Fuller    (class boottestModel scalar M, real scalar Fuller) {
	M.Fuller = Fuller; M.set_dirty(1)
}
void boottest_set_k    (class boottestModel scalar M, real scalar K) {
	M.K = K; M.set_dirty(1)
}
void boottest_set_quietly (class boottestModel scalar M, real scalar quietly )
	M.quietly = quietly
void boottest_set_beta    (class boottestModel scalar M, real colvector beta) {
	M.beta = beta; M.set_dirty(1)
}
void boottest_set_V    (class boottestModel scalar M, real matrix V ) {
	M.pV = &V; M.set_dirty(1)
}
void boottest_set_W    (class boottestModel scalar M, real matrix W ) {
	M.pW = &W; M.set_dirty(1)
}
void boottest_set_small   (class boottestModel scalar M, real scalar small   ) {
	M.small = small; M.set_dirty(1)
}
void boottest_set_cons(class boottestModel scalar M, real scalar cons  ) {
	M.cons = cons; M.set_dirty(1)
}
void boottest_set_scoreBS (class boottestModel scalar M, real scalar scoreBS ) {
	M.scoreBS = scoreBS; M.set_dirty(1)
}
void boottest_set_reps    (class boottestModel scalar M, real scalar reps    ) {
	M.reps = reps; M.set_dirty(1)
}
void boottest_set_null    (class boottestModel scalar M, real scalar null    ) {
	M.null = null; M.set_dirty(1)
}
void boottest_set_wttype  (class boottestModel scalar M, string scalar wttype) {
	M.wttype = wttype; M.set_dirty(1)
}
void boottest_set_ID      (class boottestModel scalar M, real matrix ID      ) {
	M.pID = &ID; M.set_dirty(1)
	if (cols(ID)) M.robust = 1
}
void boottest_set_robust  (class boottestModel scalar M, real scalar robust  ) {
	M.robust = robust
	if (robust==0) boottest_set_ID(M, J(0,0,0))
	M.set_dirty(1)
}
void boottest_set_R (class boottestModel scalar M, real matrix R , real colvector r ) {
	M.pR = &R; 	M.pr  = &r; M.set_dirty(1)
}
void boottest_set_R0(class boottestModel scalar M, real matrix R0, real colvector r0) {
	M.pR0 = &R0; M.pr0 = &r0; M.set_dirty(1)
}
void boottest_set_willplot(class boottestModel scalar M, real scalar willplot) {
	M.willplot = willplot
}
void boottest_set_grid(class boottestModel scalar M, real scalar _gridstart, real scalar _gridstop, real scalar _gridpoints) {
	M.gridstart = _gridstart; M.gridstop = _gridstop; M.gridpoints = _gridpoints
}
void boottest_set_madjust(class boottestModel scalar M, string scalar madjtype, real scalar NumH0s) {
	M.madjtype = strlower(madjtype)
	M.NumH0s = NumH0s
	if (M.madjtype != "bonferroni" & M.madjtype != "sidak" & M.madjtype != "")
		_error(198, `"Multiple-hypothesis adjustment type must be "Bonferroni" or "Sidak"."')
}
void boottest_set_wildtype(class boottestModel scalar M, string scalar wildtype) {
	wildtype = strlower(wildtype)
	if (.==(M.wildtype = wildtype=="rademacher" ? 0 : (wildtype=="mammen" ? 1 : (wildtype=="webb" ? 2 : (wildtype=="normal" ? 3 : .)))))
		_error(198, `"Wild type must be "Rademacher" or "Mammen" or "Webb" or "Normal"."')
	M.set_dirty(1)
}

real colvector boottestModel::get_dist() {
	if (dirty) boottest()
	make_DistCDR()
	return(DistCDR)
}
void boottestModel::make_DistCDR() {
	if (!rows(DistCDR))
		if (rows(Dist)>1)
			_sort( DistCDR=Dist[|2\.|] , 1)
		else
			DistCDR = J(0,1,0)
}

// Ties count half. Robust to missing bootstrapped values being interpreted as +infinity.
real scalar boottestModel::get_p(|real scalar analytical) {
	real scalar t; real colvector _Dist
	if (dirty) boottest()
	t = Dist[1]
	if (reps & analytical==.) {
		if (t == .) return (.)
		if (sqrt & ptype != 3) {
			if (ptype==0) { // symmetric p value
				_Dist = abs(Dist); t = abs(t)
				p = 1 -(       colsum(t:>_Dist)                       + (colsum(t:==_Dist) - 1)*.5) / (colnonmissing(Dist) - 1)
			} else if (ptype==1) // equal-tail p value
				p =    (2*min((colsum(t:> Dist) , colsum(-t:>-Dist))) + (colsum(t:== Dist) - 1)   ) / (colnonmissing(Dist) - 1)
			else // upper-tailed p value
				p = 1 -(       colsum(t:< Dist)                       + (colsum(t:==_Dist) - 1)*.5) / (colnonmissing(Dist) - 1)
		} else // upper-tailed p value or p value based on squared stats 
				p = 1 -(       colsum(t:> Dist)                       + (colsum(t:== Dist) - 1)*.5) / (colnonmissing(Dist) - 1)
	} else {
		p = small? Ftail(df, df_r, sqrt? t*t : t) : chi2tail(df, sqrt? t*t : t)
		if (sqrt & !twotailed) {
			p = p / 2
			if ((ptype==3) == (t<0))
				p = 1 - p
		}
	}
	return (p)
}

real scalar boottestModel::get_padj(|real scalar analytical) {
	(void) get_p(analytical)
	if (madjtype=="bonferroni") return(min((1, NumH0s*p)))
	if (madjtype=="sidak"     ) return(1 - (1 - p)^NumH0s)
	return (p)
}

real scalar boottestModel::get_stat() {
	if (dirty) boottest()
	return (Dist[1])
}

real scalar boottestModel::get_df() {
	if (dirty) boottest()
	return (df)
}

real scalar boottestModel::get_df_r() {
	if (dirty) boottest()
	return (df_r)
}

void _boottest_st_view(real matrix V, real scalar i, string rowvector j, string scalar selectvar) {
	if (favorspeed() | 1) {
		V = length(tokens(j))? st_data(i, j, selectvar) : st_data(i, J(1,0,0), selectvar)
	} else
		st_view(V, i, j, selectvar)
}






















































// main routine
void boottestModel::boottest() {
	real colvector rAll, numer_i, _e, eUZVR0wt, ID1, Ystar, _beta, betaEnd
	real rowvector val, YstarYstar
	real matrix betadevEx, betadevEnd, betanumer, RAll, L, LAll, vec, denom, combs, XExZVR0wt, XEndZVR0wt, ZVR0wt, t, ZExclYstar, XExYstar, Subscripts, Zi, AVR0, betadenom, eZVR0, SeuZVR0, SeZVR0, eu, VR0
	real scalar i, j, c
	pointer (real matrix) scalar _pR0, pewt, pSeZVR0, pXEndstar, pXExXEndstar, pZExclXEndstar, pu, pVR0
	pointer (real colvector) scalar peZVR0wt
	struct smatrix colvector denoms, XEndstar, XExXEndstar, ZExclXEndstar, XZi, eZi
	class AnalyticalModel scalar M_WRE
	pragma unset vec; pragma unset val; pragma unset denom; pragma unset ID1; pragma unused M_WRE

	if (!initialized) {  // for efficiency when varying r0 repeatedly to make CI, do stuff once that doesn't depend on r0
		kEx = cols(*pXEx)
		Nobs = rows(*pXEx)
		if (!cols(*pZExcl)) pZExcl = &J(Nobs,0,0)
		if (!cols(*pXEnd)) pXEnd = &J(Nobs,0,0)
		d = cols(*pXEnd) + 1
		k  = cols(*pR0)
		REst = rows(*pR) // base model contains restrictions?
		if (pZExcl != NULL) el = cols(*pZExcl) + kEx
		if (K==.) K = cols(*pZExcl)>0
		IV = K & pW==NULL
		WRE = (IV & !scoreBS) | AR
		WREnonAR = WRE & !AR

		if (weights = rows(*pwt)>1)
			sumwt = sum(*pwt)
		else
			pwt = &(sumwt = 1)
		_Nobs = weights & wttype=="fweight"? sumwt : Nobs

		if (Nclust = cols(*pID)) {
			combs = combs(1..Nclust)
			clust = boottest_clust(rows(combs)-1) // leave out no-cluster combination
			for (c=1; c<=length(clust); c++) {
				clust[c].cols         = boottest_selectindex(combs[c,]:<.)

				if (c > 1) // if this sort ordering same as last to some point and missing thereafter, no need to resort
				  if (!allof( combs[|c, min(boottest_selectindex(combs[c,] :!= combs[c-1,])) \ c,.|], .))
						_collate(*pID, clust[c].order = order(*pID, clust[c].cols))

				clust[c].info         = c==1 & (Nclust > 1 | WREnonAR)? _panelsetup(*pID, clust[c].cols, ID1) : _panelsetup(*pID, clust[c].cols)  // in some cases, save originally ordered ID marker for bootstrapping cluster var
				clust[c].N            = rows(clust[c].info)
				clust[c].multiplier   = mod(cols(clust[c].cols),2)? 1 : -1
				if (small)
				  clust[c].multiplier = clust[c].multiplier * clust[c].N/(clust[c].N-1)
			}
		} else {
			clust = boottest_clust()
			clust.multiplier = small? _Nobs / (_Nobs - 1) : 1
			clust.N = Nobs
		}
		
		if (scoreBS)
			if (Nclust)
				for (c=length(clust); c; c--)
					clust[c].ClustShare = weights? _panelsum(*pwt, clust[c].info)/sumwt : (clust[c].info[,2]-clust[c].info[,1]:+ 1)/Nobs // share of observations by group 
			else
				clust.ClustShare = weights? *pwt/sumwt : 1/_Nobs

		if (reps & wildtype==0 & clust.N*ln(2) < ln(reps)+1e-6) {
			if (!quietly) printf("\nWarning: with %g clusters, number of replications, %g, exceeds the universe of Rademacher draws, 2^%g = %g. Sampling each once. \nConsider Webb weights instead, using {cmd:weight(webb)}.\n", clust.N, reps, clust.N, 2^clust.N)
			u = J(clust.N,1,1), count_binary(clust.N, -1-WREnonAR, 1-WREnonAR) // complete Rademacher set
		} else {
			if (wildtype==3)
				u = rnormal(clust.N, reps+1, -WREnonAR, 1) // normal weights
			else if (wildtype==2) {
				u = rdiscrete(clust.N, reps+1, (1\1\1\0\1\1\1)/6) * .5 :- 2
				u = sqrt(abs(u)) :* sign(u); if (WREnonAR) u = u :- 1 // Webb weights
			}	else if (wildtype) {
				u = ( rdiscrete(clust.N, reps+1,(.5+sqrt(.05)\.5-sqrt(.05))) :- 1.5 ) * sqrt(5) :+ (.5 - WREnonAR) // Mammen
				if (!quietly & clust.N*ln(2) < ln(reps)+1e-6) printf("\nWarning: with %g clusters, number of replications, %g, exceeds the universe of Mammen draws, 2^%g = %g. \nConsider Webb weights instead, using {cmd:weight(webb)}.\n", clust.N, reps, clust.N, 2^clust.N) 
			}	else {
				u = runiform(clust.N, reps+1) :>= .5; u = u + u :- (1 + WREnonAR) // Rademacher
			}

			u[,1] = J(clust.N, 1, 1-WREnonAR)  // keep original residuals in first entry to compute base model stat
		}

		if (Nclust + WREnonAR > 1) U = u[ID1,] // for multi-way clustering and clustered WRE, also explode to one row per observation instead of per bootstrapping cluster

		if (!ML) {
			if (REst) {
				symeigensystem(*pR ' invsym(*pR * *pR') * (*pR), vec, val) // make "inverse" S,s of constraint matrices; formulas adapted from [P] makecns
				L = vec[|.,.\.,rows(*pR)|] // eigenvectors not in kernel of projection onto R
				S = vec[|.,rows(*pR)+1\.,.|] // eigenvectors in kernel
				s = L * luinv(*pR * L) * *pr
				if (AR) {
					SAR = blockdiag(S[|.,. \ rows(S)-cols(*pXEnd) , cols(S)-cols(*pXEnd)|] , I(cols(*pZExcl))) // adapt S,s from XExog, XEndog to XExog, ZEXcl. Assumes no constraints link XExog and XEndog
					sAR = s[|.\rows(s)-cols(*pXEnd)|] \ J(cols(*pZExcl),1,0)
				}
			}

			// Estimation with null imposed along with any model constraints; in IV, Z is unconstrained regardless of overlap with potentially constrained X
			_pR0 = null? pR0 : &J(0, k, 0)
			RAll = REst? *pR \ *_pR0 : *_pR0 // combine model and hypothesis constraints to prepare to "invert" them as a group too
			if (rows(RAll)) {
				LAll = invsym(RAll * RAll')
				if (!all(diagonal(LAll)))
					_error(111, "A null hypothesis constraint is inconsistent or redundant.")
				symeigensystem(RAll ' LAll * RAll, vec, val)
				LAll  = vec[|.,. \ .,rows(RAll)|]
				SAll = rows(RAll) < cols(vec)? vec[|.,rows(RAll)+1 \ .,.|] : J(rows(vec), 0, 0)
				LAll_invRAllLAll = LAll * luinv(RAll * LAll)
			} else
				SAll = J(0,0,0)

			M_DGP.setParent(this)
			M_DGP.InitExog()

			if (WRE) {
				pM_Repl = &(M_WRE = M_DGP)
				pM_Repl->SetDGP(M_DGP)
				  M_DGP.SetLIMLFullerK(1, 0, 1)
				if (!AR) pM_Repl->SetLIMLFullerK(LIML, Fuller, K)
				pM_Repl->SetAR(AR)
				pM_Repl->SetS(AR? SAR : S)
			} else
				M_DGP.SetLIMLFullerK(LIML, Fuller, K)

			M_DGP.InitEndog(pY, pXEnd)

			if (AR) {
				if (willplot & rows(*pR0)==1) { // for plotting purposes get original point estimate if not normally generated
					M_DGP.SetS(S) // no-null model in DGP
					M_DGP.InitEstimate()
					M_DGP.Estimate(s)
					cuepoint = *pR0 * M_DGP.beta - *pr0 // not true CUE estimate unless classical errors, but serves same purpose as weakiv cuepoint option
				}
				pR0 = &(J(cols(*pZExcl),kEx,0), I(cols(*pZExcl))) // for AR test, picks out coefs on excluded exogenous variables
			}
			
			M_DGP.SetS(SAll) // (potentially) constrained model in DGP; SAll imposes constraints, pR0 tests hypotheses on results
			M_DGP.InitEstimate()

			if (AR) {
				k = el
				K = 0
				pM = pM_Repl
				_pXEnd = pZExcl // in AR test, instruments supplant instrumented on RHS
			} else {
				pM = &M_DGP
				_pXEnd = pXEnd
			}
		}

		df = rows(*pR0)
		df_r = Nclust? clust[Nclust].N - 1 : _Nobs - k

		if (df==1) set_sqrt(1) // work with t/z stats instead of F/chi2

		if (small)
			multiplier = (_Nobs - k) / (_Nobs - robust) / df // divide by # of constraints because F stat is so defined
		else
			multiplier = 1
		if (!(robust | ML))
			multiplier = multiplier * _Nobs // will turn sum of squared errors in denom of t/z into mean

		if (!null) M_DGP.beta = J(0,1,0) // in case model re-dirtied and we're not imposing null, Estimate() will know to recompute beta for first r0 value tried, then stop

		initialized = 1
	} // done with one-time stuff--not dependent on r0--if constructing CI or plotting confidence curve

	if (!ML) { // GMM, 2SLS, analytical LIML
		rAll = null? *pr0 : J(0, 1, 0); if (REst) rAll =  *pr \ rAll // constant terms of model + null constraints
		sAll = rows(rAll) ? LAll_invRAllLAll * rAll : J(0,1,0)
		M_DGP.Estimate(sAll)

		if (AR) {
			pM_Repl->InitEndog(&(*pY - *pXEnd * *pr0), NULL, &(*M_DGP.pZExclY - *M_DGP.pZExclXEnd * *pr0), &(*M_DGP.pXExY - *M_DGP.pXExXEnd * *pr0))
			pM_Repl->InitEstimate()
			pM_Repl->Estimate(sAR)
		}
	}

	if (WREnonAR) {
		_e = M_DGP.e + M_DGP.e2 * M_DGP.beta[|kEx+1\.|]
		Dist = J(cols(u), 1, .)
		pu = Nclust? &U : &u
		Ystar = *M_DGP.pY :+ _e :* *pu
		XExYstar   = cross(*pXEx  , *pwt, Ystar)
		ZExclYstar = cross(*pZExcl, *pwt, Ystar)
		XEndstar = XExXEndstar = ZExclXEndstar = smatrix(d-1)
		
		if (LIML | !robust)
			YstarYstar = weights? cross(*pwt, Ystar:*Ystar) : colsum(Ystar:*Ystar)

		if (d==2) {
				XEndstar.M         = *pXEnd :+ M_DGP.e2 :* *pu
				XExXEndstar.M      = cross(*pXEx  , *pwt, XEndstar.M)
				ZExclXEndstar.M    = cross(*pZExcl, *pwt, XEndstar.M)
		} else
			for (j=d-1; j; j--) {
				XEndstar[j].M      = (*pXEnd)[,j] :+ M_DGP.e2[,j] :* *pu
				XExXEndstar  [j].M = cross(*pXEx  , *pwt, XEndstar[j].M)
				ZExclXEndstar[j].M = cross(*pZExcl, *pwt, XEndstar[j].M)
			}

		if (Nclust) {
			XZi = eZi = smatrix(clust.N)
			for (i=clust.N; i; i--) {
				Subscripts = clust.info[i,]', (.\.)
				Zi = (*pXEx)[|Subscripts|] , (*pZExcl)[|Subscripts|] // inefficient?
				if (weights) Zi = Zi :* (*pwt)[|Subscripts|]
				XZi[i].M = cross((*pXEx)[|Subscripts|], Zi) \ cross((*pXEnd)[|Subscripts|], Zi) \ cross((*pY)[|Subscripts|], Zi)
				eZi[i].M =                                    cross(M_DGP.e2[|Subscripts|], Zi) \ cross(   _e[|Subscripts|], Zi)
			}
		}

		for (j=cols(u); j; j--) { // WRE bootstrap
			pXEndstar      = &( XEndstar.M  [,j])
			pXExXEndstar   = &(XExXEndstar.M  [,j])
			pZExclXEndstar = &(ZExclXEndstar.M[,j])
			for (i=2; i<d; i++) {
				pXEndstar      = &(*pXEndstar     , XEndstar     [i].M[,j])
				pXExXEndstar   = &(*pXExXEndstar  , XExXEndstar  [i].M[,j])
				pZExclXEndstar = &(*pZExclXEndstar, ZExclXEndstar[i].M[,j])
			}

			pM_Repl->InitEndog(&(Ystar[,j]), pXEndstar, &(ZExclYstar[,j]), &(XExYstar[,j]), (LIML | !robust? YstarYstar[j] : .), pZExclXEndstar, pXExXEndstar)
			pM_Repl->InitEstimate()
			pM_Repl->InitTestDenoms(S) // prepare for replication regressions, null not imposed
			pM_Repl->Estimate(s)
			numer = *pR0 * pM_Repl->beta - *pr0

			if (robust) { // Compute denominator for this WRE test stat
				if (Nclust != 1) eZVR0 = pM_Repl->e :* pM_Repl->ZVR0
				for (c=1; c<=length(clust); c++) {
					if (c==1 & Nclust) {
						AVR0 = pM_Repl->A * pM_Repl->VR0; _beta = -pM_Repl->beta \ 1; betaEnd = _beta[|kEx+1\.|]

						SeuZVR0 = (_beta'XZi[clust.N].M + betaEnd'eZi[clust.N].M * u[clust.N,j]) * AVR0 // R0 * V * Z_i'estar_i
						denom = cross(SeuZVR0, SeuZVR0)
						for (i=clust.N-1; i; i--) {
							SeuZVR0 = (_beta'XZi[i].M + betaEnd'eZi[i].M * u[i,j]) * AVR0 // R0 * V * Z_i'estar_i
							denom = cross(SeuZVR0, SeuZVR0) + denom
						}
						if (clust.multiplier!=1) denom = denom * clust[c].multiplier
					} else {
						if (rows(clust[c].order)) _collate(eZVR0, clust[c].order) // non-bootstrapping cluster
						pSeZVR0 = !Nclust | clust[c].N==Nobs? (weights? &(eZVR0 :* *pwt) : &eZVR0) : &_panelsum(eZVR0, *pwt, clust[c].info)
						t = cross(*pSeZVR0, *pSeZVR0); if (clust[c].multiplier!=1) t = t * clust[c].multiplier; denom = c==1? t : denom + t
					}
				}
			} else
				denom = (*pR0 * pM_Repl->VR0) * pM_Repl->eec

			Dist[j] = sqrt? numer/sqrt(denom) : cross(numer, invsym(denom) * numer)
		}
	} else {

		if (ML)
			eZVR0 = *pSc * (VR0 = *pV * *pR0')
		else {
			pM->InitTestDenoms(AR? SAR : S)
			if (scoreBS | robust)
				eZVR0 = pM->e :* pM->ZVR0
		}

		if (scoreBS)
			numer = cross(Nclust? _panelsum(eZVR0, *pwt, clust.info) : (weights? eZVR0:* *pwt : eZVR0), u)
		else {
			pewt = weights? &(pM->e:* *pwt) : &pM->e
			betadenom = K? (*pM->pV * pM->A ') : *pM->pV // in IV/GMM, this is actually not denominator (V) but V * X'Z(Z'Z)^-1
			betanumer = cross( Nclust? _panelsum(*pXEx  , *pewt, clust.info) : *pXEx   :* *pewt , u) \ 
			            cross( Nclust? _panelsum(*pZExcl, *pewt, clust.info) : *pZExcl :* *pewt , u)
			numer = (*pR0)[|.,.\.,kEx|] * (betadevEx = betadenom[|.,.\kEx,.|] * betanumer)

			if (K | AR)
				numer = numer + (*pR0)[|.,kEx+1\.,.|] * (betadevEnd = betadenom[|kEx+1,.\.,.|] * betanumer)
			else
				betadevEnd = J(0,cols(u),0)
		}

		if      (AR)    numer[,1] = pM->beta[|kEx+1\.|]    // coefficients on excluded instruments in AR OLS
		else if (!null) numer[,1] = *pR0 * (ML? beta : pM->beta) - *pr0 // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this.

		// Compute denominators and then test stats
		if (robust) {
			if (df == 1) {  // special, optimized for one null constraint
				if (Nclust > 1) {
					peZVR0wt  = weights? &(eZVR0 :* *pwt) : &eZVR0
					eUZVR0wt = U :* *peZVR0wt
					if (!scoreBS) {
						XExZVR0wt  = *pXEx   :* *peZVR0wt
						XEndZVR0wt = *_pXEnd :* *peZVR0wt
					}
				}

				for (c=1; c<=length(clust); c++) {
					if (!Nclust | clust[c].N==Nobs) { // het-only robust
						SeuZVR0 = (clust[c].N==Nobs? U : u) :* eZVR0 :* *pwt
						if (!scoreBS) SeuZVR0 = SeuZVR0 - (*pXEx :* (*pwt :* pM->ZVR0)) * betadevEx - (*_pXEnd :* (*pwt :* pM->ZVR0)) * betadevEnd // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
					} else if (c == 1) { // cluster we're bootstrapping on?
						SeuZVR0 = u :* _panelsum(eZVR0, *pwt, clust.info)
						if (!scoreBS) SeuZVR0 = SeuZVR0 - _panelsum(*pXEx :* pM->ZVR0, *pwt, clust.info) * betadevEx - _panelsum(*_pXEnd :* pM->ZVR0, *pwt, clust.info) * betadevEnd
					} else {
						if (rows(clust[c].order)) {
							_collate(eUZVR0wt, clust[c].order)
							if (!scoreBS) {
								_collate(XExZVR0wt , clust[c].order)
								_collate(XEndZVR0wt, clust[c].order)
							}
						}
						SeuZVR0 = _panelsum(eUZVR0wt, clust[c].info)
						if (!scoreBS) SeuZVR0 = SeuZVR0 - _panelsum(XExZVR0wt, clust[c].info) * betadevEx - _panelsum(XEndZVR0wt, clust[c].info) * betadevEnd
					}
					if (scoreBS) SeuZVR0 = SeuZVR0 :- clust[c].ClustShare*colsum(SeuZVR0) // recenter variance if not already done. Horowitz (2001), (3.29)
					t = colsum(SeuZVR0 :* SeuZVR0); if (clust[c].multiplier!=1) t = t * clust[c].multiplier; denom = c==1? t : denom + t
				}
				Dist = (sqrt? numer :/ sqrt(denom) : (numer:*numer) :/ denom)'
			} else { // more than one null constraint
				Dist = J(cols(u), 1, .); denoms = smatrix(cols(u))

				if (!scoreBS) ZVR0wt = pM->ZVR0 :* *pwt

				for (c=1; c<=length(clust); c++) {
					if (!Nclust) // het-only robust
						SeZVR0 = eZVR0 :* *pwt
					else if (c == 1) // cluster we're bootstrapping on?
						SeZVR0 = _panelsum(eZVR0, *pwt, clust.info)
					else if (rows(clust[c].order)) { // non-bootstrapping cluster
						_collate(eZVR0, clust[c].order)
						if (!scoreBS) {
							_collate(ZVR0wt , clust[c].order)
							_collate(*pXEx  , clust[c].order)
							_collate(*_pXEnd, clust[c].order)
						}
					}
					for (i=cols(u); i; i--) {
						             SeuZVR0 = c==1? SeZVR0 :* u[,i] : _panelsum(eZVR0, (weights? U[,i] :* *pwt : U[,i]), clust[c].info)
						if (!scoreBS)SeuZVR0 = SeuZVR0  - (clust[c].N==Nobs? (*pXEx*betadevEx[,i]+*_pXEnd*betadevEnd[,i]) :* ZVR0wt : _panelsum(ZVR0wt, *pXEx*betadevEx[,i]+*_pXEnd*betadevEnd[,i], clust[c].info)) // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
						if (scoreBS & reps) SeuZVR0 = SeuZVR0 :- clust[c].ClustShare*colsum(SeuZVR0) // Center variance
						t = cross(SeuZVR0, SeuZVR0); if (clust[c].multiplier!=1) t = t * clust[c].multiplier; denoms[i].M = cols(denoms[i].M)? denoms[i].M + t : t
					}
				}
				for (i=cols(u); i; i--) {
					numer_i = numer[,i]
					Dist[i] = cross(numer_i, invsym(denoms[i].M) * numer_i)
				}
			}
		} else { // non-robust
			pVR0 = ML? &VR0 : &(pM->VR0)
			if (df == 1) {  // optimize for one null constraint
				Dist = sqrt? numer / sqrt(*pR0 * *pVR0) : (numer:*numer) / (*pR0 * *pVR0)
				if (ML)
					Dist = Dist'
				else {
					             eu = u :* pM->e
					if (scoreBS) eu = eu :- (weights? cross(clust.ClustShare, eu) : colsum(eu) * clust.ClustShare)  // Center variance if needed
					  else       eu = eu  - *pXEx * betadevEx - *_pXEnd * betadevEnd // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
					t = weights? cross(*pwt, eu :* eu) : colsum(eu :* eu)
					Dist = (Dist  :/ (sqrt? sqrt(t) : t))'
				}
			} else {
				denom = invsym(*pR0 * *pVR0)
				Dist = J(cols(u), 1, .)

				for (i=cols(u); i; i--) {
					numer_i = numer[,i]
					Dist[i] = cross(numer_i, denom * numer_i) 
					if (!(ML | LIML)) {
						             eu = u[,i] :* pM->e
						if (scoreBS) eu = eu :- (weights? cross(*pwt, eu) : colsum(eu)) * clust.ClustShare // Center variance if needed
						  else       eu = eu  - *pXEx * betadevEx[,i] - *_pXEnd * betadevEnd[,i] // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
						
						Dist[i] = Dist[i] / cross(eu, *pwt, eu)
					}
				}
			}
		}
	}
	if (multiplier!=1) Dist = Dist * (sqrt? sqrt(multiplier) : multiplier)
	DistCDR = J(0,0,0)
	dirty = 0
}

// like panelsetup() but can group on multiple columns, like sort(), and faster. But doesn't take minobs, maxobs arguments.
// Does take optional third argument, a matrix in which to store standardized ID variable, starting from 1
real matrix _panelsetup(real matrix X, real rowvector cols, | real colvector ID) {
	real matrix info; real scalar i, N; real scalar p; real rowvector t, id
	N = rows(X)
	info = J(N, 2, N); if (args()>2) ID = J(N, 1, 1)
	info[1,1] = p = 1
	id = X[1, cols]
	for (i=2; i<=N; i++) {
		if ((t=X[i,cols]) != id) {
			info[  p,2] = i - 1
			info[++p,1] = i
			id = t
		}
		if (args()>2) ID[i] = p
	}
	return (info[|.,.\p,.|])
}

// stick call to panelsum() in separate function to prevent run-time error in old Stata versions
real matrix __panelsum(real matrix X, real matrix arg2, real matrix arg3)
	return (cols(arg3)? panelsum(X, arg2, arg3) : panelsum(X, arg2))

// implement Mata's panelsum() for pre-version 13. Differs in that a single missing value in X doesn't make all results missing.
real matrix _panelsum(real matrix X, real matrix arg2, | real matrix arg3) {
	if (stataversion() >= 1300)
		return (__panelsum(X, arg2, arg3))
	
	real matrix retval, Xi, Wi; pointer(real matrix) scalar pinfo; pointer(real colvector) scalar pwt; real scalar i
	pragma unset Xi; pragma unset Wi
	if (args()==2)
		pinfo = &arg2
	else {
		pinfo = &arg3
		if (arg2!=1)
			pwt = &arg2
	}
	
	retval = J(rows(*pinfo), cols(X), .)
	for (i=rows(*pinfo); i; i--) {
		panelsubview(Xi, X, i, *pinfo)
		if (pwt==NULL)
			retval[i,] = colsum(Xi)
		else if (rows(*pwt)==1)
			retval[i,] = Xi * *pwt
		else {
			panelsubview(Wi, *pwt, i, *pinfo)
			retval[i,] = cross(Xi, Wi)'
		}
	}
	return (retval)
}

// given a vector, return indices of the non-zero elements, like selectindex() function added in Stata 13
// if v = 0 (so can't tell if row or col vector), returns rowvector J(1, 0, 0) 
real vector boottest_selectindex(real vector v) {
	real scalar rows
	rows = rows(v)
	return(select(rows>1? 1::rows : 1..cols(v), v))
}

// Return matrix that counts from 0 to 2^N-1 in binary, one column for each number, one row for each binary digit
// except use provided lo and hi values for 0 and 1
real matrix boottestModel::count_binary(real scalar N, real scalar lo, real scalar hi) {
	real matrix t
	if (N<=1) return (lo , hi)
	t = count_binary(N-1, lo, hi)
	return (J(1, cols(t), lo), J(1, cols(t), hi) \ t, t)
}

// given a pre-configured boottest linear model with one-degree null imposed, compute distance from target p value of boostrapped one associated with given value of r0
// used with optimize() to construct confidence intervals
// performs no error checking
real scalar boottestModel::r0_to_p(real scalar r0) {
	pr0 = &r0
	dirty = 1
	return (get_padj())
}

real scalar boottestModel::search(real scalar alpha, real scalar p_lo, real scalar lo, real scalar p_hi, real scalar hi) {
	real scalar mid, _p
	mid = (alpha-p_lo)/(p_hi-p_lo)*(hi-lo) + lo
	if (mreldif(lo,mid)<1e-6 | mreldif(hi,mid)<1e-6)
		return (mid)
	if ( ((_p = r0_to_p(mid)) < alpha) == (p_lo < alpha) )
		return(search(alpha, _p, mid, p_hi, hi))
	return(search(alpha, p_lo, lo, _p, mid))
}
	
// derive wild bootstrap-based CI, for case of linear model with one-degree null imposed.
void boottestModel::plot(real scalar level) {
	real scalar t, alpha, _quietly, c, i, j; real colvector lo, hi; pointer (real colvector) _pr0

	_quietly = quietly
	boottest_set_quietly(this, 1)
	boottest() // run in order to get true number of replications

	alpha = 1 - level*.01
	if (alpha>0 & cols(u)-1 <= 1/alpha - 1e6) {
		boottest_set_quietly(this, _quietly)
		if (!quietly) errprintf("\nError: need at least %g replications to resolve a %g%% two-sided confidence interval.\n", ceil(1/alpha), level)
		return (.\.)
	}
	
	_pr0 = pr0
	_editmissing(gridpoints, 25)

	if (gridstart==. | gridstop==.) {
		if (reps)
			if (AR) {
				t = abs(cuepoint) * (small? invttail(df_r, alpha/2)/invttail(df_r, get_padj(1)/2) : invnormal(alpha/2)/invnormal(get_padj(1)/2))
				lo = gridstart<.? gridstart : cuepoint - t
				hi = gridstop <.? gridstop  : cuepoint + t
			} else {
				make_DistCDR()
				lo = gridstart<.? gridstart : numer[1] + *pr0 + DistCDR[floor((   alpha/2)*(rows(DistCDR)-1))+1] * abs(numer[1]/Dist[1]) // initial guess based on distribution from main test
				hi = gridstop <.? gridstop  : numer[1] + *pr0 + DistCDR[ceil (( 1-alpha/2)*(rows(DistCDR)-1))+1] * abs(numer[1]/Dist[1])
			}
		else {
			t = abs(numer[1]/Dist[1]) * (small? -invttail(df_r, alpha/2) : invnormal(alpha/2))
			lo = gridstart<.? gridstart : numer[1] + *pr0 + t
			hi = gridstop <.? gridstop  : numer[1] + *pr0 - t
		}

		if (gridstart==. & ptype!=3) // unless upper-tailed p value, try at most 10 times to bracket confidence set by doubling on low side
			for (i=10; i & -r0_to_p(lo)<-alpha; i--)
				lo = 2 * lo - hi
		if (gridstop==. & ptype!=2) // ditto for high side
			for (i=10; i & -r0_to_p(hi)<-alpha; i--)
				hi = 2 * hi - lo
	} else {
		lo = gridstart
		hi = gridstop
	}

	plotX = rangen(lo, hi, gridpoints)
	if (cuepoint == .) cuepoint = numer[1] + *pr0 // non-AR case
	if (cuepoint < lo) { // insert original point estimate into grid
		if (gridstart == .) {
			plotX = cuepoint \ plotX
			c = 1
		}
	} else if (cuepoint > hi) {
		if (gridstop == .) {
			plotX = plotX \ cuepoint
			c = gridpoints+1
		}
	} else {
		c = floor((cuepoint - lo)/(hi - lo)*(gridpoints - 1)) + 2
		plotX = plotX[|.\c-1|] \ cuepoint \ plotX[|c\.|]
	}

	plotY = J(rows(plotX), 1, .)
	printf("{txt}")
	for (i = rows(plotX); i; i--) {
		plotY[i] = r0_to_p(plotX[i])
		printf(".")
		if (mod(i-rows(plotX)-1,50)) displayflush()
			else printf("\n")
	}

	if (level<100) {
		CI = (plotY :> alpha) :/ (plotY :< .); CI = CI[|2\.|] - CI[|.\rows(plotX)-1|]
		lo = boottest_selectindex(CI:== 1)
		hi = boottest_selectindex(CI:==-1)
		if (rows(lo)==0 & rows(hi)==0)
			CI = . , .
		else {
					 if (rows(lo)==0) lo = .
			else if (rows(hi)==0) hi = .
			else {
				if ( lo[1       ] >  hi[1       ]) lo = .  \ lo // non-rejection ranges that are not bounded within grid
				if (-lo[rows(lo)] < -hi[rows(hi)]) hi = hi \ .
			}
			CI = lo, hi
			for (i=rows(lo); i; i--)
				for (j=2; j; j--)
					if (CI[i,j]<.)
						CI[i,j] = search(alpha, plotY[CI[i,j]], plotX[CI[i,j]], plotY[CI[i,j]+1], plotX[CI[i,j]+1])
		}
	}

	if (c < .) { // now that it's done helping find CI points, remove cuepoint from grid for evenness, for Bayesian sampling purposes
		peak = plotX[c], plotY[c]
		if (c==1) { // now that it's done helping find CI points, remove cuepoint from grid for evenness, for Bayesian sampling purposes
			plotX = plotX[|2\.|]; plotY = plotY[|2\.|]
		} else if (c==gridpoints+1) {
			plotX = plotX[|.\gridpoints|]; plotY = plotY[|.\gridpoints|]
		} else {
			plotX = plotX[|.\c-1|] \ plotX[|c+1\.|]; plotY = plotY[|.\c-1|] \ plotY[|c+1\.|]
		}
	}

	boottest_set_quietly(this, _quietly)
	pr0 = _pr0 // inconsistent in bypassing function interface to class members
}

// return matrix whose rows are all the subsets of a row of numbers. Nil is at bottom.
real matrix boottestModel::combs(real rowvector X) {
	real matrix t
	if (cols(X)==1) return (X \ .)
	t = combs(X[|2\.|])
	return ((J(rows(t),1,X[1]), t) \ (J(rows(t),1,.), t))
}

// Stata interface
void boottest_stata(string scalar statname, string scalar dfname, string scalar dfrname, string scalar pname, string scalar padjname, string scalar ciname, 
	string scalar plotname, string scalar peakname, real scalar level, real scalar ML, real scalar LIML, real scalar Fuller, 
	real scalar K, real scalar AR, real scalar null, real scalar scoreBS, string scalar wildtype, string scalar ptype, string scalar madjtype, real scalar NumH0s,
	string scalar XExnames, string scalar XEndnames, real scalar cons, string scalar Ynames, string scalar bname, string scalar Vname, string scalar Wname, 
	string scalar ZExclnames, string scalar samplename, string scalar scnames, real scalar robust, string scalar IDnames, 
	string scalar wtname, string scalar wttype, string scalar Cname, string scalar C0name, real scalar reps, real scalar small, string scalar distname, ///
	real scalar gridmin, real scalar gridmax, real scalar gridpoints) {

	real matrix C, R, C0, R0, ZExcl, ID, sc, XEnd, XEx
	real colvector r, wt, r0, Y
	class boottestModel scalar M
	pragma unset ID; pragma unset wt; pragma unset XEnd; pragma unset XEx; pragma unset Y; pragma unset ZExcl; pragma unset sc

	C0 = st_matrix(C0name)
	R0 = C0[|.,.\.,cols(C0)-1|]
	r0 = C0[,cols(C0)]
	C = st_matrix(Cname)
	if (rows(C)) { // restricted OLS?
		R = C[|.,.\.,cols(C)-1|]
		r = C[,cols(C)]
	}

	_boottest_st_view(sc, ., scnames, samplename)
	_boottest_st_view(Y , ., Ynames , samplename)
	_boottest_st_view(ZExcl, ., ZExclnames , samplename)
	if (IDnames != "") ID  = st_data(., IDnames, samplename)
	if (wtname  != "") wt  = st_data(., wtname, samplename) // panelsum() doesn't like views as weights

	boottest_set_cons(M, cons)
	boottest_set_sc(M, sc)
	boottest_set_ML(M, ML)
	boottest_set_Y (M, Y)
	boottest_set_ZExcl(M, ZExcl)
	boottest_set_wt (M, wt)
	boottest_set_ID(M, ID)
	boottest_set_R (M, R , r )
	boottest_set_R0(M, R0, r0)
	boottest_set_null(M, null)
	boottest_set_small(M, small)
	boottest_set_robust(M, robust)
	boottest_set_scoreBS(M, scoreBS)
	boottest_set_wildtype(M, wildtype)
	boottest_set_ptype(M, ptype)
	boottest_set_wttype(M, wttype)
	boottest_set_reps (M, reps)
	boottest_set_LIML(M, LIML)
	boottest_set_Fuller(M, Fuller)
	boottest_set_k(M, K)
	boottest_set_AR(M, AR)
	boottest_set_grid(M, gridmin, gridmax, gridpoints)
	boottest_set_madjust(M, madjtype, NumH0s)

	_boottest_st_view(XEnd, ., XEndnames, samplename)
	boottest_set_XEnd(M, XEnd)
	_boottest_st_view(XEx, ., XExnames, samplename)
	boottest_set_XEx(M, XEx)
	if (bname != "") boottest_set_beta(M, st_matrix(bname)')
	if (Vname != "") boottest_set_V   (M, st_matrix(Vname) )
	if (Wname != "") boottest_set_W   (M, st_matrix(Wname) )
	boottest_set_willplot(M, plotname != "")

	st_numscalar(statname, M.get_stat())
	st_numscalar(pname   , M.get_p   ())
	st_numscalar(padjname, M.get_padj())
	st_numscalar(dfname  , M.get_df  ())
	st_numscalar(dfrname , M.get_df_r())
	if (distname != "") st_matrix(distname, M.get_dist())
	if (plotname != "") {
		M.plot(level)
		st_matrix(plotname, (M.plotX,M.plotY))
		if (cols(M.peak)) st_matrix(peakname, M.peak)
		if (level<100 & ciname != "") st_matrix(ciname, M.CI) // also makes plotX & plotY
	}

	M.M_DGP.setParent(NULL) // actually sets the pointer to &NULL, but that suffices to break loop in the data structure topology and avoid Mata garbage-cleaning leak
}

mata mlib create lboottest, dir("`c(sysdir_plus)'l") replace
mata mlib add lboottest *(), dir("`c(sysdir_plus)'l")
mata mlib index
end
