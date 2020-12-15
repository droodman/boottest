*! boottest 3.0.0 14 December 2020
*! Copyright (C) 2015-20 David Roodman

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

struct smatrix {
	real matrix M
}

struct ssmatrix {
	struct smatrix matrix M
}

struct structboottestClust {
	real scalar N, multiplier, even
	real colvector order
	real matrix info
}

struct structFE {
	real colvector is, wt
}

class AnalyticalModel {  // class for analyitcal OLS, 2SLS, LIML, GMM estimation--everything but iterative ML
	real scalar LIML, YY, eec, Fuller, AR, K, k, isDGP
	real matrix ZX, ZXEnd, ZZ, H_2SLS, invH, A, VR0, ZVR0, e2, Splus, pi, XXEnd, dbetads, ZExclXEnd, XExXEnd
	real colvector sAll, e, beta, beta0, ZY
	real rowvector YXEnd
	pointer(real colvector) scalar pY
	pointer(real matrix) scalar pXEnd, pXX, pXY, pZExclY, pXExY, pV, pW, pinvZZ, pXExXEx, pZXEx, pH, pR0, pS
	pointer (class boottestModel scalar) scalar parent
	struct smatrix matrix CT_ZVR0
	struct smatrix colvector WZVR0

	void new(), InitExog(), InitEndog(), InitTestDenoms(), SetS(), InitEstimate(), Estimate(), SetAR()
	pointer(real matrix) scalar partialFE()
}

class boottestModel {
	real scalar scoreBS, reps, small, weighttype, null, dirty, initialized, Neq, ML, GMM, Nobs, _Nobs, k, kEx, el, sumwt, NClustVar, robust, weights, REst, multiplier, smallsample, quietly, FEboot, NErrClustCombs, ///
		sqrt, hascons, LIML, Fuller, K, IV, WRE, WREnonAR, ptype, twotailed, df, df_r, AR, D, confpeak, willplot, notplotted, NumH0s, p, NBootClustVar, NErrClust, dH0, ///
		NFE, doKK, granular, purerobust, subcluster, NBootClust, repsFeas, u_sd, level, ptol, MaxMatSize, NWeightGrps, enumerate, bootstrapt, q, q0, interpolative, interpolating, interpolate_e, interpolate_denom
	real matrix VR0, numer, u, eu, S, SAR, SAll, LAll_invRAllLAll, plot, CI, CT_WE, infoBootData, infoBootAll, infoErrAll, J_ClustN_NBootClust, statDenom, eZVR0, SewtXV, numer0
	real colvector DistCDR, s, sAR, plotX, plotY, sAll, beta, wtFE, ClustShare, IDBootData, IDBootAll, WeightGrpStart, WeightGrpStop, gridmin, gridmax, gridpoints, numersum, r00, e0, anchor, poles
	real rowvector peak
	string scalar wttype, madjtype, seed
	pointer (real matrix) scalar pZExcl, pR, pR0, pID, pFEID, pXEnd, pXEx, pG, pX, pinfoAllData, pinfoErrData
	pointer (real colvector) scalar pr, pr0, pY, pSc, pwt, pW, pV, pe, pDist
	class AnalyticalModel scalar M_DGP
	pointer (class AnalyticalModel scalar) scalar pM_Repl, pM
	struct structboottestClust colvector Clust
	pointer (struct structboottestClust scalar) scalar pBootClust
	struct smatrix matrix KK, denom, Kcd, denom0, Jcd0
	struct smatrix colvector Kd, XZi, eZi, euZVR0, betadev, dedr, dnumerdr
  struct ssmatrix colvector ddenomdr, dJcddr
  struct ssmatrix matrix ddenomdr2
	pointer(struct smatrix matrix) scalar pJcd /*, pKJcd*/
	struct structFE rowvector FEs
  /*pointer (real matrix function) pfncross, pfndouble*/

	void new(), setsqrt(), boottest(), plot(), setXEx(), setptype(), setdirty(), setXEnd(), setY(), setZExcl(), setwt(), setsc(), setML(), setLIML(), setAR(),
		setFuller(), setk(), setquietly(), setbeta(), setV(), setW(), setsmall(), sethascons(), setscoreBS(), setreps(), setnull(), setWald(), setRao(), setwttype(), setID(), setFEID(), setlevel(), setptol(), 
		setrobust(), setR(), setR0(), setwillplot(), setgrid(), setmadjust(), setweighttype(), makeWildWeights(), makeStuffLinearInr0(), _makeStuffLinearInr0(), makeNonWREStats(), makeBootstrapcDenom(), setMaxMatSize(), setstattype(), Initialize()
  private void makeKK(), makeJ(), _clustAccum()
	real scalar r0_to_p(), search(), getp(), getpadj(), getstat(), getdf(), getdf_r(), getreps(), getrepsFeas(), makeNonWRENumers(), makeWREStats()
	real matrix getplot(), getCI(), getV(), getv()
  private real matrix count_binary(), crosstab() /*, QuarticRoots()*/
  private static real matrix combs()
  private static real colvector stableorder()
  private static void _fold()
  static void _st_view()
	real rowvector getpeak()
	real colvector getdist(), getb()
}

void AnalyticalModel::new() {
	AR = 0
  isDGP = 1  // by default, created object is for DGP rather than replication regression
}

void AnalyticalModel::SetAR(real scalar _AR) {
	if (AR = _AR) LIML = Fuller = K = 0
}

// stuff that can be done before r0 set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
void AnalyticalModel::InitExog() {
	real matrix ZExclXEx

	parent->pXEx = partialFE(parent->pXEx)
	pXExXEx = &cross(*parent->pXEx, *parent->pwt, *parent->pXEx)
	if (cols(*parent->pZExcl)) {  // GMM, 2SLS, LIML
		parent->pZExcl = partialFE(parent->pZExcl)
		ZExclXEx = cross(*parent->pZExcl, *parent->pwt, *parent->pXEx)
		pZXEx = &(*pXExXEx \ ZExclXEx)
		if (parent->IV)  // non-GMM
			pinvZZ = &invsym((ZZ = *pZXEx, (ZExclXEx' \ cross(*parent->pZExcl, *parent->pwt, *parent->pZExcl))))
	} else
		pXX = pXExXEx
}

// stuff that can be done before S & r0 set, but depend on endogenous variables, which are bootstrapped in WRE
void AnalyticalModel::InitEndog(pointer (real colvector) scalar _pY, pointer (real matrix) scalar _pXEnd, | ///
		pointer (real colvector) scalar _pZExclY, pointer (real rowvector) scalar _pXExY) {

	pY = partialFE(_pY); pXEnd = partialFE(_pXEnd)

	pXExY = _pXExY==NULL? &cross(*parent->pXEx, *parent->pwt, *pY) : _pXExY
	if (K | AR)
		pZExclY = _pZExclY==NULL? &cross(*parent->pZExcl, *parent->pwt, *pY) : _pZExclY
	if (K) {
		XExXEnd   = cross(*parent->pXEx  , *parent->pwt, *pXEnd)
		ZExclXEnd = cross(*parent->pZExcl, *parent->pwt, *pXEnd)
		ZXEnd = XExXEnd \ ZExclXEnd
		ZX = *pZXEx, ZXEnd
		XXEnd = XExXEnd \ cross(*pXEnd, *parent->pwt, *pXEnd)
		pXX = &(*pXExXEx,  XExXEnd \ XXEnd')
		YXEnd = cross(*pY, *parent->pwt, *pXEnd)
		ZY = *pXExY \ *pZExclY
		pXY = &(*pXExY \ YXEnd')
		if (LIML | (parent->robust | parent->scoreBS)==0)
			YY = cross(*pY, *parent->pwt, *pY)
		if (parent->IV)  // if GMM weight matrix not provided, prepare 2SLS one
			A = (I(parent->kEx) \ J(parent->el-parent->kEx, parent->kEx, 0)), *pinvZZ * ZXEnd // 2SLS is (A' ZX)^-1 * (A'ZY). Also apparently used in k-class and LIML robust VCV by Stata convention
		else
			A = *parent->pW * ZX
		H_2SLS = A ' ZX  // Hessian
	} else {  // OLS / AR
		pXY = pXExY
		if (AR) {
			pXY = &(*pXY \ *pZExclY)
			pXX = &ZZ
		}
	}
	k = cols(*pXX)
}

void AnalyticalModel::SetS(real matrix S) {  // set model constraint matrix (not null constraints)
	pS = &S
	if (LIML) Splus = blockdiag(S, 1)  // add an entry to S for the dep var
}

// stuff that can be done before r0 set but depends on S and endogenous variables
void AnalyticalModel::InitEstimate() {
	real rowvector val
	real matrix _ZY, TT, TPZT, vec
	pointer (real matrix) scalar pbetadenom
	pragma unset vec; pragma unset val

 	if (LIML)
		if (parent->el == k)  // exactly identified LIML = 2SLS
			K = 1
		else {
			_ZY = ZXEnd, ZY
			TT = *pXX, (*pXExY \ YXEnd') \ *pXExY', YXEnd, YY  // where T = *pXEx, *pXEnd, *pY
			(TPZT = TT)[|parent->kEx+1,parent->kEx+1\.,.|] = _ZY ' (*pinvZZ) * _ZY

			A = *pinvZZ * ZX
			H_2SLS = A ' ZX  // Hessian

			if (rows(*pS)) {  // includes H0 if pSAll really points to SAll rather than S
				TT   = Splus '   TT * Splus
				TPZT = Splus ' TPZT * Splus
			}
			eigensystemselecti( I(rows(TT)) - boottest_lusolve(TT, TPZT), 1\1, vec, val)  // eigensystemselecti(lusolve(TT, TPZT), rows(TT)\rows(TT), ... gives 1 - the eigenvalue, but can cause eigensystem() to return all missing
			K = 1/Re(val) - Fuller / (parent->_Nobs - parent->el)   // sometimes a tiny imaginary component sneaks into val
		}

  pH = K? (K==1? &H_2SLS : &((1-K)* *pXX + K*H_2SLS)) : pXX
	if (rows(*pS)) {
		pbetadenom = &(*pS * boottest_lusolve(*pS ' (*pH) * *pS, *pS'))
		invH = J(0,0,0)
	} else
		pbetadenom = &(invH = invsym(*pH))

  if (K)
		if (K==1) {  // 2SLS, GMM
			beta0 = *pbetadenom * A ' ZY
			dbetads = I(rows(beta0)) - *pbetadenom * A ' ZX
		} else {  // k-class, LIML
			beta0 = *pbetadenom * (K * A ' ZY + (1-K) * *pXY)
			dbetads = I(rows(beta0)) - *pbetadenom * (K * A ' ZX + (1-K) * *pXX)
		}
	else {  // OLS / AR
		beta0 = *pbetadenom * *pXY
		dbetads = I(rows(beta0)) - *pbetadenom * *pXX
	}
}

// stuff that doesn't depend on r0, for test stat denominators in replication regressions
void AnalyticalModel::InitTestDenoms(real matrix S) {
	real matrix AVR0; real scalar d, c; struct smatrix rowvector _CT_ZVR0; pointer (real matrix) scalar pWZVR0

	pV = rows(S)? &(S * boottest_lusolve(S ' (*pH) * S, S')) : ///
	              &(rows(invH)? invH : invsym(*pH))
	VR0 = *pV * *parent->pR0'

	if (parent->scoreBS | (parent->robust & (parent->WREnonAR & parent->NClustVar==1 & parent->NFE==0)==0)) {
		if (K) {
			AVR0 = A * VR0
			ZVR0 = *parent->pZExcl * AVR0[|parent->kEx+1,.\.,.|]; if (parent->kEx) ZVR0 = ZVR0 + *parent->pXEx * AVR0[|.,.\parent->kEx,.|]
		} else if (AR) {
			ZVR0 = *parent->pZExcl * VR0[|parent->kEx+1,.\.,.|]
			if (cols(*parent->pXEx))
				ZVR0 = ZVR0 + *parent->pXEx * VR0[|.,.\parent->kEx,.|]
		} else
			ZVR0 = *parent->pXEx * VR0

			if (parent->bootstrapt == 0) return

			if (parent->granular) {
				pWZVR0 = &(parent->weights? ZVR0 :* *parent->pwt : ZVR0)
				WZVR0 = smatrix(parent->df)
				for (d=parent->df;d;d--)
					WZVR0[d].M = (*pWZVR0)[,d]
			}
			
      if (parent->NFE & parent->robust & (parent->WREnonAR | parent->FEboot | parent->scoreBS)==0 & parent->granular < parent->NErrClustCombs) {
				if (pWZVR0==NULL) pWZVR0 = &(parent->weights? ZVR0 :* *parent->pwt : ZVR0)
				CT_ZVR0 = smatrix(parent->NErrClustCombs, parent->df)
				for (d=parent->df;d;d--)
					CT_ZVR0[1,d].M = parent->crosstab((*pWZVR0)[,d])
				if (parent->NClustVar > 1) {
					_CT_ZVR0 = CT_ZVR0[1,]
					for (c=2; c<=parent->NErrClustCombs; c++)
						for (d=parent->df;d;d--) {
							if (rows(parent->Clust[c].order))
								_CT_ZVR0[d].M = _CT_ZVR0[d].M[parent->Clust[c].order,]
							CT_ZVR0[c,d].M = _panelsum(_CT_ZVR0[d].M, parent->Clust[c].info)
						}
				}
			}
	}
}

// stuff that depends on r0 and endogenous variables: compute beta and residuals
void AnalyticalModel::Estimate(real colvector s) {
	real colvector negZeinvee; real matrix Ze; real scalar ee

	beta = rows(s)? beta0 + dbetads * s : beta0
	if (isDGP | parent->bootstrapt | parent->WREnonAR==0) {  // don't need residuals in replication regressions in bootstrap-c on WRE/non-AR
		if (AR) {
			e = *pY - *parent->pZExcl * beta[|cols(*parent->pXEx)+1\.|]
			if (cols(*parent->pXEx))
				e =  e - *parent->pXEx * beta[|.\cols(*parent->pXEx)|]
		} else if (K)
			if (parent->kEx == 0)
				e = *pY - *pXEnd * beta[|parent->kEx+1\.|]
			else
				e = *pY - *pXEnd * beta[|parent->kEx+1\.|] - *parent->pXEx * beta[|.\parent->kEx|]
		else
				e = *pY                                    - *parent->pXEx * beta[|.\parent->kEx|]

    if ((parent->robust | parent->scoreBS)==0 | (isDGP & LIML))  // useful in non-robust, residual-based bootstrap, and in computing e2 in LIML (just below)
			ee = YY - 2 * *pXY ' beta + beta ' (*pXX) * beta
		if ((parent->robust | parent->scoreBS)==0 & parent->bootstrapt==0)
			eec = parent->hascons? ee : ee - (parent->weights? cross(e, *parent->pwt) : sum(e))^2 / parent->_Nobs  // sum of squares after centering, N * Var
	}

	if (isDGP & LIML) {
		Ze = ZY - ZX * beta
		negZeinvee = Ze / -ee
		pi = boottest_lusolve(ZZ + negZeinvee * Ze', negZeinvee * (YXEnd - beta ' XXEnd) + ZXEnd)  // coefficients in reduced-form equations; Davidson & MacKinnon (2010), eq 15
		e2 = *pXEnd - *parent->pZExcl * pi[|parent->kEx+1,.\.,.|]; if (parent->kEx) e2 = e2 - *parent->pXEx * pi[|.,.\parent->kEx,.|]
	}
}

// partial fixed effects out of a data matrix
pointer(real matrix) scalar AnalyticalModel::partialFE(pointer(real matrix) scalar pIn) {
	real matrix Out, t; real scalar i
	if (parent->NFE & pIn!=NULL) {
		Out = *pIn
		for (i=parent->NFE;i;i--) {
			t = Out[parent->FEs[i].is,]
			Out[parent->FEs[i].is,] = t :- cross(parent->FEs[i].wt, t)
		}
		return(&Out)
	}
	return (pIn)
}

void boottestModel::new() {
	AR = LIML = Fuller = WRE = small = scoreBS = weighttype = Neq = ML = initialized = quietly = sqrt = hascons = IV = ptype = robust = NFE = FEboot = granular = NErrClustCombs = subcluster = reps = repsFeas = interpolating = 0
	twotailed = null = dirty = willplot = u_sd = bootstrapt = notplotted = 1
	level = 95
  ptol = 1e-6
	confpeak = MaxMatSize = .
	pXEnd = pXEx = pZExcl = pY = pSc = pID = pFEID = pR = pR0 = pwt = &J(0,0,0)
	pr = pr0 = &J(0,1,0)
	IDBootData = .
}

void boottestModel::setdirty(real scalar _dirty, | real scalar noinitialize) {
	dirty = _dirty
	if (_dirty & noinitialize!=1)
		initialized = 0
}

void boottestModel::setsqrt(real scalar _sqrt) {
	if (_sqrt < sqrt) {
		if (dirty==0) {
    	pDist = &(*pDist :* *pDist)
      multiplier = multiplier * multiplier
    }
	} else
		setdirty(1)
	sqrt = _sqrt
}

void boottestModel::setptype(string scalar ptype) {
	real scalar p
	p = cross( (strtrim(strlower(ptype)) :== ("symmetric"\"equaltail"\"lower"\"upper")), 1::4 ) - 1
	if (p<0) 
		_error(198, `"p-value type must be "symmetric", "equaltail", "lower", or "upper"."')
	this.ptype = p
	this.twotailed = p<=1
}
void boottestModel::setstattype(string scalar stattype) {
	real scalar p
	p = cross( (strtrim(strlower(stattype)) :== ("c"\"t")), 1::2 ) - 1
	if (p<0) 
		_error(198, `"statistic type must be "t" or "c"."')
	this.bootstrapt = p
	setdirty(1)
}	
void boottestModel::setXEnd    (real matrix X ) {
	this.pXEnd  = &X; setdirty(1)
}
void boottestModel::setXEx     (real matrix X ) {
	this.pXEx  = &X; setdirty(1)
}
void boottestModel::setY       (real matrix Y ) {
	this.pY  = &Y; setdirty(1)
}
void boottestModel::setZExcl   (real matrix Z ) {
	this.pZExcl  = &Z; setdirty(1)
}
void boottestModel::setwt      (real matrix wt) {
	this.pwt  = &wt; setdirty(1)
}
void boottestModel::setsc      (real matrix Sc) {
	this.pSc  = &Sc
	setdirty(1)
}
void boottestModel::setML      (real scalar ML) {
	this.ML  = ML; setdirty(1)
	if (ML) setscoreBS(1)
}
void boottestModel::setLIML    (real scalar LIML) {
	this.LIML = LIML; setdirty(1)
}
void boottestModel::setAR    (real scalar AR) {
	this.AR = AR; setdirty(1)
}
void boottestModel::setFuller    (real scalar Fuller) {
	this.Fuller = Fuller; setdirty(1)
}
void boottestModel::setk    (real scalar K) {
	this.K = K; setdirty(1)
}
void boottestModel::setquietly (real scalar quietly )
	this.quietly = quietly
void boottestModel::setbeta    (real colvector beta) {
	this.beta = beta; setdirty(1)
}
void boottestModel::setV    (real matrix V ) {
	this.pV = &V; setdirty(1)
}
void boottestModel::setW    (real matrix W ) {
	this.pW = &W; setdirty(1)
}
void boottestModel::setsmall   (real scalar small   ) {
	this.small = small; setdirty(1)
}
void boottestModel::sethascons(real scalar hascons) {
	this.hascons = hascons; setdirty(1)
}
void boottestModel::setscoreBS (real scalar scoreBS ) {
	this.scoreBS = scoreBS; setdirty(1)
}
void boottestModel::setreps    (real scalar reps    ) {
	this.reps = reps
	if (reps==0)
		setscoreBS(1)
	setdirty(1)
}
void boottestModel::setnull    (real scalar null    ) {
	this.null = null; setdirty(1)
}
void boottestModel::setWald() { // set-up for classical Wald test
	this.scoreBS = 1; this.reps = 0; this.null = 0; setdirty(1)
}
void boottestModel::setRao() { // set-up for classical Rao test
	this.scoreBS = 1; this.reps = 0; this.null = 1; setdirty(1)
}
void boottestModel::setwttype  (string scalar wttype) {
	this.wttype = wttype; setdirty(1)
}
void boottestModel::setID      (real matrix ID, | real scalar NBootClustVar, real scalar NErrClust) {
	this.pID = &ID; this.NBootClustVar = editmissing(NBootClustVar,1); this.NErrClust=editmissing(NErrClust,1); setdirty(1)
	if (cols(ID)) this.robust = 1
}
void boottestModel::setFEID(real matrix ID, real scalar NFE) {
	this.pFEID = &ID; this.NFE = NFE; setdirty(1)
}
void boottestModel::setlevel  (real scalar level  )
	this.level = level
void boottestModel::setptol  (real scalar ptol  )
	this.ptol = ptol
void boottestModel::setrobust  (real scalar robust  ) {
	this.robust = robust
	if (robust==0) setID(J(0,0,0), 1, 1)
	setdirty(1)
}
void boottestModel::setR (real matrix R, real colvector r ) {
	this.pR = &R; 	this.pr  = &r; q = rows(R); setdirty(1)
}
void boottestModel::setR0(real matrix R0, real colvector r0) {
	this.pR0 = &R0; this.pr0 = &r0; q0 = rows(R0); setdirty(1); dH0 = rows(r0)  // dH0 can differ from df in AR test
}
void boottestModel::setwillplot(real scalar willplot) {
	this.willplot = willplot
}
void boottestModel::setgrid(real rowvector gridmin, real rowvector gridmax, real rowvector gridpoints) {
	this.gridmin = gridmin; this.gridmax = gridmax; this.gridpoints = gridpoints
}
void boottestModel::setmadjust(string scalar madjtype, real scalar NumH0s) {
	this.madjtype = strlower(madjtype)
	this.NumH0s = NumH0s
	if (this.madjtype != "bonferroni" & this.madjtype != "sidak" & this.madjtype != "")
		_error(198, `"Multiple-hypothesis adjustment type must be "Bonferroni" or "Sidak"."')
}
void boottestModel::setweighttype(string scalar weighttype) {
	weighttype = strlower(weighttype)
	if (.==(this.weighttype = weighttype=="rademacher" ? 0 : (weighttype=="mammen" ? 1 : (weighttype=="webb" ? 2 : (weighttype=="normal" ? 3 : (weighttype=="gamma" ? 4 : .))))))
		_error(198, `"Wild type must be "Rademacher", "Mammen", "Webb", "Normal", or "Gamma"."')
	setdirty(1)
}
void boottestModel::setMaxMatSize(real scalar MaxMatSize) {
	this.MaxMatSize = MaxMatSize; setdirty(1)
}

real colvector boottestModel::getdist(| string scalar diststat) {
	pointer (real rowvector) scalar pnumer
	if (dirty) boottest()
	if (diststat == "numer") {
		pnumer = u_sd==1? &numer : &(numer / u_sd)
		_sort( DistCDR = (*pnumer)[|2\.|]' :+ *pr0 , 1)
	} else if (rows(DistCDR)==0)
		if (rows(*pDist)>1)
			_sort( DistCDR=(*pDist)[|2\.|] , 1)
		else
			DistCDR = J(0,1,0)
	return(DistCDR)
}

// Robust to missing bootstrapped values interpreted as +infinity.
real scalar boottestModel::getp(|real scalar classical) {
	real scalar t
	if (dirty) boottest()
	t = (*pDist)[1]
	if (t == .) return (.)
	if (reps & classical==.)
		if (sqrt & ptype != 3) {
			if (ptype==0)
					p = colsum(-abs(t) :> -abs(*pDist)) / repsFeas  // symmetric p value; do so as not to count missing entries in Dist
			else if (ptype==1)  // equal-tail p value
				p = 2 * min((colsum(t :> *pDist) , colsum(-t:>- *pDist))) / repsFeas
			else
				p = colsum( t :>   *pDist) / repsFeas  // lower-tailed p value
		} else
				p = colsum(-t :> - *pDist) / repsFeas  // upper-tailed p value or p value based on squared stats
	else {
		t = t * multiplier
    p = small? Ftail(df, df_r, sqrt? t*t : t) : chi2tail(df, sqrt? t*t : t)
		if (sqrt & twotailed==0) {
			p = p / 2
			if ((ptype==3) == (t<0))
				p = 1 - p
		}
	}
	return(p)
}

// numerator for full-sample test stat
real colvector boottestModel::getb() {
	if (dirty) boottest()
	return(u_sd == 1? numer[,1] : numer[,1] / u_sd)
}

// denominator for full-sample test stat
real matrix boottestModel::getV() {
	if (dirty) boottest()
	return (statDenom / ((u_sd == 1? smallsample : u_sd * u_sd * smallsample)  * (sqrt? multiplier*multiplier : multiplier) * df))
}

// wild weights
real matrix boottestModel::getv()
	return(u_sd==1? u[|.,2\.,.|] : u[|.,2\.,.|] / u_sd)

// Return number of bootstrap replications with feasible results
// Returns 0 if getp() not yet accessed, or doing non-bootstrapping tests
real scalar boottestModel::getrepsFeas()
	return (repsFeas)

// return number of replications, possibly reduced to 2^G
real scalar boottestModel::getreps()
	return (reps)

real scalar boottestModel::getpadj(|real scalar classical) {
	(void) getp(classical)
	if (madjtype=="bonferroni") return(min((1, NumH0s*p)))
	if (madjtype=="sidak"     ) return(1 - (1 - p)^NumH0s)
	return(p)
}

real scalar boottestModel::getstat() {
	if (dirty) boottest()
	return(multiplier * (*pDist)[1])
}
real scalar boottestModel::getdf() {
	if (dirty) boottest()
	return(df)
}
real scalar boottestModel::getdf_r() {
	if (dirty) boottest()
	return(df_r)
}
real matrix boottestModel::getplot() {
	if (notplotted) plot()
	return((plotX,plotY))
}
real rowvector boottestModel::getpeak() {
	if (notplotted) plot()
	return(peak)
}
real matrix boottestModel::getCI() {
	if (notplotted) plot()
	return(CI)
}

void boottestModel::_st_view(real matrix V, real scalar i, string rowvector j, string scalar selectvar) {
	if (favorspeed() | 1)
		V = length(tokens(j))? st_data(i, j, selectvar) : st_data(i, J(1,0,0), selectvar)
	else
		st_view(V, i, j, selectvar)
}

// helper for summing over clusters while factoring in cluster-specific parity and small-sample adjustments
// replace X with Y if c=1; otherwise add it
void boottestModel::_clustAccum(real matrix X, real scalar c, real matrix Y)
  X = c == 1?
          (Clust.   even?
            (Clust.   multiplier != 1?   Clust.   multiplier  * Y :  Y) :
            (Clust.   multiplier != 1? (-Clust.   multiplier) * Y : -Y)) :
      X + (Clust[c].even?                
            (Clust[c].multiplier != 1?   Clust[c].multiplier  * Y :  Y) :
            (Clust[c].multiplier != 1? (-Clust[c].multiplier) * Y : -Y))

void boottestModel::Initialize() {  // for efficiency when varying r0 repeatedly to make CI, do stuff once that doesn't depend on r0
  real colvector sortID, o, _FEID
	real rowvector val, ClustCols
	real matrix RAll, L, LAll, vec, Combs, t, IDErr
	real scalar i, j, c, minN, sumN, _reps, i_FE, qq0
	pointer (real matrix) scalar _pR0, pIDAll
	class AnalyticalModel scalar M_WRE
	pragma unset vec; pragma unset val; pragma unused M_WRE

  Nobs = rows(*pXEx)
  NClustVar = cols(*pID)
  k  = cols(*pR0)
  kEx = cols(*pXEx)
  if (cols(*pZExcl)==0) pZExcl = &J(Nobs,0,0)
  if (cols(*pXEnd)==0) pXEnd = &J(Nobs,0,0)
  D = cols(*pXEnd) + 1
  REst = rows(*pR)  // base model contains restrictions?
  if (pZExcl != NULL) el = cols(*pZExcl) + kEx  // l = # of exogenous variables
  if (K==.) K = cols(*pZExcl)>0  // if K in K-class estimation not specified, it's 0 or 1 for OLS or 2SLS
  GMM = pW != NULL
  IV = K & GMM==0
  WRE = (K & scoreBS==0) | AR
  WREnonAR = WRE & AR==0

  if (weights = rows(*pwt)>1)
    sumwt = sum(*pwt)
  else
    pwt = &(sumwt = 1)
  _Nobs = weights & wttype=="fweight"? sumwt : Nobs

  if (WREnonAR)
    if (NClustVar)
      infoBootData = _panelsetup(*pID, 1..NBootClustVar, IDBootData)
    else
      pinfoErrData = &(infoBootData = J(Nobs,0,0))
  else if (NClustVar) {
    if (NClustVar > NBootClustVar)  // bootstrap Cluster grouping defs rel to original data
      infoBootData = _panelsetup(*pID, 1..NBootClustVar)
    else
      infoBootData = _panelsetup(*pID, 1..NClustVar)
  } else
    pinfoErrData = pinfoAllData = &(infoBootData = J(Nobs,0,0))  // causes no collapsing of data in _panelsum() calls, only multiplying by weights if any
  NBootClust = rows(infoBootData)

  if (bootstrapt) {
    if (NClustVar) {
      minN = .; sumN = 0

      Combs = combs(NErrClust)  // represent all error clustering combinations. First is intersection of all error clustering vars
      Clust = structboottestClust(rows(Combs)-1)  // leave out no-cluster combination
      NErrClustCombs = length(Clust)
      subcluster = NClustVar - NErrClust

      pinfoAllData = NClustVar > NBootClustVar? &_panelsetup(*pID,            1..NClustVar) : &infoBootData  // info for grouping by intersections of all bootstrap & clustering vars wrt data; used to speed crosstab EZVR0 wrt bootstrapping cluster & intersection of all error clusterings
      pinfoErrData = NClustVar > NErrClust    ? &_panelsetup(*pID, subcluster+1..NClustVar) : pinfoAllData  // info for intersections of error clustering wrt data
       IDErr = rows(*pinfoErrData)==Nobs? *pID :   (*pID)[(*pinfoErrData)[,1],]   // version of ID matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
      pIDAll = rows(*pinfoAllData)==Nobs?  pID : &((*pID)[(*pinfoAllData)[,1],])  // version of ID matrix with one row for each all-bootstrap & error cluster-var intersection instead of 1 row for each obs

      if (subcluster) {  // for subcluster bootstrap, bootstrapping cluster is not among error clustering combinations
        pBootClust = &(structboottestClust())
        pBootClust->info = _panelsetup(*pIDAll, 1..NBootClustVar)  // bootstrapping cluster info w.r.t. all-bootstrap & error-cluster intersections
        NBootClust = rows(pBootClust->info)
      } else
        pBootClust = &(Clust[2^(NClustVar - NBootClustVar)])  // location of bootstrap clustering within list of cluster combinations

      for (c=1; c<=NErrClustCombs; c++) {  // for each error clustering combination
        ClustCols = subcluster :+ boottest_selectindex(Combs[c,])
        Clust[c].even = mod(cols(ClustCols),2)

        if (c == 1)
          if (subcluster) {
            IDErr = IDErr[ Clust.order = stableorder(IDErr, ClustCols), ]
            Clust.info       = _panelsetup(IDErr, ClustCols)
          } else
            Clust.info       = J(rows(*pinfoAllData),0,0)  // causes no collapsing of data in _panelsum() calls
        else {
          if (any( Combs[|c, min(boottest_selectindex(Combs[c,] :!= Combs[c-1,])) \ c,.|])) // if this sort ordering same as last to some point and missing thereafter, no need to re-sort
            IDErr = IDErr[ Clust[c].order = stableorder(IDErr, ClustCols), ]

          Clust[c].info      = _panelsetup(IDErr, ClustCols)
        }

        Clust[c].N           = rows(Clust[c].info)
        sumN = sumN + Clust[c].N

        if (small) {
          Clust[c].multiplier = Clust[c].N/(Clust[c].N-1)
          if (Clust[c].N < minN) minN = Clust[c].N
        } else
          Clust[c].multiplier = 1
      }

      if (scoreBS | WREnonAR==0)
        ClustShare = weights? _panelsum(*pwt, *pinfoErrData)/sumwt : ((*pinfoErrData)[,2]-(*pinfoErrData)[,1]:+ 1)/Nobs // share of observations by group 

    } else {  // if no clustering, cast "robust" as clustering by observation
      pBootClust = &(Clust = structboottestClust())
      Clust.multiplier = small? _Nobs / (_Nobs - 1) : 1
      Clust.even = 1
      sumN = Clust.N = Nobs
      Clust.info = J(Nobs, 0, 0)  // signals _panelsum not to aggregate
      NErrClustCombs = 1
      if (scoreBS | WREnonAR==0)
        ClustShare = weights? *pwt/sumwt : 1/_Nobs
    }

    purerobust = NClustVar & (scoreBS | subcluster)==0 & NBootClust==Nobs  // do we ever error-cluster *and* bootstrap-cluster by individual?
    granular   = NClustVar & scoreBS==0 & (purerobust | 5*Nobs*k+1/25*2*NBootClust*k^2+1/25*2*NBootClust^2*k+NBootClust+1/25*2*NBootClust^2*reps+2*NBootClust*reps > 3*Nobs*reps+Nobs*k+1/25*(2*NBootClust*k*reps+2*k*reps-2*k*NBootClust-2*NBootClust*reps+2*NBootClust*k*reps)+2*NBootClust*reps)
    doKK = MaxMatSize<. | (willplot==0 & granular==0 & NBootClust * (sumN + NBootClust*reps +.5*(NErrClustCombs)) < 2*reps*(2*sumN + NErrClustCombs)) // estimate compute time for U :* (sum_c KK) * U vs. sum_c(colsum((J_c*U):*(J_c*U)))

    if (robust & purerobust==0) {
      if (subcluster | granular)
        infoErrAll = _panelsetup(*pIDAll, subcluster+1..NClustVar)  // info for error clusters wrt data collapsed to intersections of all bootstrapping & error clusters; used to speed crosstab EZVR0 wrt bootstrapping cluster & intersection of all error clusterings
      if (scoreBS & reps)
        J_ClustN_NBootClust = J(Clust.N, NBootClust, 0)
    }
  } else
    minN = rows(infoBootData)

  if (NFE) {
    sortID = (*pFEID)[o = stableorder(*pFEID, 1)]
    i_FE = 1; FEboot = reps>0; j = Nobs; _FEID = wtFE = J(Nobs, 1, 1)
    FEs = structFE(NFE)
    for (i=Nobs-1;i;i--) {
      if (sortID[i] != sortID[i+1]) {
        FEs[i_FE].is = o[|i+1\j|]
        if (weights) {
          t  = (*pwt)[FEs[i_FE].is]
          FEs[i_FE].wt = t / colsum(t)
        } else
          FEs[i_FE].wt = J(j-i, 1, 1/(j-i))
        if (reps & robust & granular < NErrClust)
          wtFE[FEs[i_FE].is] = FEs[i_FE].wt

        j = i
        
        if (reps & FEboot & NClustVar) {  // are all of this FE's obs in same bootstrapping cluster? (But no need to check if reps=0 for then CT_WE in 2nd term of (62) orthogonal to u = col of 1's)
          t = (*pID)[FEs[i_FE].is, 1..NBootClustVar]
          FEboot = all(t :== t[1,])
        }
        ++i_FE
      }
      _FEID[o[i]] = i_FE
    }
    FEs[NFE].is = FEs[NFE].is = o[|.\j|]
    if (weights) {
      t  = (*pwt)[FEs[NFE].is]
      FEs[NFE].wt = t / colsum(t)
    } else
      FEs[NFE].wt = J(j-i,1,1/(j-i))
    if (reps & robust & granular < NErrClust)
      wtFE[FEs[NFE].is] = FEs[NFE].wt

    if (reps & FEboot & NClustVar) {  // are all of this FE's obs in same bootstrapping cluster?
      t = (*pID)[FEs[NFE].is, 1..NBootClustVar]
      FEboot = all(t :== t[1,])
    }

    pFEID = &_FEID  // ordinal fixed effect ID

    if (robust & FEboot==0 & granular < NErrClust & reps & FEboot==0 & bootstrapt)
      infoBootAll = _panelsetup(*pIDAll, 1..NBootClustVar)  // info for bootstrapping clusters wrt data collapsed to intersections of all bootstrapping & error clusters
  }

  if (granular & (WREnonAR & purerobust)==0)
    if (NFE)
      (void) _panelsetup(*pID   , 1..NBootClustVar, IDBootData)
    else
      (void) _panelsetup(*pIDAll, 1..NBootClustVar, IDBootAll )

  if (ML)
    df = rows(*pR0)
  else {
    if (REst) {  // restricted estimation, e.g., cnsreg?
      symeigensystem(*pR ' boottest_lusolve(*pR * *pR', *pR), vec, val)  // make "inverse" S,s of constraint matrices; formulas adapted from [P] makecns
      L = vec[|.,.\.,q|]  // eigenvectors not in kernel of projection onto R
      S = q < k? vec[|.,q+1\.,.|] : J(k,0,0)  // eigenvectors in kernel
      s = L * luinv(*pR * L) * *pr
      if (AR) {
        SAR = blockdiag(S[|.,. \ k-cols(*pXEnd) , k-q-cols(*pXEnd)|] , I(cols(*pZExcl))) // adapt S,s from XExog, XEndog to XExog, ZEXcl. Assumes no constraints link XExog and XEndog
        sAR = s[|.\k-cols(*pXEnd)|] \ J(cols(*pZExcl),1,0)
      }
    }

    // Estimation with null imposed along with any model constraints; in IV, Z is unconstrained regardless of overlap with potentially constrained X
    _pR0 = null? pR0 : &J(0, k, 0)
    RAll = REst? *pR \ *_pR0 : *_pR0  // combine model and hypothesis constraints to prepare to "invert" them as a group too
    if (qq0 = rows(RAll)) {
      LAll = invsym(RAll * RAll')
      if (!all(diagonal(LAll)))
        _error(111, "A null hypothesis constraint is inconsistent or redundant.")
      symeigensystem(RAll ' LAll * RAll, vec, val)
      LAll  = vec[|.,. \ .,qq0|]
      SAll = qq0 < cols(vec)? vec[|.,qq0+1 \ .,.|] : J(rows(vec), 0, 0)
      LAll_invRAllLAll = LAll * luinv(RAll * LAll)
    } else
      SAll = J(0,0,0)

    M_DGP.parent = &this
    M_DGP.InitExog()

    if (WRE) {
      pM_Repl = &(M_WRE = M_DGP)
      pM_Repl->isDGP = 0
      M_DGP.LIML = 1; M_DGP.Fuller = 0; M_DGP.K = 1
      if (AR==0) { 
        pM_Repl->LIML = this.LIML; pM_Repl->Fuller = this.Fuller; pM_Repl->K = this.K
      }
      pM_Repl->SetAR(AR)
      pM_Repl->SetS(AR? SAR : S)
    } else {
      M_DGP.LIML = this.LIML; M_DGP.Fuller = this.Fuller; M_DGP.K = this.K
    }

    M_DGP.InitEndog(pY, pXEnd)

    if (AR) {
      if (willplot) {  // for plotting/CI purposes get original point estimate if not normally generated
        M_DGP.SetS(S)  // no-null model in DGP
        M_DGP.InitEstimate()
        M_DGP.Estimate(s)
        confpeak = *pR0 * M_DGP.beta  // estimated coordinate of confidence peak
      }
      pR0 = &(J(cols(*pZExcl),kEx,0), I(cols(*pZExcl)))  // for AR test, picks out coefs on excluded exogenous variables
    }
    df = rows(*pR0)
    
    if (granular) {
      euZVR0 = smatrix(df)
      pX = AR? &(*pXEx, *pZExcl) : pXEx
    }

    M_DGP.SetS(SAll)  // (potentially) constrained model in DGP; SAll imposes constraints, pR0 tests hypotheses on results
    M_DGP.InitEstimate()

    if (AR) {
      k = el
      K = 0
      pM = pM_Repl
    } else {
      M_DGP.InitTestDenoms(S)
      pM = &M_DGP
    }
  }

  if (bootstrapt) {
    denom = smatrix(df,df)
    if (WREnonAR==0 & robust) {
      if (reps) Kd = smatrix(df)
      Kcd = smatrix(NErrClustCombs, df)
      if (doKK)
        KK = smatrix(df, df)
      else
        pJcd = reps? &smatrix(NErrClustCombs, df) : &Kcd  // if reps = 0, Kcd will be multiplied by u, which is all 1's, and will constitute Jcd
    }
  }

  if (WREnonAR & NClustVar)
    XZi = eZi = smatrix(NBootClust)

  if (small) df_r = NClustVar? minN - 1 : _Nobs - k - NFE

  if (df==1) setsqrt(1)  // work with t/z stats instead of F/chi2

  if (small)
    multiplier = (smallsample = (_Nobs - k - NFE) / (_Nobs - robust)) / df  // divide by # of constraints because F stat is so defined
  else
    multiplier = smallsample = 1

  if (GMM)
    multiplier = multiplier / _Nobs  // variance of test stat uses GMM weight matrix from ivreg2, ivregress, divided by N
  else if ((robust | ML)==0)
    multiplier = multiplier * _Nobs  // will turn sum of squared errors in denom of t/z into mean; but division by N already in GMM weight matrix so skip for GMM
  
  if (sqrt) multiplier = sqrt(multiplier)

  enumerate = 0
  if (weighttype<=1 & reps)
    if (NBootClust*ln(2) < ln(reps)+1e-6)
      if (weighttype==1) {
        if (quietly==0)
          printf("\nWarning: with %g Clusters, the number of replications, %g, exceeds the universe of Mammen draws, 2^%g = %g. \nConsider Webb weights instead, using {cmd:weight(webb)}.\n", NBootClust, reps, NBootClust, 2^NBootClust) 
      } else {
        if (quietly==0)
          printf("\nWarning: with %g Clusters, the number of replications, %g, exceeds the universe of Rademacher draws, 2^%g = %g. Sampling each once. \nConsider Webb weights instead, using {cmd:weight(webb)}.\n", NBootClust, reps, NBootClust, 2^NBootClust)
        enumerate = 1
        MaxMatSize = .
      }
    else if (quietly==0 & (NBootClust>12 | weighttype) & floor(level/100 * (reps+1)) != level/100 * (reps+1)) {
        printf("\nNote: The bootstrap usually performs best when the confidence level (here, %g%%)\n", level)
        printf("      times the number of replications plus 1 (%g+1=%g) is an integer.\n", reps, reps+1)
      }

  NWeightGrps = MaxMatSize == .? 1 : ceil((reps+1) * max((rows(IDBootData), rows(IDBootAll), NBootClust)) * 8 / MaxMatSize / 1.0X+1E) // 1.0X+1E = giga(byte)
  betadev = smatrix(NWeightGrps)
  if (NWeightGrps == 1) {
    makeWildWeights(reps, 1)  // make all wild weights, once
    if (enumerate) reps = cols(u) - 1  // replications reduced to 2^G
    NWeightGrps    = 1
    WeightGrpStart = 1
    WeightGrpStop  = reps+1
  } else {
    seed = rseed()
    _reps = ceil((reps+1) / NWeightGrps)
     WeightGrpStart = (0::NWeightGrps-1) * _reps :+ 1
    (WeightGrpStop  = (1::NWeightGrps  ) * _reps     )[NWeightGrps] = reps+1
  }

  if (bootstrapt & (WREnonAR | df>1 | MaxMatSize<.)) // unless nonWRE or df=1 or splitting weight matrix, code will create Dist element-by-element, so pre-allocate vector now
    pDist = &J(reps+1, 1, .)
  if (NWeightGrps>1 | WREnonAR | (null==0 & df<=2))
    numer = J(df, reps+1, .)
    
  interpolative = reps & WREnonAR==0 & null & scoreBS==0
  interpolate_denom = interpolative & robust
  interpolate_e = interpolative & ((interpolate_denom & doKK==0 & reps & granular) | ((robust==0 | GMM) & (ML | GMM)==0))  // doesn't look right cases where makeNonWREStats refers directly to residuals, which may thus need interpolation
  if (interpolative) {
    /*if (doKK) {
      pKJcd = &Kcd  // interpolate K matrices
      pfncross  = &boottestCrossK()
      pfndouble = &boottestDoubleK()
    } else {
      pKJcd = pJcd  // interpolate J matrices
      pfncross  = &boottestCrossJ()
      pfndouble = &boottestDoubleJ()
    }*/
    dnumerdr = smatrix(dH0)
    if (interpolate_e) dedr = dnumerdr
    if (interpolate_denom) {
      ddenomdr = dJcddr = ssmatrix(dH0)
      ddenomdr2 = ssmatrix(dH0, dH0)
      for (i=dH0;i;i--) {
        ddenomdr[i].M = smatrix(df,df)
        for (j=i;j;j--)
          ddenomdr2[i,j].M = ddenomdr[i].M
        dJcddr[i].M = smatrix(NErrClustCombs, df)
      }
    }
  }
}

// main routine
void boottestModel::boottest() {
	real scalar g

	if (initialized==0) Initialize()

	makeStuffLinearInr0()  // make stuff that depends linearly on r0, possibly by interpolating

  if (doKK & interpolating==0 & WREnonAR==0)
    makeKK()

	if (NWeightGrps > 1)
		rseed(seed)

  for (g=1; g<=NWeightGrps; g++) {  // do group 1 first because it includes col 1, which is all that might need updating in constructing CI in WCU
		if (NWeightGrps > 1)
			makeWildWeights(WeightGrpStop[g] - WeightGrpStart[g] + (g>1), g==1)

		if (WREnonAR) {
			if (makeWREStats(WeightGrpStart[g], WeightGrpStop[g]))
				return
			if (bootstrapt==0)
				makeBootstrapcDenom(WeightGrpStart[g], WeightGrpStop[g])
		} else {
			if (makeNonWRENumers())
				return
			if (bootstrapt)
				makeNonWREStats(g, WeightGrpStart[g], WeightGrpStop[g])
			else
				makeBootstrapcDenom(WeightGrpStart[g], WeightGrpStop[g])
		}
	}

	repsFeas = (*pDist)[1]==.? 0 : colnonmissing(*pDist) - 1
	DistCDR = J(0,0,0)
	setdirty(0)
	initialized = 1
}

// compute bootstrap-c denominator from all bootstrap numerators
void boottestModel::makeBootstrapcDenom(real scalar thisWeightGrpStart, real scalar thisWeightGrpStop) {
	real colvector t
  if (thisWeightGrpStart == 1) {
		t = numer[,1]
    statDenom = numer * numer' - t * t'
		numersum = rowsum(numer) - t
	} else {
		statDenom = statDenom + numer * numer'
		numersum = numersum + rowsum(numer)
	}
	if (thisWeightGrpStop > reps) {  // last weight group?
		statDenom = (statDenom - numersum * numersum' / reps) / reps
		pDist = &((sqrt? numer:/sqrt(statDenom) : colsum(numer :* boottest_lusolve(statDenom, numer)))')
	}
}

// draw wild weight matrix of width _reps. If first=1, insert column of 1s at front
void boottestModel::makeWildWeights(real scalar _reps, real scalar first) {
	if (_reps) {
		if (enumerate)
			u = J(NBootClust,1,1), count_binary(NBootClust, -1-WREnonAR, 1-WREnonAR)  // complete Rademacher set
		else if (weighttype==3)
			u = rnormal(NBootClust, _reps+1, -WREnonAR, 1)  // normal weights
		else if (weighttype==4)
			u = rgamma(NBootClust, _reps+1, 4, .5) :- (2 + WREnonAR)  // Gamma weights
		else if (weighttype==2)
			if (WREnonAR)
				u = sqrt(2 * ceil(runiform(NBootClust, _reps+first) * 3)) :* ((runiform(NBootClust, _reps+1):>=.5):-.5) :- 1  // Webb weights, minus 1 for convenience in WRE
			else {
				u = sqrt(    ceil(runiform(NBootClust, _reps+first) * 3)) :* ((runiform(NBootClust, _reps+1):>=.5):-.5)       // Webb weights, divided by sqrt(2)
				u_sd = 1.6a09e667f3bcdX-001 /*sqrt(.5)*/
			}
		else if (weighttype)
			if (WREnonAR)
				u = ( rdiscrete(NBootClust, _reps+first, 1.727c9716ffb76X-001\1.1b06d1d200914X-002 /*.5+sqrt(.05)\.5-sqrt(.05)*/) :- 1.5 ) * 1.1e3779b97f4a8X+001 /*sqrt(5)*/ :- .5  // Mammen weights, minus 1 for convenience in WRE
			else {
				u = ( rdiscrete(NBootClust, _reps+first, 1.727c9716ffb76X-001\1.1b06d1d200914X-002 /*.5+sqrt(.05)\.5-sqrt(.05)*/) :- 1.5 ) :+ 1.c9f25c5bfedd9X-003 /*.5/sqrt(5)*/  // Mammen weights, divided by sqrt(5)
				u_sd = 1.c9f25c5bfedd9X-002 /*sqrt(.2)*/
			}
		else if (WREnonAR) {
      u = runiform(NBootClust, _reps+first) :<  .5; u = (-2) * u  // Rademacher weights, minus 1 for convenience in WRE
    } else {
      u = runiform(NBootClust, _reps+first) :>= .5; u = u :- .5   // Rademacher weights, divided by 2
      u_sd = .5
    }

		if (first)
			u[,1] = J(NBootClust, 1, WREnonAR? 0 : u_sd)  // keep original residuals in first entry to compute base model stat
	} else
		u = J(0,1,0)  // in places, cols(u) indicates number of reps -- 1 for classical tests
}


real scalar boottestModel::makeWREStats(real scalar thisWeightGrpStart, real scalar thisWeightGrpStop) {
	pragma unused thisWeightGrpStop
	real scalar c, j, i
	real colvector _e, _beta, betaEnd, _u, numer_j
	real matrix Subscripts, Zi, AVR0, t, XExi
	pointer (real matrix) scalar peZVR0

	if (initialized & null==0) {  // if not imposing null and we have returned, then df=1 or 2; we're plotting and only test stat, not distribution, changes with r0
		numer[,1] = *pR0 * pM_Repl->beta - *pr0
		(*pDist)[1] = df==1? numer[1] / sqrt(denom.M) : numer[,1] ' boottest_lusolve(denom.M, numer[,1])
		return(1)  // no need to continue calc
	}

	if (thisWeightGrpStart == 1) {  // first/only weight group? initialize a couple of things
		_e = M_DGP.e + M_DGP.e2 * M_DGP.beta[|kEx+1\.|]

		if (NClustVar & NFE==0)  // prep for optimized computation for bootstrapping cluster when no FE
			for (i=NBootClust; i; i--) {
				Subscripts = (infoBootData)[i,]', (.\.)
				XExi = kEx? (*pXEx)[|Subscripts|] : J(Subscripts[2,1]-Subscripts[1,1]+1,0,0)
				Zi = XExi , (*pZExcl)[|Subscripts|]  // inefficient?
				if (weights) Zi = Zi :* (*pwt)[|Subscripts|]
				XZi[i].M = cross(XExi, Zi) \ cross((*pXEnd)[|Subscripts|], Zi) \ cross((*pY)[|Subscripts|], Zi)
				eZi[i].M =                   cross(M_DGP.e2[|Subscripts|], Zi) \ cross(   _e[|Subscripts|], Zi)
			}
	}

	for (j=cols(u); j; j--) {  // WRE bootstrap
		_u = u[IDBootData,j]

		pM_Repl->InitEndog(&(*M_DGP.pY:+_e:*_u) , &(*pXEnd:+M_DGP.e2:*_u))
		pM_Repl->InitEstimate()
		pM_Repl->InitTestDenoms(S)  // prepare for replication regressions, null not imposed
		pM_Repl->Estimate(s)

		numer_j = null | thisWeightGrpStart==1 & j==1? *pR0 * pM_Repl->beta - *pr0 : *pR0 * (pM_Repl->beta - M_DGP.beta0)

		if (bootstrapt) {
			if (robust) {  // Compute denominator for this WRE test stat
				denom = smatrix()
				if (NClustVar != 1 | NFE)  // collapse meat+sandwich  to all-Cluster-var intersections. If no collapsing needed, _panelsum() will still fold in any weights
					peZVR0 = &_panelsum(pM_Repl->ZVR0, weights? *pwt :* pM_Repl->e : pM_Repl->e, *pinfoErrData)  // really eZAVR0, where e is wildized residual, not residual from replication fit (estar)
				for (c=1; c<=NErrClustCombs; c++) {
					if (NClustVar != 1 & rows(Clust[c].order))
						peZVR0 = &((*peZVR0)[Clust[c].order,])
					if (*pBootClust==Clust[c] & NClustVar & NFE==0) {  // optimized computation for bootstrapping Cluster when no FE
						AVR0 = pM_Repl->A * pM_Repl->VR0; _beta = -pM_Repl->beta \ 1; betaEnd = _beta[|kEx+1\.|]
						pragma unset t
						for (i=1; i<=NBootClust; i++) {
							pG = &((_beta'XZi[i].M + betaEnd'eZi[i].M * u[i,j]) * AVR0) // R0 * V * Z_i'estar_i
							t = i==1? cross(*pG,*pG) : t + cross(*pG,*pG)
						}
            _clustAccum(denom.M, c, t)
					} else {
						pG = &_panelsum(*peZVR0, Clust[c].info)
            _clustAccum(denom.M, c, cross(*pG,*pG))
					}
				}
			} else
				denom.M = (*pR0 * pM_Repl->VR0) * pM_Repl->eec

 			(*pDist)[j+thisWeightGrpStart-1] = sqrt? numer_j/sqrt(denom.M) : cross(numer_j, boottest_lusolve(denom.M, numer_j))
		}
		numer[,j+thisWeightGrpStart-1] = numer_j  // slight inefficiency: in usual bootstrap-t case, only need to save numerators in numer if getdist("numer") is coming because of svmat(numer)
	}

	if (thisWeightGrpStart==1 & bootstrapt) statDenom = denom.M  // original-sample denominator

	return(0)  // don't skip remaining calcs in calling procedure
}


// Construct stuff that depends linearly on r0, possibly by interpolation
void boottestModel::makeStuffLinearInr0() {
	real scalar h1, h2, d1, d2, c; real matrix t; real colvector Delta, newAnchor

  if (interpolative) {
    if (rows(anchor)==0) {  // first call? save current r0 as anchor for interpolation
      _makeStuffLinearInr0(anchor = *pr0)
      numer0 = numer
      if (interpolate_e) e0 = *pe
      if (interpolate_denom) Jcd0 = *pJcd
      return
    }

    if (rows(poles) == 0) {  // not enough info yet to interpolate
      poles = *pr0 - anchor  // second call: from anchor make set of orthogonal poles, which equal anchor except in one dimension
      if (interpolate_denom)  // grab quadratic denominator from *previous* (1st) evaluation
        denom0 = /*doKK? KK :*/ denom
      newAnchor = J(dH0,1,1)  // all anchors new
    } else
    	newAnchor = abs(*pr0 - anchor) :> 2 * abs(poles)  //  been here at least twice: interpolate unless current r0 stretches range > 2X in some dimension(s)

    if (any(newAnchor)) {  // prep interpolation
      for (h1=1;h1<=dH0;h1++)
        if (newAnchor[h1]) {
        	poles[h1] = (*pr0)[h1] - anchor[h1]
        	(t = anchor)[h1] = (*pr0)[h1]  // if dH0>1 this creates anchor points that are not graphed, an inefficiency. But simpler to make the deviations from 1st point orthogonal

          _makeStuffLinearInr0(t)  // calculate linear stuff at new anchor

          dnumerdr[h1].M = (numer - numer0) / poles[h1]
          if (interpolate_e)
            dedr[h1].M = (*pe - e0) / poles[h1]
          if (interpolate_denom) {  // prepare to interpolate denominators. df > 1 for an AR test with >1 instruments. Quadratic interactions mean changing anchor forces recalc 
            for (c=rows(Kcd);c;c--)
              for (d1=df;d1;d1--)
                dJcddr[h1].M[c,d1].M = ((*pJcd)[c,d1].M - Jcd0[c,d1].M) / poles[h1]
            for (d1=df;d1;d1--)
              for (d2=1;d2<=d1;d2++) {
                for (c=1;c<=NErrClustCombs;c++) {
                                t =     /*(*pfncross)*/colsum( Jcd0       [c,d1].M :* dJcddr[h1].M[c,d2].M)
                  if (d1 != d2) t = t + /*(*pfncross)*/colsum(dJcddr[h1].M[c,d1].M :*  Jcd0       [c,d2].M)  // for diagonal items, faster to just double after the c loop
                  _clustAccum(ddenomdr[h1].M[d1,d2].M, c, t)
                }
                if (d1==d2) ddenomdr[h1].M[d1,d2].M = /*(*pfndouble)*/ ddenomdr[h1].M[d1,d2].M + ddenomdr[h1].M[d1,d2].M
                /*if (doKK)
                  _fold(ddenomdr[h1].M[d1,d2].M)*/
              }
          }
        }
      if (interpolate_denom)
        for (h1=1;h1<=dH0;h1++)
          for (h2=h1;h2;h2--)
            if (newAnchor[h1] | newAnchor[h2])
              for (d1=df;d1;d1--)
                for (d2=d1;d2;d2--) {
                  for (c=1;c<=NErrClustCombs;c++)
                    _clustAccum(ddenomdr2[h1,h2].M[d1,d2].M, c, /*(*pfncross)*/colsum(dJcddr[h1].M[c,d1].M :* dJcddr[h2].M[c,d2].M))
                  /*if (doKK) _fold(ddenomdr2[h1,h2].M[d1,d2].M)*/
                }

      Delta = poles
      interpolating = 1

    } else {  // routine linear interpolation if the anchors not moved

      Delta = *pr0 - anchor  // interpolate!
      numer = numer0 + dnumerdr.M * Delta[1]; if (dH0 > 1) numer = numer + dnumerdr[2].M * Delta[2]
      if (interpolate_e) {
        pe = &(e0 + dedr.M * Delta[1]); if (dH0 > 1) pe = &(*pe + dedr[2].M * Delta[2])
      }
    }

    if (interpolate_denom)  // even if an anchor was just moved, and linear components just compued from scratch, we do the quadratic interpolation now, from the updated linear factors
      if (dH0==1)
        for (d1=df;d1;d1--)
          for (d2=d1;d2;d2--)
                (/*doKK? KK :*/ denom)[d1,d2].M = denom0[d1,d2].M + ddenomdr.M[d1,d2].M * Delta + ddenomdr2.M[d1,d2].M * (Delta * Delta)
      else  // dH0==2
        for (d1=df;d1;d1--)
          for (d2=d1;d2;d2--)
            (/*doKK? KK :*/ denom)[d1,d2].M = denom0[d1,d2].M + 
                                          ddenomdr[1].M[d1,d2].M * Delta[1] + 
                                          ddenomdr[2].M[d1,d2].M * Delta[2] + 
                                          ddenomdr2[1,1].M[d1,d2].M * (Delta[1] * Delta[1]) + 
                                          ddenomdr2[2,1].M[d1,d2].M * (Delta[1] * Delta[2]) + 
                                          ddenomdr2[2,2].M[d1,d2].M * (Delta[2] * Delta[2])
    return
  }

  _makeStuffLinearInr0(*pr0)  // still here? non-interpolative case
}

// Construct stuff that depends linearly on r0 and doesn't depend on u. Only does one bit of WRE prep. No interpolation.
void boottestModel::_makeStuffLinearInr0(real colvector r0) {
  pointer (real matrix) scalar pewt; real scalar g, d, i, c; real matrix t; pointer (real matrix) scalar peZVR0, pt; real colvector rAll

  if (ML==0) {
    if (AR) {
      pM_Repl->InitEndog(&(*pY - *pXEnd * r0), NULL, &(*M_DGP.pZExclY - M_DGP.ZExclXEnd * r0), &(*M_DGP.pXExY - M_DGP.XExXEnd * r0))
      pM_Repl->InitEstimate()
      pM_Repl->Estimate(sAR)
      pM_Repl->InitTestDenoms(SAR)
    } else if (initialized==0 | null) {  // don't need to recompute if we're not imposing the null  -- *** THIS BIT RUNS FOR WRE too
      rAll = null? r0 : J(0, 1, 0); if (REst) rAll =  *pr \ rAll  // constant terms of model + null constraints
      sAll = rows(rAll) ? LAll_invRAllLAll * rAll : J(0,1,0)
      M_DGP.Estimate(sAll)
    }
    pe = &(pM->e)
  }
	if (WREnonAR) return

  if (initialized==0 | null) {  // if are imposing null or we are not, but this is first call, then build stuff
    if (ML)
      eZVR0 = *pSc * (VR0 = *pV * *pR0')
    else if (scoreBS | (robust & granular < NErrClustCombs))
      eZVR0 = *pe :* pM->ZVR0
		if (scoreBS)
			SewtXV = reps? (NClustVar? _panelsum(eZVR0, *pwt, infoBootData) : (weights?       eZVR0 :* *pwt  :        eZVR0) ) :
                                                                        (weights? cross(eZVR0,   *pwt) : colsum(eZVR0)')
    else {  // same calc as in score BS but broken apart to grab intermediate stuff, and assuming residuals defined
      pewt = weights? &(*pe :* *pwt) : pe
      pt = &_panelsum(*pXEx, *pewt, infoBootData)
      if (AR)
        pt = &(*pt, _panelsum(*pZExcl, *pewt, infoBootData))
      SewtXV = *pM->pV * (*pt)'
    }

    if (NWeightGrps == 1)
      if (scoreBS)
        numer = reps? cross(SewtXV, u) : SewtXV * u_sd
      else if (robust==0 | granular)
			  numer = *pR0 * (betadev.M = SewtXV * u)
      else
        numer = (*pR0 * SewtXV) * u
    else
      if (scoreBS)
        for (g=NWeightGrps; g; g--)
          numer[|WeightGrpStart[g] \ WeightGrpStop[g]|] = reps? cross(SewtXV, u) : SewtXV * u_sd
      else if (robust==0 | granular)
        for (g=NWeightGrps; g; g--)
		      numer[|WeightGrpStart[g] \ WeightGrpStop[g]|] = *pR0 * (betadev[g].M = SewtXV * u)
      else
        for (g=NWeightGrps; g; g--)
	        numer[|WeightGrpStart[g] \ WeightGrpStop[g]|] = (*pR0 * SewtXV) * u
  }

  if      ( AR  ) numer[,1] = u_sd * pM->beta[|kEx+1\.|] // coefficients on excluded instruments in AR OLS
  else if (!null) numer[,1] = u_sd * (*pR0 * (ML? beta : pM->beta) - r0) // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this.

  if (initialized & null==0) return

  if (robust & GMM==0 & granular < NErrClustCombs) {
    peZVR0 = &_panelsum(eZVR0, *pwt, *pinfoAllData)  // collapse data to all-boot & error-cluster-var intersections. If no collapsing interpolate_ed, _panelsum() will still fold in any weights
    if (reps) {
      if (scoreBS)
        for (d=df;d;d--)
          Kd[d].M = J_ClustN_NBootClust  // inefficient, but we're not optimizing for the score bootstrap
      else
        for (d=df;d;d--) {
          t = pM->ZVR0[,d]; if (weights) t = t :* *pwt
          if (AR)  // final term in (64) in paper, for c=intersection of all error clusters
            Kd[d].M = (_panelsum(*pXEx, t, *pinfoErrData), _panelsum(*pZExcl, t, *pinfoErrData)) * SewtXV
          else
            Kd[d].M =  _panelsum(*pXEx, t, *pinfoErrData                                       ) * SewtXV
        }

      if (NFE & FEboot==0)
        CT_WE = _panelsum(crosstab(wtFE :* *pe), infoBootAll)'

      // subtract crosstab of E:*ZVR0 wrt bootstrapping cluster combo and all-cluster-var intersections
      if (*pBootClust == Clust[1])  // crosstab c,c* is square
        for (d=df;d;d--) {  // if bootstrapping on all-cluster-var intersections (including one-way clustering), the base crosstab is diagonal
          for (i=Clust.N;i;i--)
            Kd[d].M[i,i] = Kd[d].M[i,i] - (*peZVR0)[i,d]
          if (scoreBS)
            Kd[d].M = Kd[d].M :+ ClustShare * (*peZVR0)[,d]' // for score bootstrap, recenter; "+" because we subtracted *peZVR0

        }
      else
        if (subcluster) // crosstab c,c* is wide
          for (d=df;d;d--) {
            for (i=Clust.N;i;i--) {
              t = infoErrAll[i,]'
              Kd.M[|(i\i), t|] = Kd.M[|(i\i), t|] - (*peZVR0)[|t, (d\d)|]'
            }
            if (scoreBS)
	            Kd[d].M = Kd[d].M - ClustShare * colsum(Kd[d].M) // for score bootstrap, recenter
          }
        else // crosstab c,c* is tall
          for (d=df;d;d--) {
            for (i=NBootClust;i;i--) {
              t = pBootClust->info[i,]'
              Kd[d].M[|t, (i\i)|] = Kd[d].M[|t, (i\i)|] - (*peZVR0)[|t, (d\d)|]
            }
            if (scoreBS)
	            Kd[d].M = Kd[d].M - ClustShare * colsum(Kd[d].M) // for score bootstrap, recenter
          }

      for (c=1+granular; c<=NErrClustCombs; c++) {
        if (rows(Clust[c].order))
          for (d=df;d;d--)
            Kd[d].M = Kd[d].M[Clust[c].order,]
        for (d=df;d;d--) {
          Kcd[c,d].M = _panelsum(Kd[d].M, Clust[c].info)
          if (NFE & FEboot==0)
            Kcd[c,d].M = Kcd[c,d].M + _panelsum(pM->CT_ZVR0[c,d].M, Clust[c].info) * CT_WE
        }
      }
    } else {  // reps = 0. In this case, only 1st term of (64) is non-zero after multiplying by v* (= all 1's), and it is then a one-way sum by c
      if (scoreBS)
        peZVR0 = &(*peZVR0 :- ClustShare * colsum(*peZVR0))  // recenter if OLS

      for (c=1; c<=NErrClustCombs; c++) {
        if (rows(Clust[c].order))
          peZVR0 = &((*peZVR0)[Clust[c].order,])
        pt = &_panelsum(*peZVR0, Clust[c].info)  // when c=1 (unless subcluster bootstrap), two args have same # of rows, &_panelsum() returns 1st arg by reference. Using & then prevents unnecessary cloning.
        for (d=df;d;d--)
          Kcd[c,d].M = (*pt)[,d]
      }
    }
  }

  makeJ(1)  // compute J = K * u; if NWeightGrps > 1, then this is for 1st group
}

// prep K'K-based denominator computation loop
void boottestModel::makeKK() {
  real scalar i, j, c
  for (i=df;i;i--)
    for (j=i;j;j--) {
      for (c=1;c<=NErrClustCombs;c++)
        _clustAccum(KK[i,j].M, c, cross(Kcd[c,i].M, Kcd[c,j].M))
      _fold(KK[i,j].M)
    }
}

void boottestModel::makeJ(real scalar g) {
  real scalar c, d

  if (doKK | reps==0 | robust==0)
    return

  if (granular)  // prep optimized treatment when bootstrapping by many/small groups
    if (purerobust)
      eu = *M_DGP.partialFE(&(*pe :* u)) - *pX * betadev[g].M
    else {  // clusters small but not all singletons
      if (NFE) {
        eu = *M_DGP.partialFE(&(*pe :* u[IDBootData,]))
        for (d=df;d;d--)
          (*pJcd)[1,d].M = _panelsum(eu,            pM->WZVR0[d].M, *pinfoErrData)                               - _panelsum(*pX, pM->WZVR0[d].M, *pinfoErrData) * betadev[g].M
      } else
        for (d=df;d;d--)
          (*pJcd)[1,d].M = _panelsum(_panelsum(*pe, pM->WZVR0[d].M, *pinfoAllData) :* u[IDBootAll,], infoErrAll) - _panelsum(*pX, pM->WZVR0[d].M, *pinfoErrData) * betadev[g].M
    }

  for (c=NErrClustCombs; c>granular; c--)
    for (d=df;d;d--)
	      (*pJcd)[c,d].M = Kcd[c,d].M * u
}

real scalar boottestModel::makeNonWRENumers() {
  if (initialized & null==0) {  // if not imposing null and have returned, then thisWeightGrpStart=1, df=1 or 2; and distribution doesn't change with r0, only test stat
		(*pDist)[1] = df==1? numer[1] / sqrt(statDenom) : numer[,1] ' boottest_lusolve(statDenom, numer[,1])
		return(1)  // skip rest of calc
	}
  return(0)
}

void boottestModel::makeNonWREStats(real scalar g, real scalar thisWeightGrpStart, real scalar thisWeightGrpStop) {
	real scalar i, c, j, l; real matrix eueu, t; real colvector numer_l; pointer (real matrix) scalar pVR0; real rowvector t1, t2, t12

	if (robust & GMM==0) {
    if (doKK)  // core denominator computation loop
			for (i=df;i;i--)
				for (j=i;j;j--)
					denom[i,j].M = colsum(u :* KK[i,j].M * u)  // (60), 2nd version GGB/2+GB + 3GG < 3GB
		else if (!interpolating) {  // alternative core computational loop, avoiding computing K'K which has cubic time cost in numbers of bootstrapping clusters
			if (g > 1)
        makeJ(g)

      if (purerobust)
        eueu = eu :* eu
      for (i=df;i;i--)
        for (j=i;j;j--) {
          _clustAccum(denom[i,j].M, 1, purerobust? cross(pM->WZVR0[i].M, pM->WZVR0[j].M, eueu) : colsum((*pJcd)[1,i].M :* (*pJcd)[1,j].M))  // (60), 1st version
          for (c=2;c<=NErrClustCombs;c++)
            _clustAccum(denom[i,j].M, c, colsum((*pJcd)[c,i].M :* (*pJcd)[c,j].M))
        }
		}

		if (df == 1) {
			if (NWeightGrps > 1)
				(*pDist)[|thisWeightGrpStart \ thisWeightGrpStop|] =   (numer :/ sqrt(denom.M))'
			else
          pDist                                            = &((numer :/ sqrt(denom.M))')
			if (thisWeightGrpStart==1)
				statDenom = denom.M[1]  // original-sample denominator
		} else if (df==2) {  // hand-code 2D numer'inv(denom)*numer
    	t1 = numer[1,]; t2 = numer[2,]; t12 = t1:*t2
			if (NWeightGrps > 1)
				(*pDist)[|thisWeightGrpStart \ thisWeightGrpStop|] =   ( (t1:*t1:*denom[2,2].M - (t12+t12):*denom[2,1].M + t2:*t2:*denom[1,1].M) :/ (denom[1,1].M:*denom[2,2].M - denom[2,1].M:*denom[2,1].M) )'
			else                                                                                                                         
          pDist                                            = &(( (t1:*t1:*denom[2,2].M - (t12+t12):*denom[2,1].M + t2:*t2:*denom[1,1].M) :/ (denom[1,1].M:*denom[2,2].M - denom[2,1].M:*denom[2,1].M) )')
			if (thisWeightGrpStart==1)
				statDenom = denom[1,1].M[1], denom[2,1].M[1] \ denom[2,1].M[1], denom[2,2].M[1]  // original-sample denominator
    } else {  // build each replication's denominator from vectors that hold values for each position in denominator, all replications
			t = J(df,df,.)
			for (l=cols(u); l; l--) {
				for (i=df;i;i--)
					for (j=i;j;j--)
						t[i,j] = denom[i,j].M[l]
				_makesymmetric(t)
				numer_l = numer[,l]
				(*pDist)[l+thisWeightGrpStart-1] = numer_l ' boottest_lusolve(t, numer_l)  // in degenerate cases, cross() would turn cross(.,.) into 0
			}
			if (thisWeightGrpStart==1)
				statDenom = t  // original-sample denominator
		}

	} else { // non-robust or GMM

		pVR0 = ML? &VR0 : &(pM->VR0)
		if (df == 1) {  // optimize for one null constraint
			denom.M = *pR0 * *pVR0

			if ((ML | GMM)==0) {
                     eu = reps? u :* *pe : *pe
				if (scoreBS) eu = eu :- (weights? cross(ClustShare, eu) : colsum(eu) * ClustShare)  // Center variance if interpolate_ed
				        else eu = eu  - (*pXEx, *pZExcl) * betadev[g].M  // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
				denom.M = denom.M :* (weights? cross(*pwt, eu :* eu) : colsum(eu :* eu))
			}
			if (NWeightGrps > 1)
				(*pDist)[|thisWeightGrpStart \ thisWeightGrpStop|] =   (numer :/ sqrt(denom.M))'
			else
				  pDist                                            = &((numer :/ sqrt(denom.M))')
			if (thisWeightGrpStart==1)
				statDenom = denom.M[1]  // original-sample denominator
		} else {
			denom.M = *pR0 * *pVR0

			if (ML | GMM) {
				for (l=cols(u); l; l--) {
					numer_l = numer[,l]
					(*pDist)[l+thisWeightGrpStart-1] = cross(numer_l, boottest_lusolve(denom.M, numer_l)) 
				}
				if (thisWeightGrpStart==1)
					statDenom = denom.M  // original-sample denominator
			} else {
				for (l=cols(u); l; l--) {
					numer_l = numer[,l]
					(*pDist)[l+thisWeightGrpStart-1] = cross(numer_l, boottest_lusolve(denom.M, numer_l)) 
                       eu = reps? u[,l] :* *pe : *pe
					if (scoreBS) eu = eu :- (weights? cross(*pwt, eu) : colsum(eu)) * ClustShare  // Center variance if interpolate_ed
					        else eu = eu  - (*pXEx, *pZExcl) * betadev[g].M[,l]  // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
					(*pDist)[l+thisWeightGrpStart-1] = (*pDist)[l+thisWeightGrpStart-1] / (t = cross(eu, *pwt, eu))
				}
				if (thisWeightGrpStart==1)
					statDenom = denom.M * t  // original-sample denominator
			}
		}
	}
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

// Do panelsum() except that a single missing value in X doesn't make all results missing and
// efficiently handles case when all groups have one row.
real matrix _panelsum(real matrix X, real matrix arg2, | real matrix arg3) {
	if (args()==2) {
		if (rows(arg2)==0 | rows(arg2)==rows(X))
			return(X)
	} else if (rows(arg3)==0 | rows(arg3)==rows(X))
		return(arg2==1? X : X :* arg2) // if no collapsing called for, still fold in provided weights
	return (cols(arg3)? panelsum(X, arg2, arg3) : panelsum(X, arg2))
}

// given a vector, return indices of the non-zero elements, like selectindex() function added in Stata 13
// if v = 0 (so can't tell if row or col vector), returns rowvector J(1, 0, 0) 
real vector boottest_selectindex(real vector v) {
	real scalar rows
	if (v==0) return(J(1,0,0))
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

// cross-tab sum of a column vector w.r.t. intersection-of-error & bootstrap-clustering-vars and fixed-effect var
real matrix boottestModel::crosstab(real colvector v) {
	real matrix retval; real scalar i, j, t; real colvector _FEID, _v
	retval = J(rows(*pinfoAllData), NFE, 0)
	for (i=rows(*pinfoAllData);i;i--) {
		_FEID = panelsubmatrix(*pFEID, i, *pinfoAllData)
		_v    = panelsubmatrix(v     , i, *pinfoAllData)
		for (j=rows(_FEID);j;j--) {
			t = _FEID[j] 
			retval[i,t] = retval[i,t] + _v[j]
		}
	}
	return(retval)
}

// given a pre-configured boottest linear model with one-degree null imposed, compute distance from target p value of boostrapped one associated with given value of r0
// used with optimize() to construct confidence intervals
// performs no error checking
real scalar boottestModel::r0_to_p(real colvector r0) {
	pr0 = &r0
	setdirty(1, 1) // set dirty = 1, but leave initialized=0, which we want when only changing r0
	return (getpadj())
}

real scalar boottestModel::search(real scalar alpha, real scalar p_lo, real scalar lo, real scalar p_hi, real scalar hi) {
	real scalar mid, _p
	mid = lo + (alpha-p_lo)/(p_hi-p_lo)*(hi-lo)  // interpolate linearly
//	mid = lo + ( reldif(p_hi-p_lo, 1/repsFeas)<1e-6 & repsFeas? (hi - lo)/2 : (alpha-p_lo)/(p_hi-p_lo)*(hi-lo) ) // interpolate linearly until bracketing a single bump-up; then switch to binary search
	if (reldif(lo,mid)<ptol | reldif(hi,mid)<ptol | (repsFeas & abs(p_hi-p_lo)<(1+(ptype==1))/repsFeas*1.000001))
		return (mid)
	if ( ((_p = r0_to_p(mid)) < alpha) == (p_lo < alpha) )
		return(search(alpha, _p, mid, p_hi, hi))
	return(search(alpha, p_lo, lo, _p, mid))
}


// derive wild bootstrap-based CI, for case of linear model with one-degree null imposed.
void boottestModel::plot() {
	real scalar t, alpha, _quietly, c, d, i, j, ARTrialHalfWidth, p_lo, p_hi; real colvector lo, hi; pointer (real colvector) scalar _pr0; real rowvector b; real matrix V

	_quietly = quietly; _pr0 = pr0
	setquietly(1)

  _editmissing(gridpoints, 25)

  if (AR==0) {
    b = getb()
    V = getV()
  }

  if (dH0==2) {
    lo = hi = J(2, 1, .)

    if (AR)
      ARTrialHalfWidth = abs(confpeak) / (-2/invnormal(getpadj(1)/2))
    for(d=df;d;d--) {
      if (AR) {
        lo[d] = editmissing(gridmin[d], confpeak[d] - ARTrialHalfWidth[d])
        hi[d] = editmissing(gridmin[d], confpeak[d] + ARTrialHalfWidth[d])
      } else {
        lo[d] = editmissing(gridmin[d], b[d] + (*pr0)[d] - 3 * sqrt(V[d,d]))
        hi[d] = editmissing(gridmax[d], b[d] + (*pr0)[d] + 3 * sqrt(V[d,d]))
      }

      stata("_natscale " + strofreal(lo[d]) + " " + strofreal(hi[d]) + " 4")
      if (gridmin[d]==.) {
        stata("local min = r(min)")  // for some reason st_global("r(min)") doesn't work
        lo[d] = strtoreal(st_local("min"))
      }
      if (gridmax[d]==.) {
        stata("local max = r(max)")
        hi[d] = strtoreal(st_local("max"))
      }
    }
    plotX = (rangen(lo[1], hi[1], gridpoints[1]) # J(gridpoints[2],1,1)), (J(gridpoints[1],1,1) # rangen(lo[2], hi[2], gridpoints[2]))
    plotY = J(rows(plotX), 1, .)

  } else {  // 1D plot

    alpha = 1 - level*.01
    if (alpha > 0 & cols(u)-1 <= 1/alpha-1e6) {
      setquietly(_quietly)
      if (quietly==0) errprintf("\nError: need at least %g replications to resolve a %g%% two-sided confidence interval.\n", ceil(1/alpha), level)
      return
    }
    
    if (AR)
      ARTrialHalfWidth = abs(confpeak) / (small? invttail(df_r, alpha/2)/invttail(df_r, getpadj(1)/2) : invnormal(alpha/2)/invnormal(getpadj(1)/2))
    else
      confpeak = b + *pr0
    if (gridmin[1]==. | gridmax[1]==.) {
      if (reps)
        if (AR) {
          lo = editmissing(gridmin[1], confpeak - ARTrialHalfWidth)
          hi = editmissing(gridmax[1], confpeak + ARTrialHalfWidth)
        } else {
          (void) getdist()
          lo = editmissing(gridmin[1], confpeak + 5 * DistCDR[floor((   alpha/2)*(repsFeas-1))+1] * sqrt(V)) // initial guess based on distribution from main test
          hi = editmissing(gridmax[1], confpeak + 5 * DistCDR[ceil (( 1-alpha/2)*(repsFeas-1))+1] * sqrt(V))
        }
      else {
        t = sqrt(statDenom) * (small? invttail(df_r, alpha/2) : -invnormal(alpha/2))
        lo = editmissing(gridmin[1], confpeak - t)
        hi = editmissing(gridmax[1], confpeak + t)
        if (scoreBS & (null | willplot)==0) {  // if doing simple Wald test with no graph, we're done
          CI = lo, hi
          return
        }
      }

      if (gridmin[1]==. & ptype!=2)  // unless upper-tailed p value, try at most 10 times to bracket confidence set by symmetrically widening
        for (i=10; i & -(p_lo=r0_to_p(lo)) < -alpha; i--) {
          t = hi - lo
          lo = lo - t
          if (gridmax[1]==. & twotailed) hi = hi + t  // maintain rough symmetry unless user specified upper bound
        }
      if (gridmax[1]==. & ptype!=3)  // ditto for high side
	      for (i=10; i & -(p_hi=r0_to_p(hi)) < -alpha; i--) {
          t = hi - lo
          if (gridmin[1]==. & twotailed) lo = lo - t
          hi = hi + t
        }
    } else {
      lo = gridmin[1]
      hi = gridmax[1]
    }

    plotX = rangen(lo, hi, gridpoints[1])
    plotY = J(rows(plotX), 1, .)
    plotY[1]  = p_lo; plotY[rows(plotX)]  = p_hi
    if (confpeak < lo) { // insert original point estimate into grid
      if (gridmin[1] == .) {
        plotX = confpeak            \ plotX
        plotY = (twotailed? 1 : .5) \ plotY
        c = 1
      }
    } else if (confpeak > hi) {
      if (gridmax[1] == .) {
        plotX = plotX \ confpeak
        plotY = plotY \ (twotailed? 1 : .5)
        c = gridpoints[1] + 1
      }
    } else {
      c = floor((confpeak - lo)/(hi - lo)*(gridpoints[1] - 1)) + 2
      plotX = plotX[|.\c-1|] \ confpeak            \ plotX[|c\.|]
      plotY = plotY[|.\c-1|] \ (twotailed? 1 : .5) \ plotY[|c\.|]
    }
  }  // end 1D plot

  i = 1
  if (_quietly==0 & WREnonAR) {
    printf("{txt}")
    do {  // loop so that 1st 2 points are extremes, for efficient interpolation in case when these are first calls--gridmin() and gridmax() manually set
      if (plotY[i] == .) plotY[i] = r0_to_p(plotX[i,]')
      printf(".")
      if (mod(i-rows(plotX)-2, 50)) displayflush()
        else printf("\n")
    } while (1 < (i = mod(i-2,rows(plotX))+1))
    printf("\n")
  } else
    do {
      if (plotY[i] == .) plotY[i] = r0_to_p(plotX[i,]')
    } while (1 < (i = mod(i-2,rows(plotX))+1))

/* determine crossovers analytically, by solving quartic equations 
 appealing in principle, but doesn't seem faster in practice and runs into complications such as handling cubic subcase of quartic, and judging whether a complex root
 is actually real root with imprecision, and requiring doKK==0 approach

 if (dH0==1 & level<100)
    if (interpolating & doKK==0 & interpolate_denom) {
real matrix roots, t0; real scalar numer0sq, numer0tw
      numer0sq = numer0[1] * numer0[1]; numer0tw = numer0[1] + numer0[1]
"ddenomdr2.M.M[404]"
ddenomdr2.M.M[404]
      roots = QuarticRoots((dnumerdr.M[1]*dnumerdr.M[1]) * ddenomdr2.M.M,
                           (dnumerdr.M[1]*dnumerdr.M[1]) * ddenomdr.M.M + (dnumerdr.M[1] * numer0tw) * ddenomdr2.M.M,
                           ddenomdr2.M.M * numer0sq + (dnumerdr.M[1] * numer0tw) * ddenomdr.M.M + (dnumerdr.M[1]*dnumerdr.M[1]) * denom0.M - dnumerdr.M :* dnumerdr.M * denom0.M[1],
                           ddenomdr.M.M  * numer0sq + (dnumerdr.M[1] * numer0tw) * denom0.M - dnumerdr.M :* numer0 * (denom0.M[1]+denom0.M[1]),
                           numer0sq * denom0.M - numer0 :* numer0 * denom0.M[1])
      t = (numer0 :+ dnumerdr.M :* roots) :> 0; t0 = (numer0[1] :+ dnumerdr.M[1] :* roots) :> 0
      roots = roots :/ (t :& t0 :| (t :| t0):==0)  // annihilate quartic solutions corresponding to opposite-signed t stats

      if (allof(colnonmissing(roots)[|2\.|], 1)) { // simple case of one crossing for each replication
        roots = colmin(roots)[|2\.|]'
/*external matrix _roots
_roots=roots*/
        if (ptype) {
          _sort(roots, 1)
          if (ptype==1)  // equal-tail p value
            CI = roots[round(alpha*.5*reps)], roots[round((1-alpha*.5)*reps)]
          else if (ptype==2) // lower-tailed p value
            CI = ., roots[round((1-alpha)*reps)]
          else  // lower-tailed p value
            CI = roots[round(alpha*reps)], .
        } else {  // symmetric
          roots = sort(abs(roots :- confpeak), -1)
          t = roots[round(alpha*reps)]
"alpha,confpeak, mean(roots),t"
alpha,confpeak, mean(roots),t
          CI = confpeak - t, confpeak + t
        }
"CI"
CI
      }

"selectindex(colnonmissing(roots):!=1)\roots[,selectindex(colnonmissing(roots):!=1)]"
selectindex(colnonmissing(roots):!=1)\roots[,selectindex(colnonmissing(roots):!=1)]
external real matrix _numer0,_dnumerdr,_denom0,_ddenomdr,_ddenomdr2
_numer0=numer0; _dnumerdr=dnumerdr.M;_denom0=denom0.M;_ddenomdr=ddenomdr.M.M; _ddenomdr2=ddenomdr2.M.M
    } else*/ {
      CI = (plotY :> alpha) :/ (plotY :< .); CI = CI[|2\.|] - CI[|.\rows(plotX)-1|]
      lo = boottest_selectindex(CI:== 1)
      hi = boottest_selectindex(CI:==-1)
      if (rows(lo)==0 & rows(hi)==0)
        CI = . , .
      else {
             if (rows(lo)==0) lo = .
        else if (rows(hi)==0) hi = .
        else {
          if ( lo[1       ] >  hi[1       ]) lo = .  \ lo // non-rejection ranges that are not within grid range
          if (-lo[rows(lo)] < -hi[rows(hi)]) hi = hi \ .
        }
        CI = lo, hi
        for (i=rows(lo); i; i--)
          for (j=2; j; j--)
            if (CI[i,j]<.)
              CI[i,j] = search(alpha, plotY[CI[i,j]], plotX[CI[i,j]], plotY[CI[i,j]+1], plotX[CI[i,j]+1])
      }
    }

    if (c < .) {  // now that it's done helping graph look good, remove peak point from returned grid for evenness, for Bayesian sampling purposes
      peak = plotX[c], plotY[c]
      if (c==1) {
        plotX = plotX[|2\.|]; plotY = plotY[|2\.|]
      } else if (c==gridpoints[1]+1) {
        plotX = plotX[|.\gridpoints[1]|]; plotY = plotY[|.\gridpoints[1]|]
      } else {
        plotX = plotX[|.\c-1|] \ plotX[|c+1\.|]; plotY = plotY[|.\c-1|] \ plotY[|c+1\.|]
      }
    }

	setquietly(_quietly); pr0 = _pr0; dirty = 1  // restore backups
	notplotted = 0
}

// Return real roots of quartic equation using Euler's method, mathforum.org/dr.math/faq/faq.cubic.equations.html
// . = complex
// modifies arguments!
/*real matrix boottestModel::QuarticRoots(real rowvector _a, real rowvector a, real rowvector b, real rowvector c, real rowvector d) {
  real rowvector e, f, fourg, h, i, j, a2, h2, h3; complex rowvector alpha, z1, z2, p, q, r; complex matrix x

  a = a :/ _a; b = b :/ _a; c = c :/ _a; d = d :/ _a  // need to handle _a=0?

	    e = b - .375*(a2=a:*a)  // substitute with y = x - b/4
	    f = c + (.125*a2 - .5*b):*a
	fourg = 4*d - a2:*(.046875*a2 - .25*b) - a:*c

	h = .5*e  // auxilliary cubic equation
	i = .0625*(e:*e-fourg)
	j = -.015625*f:*f

	p = (i-(h2=h:*h)/3) / 3  // substite with z = y - h/3
	q = (((2/27)*h2-i/3):*h+j) * .5

	alpha = (sqrt(C(q:*q + p:*p:*p)) - q) :^ (1/3)
	z1 = alpha - p :/ alpha  // roots of p,q equation
	alpha = alpha * (-.5+1.bb67ae8584cabX-001i /*exp(2i*pi()/3)*/)
	z2 = alpha - p :/ alpha

	h3 = h/3
	p = sqrt(z1 - h3)  // square roots of roots of auxilliary cubic
	q = sqrt(z2 - h3)
	r = f:/(p:*q) * -.125

	x = p+q+r\p-q-r\-p+q-r\-p-q+r
	return (Re(x) :/ (abs(Im(x)):<1e-5) :- 0.25*a)
}*/


/*real matrix boottestCrossK(real matrix A, real matrix B) return(cross(A,B))
real matrix boottestCrossJ(real matrix A, real matrix B) return(colsum(A:*B))
real matrix boottestDoubleK(real matrix A) return(A+A')
real matrix boottestDoubleJ(real matrix A) return(A+A)*/


// fold a square matrix in half diagonally such that v'Xv is unchanged, but will be computed faster because folded X is half 0's
void boottestModel::_fold(real matrix X)
  X = lowertriangle(X) + uppertriangle(X,0)'

// return matrix whose rows are all the subsets of a row of numbers. Nil is at bottom.
real matrix boottestModel::combs(real scalar d) {
	real matrix retval; real scalar i
	retval = J(2^d, 0, 0)
	for (i=d;i;i--)
		retval = retval, J(2^(d-i),1,1) # (1\0) # J(2^(i-1),1,1) 
	return (retval)
}

// Like Mata's order() but does a stable sort
real colvector boottestModel::stableorder(real matrix X, real rowvector idx)
	return (order((X, (1::rows(X))), (idx,cols(X)+1)))
	
// do lusolve() for precision, but if that fails, use generalized inverse from invsym() NOT
real matrix boottest_lusolve(real matrix A, real matrix B) {
	real matrix retval
return (invsym(A) * B)  // NOT--this is faster and in general precise enough for bootstrap
	if (hasmissing(retval = lusolve(A, B)))
		return (invsym(A) * B)
	return(retval)
}

// Stata interface
void boottest_stata(string scalar statname, string scalar dfname, string scalar dfrname, string scalar pname, string scalar padjname, string scalar ciname, 
	string scalar plotname, string scalar peakname, real scalar level, real scalar ptol, real scalar ML, real scalar LIML, real scalar Fuller, 
	real scalar K, real scalar AR, real scalar null, real scalar scoreBS, string scalar weighttype, string scalar ptype, string scalar statistic, string scalar madjtype, real scalar NumH0s,
	string scalar XExnames, string scalar XEndnames, real scalar hascons, string scalar Ynames, string scalar bname, string scalar Vname, string scalar Wname, 
	string scalar ZExclnames, string scalar samplename, string scalar scnames, real scalar robust, string scalar IDnames, real scalar NBootClustVar, real scalar NErrClust, 
	string scalar FEname, real scalar NFE, string scalar wtname, string scalar wttype, string scalar Cname, string scalar C0name, real scalar reps, string scalar repsname, string scalar repsFeasname, 
	real scalar small, string scalar diststat, string scalar distname, string scalar gridmin, string scalar gridmax, string scalar gridpoints, real scalar MaxMatSize, real scalar quietly,
	string scalar b0name, string scalar V0name, string scalar vname) {

	real matrix C, R, C0, R0, ZExcl, ID, FEID, sc, XEnd, XEx
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

	M._st_view(sc, ., scnames, samplename)
	M._st_view(Y , ., Ynames , samplename)
	M._st_view(ZExcl, ., ZExclnames , samplename)
	if (FEname != "" ) FEID = st_data(., FEname , samplename)
	if (IDnames != "") ID   = st_data(., IDnames, samplename)
	if (wtname  != "") wt   = st_data(., wtname , samplename) // panelsum() doesn't like views as weights

	M.setMaxMatSize(MaxMatSize)
	M.sethascons(hascons)
	M.setsc(sc)
	M.setML(ML)
	M.setY (Y)
	M.setZExcl(ZExcl)
	M.setwt (wt)
	M.setID(ID, NBootClustVar, NErrClust)
	M.setFEID(FEID, NFE)
	M.setR (R, r)
	M.setR0(R0, r0)
	M.setnull(null)
	M.setsmall(small)
	M.setrobust(robust)
	M.setscoreBS(scoreBS)
	M.setweighttype(weighttype)
	M.setptype(ptype)
	M.setstattype(statistic)
	M.setwttype(wttype)
	M.setreps(reps)
	M.setLIML(LIML)
	M.setFuller(Fuller)
	M.setk(K)
	M.setAR(AR)
	M.setgrid(strtoreal(tokens(gridmin)), strtoreal(tokens(gridmax)), strtoreal(tokens(gridpoints)))
	M.setmadjust(madjtype, NumH0s)
	M.setlevel(level)
	M.setptol(ptol)
	M.setquietly(quietly)

	M._st_view(XEnd, ., XEndnames, samplename)
	M.setXEnd(XEnd)
	M._st_view(XEx, ., XExnames, samplename)
	M.setXEx(XEx)
	if (bname != "") M.setbeta(st_matrix(bname)')
	if (Vname != "") M.setV   (st_matrix(Vname) )
	if (Wname != "") M.setW   (st_matrix(Wname) )
	M.setwillplot(plotname != "" | ciname != "")  // indicates whether to optimize for many p evaluations
	st_numscalar(statname, M.getstat())
	st_numscalar(pname   , M.getp   ())
	st_numscalar(repsname, M.getreps())
	st_numscalar(repsFeasname, M.getrepsFeas())
	st_numscalar(padjname, M.getpadj())
	st_numscalar(dfname  , M.getdf  ())
	st_numscalar(dfrname , M.getdf_r())
	st_matrix(b0name, M.getb()')
	st_matrix(V0name, M.getV())
	if (plotname != "" | (level<100 & ciname != "")) {
		if (plotname != "") st_matrix(plotname, M.getplot())
		if (cols(M.peak)) st_matrix(peakname, M.getpeak())
		if (level<100 & ciname != "") st_matrix(ciname, M.getCI())
	}
	if (distname != "") st_matrix(distname, M.getdist(diststat))
	if (vname != "" & reps) st_matrix(vname, M.getv())
	M.M_DGP.parent = NULL // break loop in data structure topology to avoid Mata garbage-cleaning leak
}

mata mlib create lboottest, dir("`c(sysdir_plus)'l") replace
mata mlib add lboottest *(), dir("`c(sysdir_plus)'l")
mata mlib index
end
