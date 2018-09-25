*! boottest 2.2.2 25 September 2018
*! Copyright (C) 2015-18 David Roodman

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
string scalar      boottestVersion() return("02.02.02")

struct smatrix {
	real matrix M
}

struct boottest_clust {
	real scalar N, multiplier
	real colvector order
	real matrix info
}

struct structFE {
	pointer (struct structFE scalar) scalar next
	real colvector is, wt
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
	struct smatrix matrix CT_ZVR0
	struct smatrix colvector WZVR0

	void new(), InitExog(), InitEndog(), InitTestDenoms(), SetDGP(), SetS(), InitEstimate(), Estimate(), setParent(), SetLIMLFullerK(), SetAR()
	pointer(real matrix) scalar partialFE()
}

class boottestModel {
	real scalar scoreBS, reps, small, weighttype, null, dirty, initialized, Neq, ML, Nobs, _Nobs, k, kEx, el, sumwt, NClustVar, robust, weights, REst, multiplier, quietly, FEboot, NErrClustCombs, ///
		sqrt, hascons, LIML, Fuller, K, IV, WRE, WREnonAR, ptype, twotailed, df, df_r, AR, D, cuepoint, willplot, plotted, NumH0s, p, NBootClustVar, NErrClust, ///
		NFE, doQQ, granular, purerobust, subcluster, NBootClust, repsFeas, u_sd, level, MaxMatSize, NWeightGrps, enumerate
	real matrix QQ, eZVR0, SewtXV, VR0, betadev, numer, u, U, S, SAR, SAll, LAll_invRAllLAll, plot, CI, CT_WE, infoBootAll, infoErrAll, infoAllData, J_ClustN_NBootClust, denom0
	pointer (real matrix) scalar pZExcl, pR, pR0, pID, pFEID, pXEnd, pXEx, pG, pX, pinfoBootData, pinfoErrData
	pointer (real colvector) scalar pr, pr0, pY, pSc, pwt, pW, pV
	string scalar wttype, madjtype, seed
	real colvector Dist, DistCDR, s, sAR, plotX, plotY, sAll, beta, wtFE, ClustShare, IDBootData, IDBootAll, WeightGrpStart, WeightGrpStop, gridstart, gridstop, gridpoints
	real rowvector peak
	struct boottest_clust colvector Clust
	class AnalyticalModel scalar M_DGP
	pointer (class AnalyticalModel scalar) scalar pM_Repl, pM
	struct smatrix matrix denom
	struct smatrix colvector Kcd, XEndstar, XExXEndstar, ZExclXEndstar, XZi, eZi, euZVR0
	struct structFE rowvector FEs
	pointer (real matrix) matrix pQ
	pointer (struct boottest_clust scalar) scalar pBootClust

	void new(), set_sqrt(), boottest(), make_DistCDR(), plot(), contourplot(), setXEx(), setptype(), setdirty(), setXEnd(), setY(), setZExcl(), setwt(), setsc(), setML(), setLIML(), setAR(), 
		setFuller(), setk(), setquietly(), setbeta(), setV(), setW(), setsmall(), sethascons(), setscoreBS(), setreps(), setnull(), setWald(), setRao(), setwttype(), setID(), setFEID(), setlevel(), 
		setrobust(), setR(), setR0(), setwillplot(), setgrid(), setmadjust(), setweighttype(), MakeWildWeights(), MakeNonWREStats(), MaxMatSize(), _st_view()
	real scalar r0_to_p(), search(), getp(), getpadj(), getstat(), getdf(), getdf_r(), getreps(), getrepsFeas(), MakeNonWRENumers(), MakeWREStats()
	real matrix combs(), count_binary(), crosstab(), getplot(), getCI()
	real rowvector getpeak()
	real colvector getdist(), stableorder()
}
void AnalyticalModel::new()
	AR = 0

void AnalyticalModel::setParent(class boottestModel scalar B) parent = &B

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

	parent->pXEx = partialFE(parent->pXEx)
	pXExXEx = &cross(*parent->pXEx, *parent->pwt, *parent->pXEx)
	if (cols(*parent->pZExcl)) { // GMM, 2SLS, LIML
		parent->pZExcl = partialFE(parent->pZExcl)
		ZExclXEx = cross(*parent->pZExcl, *parent->pwt, *parent->pXEx)
		pZXEx = &(*pXExXEx \ ZExclXEx)
		if (parent->IV)
			pinvZZ = &invsym((ZZ = *pZXEx, (ZExclXEx' \ cross(*parent->pZExcl, *parent->pwt, *parent->pZExcl))))
	} else
		pXX = pXExXEx
}

void AnalyticalModel::SetDGP(class AnalyticalModel scalar _DGP)
	DGP = &_DGP

void AnalyticalModel::SetS(real matrix S) {
	pS = &S // DGP==NULL means this is the DGP
	if (LIML) Splus = blockdiag(S, 1)  // add an entry to S for the dep var
}

// stuff that can be done before S & r0 set, but depend on endogenous variables, which are bootstrapped in WRE
void AnalyticalModel::InitEndog(pointer (real colvector) scalar _pY, pointer (real matrix) scalar _pXEnd, | ///
		pointer (real colvector) scalar _pZExclY, pointer (real rowvector) scalar _pXExY, real scalar _YY, pointer (real matrix) scalar _pZExclXEnd, pointer (real matrix) scalar _pXExXEnd) {

	pY = partialFE(_pY); pXEnd = partialFE(_pXEnd)
	
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
			eigensystemselecti(luinv(TT) * TPZT, rows(TT)\rows(TT), vec, val)
			K = 1/Re(1-val) - Fuller / (parent->_Nobs - parent->el)   // sometimes a tiny imaginary component sneaks in
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
	real matrix AVR0; real scalar d, c; struct smatrix rowvector _CT_ZVR0; pointer (real matrix) scalar pWZVR0

	pV = rows(S)? &(S * invsym(S ' (*pH) * S) * S') : ///
	              &(rows(invH)? invH : invsym(*pH))
	VR0 = *pV * *parent->pR0'

	if (parent->scoreBS | (parent->robust & !(parent->WREnonAR & parent->NClustVar==1 & !parent->NFE))) {
		if (K) {
			AVR0 = A * VR0
			ZVR0 = *parent->pZExcl * AVR0[|parent->kEx+1,.\.,.|]; if (parent->kEx) ZVR0 = ZVR0 + *parent->pXEx * AVR0[|.,.\parent->kEx,.|]
		} else if (AR) {
			ZVR0 = *parent->pZExcl * VR0[|parent->kEx+1,.\.,.|]
			if (cols(*parent->pXEx))
				ZVR0 = ZVR0 + *parent->pXEx * VR0[|.,.\parent->kEx,.|]
		} else
			ZVR0 = *parent->pXEx * VR0

			if (parent->granular) {
				pWZVR0 = &(parent->weights? ZVR0 :* *parent->pwt : ZVR0)
				WZVR0 = smatrix(parent->df)
				for (d=parent->df;d;d--)
					WZVR0[d].M = (*pWZVR0)[,d]
			}

			if (parent->NFE & parent->robust & !parent->WREnonAR & !parent->FEboot & !parent->scoreBS & parent->granular<parent->NErrClustCombs) {
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
	real matrix invZMeZ; real colvector negZeinvee

	beta = rows(s)? beta0 + dbetads * s : beta0

	if (AR) {
		e = *pY - *parent->pZExcl * beta[|cols(*parent->pXEx)+1\.|]
		if (cols(*parent->pXEx))
			e =  e - *parent->pXEx * beta[|.\cols(*parent->pXEx)|]
	} else if (parent->IV)
		if (parent->kEx == 0)
			e = *pY - *pXEnd * beta[|parent->kEx+1\.|]
		else
			e = *pY - *pXEnd * beta[|parent->kEx+1\.|] - *parent->pXEx * beta[|.\parent->kEx|]
	else
			e = *pY                                    - *parent->pXEx * beta[|.\parent->kEx|]

	if (!(parent->robust | parent->scoreBS) | (DGP==NULL & LIML)) // useful in non-robust, residual-based bootstrap, and in computing e2 in LIML (just below)
		ee = YY - 2 * *pXY ' beta + beta ' (*pXX) * beta
	if (!(parent->robust | parent->scoreBS))
		eec = parent->hascons? ee : ee - (parent->weights? cross(e, *parent->pwt) : sum(e))^2 / parent->_Nobs // sum of squares after centering, N * Var

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

void boottestModel::new() {
	AR = LIML = Fuller = WRE = small = scoreBS = weighttype = Neq = ML = initialized = quietly = sqrt = hascons = IV = ptype = robust = plotted = NFE = FEboot = granular = NErrClustCombs = subcluster = reps = repsFeas = 0
	twotailed = null = dirty = willplot = u_sd = 1
	level = 95
	cuepoint = MaxMatSize = .
	pXEnd = pXEx = pZExcl = pY = pSc = pID = pFEID = pR = pR0 = pwt = &J(0,0,0)
	pr = pr0 = &J(0,1,0)
}

void boottestModel::setdirty(real scalar _dirty, | real scalar noinitialize) {
	dirty = _dirty
	if (_dirty & noinitialize!=1)
		initialized = 0
}
void boottestModel::set_sqrt(real scalar _sqrt) {
	if (_sqrt < sqrt)
		if (!dirty) Dist = Dist :* Dist
	else
		setdirty(1)
	sqrt = _sqrt
}
void boottestModel::setptype(string scalar ptype)     {
	real scalar p
	p = cross( (strtrim(strlower(ptype)) :== ("symmetric"\"equaltail"\"lower"\"upper")), 1::4 ) - 1
	if (p<0) 
		_error(198, `"p-value type must be "symmetric", "equaltail", "lower", or "upper.""')
	this.ptype = p
	this.twotailed = p<=1
}
void boottestModel::setXEnd    (real matrix X       ) {
	this.pXEnd  = &X; setdirty(1)
}
void boottestModel::setXEx    (real matrix X        ) {
	this.pXEx  = &X; setdirty(1)
}
void boottestModel::setY       (real matrix Y       ) {
	this.pY  = &Y; setdirty(1)
}
void boottestModel::setZExcl   (real matrix Z       ) {
	this.pZExcl  = &Z; setdirty(1)
}
void boottestModel::setwt       (real matrix wt     ) {
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
	this.reps = reps; setdirty(1)
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
void boottestModel::setrobust  (real scalar robust  ) {
	this.robust = robust
	if (robust==0) setID(J(0,0,0), 1, 1)
	setdirty(1)
}
void boottestModel::setR (real matrix R, real colvector r ) {
	this.pR = &R; 	this.pr  = &r; setdirty(1)
}
void boottestModel::setR0(real matrix R0, real colvector r0) {
	this.pR0 = &R0; this.pr0 = &r0; setdirty(1)
}
void boottestModel::setwillplot(real scalar willplot) {
	this.willplot = willplot
}
void boottestModel::setgrid(real rowvector gridstart, real rowvector gridstop, real rowvector gridpoints) {
	this.gridstart = gridstart; this.gridstop = gridstop; this.gridpoints = gridpoints
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
void boottestModel::MaxMatSize (real scalar MaxMatSize) {
	this.MaxMatSize = MaxMatSize; setdirty(1)
}

real colvector boottestModel::getdist(| string scalar diststat) {
	if (dirty) boottest()
	make_DistCDR(diststat)
	return(DistCDR)
}
void boottestModel::make_DistCDR(| string scalar diststat) {
	pointer (real rowvector) scalar pnumer
	if (diststat == "numer") {
		pnumer = u_sd==1? &numer : &(numer / u_sd)
		_sort( DistCDR = (*pnumer)[|2\.|]' :+ ((*pnumer)[1] + *pr0) , 1)
	} else if (!rows(DistCDR))
		if (rows(Dist)>1)
			_sort( DistCDR=Dist[|2\.|] , 1)
		else
			DistCDR = J(0,1,0)
}

// Ties count half. Robust to missing bootstrapped values being interpreted as +infinity.
real scalar boottestModel::getp(|real scalar analytical) {
	real scalar t
	if (dirty) boottest()
	t = Dist[1]
	if (t == .) return (.)
	if (reps & analytical==.) {
		repsFeas = colnonmissing(Dist) - 1
		if (sqrt & ptype != 3) {
			if (ptype==0)
				p = colsum(-abs(t):>-abs(Dist)) / repsFeas  // symmetric p value; do so as not to count missing entries in Dist
			else if (ptype==1) // equal-tail p value
				p = 2*min((colsum(t:>Dist) , colsum(-t:>-Dist))) / repsFeas
			else
				p = colsum( t:> Dist) / repsFeas  // lower-tailed p value
		} else
				p = colsum(-t:>-Dist) / repsFeas  // upper-tailed p value or p value based on squared stats
	} else {
		p = small? Ftail(df, df_r, sqrt? t*t : t) : chi2tail(df, sqrt? t*t : t)
		if (sqrt & !twotailed) {
			p = p / 2
			if ((ptype==3) == (t<0))
				p = 1 - p
		}
	}
	return(p)
}
// Return number of bootstrap replications with feasible results
// Returns 0 if getp() not yet accessed, or doing non-bootstrapping tests
real scalar boottestModel::getrepsFeas()
	return (repsFeas)
// return number of replications, possibly reduced to 2^G
real scalar boottestModel::getreps()
	return (reps)

real scalar boottestModel::getpadj(|real scalar analytical) {
	(void) getp(analytical)
	if (madjtype=="bonferroni") return(min((1, NumH0s*p)))
	if (madjtype=="sidak"     ) return(1 - (1 - p)^NumH0s)
	return(p)
}

real scalar boottestModel::getstat() {
	if (dirty) boottest()
	return(Dist[1])
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
	if (!plotted)
		if (rows(*pr0)==1)
			plot()
		else if (rows(*pr0)==2)
			contourplot()
	return((plotX,plotY))
}
real rowvector boottestModel::getpeak() {
	if (!plotted) plot()
	return(peak)
}
real matrix boottestModel::getCI() {
	if (!plotted) plot()
	return(CI)
}

void boottestModel::_st_view(real matrix V, real scalar i, string rowvector j, string scalar selectvar) {
	if (favorspeed() | 1)
		V = length(tokens(j))? st_data(i, j, selectvar) : st_data(i, J(1,0,0), selectvar)
	else
		st_view(V, i, j, selectvar)
}















































































// main routine
void boottestModel::boottest() {
	real colvector rAll, sortID, o, _FEID
	real rowvector val, ClustCols
	real matrix RAll, L, LAll, vec, Combs, t, IDErr
	real scalar i, j, c, minN, sumN, _reps, g, i_FE
	pointer (real matrix) scalar _pR0, pIDAll
	class AnalyticalModel scalar M_WRE
	pragma unset vec; pragma unset val; pragma unused M_WRE

	if (!initialized) {  // for efficiency when varying r0 repeatedly to make CI, do stuff once that doesn't depend on r0
		Nobs = rows(*pXEx)
		kEx = cols(*pXEx)
		if (!cols(*pZExcl)) pZExcl = &J(Nobs,0,0)
		if (!cols(*pXEnd)) pXEnd = &J(Nobs,0,0)
		D = cols(*pXEnd) + 1
		k  = cols(*pR0)
		REst = rows(*pR) // base model contains restrictions?
		if (pZExcl != NULL) el = cols(*pZExcl) + kEx
		if (K==.) K = cols(*pZExcl)>0
		IV = K & 1 // & pW==NULL
		WRE = (IV & !scoreBS) | AR
		WREnonAR = WRE & !AR

		if (weights = rows(*pwt)>1)
			sumwt = sum(*pwt)
		else
			pwt = &(sumwt = 1)
		_Nobs = weights & wttype=="fweight"? sumwt : Nobs

		IDBootData = J(0,1,0)
		if (NClustVar = cols(*pID)) {
			minN = .; sumN = 0

			Combs = combs(NErrClust) // represent all error clustering combinations. First is intersection of all error clustering vars
			Clust = boottest_clust(rows(Combs)-1) // leave out no-cluster combination
			NErrClustCombs = length(Clust)
			subcluster = NClustVar - NErrClust

			 infoAllData = _panelsetup(*pID, 1..NClustVar) // info for grouping by intersections of all bootstrap & clustering vars wrt data; used to speed crosstab EZVR0 wrt bootstrapping cluster & intersection of all error clusterings
			pinfoErrData = NClustVar > NErrClust ? &_panelsetup(*pID, subcluster+1..NClustVar) : &infoAllData // info for intersections of error clustering wrt data
			 IDErr = rows(*pinfoErrData)==Nobs? *pID :   (*pID)[(*pinfoErrData)[,1],] // version of ID matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
			pIDAll = rows(  infoAllData)==Nobs?  pID : &((*pID)[   infoAllData [,1],]) // version of ID matrix with one row for each all-bootstrap & error cluster-var intersection instead of 1 row for each obs

			if (subcluster) { // for subcluster bootstrap, bootstrapping cluster is not among error clustering combinations
				pBootClust = &(boottest_clust())
				pBootClust->info = _panelsetup(*pIDAll, 1..NBootClustVar) // bootstrapping cluster info w.r.t. all-bootstrap & error-cluster intersections
				NBootClust = rows(pBootClust->info)
			} else
				pBootClust = &(Clust[2^(NClustVar - NBootClustVar)]) // location of bootstrap clustering within list of cluster combinations

			for (c=1; c<=NErrClustCombs; c++) { // for each error clustering combination
				ClustCols             = subcluster :+ boottest_selectindex(Combs[c,])
				Clust[c].multiplier   = 2 * mod(cols(ClustCols),2) - 1

				if (c == 1)
				  if (subcluster) {
						IDErr = IDErr[ Clust.order = stableorder(IDErr, ClustCols), ]
						Clust.info       = _panelsetup(IDErr, ClustCols)
					} else
						Clust.info       = J(rows(infoAllData),0,0)  // causes no collapsing of data in _panelsum() calls
				else {
					if (any( Combs[|c, min(boottest_selectindex(Combs[c,] :!= Combs[c-1,])) \ c,.|])) // if this sort ordering same as last to some point and missing thereafter, no need to re-sort
						IDErr = IDErr[ Clust[c].order = stableorder(IDErr, ClustCols), ]

					Clust[c].info       = _panelsetup(IDErr, ClustCols)
				}

				Clust[c].N            = rows(Clust[c].info)
				sumN = sumN + Clust[c].N

				if (small) {
				  Clust[c].multiplier = Clust[c].multiplier * Clust[c].N/(Clust[c].N-1)
					if (minN > Clust[c].N) minN = Clust[c].N
				}
			}

			if (scoreBS | !(WREnonAR | (1 | hascons)))
				ClustShare = weights? _panelsum(*pwt, *pinfoErrData)/sumwt : ((*pinfoErrData)[,2]-(*pinfoErrData)[,1]:+ 1)/Nobs // share of observations by group 

		} else { // if no clustering, cast "robust" as Clustering by observation
			pBootClust = &(Clust = boottest_clust())
			Clust.multiplier = small? _Nobs / (_Nobs - 1) : 1
			sumN = Clust.N = Nobs
			Clust.info = J(Nobs, 0, 0) // signals _panelsum not to aggregate
			NErrClustCombs = 1
			if (scoreBS | !(WREnonAR | (1 |hascons)))
				ClustShare = weights? *pwt/sumwt : 1/_Nobs
		}

		if (WREnonAR)
			if (NClustVar)
				pinfoBootData = &_panelsetup(*pID, 1..NBootClustVar, IDBootData)
			else
				pinfoErrData = pinfoBootData = &J(Nobs,0,0)
		else if (NClustVar) {
			if (NClustVar > NBootClustVar) // bootstrap Cluster grouping defs rel to original data
				pinfoBootData = &_panelsetup(*pID, 1..NBootClustVar)
			else
				pinfoBootData = &infoAllData
		} else
			pinfoBootData = &(infoAllData = J(Nobs,0,0)) // causes no collapsing of data in _panelsum() calls, only multiplying by weights if any
		NBootClust  = rows(*pinfoBootData)

		purerobust = NClustVar & !scoreBS & NBootClust==Nobs & !subcluster // do we ever error-cluster *and* bootstrap-cluster by individual?
		granular   = NClustVar & !scoreBS & (purerobust | 5*Nobs*k+1/25*2*NBootClust*k^2+1/25*2*NBootClust^2*k+NBootClust+1/25*2*NBootClust^2*reps+2*NBootClust*reps > 3*Nobs*reps+Nobs*k+1/25*(2*NBootClust*k*reps+2*k*reps-2*k*NBootClust-2*NBootClust*reps+2*NBootClust*k*reps)+2*NBootClust*reps)
		doQQ = !granular & NBootClust * (sumN + NBootClust*reps +.5*(NErrClustCombs)) < 2*reps*(2*sumN + NErrClustCombs) // estimate compute time for U :* (sum_c QQ) * U vs. sum_c(colsum((Q_c*U):*(Q_c*U)))

		if (robust & !purerobust) {
			if (subcluster | granular)
				infoErrAll = _panelsetup(*pIDAll, subcluster+1..NClustVar) // info for error clusters wrt data collapsed to intersections of all bootstrapping & error clusters; used to speed crosstab EZVR0 wrt bootstrapping cluster & intersection of all error clusterings
			if (scoreBS)
				J_ClustN_NBootClust = J(Clust.N, NBootClust, 0)
		}

		if (NFE) {
			sortID = (*pFEID)[o = stableorder(*pFEID, 1)]
			i_FE = 1; FEboot = reps>0; j = Nobs; _FEID = wtFE = J(Nobs, 1, 1)
			FEs = structFE(NFE)
			for (i=Nobs-1;i;i--) {
				if (sortID[i] != sortID[i+1]) {
					FEs[i_FE].is = o[|i+1\j|]
					if (weights) {
						t  = (*pwt)[FEs[i_FE].is]
						wtFE[FEs[i_FE].is] = FEs[i_FE].wt = t / colsum(t)
					} else
						wtFE[FEs[i_FE].is] = FEs[i_FE].wt = J(j-i, 1, 1/(j-i))
					j = i
					
					if (FEboot & NClustVar) { // are all of this FE's obs in same bootstrapping cluster?
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
				wtFE[FEs[NFE].is] = FEs[NFE].wt = t / colsum(t)
			} else
				wtFE[FEs[NFE].is] = FEs[NFE].wt = J(j-i,1,1/(j-i))
			if (FEboot & NClustVar) { // are all of this FE's obs in same bootstrapping cluster?
				t = (*pID)[FEs[NFE].is, 1..NBootClustVar]
				FEboot = all(t :== t[1,])
			}

			pFEID = &_FEID // ordinal fixed effect ID

			if (robust & !(scoreBS | FEboot) & granular < NErrClustCombs)
				infoBootAll = _panelsetup(*pIDAll, 1..NBootClustVar) // info for bootstrapping clusters wrt data collapsed to intersections of all bootstrapping & error clusters
		}

		if (granular & !(WREnonAR & purerobust))
			if (NFE)
				(void) _panelsetup(*pID   , 1..NBootClustVar, IDBootData)
			else
				(void) _panelsetup(*pIDAll, 1..NBootClustVar, IDBootAll )

		if (ML)
			df = rows(*pR0)
		else {
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
				if (willplot) { // for plotting purposes get original point estimate if not normally generated
					M_DGP.SetS(S) // no-null model in DGP
					M_DGP.InitEstimate()
					M_DGP.Estimate(s)
					cuepoint = *pR0 * M_DGP.beta - *pr0 // not true CUE estimate unless classical errors, but serves same purpose as weakiv cuepoint option
				}
				pR0 = &(J(cols(*pZExcl),kEx,0), I(cols(*pZExcl))) // for AR test, picks out coefs on excluded exogenous variables
			}
			df = rows(*pR0)
			
			if (granular) {
				euZVR0 = smatrix(df)
				pX = AR? &(*pXEx, *pZExcl) : pXEx
			}

			M_DGP.SetS(SAll) // (potentially) constrained model in DGP; SAll imposes constraints, pR0 tests hypotheses on results
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

		if (small) df_r = NClustVar? minN - 1 : _Nobs - k - NFE

		if (df==1) set_sqrt(1) // work with t/z stats instead of F/chi2

		if (small)
			multiplier = (_Nobs - k - NFE) / (_Nobs - robust) / df // divide by # of constraints because F stat is so defined
		else
			multiplier = 1
		if (!(robust | ML))
			multiplier = multiplier * _Nobs // will turn sum of squared errors in denom of t/z into mean
		if (sqrt) multiplier = sqrt(multiplier)
			
		if (!null) M_DGP.beta = J(0,1,0) // in case model re-dirtied and we're not imposing null, Estimate() will know to recompute beta for first r0 value tried, then stop

		denom = smatrix(df,df)
		if (WREnonAR) {
			XEndstar = XExXEndstar = ZExclXEndstar = smatrix(D-1)
			if (NClustVar)
				XZi = eZi = smatrix(NBootClust)
		} else if (robust) {
			Kcd = smatrix(df)
			pQ = J(NErrClustCombs, df, NULL)
		}

		enumerate = 0
		if (weighttype<=1 & reps)
			if (NBootClust*ln(2) < ln(reps)+1e-6)
				if (weighttype==1) {
					if (!quietly)
						printf("\nWarning: with %g Clusters, the number of replications, %g, exceeds the universe of Mammen draws, 2^%g = %g. \nConsider Webb weights instead, using {cmd:weight(webb)}.\n", NBootClust, reps, NBootClust, 2^NBootClust) 
				} else {
					if (!quietly)
						printf("\nWarning: with %g Clusters, the number of replications, %g, exceeds the universe of Rademacher draws, 2^%g = %g. Sampling each once. \nConsider Webb weights instead, using {cmd:weight(webb)}.\n", NBootClust, reps, NBootClust, 2^NBootClust)
					enumerate = 1
					MaxMatSize = .
				}
			else if (!quietly & (NBootClust>12 | weighttype) & floor(level/100 * (reps+1)) != level/100 * (reps+1)) {
					printf("\nNote: The bootstrap usually performs best when the confidence level (here, %g%%)\n", level)
					printf("      times the number of replications plus 1 (%g+1=%g) is an integer.\n", reps, reps+1)
				}

		NWeightGrps = MaxMatSize == .? 1 : ceil((reps+1) * max((rows(IDBootData), rows(IDBootAll), NBootClust)) * 8 / MaxMatSize / 1.0X+1E) // 1.0X+1E = giga(byte)
		if (NWeightGrps == 1) {
			MakeWildWeights(reps, 1) // make all wild weights, once
			if (enumerate) reps = cols(u) - 1 // replications reduced to 2^G
			NWeightGrps    = 1
			WeightGrpStart = 1
			WeightGrpStop  = reps+1
		} else {
			seed = rseed()
			_reps = ceil((reps+1) / NWeightGrps)
			 WeightGrpStart = (0::NWeightGrps-1) * _reps :+ 1
			(WeightGrpStop  = (1::NWeightGrps  ) * _reps     )[NWeightGrps] = reps+1
		}

		if (WREnonAR | df>1 | MaxMatSize<.) // unless nonWRE or df=1 or splitting weight matrix, code will create Dist element-by-element, so pre-allocate vector now
			Dist = J(reps+1, 1, .)
	} // done with one-time stuff--not dependent on r0--if constructing CI or plotting confidence curve

	if (!ML) // GMM, 2SLS, analytical LIML
		if (AR) {
			pM_Repl->InitEndog(&(*pY - *pXEnd * *pr0), NULL, &(*M_DGP.pZExclY - *M_DGP.pZExclXEnd * *pr0), &(*M_DGP.pXExY - *M_DGP.pXExXEnd * *pr0))
			pM_Repl->InitEstimate()
			pM_Repl->Estimate(sAR)
			pM_Repl->InitTestDenoms(SAR)
		} else if (null | !rows(M_DGP.beta)) { // don't need to recompute if we're not actually imposing the null
			rAll = null? *pr0 : J(0, 1, 0); if (REst) rAll =  *pr \ rAll // constant terms of model + null constraints
			sAll = rows(rAll) ? LAll_invRAllLAll * rAll : J(0,1,0)
			M_DGP.Estimate(sAll)
		}

	if (NWeightGrps > 1)
		rseed(seed)

	for (g=1; g<=NWeightGrps; g++) { // do group 1 first because it includes col 1, which is all that might need updating in constructing CI in WCU
		if (NWeightGrps > 1)
			MakeWildWeights(WeightGrpStop[g] - WeightGrpStart[g] + (g>1), g==1)

		if (WREnonAR) {
			if (MakeWREStats(WeightGrpStart[g], WeightGrpStop[g]))
				return
		} else {
			if (MakeNonWRENumers(WeightGrpStart[g], WeightGrpStop[g]))
				return
			MakeNonWREStats(WeightGrpStart[g], WeightGrpStop[g])
		}
	}

	if (multiplier!=1) Dist = Dist * multiplier
	DistCDR = J(0,0,0)
	setdirty(0)
	initialized = 1
}


// draw wild weight matrix of width _reps. If first=1, insert column of 1s at front
// under Rademacher, if warning recommending Webb weights is triggered, then 
void boottestModel::MakeWildWeights(real scalar _reps, real scalar first) {
	if (enumerate)
		u = J(NBootClust,1,1), count_binary(NBootClust, -1-WREnonAR, 1-WREnonAR) // complete Rademacher set
	else if (weighttype==3)
		u = rnormal(NBootClust, _reps+1, -WREnonAR, 1) // normal weights
	else if (weighttype==4)
		u = rgamma(NBootClust, _reps+1, 4, .5) :- (2 + WREnonAR) // Gamma weights
	else if (weighttype==2)
		if (WREnonAR) {
			u = sqrt(2 * ceil(runiform(NBootClust, _reps+first) * 3)) :* ((runiform(NBootClust, _reps+1):>=.5):-.5) :- 1 // Webb weights, minus 1 for convenience in WRE
		} else {
			u = sqrt(    ceil(runiform(NBootClust, _reps+first) * 3)) :* ((runiform(NBootClust, _reps+1):>=.5):-.5)      // Webb weights, divided by sqrt(2)
			u_sd = sqrt(.5)
		}
	else if (weighttype)
		if (WREnonAR)
			u = ( rdiscrete(NBootClust, _reps+first,(.5+sqrt(.05)\.5-sqrt(.05))) :- 1.5 ) *sqrt(5) :- .5 // Mammen weights, minus 1 for convenience in WRE
		else {
			u = ( rdiscrete(NBootClust, _reps+first,(.5+sqrt(.05)\.5-sqrt(.05))) :- 1.5 )          :+ .5/sqrt(5) // Mammen weights, divided by sqrt(5)
			u_sd = sqrt(.2)
		}
	else
		if (WREnonAR) {
			u = runiform(NBootClust, _reps+first) :<  .5; u = (-2) * u // Rademacher weights, minus 1 for convenience in WRE
		} else {
			u = runiform(NBootClust, _reps+first) :>= .5; u = u :- .5  // Rademacher weights, divided by 2
			u_sd = .5
		}
	if (first)
		u[,1] = J(NBootClust, 1, WREnonAR? 0 : u_sd)  // keep original residuals in first entry to compute base model stat

	U = WREnonAR? u[IDBootData,] : J(0,0,0)
}


real scalar boottestModel::MakeWREStats(real scalar thisWeightGrpStart, real scalar thisWeightGrpStop) {
	real scalar c, j, i
	real colvector _e, Ystar, _beta, betaEnd
	real rowvector YstarYstar
	real matrix ZExclYstar, XExYstar, Subscripts, Zi, AVR0, t, XExi
	pointer (real matrix) scalar peZVR0, pu, pXEndstar, pXExXEndstar, pZExclXEndstar
	pragma unused thisWeightGrpStop
	
	if (initialized & !null) {  // if not imposing null and we have returned, then df=1; and distribution doesn't change with r0, only test stat
		numer[1] = *pR0 * pM_Repl->beta - *pr0
		Dist[1] = numer[1] / sqrt(denom.M[1]) * multiplier
		return(1) // no need to continue calc
	}

	if (thisWeightGrpStart == 1) { // first/only weight group? initialize a couple of things
		_e = M_DGP.e + M_DGP.e2 * M_DGP.beta[|kEx+1\.|]
		pu = NClustVar? &U : &u
	}

	Ystar = *M_DGP.pY :+ _e :* *pu
	XExYstar   = cross(*pXEx  , *pwt, Ystar)
	ZExclYstar = cross(*pZExcl, *pwt, Ystar)
	
	if ((LIML | !robust) & !NFE)
		YstarYstar = weights? cross(*pwt, Ystar:*Ystar) : colsum(Ystar:*Ystar)

	if (D==2) {
			XEndstar.M         =  *pXEnd      :+ M_DGP.e2     :* *pu
			XExXEndstar.M      = cross(*pXEx  , *pwt, XEndstar.M)
			ZExclXEndstar.M    = cross(*pZExcl, *pwt, XEndstar.M)
	} else
		for (j=D-1; j; j--) {
			XEndstar[j].M      = (*pXEnd)[,j] :+ M_DGP.e2[,j] :* *pu
			XExXEndstar  [j].M = cross(*pXEx  , *pwt, XEndstar[j].M)
			ZExclXEndstar[j].M = cross(*pZExcl, *pwt, XEndstar[j].M)
		}

	if (NClustVar & !NFE) // prep for optimized computation for bootstrapping cluster when no FE
		for (i=NBootClust; i; i--) {
			Subscripts = (*pinfoBootData)[i,]', (.\.)
			XExi = kEx? (*pXEx)[|Subscripts|] : J((*pinfoBootData)[i,2]-(*pinfoBootData)[i,1]+1,0,0)
			Zi = XExi , (*pZExcl)[|Subscripts|] // inefficient?
			if (weights) Zi = Zi :* (*pwt)[|Subscripts|]
			XZi[i].M = cross(XExi, Zi) \ cross((*pXEnd)[|Subscripts|], Zi) \ cross((*pY)[|Subscripts|], Zi)
			eZi[i].M =                   cross(M_DGP.e2[|Subscripts|], Zi) \ cross(   _e[|Subscripts|], Zi)
		}

	for (j=cols(u); j; j--) { // WRE bootstrap
		pXEndstar      = &(   XEndstar.M  [,j])
		pXExXEndstar   = &(XExXEndstar.M  [,j])
		pZExclXEndstar = &(ZExclXEndstar.M[,j])
		for (i=2; i<D; i++) {
			pXEndstar      = &(*pXEndstar     , XEndstar     [i].M[,j])
			pXExXEndstar   = &(*pXExXEndstar  , XExXEndstar  [i].M[,j])
			pZExclXEndstar = &(*pZExclXEndstar, ZExclXEndstar[i].M[,j])
		}

		pM_Repl->InitEndog(&(Ystar[,j]), pXEndstar, &(ZExclYstar[,j]), &(XExYstar[,j]), (cols(YstarYstar)? YstarYstar[j] : .), pZExclXEndstar, pXExXEndstar)
		pM_Repl->InitEstimate()
		pM_Repl->InitTestDenoms(S) // prepare for replication regressions, null not imposed
		pM_Repl->Estimate(s)
		numer = null | j==1? *pR0 * pM_Repl->beta - *pr0 : *pR0 * (pM_Repl->beta - M_DGP.beta0)

		if (robust) { // Compute denominator for this WRE test stat
			denom = smatrix()
			if (NClustVar != 1 | NFE) // collapse meat+sandwich  to all-Cluster-var intersections. If no collapsing needed, _panelsum() will still fold in any weights
				peZVR0 = &_panelsum(pM_Repl->ZVR0, weights? *pwt :* pM_Repl->e : pM_Repl->e, *pinfoErrData)  // really eZAVR0, where e is wildized residual, not residual from replication fit (estar)
			for (c=1; c<=NErrClustCombs; c++) {
				if (NClustVar != 1 & rows(Clust[c].order))
					peZVR0 = &((*peZVR0)[Clust[c].order,])
				if (*pBootClust==Clust[c] & NClustVar & !NFE) { // optimized computation for bootstrapping Cluster when no FE
					AVR0 = pM_Repl->A * pM_Repl->VR0; _beta = -pM_Repl->beta \ 1; betaEnd = _beta[|kEx+1\.|]
					pragma unset t
					for (i=1; i<=NBootClust; i++) {
						pG = &((_beta'XZi[i].M + betaEnd'eZi[i].M * u[i,j]) * AVR0) // R0 * V * Z_i'estar_i
						t = i==1? cross(*pG,*pG) : t + cross(*pG,*pG)
					}
				} else {
					pG = &_panelsum(*peZVR0, Clust[c].info)
					t = cross(*pG,*pG)
				}
				if (Clust[c].multiplier!=1) t = t * Clust[c].multiplier; denom.M = c==1? t : denom.M + t
			}
		} else
			denom.M = (*pR0 * pM_Repl->VR0) * pM_Repl->eec

		Dist[j+thisWeightGrpStart-1] = sqrt? numer/sqrt(denom.M) : cross(numer, invsym(denom.M) * numer)

	}
	if (thisWeightGrpStart==1 & df==2) denom0 = denom.M // original-sample denominator

	return(0) // don't skip remaining calcs in calling procedure
}


real scalar boottestModel::MakeNonWRENumers(real scalar thisWeightGrpStart, real scalar thisWeightGrpStop) {
	pointer (real matrix) scalar pewt; pointer (real colvector) scalar pt
	pragma unused thisWeightGrpStop

	if (!initialized | null) {  // if are imposing null or we are not, but this is first call, then build stuff
		if (thisWeightGrpStart == 1) // first/only weight group? initialize some things
			if (ML)
				eZVR0 = *pSc * (VR0 = *pV * *pR0')
			else if (scoreBS | (robust & granular < NErrClustCombs))
				eZVR0 = pM->e :* pM->ZVR0
		if (scoreBS)
			numer = cross(NClustVar? _panelsum(eZVR0, *pwt, *pinfoBootData) : (weights? eZVR0:* *pwt : eZVR0), u)
		else { // same calc as in score BS but broken apart to grab intermediate stuff
			if (thisWeightGrpStart == 1) { // first/only weight group? initialize some things
				pewt = weights? &(pM->e :* *pwt) : &pM->e
				pt = &_panelsum(*pXEx, *pewt, *pinfoBootData)
				if (AR)
					pt = &(*pt, _panelsum(*pZExcl, *pewt, *pinfoBootData))
				SewtXV = *pM->pV * (*pt)'
			}
			if (!robust | granular)
				numer = *pR0 * (betadev = SewtXV * u)
			else
				numer = (*pR0 * SewtXV) * u
		}
	}

	if      ( AR  ) numer[,1] = u_sd * pM->beta[|kEx+1\.|] // coefficients on excluded instruments in AR OLS
	else if (!null) numer[,1] = u_sd * (*pR0 * (ML? beta : pM->beta) - *pr0) // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this.

	if (initialized & !null) {  // if not imposing null and we have returned, then df=1; and distribution doesn't change with r0, only test stat
		Dist[1] = numer[1] / sqrt(denom.M[1]) * multiplier
		return(1) // skip rest of calc
	}
	return(0)
}









void boottestModel::MakeNonWREStats(real scalar thisWeightGrpStart, real scalar thisWeightGrpStop) {
	real scalar d, i, c, j, l; real matrix eu, eueu, t; real colvector numer_l; pointer (real matrix) scalar peZVR0, pVR0

	if (robust) {
		if (granular < NErrClustCombs & thisWeightGrpStart==1) { // if pure robust and no multi-way clustering, initialized stuff that depends on r0 but not u
			if (scoreBS)
				for (d=df;d;d--)
					Kcd[d].M = J_ClustN_NBootClust // inefficient, but we're not optimizing for the score bootstrap
			else
				for (d=df;d;d--) {
					t = pM->ZVR0[,d]; if (weights) t = t :* *pwt
					if (AR) // final term in formula for K_cd in paper, for c=intersection of all error clusters
						Kcd[d].M = (_panelsum(*pXEx, t, *pinfoErrData), _panelsum(*pZExcl, t, *pinfoErrData)) * SewtXV
					else
						Kcd[d].M =  _panelsum(*pXEx, t, *pinfoErrData                                       ) * SewtXV
				}

			if (NFE & !FEboot & !scoreBS)
				CT_WE = _panelsum(crosstab(wtFE :* pM->e), infoBootAll)'

			// subtract crosstab of E:*ZVR0 wrt bootstrapping cluster combo and all-cluster-var intersections
			peZVR0 = &_panelsum(eZVR0, *pwt, infoAllData) // collapse data to all-boot & error-cluster-var intersections. If no collapsing needed, _panelsum() will still fold in any weights
			if (*pBootClust == Clust[1]) // crosstab c,c* is square
				for (d=df;d;d--) { // if bootstrapping on all-cluster-var intersections (including one-way clustering), the base crosstab is diagonal
					for (i=Clust.N;i;i--)
						Kcd[d].M[i,i] = Kcd[d].M[i,i] - (*peZVR0)[i,d]
					if (scoreBS | !(1 | hascons))
						Kcd[d].M = Kcd[d].M :+ ClustShare * (*peZVR0)[,d]' // for score bootstrap, recenter; "+" because we subtracted *peZVR0
				}
			else
				if (subcluster) // crosstab c,c* is wide
					for (d=df;d;d--) {
						for (i=Clust.N;i;i--) {
							t = infoErrAll[i,]'
							Kcd.M[|(i\i), t|] = Kcd.M[|(i\i), t|] - (*peZVR0)[|t, (d\d)|]'
						}
						if (scoreBS | !(1 |hascons))
							Kcd[d].M = Kcd[d].M - ClustShare * colsum(Kcd[d].M) // for score bootstrap, recenter
					}
				else // crosstab c,c* is tall
					for (d=df;d;d--) {
						for (i=NBootClust;i;i--) {
							t = pBootClust->info[i,]'
							Kcd[d].M[|t, (i\i)|] = Kcd[d].M[|t, (i\i)|] - (*peZVR0)[|t, (d\d)|]
						}
						if (scoreBS | ! (1 |hascons))
							Kcd[d].M = Kcd[d].M - ClustShare * colsum(Kcd[d].M) // for score bootstrap, recenter
					}

			for (c=1+granular; c<=NErrClustCombs; c++) {
				if (rows(Clust[c].order))
					for (d=df;d;d--)
						Kcd[d].M = Kcd[d].M[Clust[c].order,]
				for (d=df;d;d--) {
					pQ[c,d] = &_panelsum(Kcd[d].M, Clust[c].info) // when c=1 (unless subcluster bootstrap), two args have same # of rows, &_panelsum() returns 1st arg by reference. Using & then prevents unnecessary cloning.
					if (!scoreBS & NFE & !FEboot)
						pQ[c,d] = &(*pQ[c,d] + _panelsum(pM->CT_ZVR0[c,d].M, Clust[c].info) * CT_WE)
				}
			}
		}

		if (doQQ) // core denominator computation loop
			for (i=df;i;i--)
				for (j=i;j;j--) {
					if (df>1 | thisWeightGrpStart==1) { // when df=1, only compute QQ in first weight group, since doesn't depend on u; could do same for all QQ[i,j]'s when df>1, but if we're grouping weights, we're worrying about space over speed
						c = NErrClustCombs
						QQ = cross(*pQ[c,i], *pQ[c,j]); if (Clust[c].multiplier!=1) QQ = QQ * Clust[c].multiplier
						for (c--;c;c--) {
							t = cross(*pQ[c,i], *pQ[c,j]); if (Clust[c].multiplier!=1) t = t * Clust[c].multiplier; QQ = QQ + t
						}
					}
					denom[i,j].M = colsum(u :* QQ * u)
				}
		else { // alternative core computational loop, avoiding computing Q'Q which has cubic time cost in numbers of bootstrapping clusters
			if (granular) // prep optimized treatment when bootstrapping by many/small groups
				if (purerobust) {
					eu = *M_DGP.partialFE(&(pM->e :* u)) - *pX * betadev
					eueu = eu:*eu
				} else { // clusters small but not all singletons
					if (NFE) {
						eu = *M_DGP.partialFE(&(pM->e :* u[IDBootData,]))
						for (d=df;d;d--)
							pQ[1,d] = &(_panelsum(eu, pM->WZVR0[d].M, *pinfoErrData)                                          - _panelsum(*pX, pM->WZVR0[d].M, *pinfoErrData) * betadev)
					} else
						for (d=df;d;d--)
							pQ[1,d] = &(_panelsum(_panelsum(pM->e, pM->WZVR0[d].M, infoAllData) :* u[IDBootAll,], infoErrAll) - _panelsum(*pX, pM->WZVR0[d].M, *pinfoErrData) * betadev)
				}

			for (c=NErrClustCombs; c>granular; c--)
				for (d=df;d;d--)
					pQ[c,d] = &(*pQ[c,d] * u)

			for (i=df;i;i--)
				for (j=i;j;j--) {
					denom[i,j].M = purerobust? cross(pM->WZVR0[i].M, pM->WZVR0[j].M, eueu) : colsum(*pQ[1,i] :* *pQ[1,j])

					if (Clust.multiplier!=1) denom[i,j].M = denom[i,j].M  * Clust.multiplier

					for (c=2;c<=NErrClustCombs;c++) {
						t = colsum(*pQ[c,i] :* *pQ[c,j])
						if (Clust[c].multiplier!=1) t = t * Clust[c].multiplier
						denom[i,j].M = denom[i,j].M + t
					}
				}
		}

		if (df == 1) {
			if (MaxMatSize < .)
				Dist[|thisWeightGrpStart \ thisWeightGrpStop|] = (numer :/ sqrt(denom.M))'
			else
				Dist                                           = (numer :/ sqrt(denom.M))'
			if (thisWeightGrpStart==1 & df==2) denom0 = denom.M[1] // original-sample denominator
		} else { // build each replication's denominator from vectors that hold values for each position in denominator, all replications
			t = J(df,df,.)
			for (l=cols(u); l; l--) {
				for (i=df;i;i--)
					for (j=i;j;j--)
						t[i,j] = denom[i,j].M[l]
				_makesymmetric(t)
				numer_l = numer[,l]
				Dist[l+thisWeightGrpStart-1] = cross(numer_l, invsym(t) * numer_l)
			}
			if (thisWeightGrpStart==1 & df==2) denom0 = t // original-sample denominator
		}

	} else { // non-robust

		pVR0 = ML? &VR0 : &(pM->VR0)
		if (df == 1) {  // optimize for one null constraint
			denom.M = *pR0 * *pVR0
			if (!ML) {
										  eu = u :* pM->e
				if (!scoreBS) eu = eu  - (*pXEx, *pZExcl) * betadev // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
				if (scoreBS | !(1 |hascons)) eu = eu :- (weights? cross(ClustShare, eu) : colsum(eu) * ClustShare)  // Center variance if needed
				denom.M = denom.M :* (weights? cross(*pwt, eu :* eu) : colsum(eu :* eu))
			}
			if (MaxMatSize < .)
				Dist[|thisWeightGrpStart \ thisWeightGrpStop|] = (numer :/ sqrt(denom.M))'
			else
				Dist                                           = (numer :/ sqrt(denom.M))'
			if (thisWeightGrpStart==1 & df==2) denom0 = denom.M[0] // original-sample denominator
		} else {
			denom.M = invsym(*pR0 * *pVR0)

			for (l=cols(u); l; l--) {
				numer_l = numer[,l]
				Dist[l+thisWeightGrpStart-1] = cross(numer_l, denom.M * numer_l) 
				if (!(ML | LIML)) {
											  eu = u[,l] :* pM->e
					if (!scoreBS) eu = eu  - (*pXEx, *pZExcl) * betadev[,l] // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
					if (scoreBS | !(1 |hascons)) eu = eu :- (weights? cross(*pwt, eu) : colsum(eu)) * ClustShare // Center variance if needed
					Dist[l+thisWeightGrpStart-1] = Dist[l+thisWeightGrpStart-1] / (t = cross(eu, *pwt, eu))
				}
			}
			if (thisWeightGrpStart==1 & df==2) denom0 = denom.M * t // original-sample denominator
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

// stick call to panelsum() in separate function to prevent run-time error in old Stata versions
real matrix __panelsum(real matrix X, real matrix arg2, real matrix arg3)
	return (cols(arg3)? panelsum(X, arg2, arg3) : panelsum(X, arg2))

// implement Mata's panelsum() for pre-version 13. Differs in that a single missing value in X doesn't make all results missing.
// efficiently handles case when all groups have one row
real matrix _panelsum(real matrix X, real matrix arg2, | real matrix arg3) {
	real matrix retval, Xi, Wi; pointer(real matrix) scalar pinfo; pointer(real colvector) scalar pwt; real scalar i
	pragma unset Xi; pragma unset Wi

	if (args()==2) {
		if (rows(arg2)==0 | rows(arg2)==rows(X))
			return(X)
	} else if (rows(arg3)==0 | rows(arg3)==rows(X))
			return(arg2==1? X : X :* arg2) // if no collapsing called for, still fold in provided weights

	if (stataversion() >= 1300)
		return (__panelsum(X, arg2, arg3))
	
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
			retval[i,] = cross(Wi, Xi)
		}
	}
	return (retval)
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

// cross-tab sum of a column vector w.r.t. intersection-of-error & bootstrap-clustering-vars and fixed-effect var
real matrix boottestModel::crosstab(real colvector v) {
	real matrix retval; real scalar i, j, t; real colvector _FEID, _v
	retval = J(rows(infoAllData), NFE, 0)
	for (i=rows(infoAllData);i;i--) {
		_FEID = panelsubmatrix(*pFEID, i, infoAllData)
		_v    = panelsubmatrix(v     , i, infoAllData)
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
//	mid = lo + ( mreldif(p_hi-p_lo, 1/repsFeas)<1e-6 & repsFeas? (hi - lo)/2 : (alpha-p_lo)/(p_hi-p_lo)*(hi-lo) ) // interpolate linearly until bracketing a single bump-up; then switch to binary search
	if (mreldif(lo,mid)<1e-6 | mreldif(hi,mid)<1e-6 | (repsFeas & abs(p_hi-p_lo)<(1+(ptype==1))/repsFeas*1.000001))
		return (mid)
	if ( ((_p = r0_to_p(mid)) < alpha) == (p_lo < alpha) )

		return(search(alpha, _p, mid, p_hi, hi))
	return(search(alpha, p_lo, lo, _p, mid))
}
	
// derive wild bootstrap-based CI, for case of linear model with one-degree null imposed.
void boottestModel::plot() {
	real scalar t, alpha, _quietly, c, i, j; real colvector lo, hi; pointer (real colvector) scalar _pr0

	_quietly = quietly
	setquietly(1)
	boottest() // run in order to get true number of replications

	alpha = 1 - level*.01
	if (alpha>0 & cols(u)-1 <= 1/alpha-1e6) {
		setquietly(_quietly)
		if (!quietly) errprintf("\nError: need at least %g replications to resolve a %g%% two-sided confidence interval.\n", ceil(1/alpha), level)
		return (.\.)
	}
	
	_pr0 = pr0
	_editmissing(gridpoints, 25)

	if (gridstart[1]==. | gridstop[1]==.) {
		if (reps)
			if (AR) {
				t = abs(cuepoint) / (small? invttail(df_r, alpha/2)/invttail(df_r, getpadj(1)/2) : invnormal(alpha/2)/invnormal(getpadj(1)/2))
				lo = editmissing(gridstart[1], cuepoint - t)
				hi = editmissing(gridstop [1], cuepoint + t)
			} else {
				make_DistCDR()
				lo = editmissing(gridstart[1], numer[1]/u_sd + *pr0 + DistCDR[floor((   alpha/2)*(repsFeas-1))+1] * abs(numer[1]/Dist[1])) // initial guess based on distribution from main test
				hi = editmissing(gridstop [1], numer[1]/u_sd + *pr0 + DistCDR[ceil (( 1-alpha/2)*(repsFeas-1))+1] * abs(numer[1]/Dist[1]))
			}
		else {
			t = abs(numer/Dist) * (small? invttail(df_r, alpha/2) : -invnormal(alpha/2))
			lo = editmissing(gridstart[1], (numer - t)/u_sd + *pr0)
			hi = editmissing(gridstop [1], (numer + t)/u_sd + *pr0)
			if (scoreBS & !null & !willplot) { // if doing simple Wald test with no graph, we're done
				CI = lo, hi
				return
			}
		}
		
		if (gridstart[1]==. & ptype!=3) // unless upper-tailed p value, try at most 10 times to bracket confidence set by symmetrically widening
			for (i=10; i & -r0_to_p(lo)<-alpha; i--) {
				t = hi - lo
				lo = lo - t
				if (gridstop[1]==.) hi = hi + t // maintain rough symmetry unless user specified upper bound
			}
		if (gridstop[1]==. & ptype!=2) // ditto for high side
			for (i=10; i & -r0_to_p(hi)<-alpha; i--) {
				t = hi - lo
				if (gridstart[1]==.) lo = lo - t // maintain rough symmetry unless user specified lower bound
				hi = hi + t
			}
	} else {
		lo = gridstart[1]
		hi = gridstop[1]
	}

	plotX = rangen(lo, hi, gridpoints[1])
	if (cuepoint == .) cuepoint = numer[1] / u_sd + *pr0 // non-AR case
	if (cuepoint < lo) { // insert original point estimate into grid
		if (gridstart[1] == .) {
			plotX = cuepoint \ plotX
			c = 1
		}
	} else if (cuepoint > hi) {
		if (gridstop[1] == .) {
			plotX = plotX \ cuepoint
			c = gridpoints[1]+1
		}
	} else {
		c = floor((cuepoint - lo)/(hi - lo)*(gridpoints[1] - 1)) + 2
		plotX = plotX[|.\c-1|] \ cuepoint \ plotX[|c\.|]
	}

	plotY = J(rows(plotX), 1, .)
	printf("{txt}")
	for (i = rows(plotX); i; i--) {
		plotY[i] = r0_to_p(plotX[i])
		if (!_quietly) {
			printf(".")
			if (mod(i-rows(plotX)-1,50)) displayflush()
				else printf("\n")
		}
	}
	printf("\n")

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
		} else if (c==gridpoints[1]+1) {
			plotX = plotX[|.\gridpoints[1]|]; plotY = plotY[|.\gridpoints[1]|]
		} else {
			plotX = plotX[|.\c-1|] \ plotX[|c+1\.|]; plotY = plotY[|.\c-1|] \ plotY[|c+1\.|]
		}
	}

	setquietly(_quietly)
	pr0 = _pr0
	plotted = 1
}

void boottestModel::contourplot() {
	real scalar t, _quietly, i, d; real colvector lo, hi; pointer (real colvector) scalar _pr0

	_quietly = quietly
	setquietly(1)
	_pr0 = pr0

	_editmissing(gridpoints, 25)
	lo = hi = J(2, 1, .)

	if (AR)
		t = abs(cuepoint) / (-2/invnormal(getpadj(1)/2))

	for(d=df;d;d--) {
		if (AR) {
			lo[d] = editmissing(gridstart[d], cuepoint[d] - t[d])
			hi[d] = editmissing(gridstart[d], cuepoint[d] + t[d])
		} else {
			lo[d] = editmissing(gridstart[d], (numer[d,1] + (*pr0)[d] - 3 * sqrt(denom0[d,d]))/u_sd)
			hi[d] = editmissing(gridstop [d], (numer[d,1] + (*pr0)[d] + 3 * sqrt(denom0[d,d]))/u_sd)
		}

		stata("_natscale " + strofreal(lo[d]) + " " + strofreal(hi[d]) + " 4")
		if (gridstart[d]==.) {
			stata("local min = r(min)")  // for some reason, st_global("r(min)") doesn't work
			lo[d] = strtoreal(st_local("min"))
		}
		if (gridstop[d]==.) {
			stata("local max = r(max)")
			hi[d] = strtoreal(st_local("max"))
		}
	}

	plotX = (rangen(lo[1], hi[1], gridpoints[1]) # J(gridpoints[2],1,1)), (J(gridpoints[1],1,1) # rangen(lo[2], hi[2], gridpoints[2]))
	plotY = J(rows(plotX), 1, .)
	printf("{txt}")
	for (i = rows(plotX); i; i--) {
		plotY[i] = r0_to_p(plotX[i,]')
		if (!_quietly) {
			printf(".")
			if (mod(i-rows(plotX)-1,50)) displayflush()
				else printf("\n")
		}
	}
	printf("\n")

	setquietly(_quietly)
	pr0 = _pr0
	plotted = 1
}

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

// Stata interface
void boottest_stata(string scalar statname, string scalar dfname, string scalar dfrname, string scalar pname, string scalar padjname, string scalar ciname, 
	string scalar plotname, string scalar peakname, real scalar level, real scalar ML, real scalar LIML, real scalar Fuller, 
	real scalar K, real scalar AR, real scalar null, real scalar scoreBS, string scalar weighttype, string scalar ptype, string scalar madjtype, real scalar NumH0s,
	string scalar XExnames, string scalar XEndnames, real scalar hascons, string scalar Ynames, string scalar bname, string scalar Vname, string scalar Wname, 
	string scalar ZExclnames, string scalar samplename, string scalar scnames, real scalar robust, string scalar IDnames, real scalar NBootClustVar, real scalar NErrClust, 
	string scalar FEname, real scalar NFE, string scalar wtname, string scalar wttype, string scalar Cname, string scalar C0name, real scalar reps, string scalar repsname, string scalar repsFeasname, 
	real scalar small, string scalar diststat, string scalar distname, string scalar gridstart, string scalar gridstop, string scalar gridpoints, real scalar MaxMatSize, real scalar quietly) {

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

	M.MaxMatSize(MaxMatSize)
	M.sethascons(hascons)
	M.setsc(sc)
	M.setML(ML)
	M.setY (Y)
	M.setZExcl(ZExcl)
	M.setwt (wt)
	M.setID(ID, NBootClustVar, NErrClust)
	M.setFEID(FEID, NFE)
	M.setR (R , r )
	M.setR0(R0, r0)
	M.setnull(null)
	M.setsmall(small)
	M.setrobust(robust)
	M.setscoreBS(scoreBS)
	M.setweighttype(weighttype)
	M.setptype(ptype)
	M.setwttype(wttype)
	M.setreps(reps)
	M.setLIML(LIML)
	M.setFuller(Fuller)
	M.setk(K)
	M.setAR(AR)
	M.setgrid(strtoreal(tokens(gridstart)), strtoreal(tokens(gridstop)), strtoreal(tokens(gridpoints)))
	M.setmadjust(madjtype, NumH0s)
	M.setlevel(level)
	M.setquietly(quietly)

	M._st_view(XEnd, ., XEndnames, samplename)
	M.setXEnd(XEnd)
	M._st_view(XEx, ., XExnames, samplename)
	M.setXEx(XEx)
	if (bname != "") M.setbeta(st_matrix(bname)')
	if (Vname != "") M.setV   (st_matrix(Vname) )
	if (Wname != "") M.setW   (st_matrix(Wname) )
	M.setwillplot(plotname != "") // can make point estimate a little faster if not going to plot

	st_numscalar(statname, M.getstat())
	st_numscalar(pname   , M.getp   ())
	st_numscalar(repsname, M.getreps())
	st_numscalar(repsFeasname, M.getrepsFeas())
	st_numscalar(padjname, M.getpadj())
	st_numscalar(dfname  , M.getdf  ())
	st_numscalar(dfrname , M.getdf_r())
	if (distname != "") st_matrix(distname, M.getdist(diststat))
	if (plotname != "" | (level<100 & ciname != "")) {
		if (plotname != "") st_matrix(plotname, M.getplot())
		if (cols(M.peak)) st_matrix(peakname, M.getpeak())
		if (level<100 & ciname != "") st_matrix(ciname, M.getCI())
	}

	M.M_DGP.setParent(NULL) // actually sets the pointer to &NULL, but suffices to break loop in the data structure topology and avoid Mata garbage-cleaning leak
}

mata mlib create lboottest, dir("`c(sysdir_plus)'l") replace
mata mlib add lboottest *(), dir("`c(sysdir_plus)'l")
mata mlib index
end
