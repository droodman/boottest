*! boottest 3.0.2 18 December 2020
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
	real scalar LIML, YY, eec, Fuller, ARubin, kappa, k, isDGP
	real matrix ZX, ZXend, ZZ, H_2SLS, invH, V, AR, ZAR, uddot2, Tplus, pi, XXend, dbetadt, ZexclXend, XexXend
	real colvector t0, uddot, beta, beta0, ZY
	real rowvector YXEnd
	pointer(real colvector) scalar pY
	pointer(real matrix) scalar pXend, pXX, pXY, pZexclY, pXexY, pA, pW, pinvZZ, pXexXex, pZXex, pH, pR, pT
	pointer (class boottestModel scalar) scalar parent
	struct smatrix matrix CT_ZAR
	struct smatrix colvector WZAR

	void new(), InitExog(), InitEndog(), InitTestDenoms(), SetT(), InitEstimate(), Estimate(), SetAR()
	pointer(real matrix) scalar partialFE()
}

class boottestModel {
	real scalar scoreBS, B, small, weighttype, null, dirty, initialized, Neq, ML, GMM, Nobs, _Nobs, k, kEx, el, sumwt, NClustVar, weights, REst, multiplier, smallsample, quietly, FEboot, NErrClustCombs, ///
		sqrt, hascons, LIML, Fuller, kappa, IV, WRE, WREnonARubin, ptype, twotailed, df, df_r, ARubin, D, confpeak, willplot, notplotted, NumH0s, p, NBootClustVar, NErrClust, ///
		NFE, granular, purerobust, subcluster, NBootClust, BFeast, u_sd, level, ptol, MaxMatSize, Ng, enumerate, bootstrapt, q1, q, interpolable, interpolating, interpolate_u, robust
	real matrix AR, numer, v, ustar, T, TARubin, T0, L0invR0L0, plot, CI, CT_WE, infoBootData, infoBootAll, infoErrAll, J_ClustN_NBootClust, statDenom, uZAR, SuwtXA, numer0, betadev
	real colvector DistCDR, t, tARubin, plotX, plotY, t0, beta, wtFE, ClustShare, IDBootData, IDBootAll, WeightGrpStart, WeightGrpStop, gridmin, gridmax, gridpoints, numersum, r00, uddot0, anchor, poles
	real rowvector peak
	string scalar wttype, madjtype, seed
	pointer (real matrix) scalar pZexcl, pR1, pR, pID, pFEID, pXend, pXex, pX, pinfoAllData, pinfoErrData
	pointer (real colvector) scalar pr1, pr, pY, pSc, pwt, pW, pA, puddot, pDist
	class AnalyticalModel scalar M_DGP
	pointer (class AnalyticalModel scalar) scalar pM_Repl, pM
	struct structboottestClust colvector Clust
	pointer (struct structboottestClust scalar) scalar pBootClust
	struct smatrix matrix denom, Kcd, denom0, Jcd0
	struct smatrix colvector Kd, XZi, eZi, euZAR, dudr, dnumerdr
  struct ssmatrix colvector ddenomdr, dJcddr
  struct ssmatrix matrix ddenomdr2
	pointer(struct smatrix matrix) scalar pJcd
	struct structFE rowvector FEs

	void new(), setsqrt(), boottest(), plot(), setXEx(), setptype(), setdirty(), setXEnd(), setY(), setZExcl(), setwt(), setsc(), setML(), setLIML(), setARubin(),
		setFuller(), setkappa(), setquietly(), setbeta(), setA(), setW(), setsmall(), sethascons(), setscoreBS(), setB(), setnull(), setWald(), setRao(), setwttype(), setID(), setFEID(), setlevel(), setptol(), 
		setrobust(), setR1(), setR(), setwillplot(), setgrid(), setmadjust(), setweighttype(), makeWildWeights(), makeInterpolables(), _makeInterpolables(), makeNonWREStats(), makeBootstrapcDenom(), setMaxMatSize(), setstattype(), Initialize()
  private void makeNumerAndJ(), _clustAccum(), makeWREStats(), PrepARubin()
	real matrix getplot(), getCI(), getV(), getv()
	real scalar getp(), getpadj(), getstat(), getdf(), getdf_r(), getreps(), getrepsFeas(), getNBootClust()
	private real scalar r_to_p(), search()
  private real matrix count_binary(), crosstab()
  private static real matrix combs()
  private static real colvector stableorder()
  static void _st_view()
	real rowvector getpeak()
	real colvector getdist(), getb()
}

void AnalyticalModel::new() {
	ARubin = 0
  isDGP = 1  // by default, created object is for DGP rather than replication regression
}

void AnalyticalModel::SetAR(real scalar _AR) {
	if (ARubin = _AR) LIML = Fuller = kappa = 0
}

// stuff that can be done before r set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
void AnalyticalModel::InitExog() {
	real matrix ZexclXex

	parent->pXex = partialFE(parent->pXex)
	pXexXex = &cross(*parent->pXex, *parent->pwt, *parent->pXex)
	if (cols(*parent->pZexcl)) {  // GMM, 2SLS, LIML
		parent->pZexcl = partialFE(parent->pZexcl)
		ZexclXex = cross(*parent->pZexcl, *parent->pwt, *parent->pXex)
		pZXex = &(*pXexXex \ ZexclXex)
		if (parent->IV)  // non-GMM
			pinvZZ = &invsym((ZZ = *pZXex, (ZexclXex' \ cross(*parent->pZexcl, *parent->pwt, *parent->pZexcl))))
	} else
		pXX = pXexXex
}

// stuff that can be done before T & r set, but depend on endogenous variables, which are bootstrapped in WRE
void AnalyticalModel::InitEndog(pointer (real colvector) scalar _pY, pointer (real matrix) scalar _pXend, | ///
		pointer (real colvector) scalar _pZexclY, pointer (real rowvector) scalar _pXexY) {

	pY = partialFE(_pY); pXend = partialFE(_pXend)
	pXexY = _pXexY==NULL? &cross(*parent->pXex, *parent->pwt, *pY) : _pXexY
	if (kappa | ARubin)
		pZexclY = _pZexclY==NULL? &cross(*parent->pZexcl, *parent->pwt, *pY) : _pZexclY
	if (kappa) {
		XexXend   = cross(*parent->pXex  , *parent->pwt, *pXend)
		ZexclXend = cross(*parent->pZexcl, *parent->pwt, *pXend)
		ZXend = XexXend \ ZexclXend
		ZX = *pZXex, ZXend
		XXend = XexXend \ cross(*pXend, *parent->pwt, *pXend)
		pXX = &(*pXexXex,  XexXend \ XXend')
		YXEnd = cross(*pY, *parent->pwt, *pXend)
		ZY = *pXexY \ *pZexclY
		pXY = &(*pXexY \ YXEnd')
		if (LIML | (parent->robust | parent->scoreBS)==0)
			YY = cross(*pY, *parent->pwt, *pY)
		if (parent->IV)  // if GMM weight matrix not provided, prepare 2SLS one
			V = (I(parent->kEx) \ J(parent->el-parent->kEx, parent->kEx, 0)), *pinvZZ * ZXend // 2SLS is (V' ZX)^-1 * (V'ZY). Also apparently used in k-class and LIML robust VCV by Stata convention
		else
			V = *parent->pW * ZX
		H_2SLS = V ' ZX  // Hessian
	} else {  // OLS / ARubin
		pXY = pXexY
		if (ARubin) {
			pXY = &(*pXY \ *pZexclY)
			pXX = &ZZ
		}
	}
	k = cols(*pXX)
}

void AnalyticalModel::SetT(real matrix T) {  // set model constraint matrix (not null constraints)
	pT = &T
	if (LIML) Tplus = blockdiag(T, 1)  // add an entry to T for the dep var
}

// stuff that can be done before r set but depends on T and endogenous variables
void AnalyticalModel::InitEstimate() {
	real rowvector val
	real matrix _ZY, TT, TPZT, vec
	pointer (real matrix) scalar pbetadenom
	pragma unset vec; pragma unset val

 	if (LIML)
		if (parent->el == k)  // exactly identified LIML = 2SLS
			kappa = 1
		else {
			_ZY = ZXend, ZY
			TT = *pXX, (*pXexY \ YXEnd') \ *pXexY', YXEnd, YY  // where T = *pXex, *pXend, *pY
			(TPZT = TT)[|parent->kEx+1,parent->kEx+1\.,.|] = _ZY ' (*pinvZZ) * _ZY

			V = *pinvZZ * ZX
			H_2SLS = V ' ZX  // Hessian

			if (rows(*pT)) {  // includes H0 if pSAll really points to T0 rather than T
				TT   = Tplus '   TT * Tplus
				TPZT = Tplus ' TPZT * Tplus
			}
			eigensystemselecti( I(rows(TT)) - invsym(TT) * TPZT, 1\1, vec, val)  // eigensystemselecti(invsym(TT) * TPZT, rows(TT)\rows(TT), ... gives 1 - the eigenvalue, but can cause eigensystem() to return all missing
			kappa = 1/Re(val) - Fuller / (parent->_Nobs - parent->el)   // sometimes a tiny imaginary component sneaks into val
		}

  pH = kappa? (kappa==1? &H_2SLS : &((1-kappa)* *pXX + kappa*H_2SLS)) : pXX
	if (rows(*pT)) {
		pbetadenom = &(*pT * invsym(*pT ' (*pH) * *pT) * *pT')
		invH = J(0,0,0)
	} else
		pbetadenom = &(invH = invsym(*pH))

  if (kappa)
		if (kappa==1) {  // 2SLS, GMM
			beta0 = *pbetadenom * V ' ZY
			dbetadt = I(rows(beta0)) - *pbetadenom * V ' ZX
		} else {  // k-class, LIML
			beta0 = *pbetadenom * (kappa * V ' ZY + (1-kappa) * *pXY)
			dbetadt = I(rows(beta0)) - *pbetadenom * (kappa * V ' ZX + (1-kappa) * *pXX)
		}
	else {  // OLS / ARubin
		beta0 = *pbetadenom * *pXY
		dbetadt = I(rows(beta0)) - *pbetadenom * *pXX
	}
}

// stuff that doesn't depend on r, for test stat denominators in replication regressions
void AnalyticalModel::InitTestDenoms(real matrix T) {
	real matrix VAR; real scalar d, c; struct smatrix rowvector _CT_ZAR; pointer (real matrix) scalar pWZAR

	pA = rows(T)? &(T * invsym(T ' (*pH) * T) * T') : ///
	              &(rows(invH)? invH : invsym(*pH))
	AR = *pA * *parent->pR'

	if (parent->scoreBS | (parent->robust & (parent->WREnonARubin & parent->NClustVar==1 & parent->NFE==0)==0)) {
		if (kappa) {
			VAR = V * AR
			ZAR = *parent->pZexcl * VAR[|parent->kEx+1,.\.,.|]; if (parent->kEx) ZAR = ZAR + *parent->pXex * VAR[|.,.\parent->kEx,.|]
		} else if (ARubin) {
			ZAR = *parent->pZexcl * AR[|parent->kEx+1,.\.,.|]
			if (cols(*parent->pXex))
				ZAR = ZAR + *parent->pXex * AR[|.,.\parent->kEx,.|]
		} else
			ZAR = *parent->pXex * AR

			if (parent->bootstrapt == 0) return

			if (parent->granular) {
				pWZAR = &(parent->weights? ZAR :* *parent->pwt : ZAR)
				WZAR = smatrix(parent->df)
				for (d=parent->df;d;d--)
					WZAR[d].M = (*pWZAR)[,d]
			}
			
      if (parent->NFE & parent->robust & (parent->WREnonARubin | parent->FEboot | parent->scoreBS)==0 & parent->granular < parent->NErrClustCombs) {
				if (pWZAR==NULL) pWZAR = &(parent->weights? ZAR :* *parent->pwt : ZAR)
				CT_ZAR = smatrix(parent->NErrClustCombs, parent->df)
				for (d=parent->df;d;d--)
					CT_ZAR[1,d].M = parent->crosstab((*pWZAR)[,d])
				if (parent->NClustVar > 1) {
					_CT_ZAR = CT_ZAR[1,]
					for (c=2; c<=parent->NErrClustCombs; c++)
						for (d=parent->df;d;d--) {
							if (rows(parent->Clust[c].order))
								_CT_ZAR[d].M = _CT_ZAR[d].M[parent->Clust[c].order,]
							CT_ZAR[c,d].M = _panelsum(_CT_ZAR[d].M, parent->Clust[c].info)
						}
				}
			}
	}
}

// stuff that depends on r and endogenous variables: compute beta and residuals
void AnalyticalModel::Estimate(real colvector t) {
	real colvector negZeinvee; real matrix Ze; real scalar ee

	beta = rows(t)? beta0 + dbetadt * t : beta0
	if (isDGP | parent->bootstrapt | parent->WREnonARubin==0) {  // don't need residuals in replication regressions in bootstrap-c on WRE/non-ARubin
		if (ARubin) {
			uddot = *pY - *parent->pZexcl * beta[|cols(*parent->pXex)+1\.|]
			if (cols(*parent->pXex))
				uddot =  uddot - *parent->pXex * beta[|.\cols(*parent->pXex)|]
		} else if (kappa)
			if (parent->kEx == 0)
				uddot = *pY - *pXend * beta[|parent->kEx+1\.|]
			else
				uddot = *pY - *pXend * beta[|parent->kEx+1\.|] - *parent->pXex * beta[|.\parent->kEx|]
		else
				uddot = *pY                                    - *parent->pXex * beta[|.\parent->kEx|]

    if ((parent->robust | parent->scoreBS)==0 | (isDGP & LIML))  // useful in non-robust, residual-based bootstrap, and in computing uddot^2 in LIML (just below)
			ee = YY - 2 * *pXY ' beta + beta ' (*pXX) * beta
		if ((parent->robust | parent->scoreBS)==0 & parent->bootstrapt==0)
			eec = parent->hascons? ee : ee - (parent->weights? cross(uddot, *parent->pwt) : sum(uddot))^2 / parent->_Nobs  // sum of squares after centering, N * Var
	}

	if (isDGP & LIML) {
		Ze = ZY - ZX * beta
		negZeinvee = Ze / -ee
		pi = invsym(ZZ + negZeinvee * Ze') * (negZeinvee * (YXEnd - beta ' XXend) + ZXend)  // coefficients in reduced-form equations; Davidson & MacKinnon (2010), eq 15
		uddot2 = *pXend - *parent->pZexcl * pi[|parent->kEx+1,.\.,.|]; if (parent->kEx) uddot2 = uddot2 - *parent->pXex * pi[|.,.\parent->kEx,.|]
	}
}

// partial fixed effects out of a data matrix
pointer(real matrix) scalar AnalyticalModel::partialFE(pointer(real matrix) scalar pIn) {
	real matrix Out, tmp; real scalar i
	if (parent->NFE & pIn!=NULL) {
		Out = *pIn
		for (i=parent->NFE;i;i--) {
			tmp = Out[parent->FEs[i].is,]
			Out[parent->FEs[i].is,] = tmp :- cross(parent->FEs[i].wt, tmp)
		}
		return(&Out)
	}
	return (pIn)
}

void boottestModel::new() {
	ARubin = LIML = Fuller = WRE = small = scoreBS = weighttype = Neq = ML = initialized = quietly = sqrt = hascons = IV = ptype = robust = NFE = FEboot = granular = NErrClustCombs = subcluster = B = BFeast = interpolating = 0
	twotailed = null = dirty = willplot = u_sd = bootstrapt = notplotted = 1
	level = 95
  ptol = 1e-6
	confpeak = MaxMatSize = .
	pXend = pXex = pZexcl = pY = pSc = pID = pFEID = pR1 = pR = pwt = &J(0,0,0)
	pr1 = pr = &J(0,1,0)
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
	this.pXend  = &X; setdirty(1)
}
void boottestModel::setXEx     (real matrix X ) {
	this.pXex  = &X; setdirty(1)
}
void boottestModel::setY       (real matrix Y ) {
	this.pY  = &Y; setdirty(1)
}
void boottestModel::setZExcl   (real matrix Z ) {
	this.pZexcl  = &Z; setdirty(1)
}
void boottestModel::setwt      (real matrix wt) {
	this.pwt  = &wt; setdirty(1)
}
void boottestModel::setsc(real matrix Sc) {
	this.pSc  = &Sc
	setdirty(1)
}
void boottestModel::setML(real scalar ML) {
	this.ML  = ML; setdirty(1)
	if (ML) setscoreBS(1)
}
void boottestModel::setLIML(real scalar LIML) {
	this.LIML = LIML; setdirty(1)
}
void boottestModel::setARubin(real scalar ARubin) {
	this.ARubin = ARubin; setdirty(1)
}
void boottestModel::setFuller    (real scalar Fuller) {
	this.Fuller = Fuller; setdirty(1)
}
void boottestModel::setkappa(real scalar kappa) {  // kappa as in k-class
	this.kappa = kappa; setdirty(1)
}
void boottestModel::setquietly(real scalar quietly )
	this.quietly = quietly
void boottestModel::setbeta(real colvector beta) {
	this.beta = beta; setdirty(1)
}
void boottestModel::setA(real matrix V) {
	this.pA = &V; setdirty(1)
}
void boottestModel::setW(real matrix W) {
	this.pW = &W; setdirty(1)
}
void boottestModel::setsmall(real scalar small) {
	this.small = small; setdirty(1)
}
void boottestModel::sethascons(real scalar hascons) {
	this.hascons = hascons; setdirty(1)
}
void boottestModel::setscoreBS (real scalar scoreBS) {
	this.scoreBS = scoreBS; setdirty(1)
}
void boottestModel::setB(real scalar B) {
	this.B = B
	if (B==0)
		setscoreBS(1)
	setdirty(1)
}
void boottestModel::setnull    (real scalar null) {
	this.null = null; setdirty(1)
}
void boottestModel::setWald() { // set-up for classical Wald test
	this.scoreBS = 1; this.B = 0; this.null = 0; setdirty(1)
}
void boottestModel::setRao() { // set-up for classical Rao test
	this.scoreBS = 1; this.B = 0; this.null = 1; setdirty(1)
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
void boottestModel::setlevel(real scalar level)
	this.level = level
void boottestModel::setptol(real scalar ptol)
	this.ptol = ptol
void boottestModel::setrobust(real scalar robust) {
	this.robust = robust
	if (robust==0) setID(J(0,0,0), 1, 1)
	setdirty(1)
}
void boottestModel::setR1(real matrix R1, real colvector r1) {
	this.pR1 = &R1; 	this.pr1  = &r1; q1 = rows(R1); setdirty(1)
}
void boottestModel::setR(real matrix R, real colvector r) {
	this.pR = &R; this.pr = &r; q = rows(R); setdirty(1)  // q can differ from df in ARubin test
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
		_sort( DistCDR = (*pnumer)[|2\.|]' :+ *pr , 1)
	} else if (rows(DistCDR)==0)
		if (rows(*pDist)>1)
			_sort( DistCDR=(*pDist)[|2\.|] , 1)
		else
			DistCDR = J(0,1,0)
	return(DistCDR)
}

// Robust to missing bootstrapped values interpreted as +infinity.
real scalar boottestModel::getp(|real scalar classical) {
	real scalar tmp
	if (dirty) boottest()
	tmp = (*pDist)[1]
	if (tmp == .) return (.)
	if (B & classical==.)
		if (sqrt & ptype != 3) {
			if (ptype==0)
					p = colsum(-abs(tmp) :> -abs(*pDist)) / BFeast  // symmetric p value; do so as not to count missing entries in Dist
			else if (ptype==1)  // equal-tail p value
				p = 2 * min((colsum(tmp :> *pDist) , colsum(-tmp:>- *pDist))) / BFeast
			else
				p = colsum( tmp :>   *pDist) / BFeast  // lower-tailed p value
		} else
				p = colsum(-tmp :> - *pDist) / BFeast  // upper-tailed p value or p value based on squared stats
	else {
		tmp = tmp * multiplier
    p = small? Ftail(df, df_r, sqrt? tmp*tmp : tmp) : chi2tail(df, sqrt? tmp*tmp : tmp)
		if (sqrt & twotailed==0) {
			p = p / 2
			if ((ptype==3) == (tmp<0))
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
	return(u_sd==1? v[|.,2\.,.|] : v[|.,2\.,.|] / u_sd)

// Return number of bootstrap replications with feasible results
// Returns 0 if getp() not yet accessed, or doing non-bootstrapping tests
real scalar boottestModel::getrepsFeas()
	return (BFeast)

real scalar boottestModel::getNBootClust()
	return (NBootClust)

// return number of replications, possibly reduced to 2^G
real scalar boottestModel::getreps()
	return (B)

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
real rowvector boottestModel::getpeak() {  // x and y values of confidence curve peak (at least in OLS & ARubin)
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


void boottestModel::Initialize() {  // for efficiency when varying r repeatedly to make CI, do stuff once that doesn't depend on r
  real colvector sortID, o, _FEID
	real rowvector val, ClustCols
	real matrix r, L, L0, vec, Combs, tmp, IDErr
	real scalar i, j, c, minN, sumN, _reps, i_FE, q1q
	pointer (real matrix) scalar _pR0, pIDAll
	class AnalyticalModel scalar M_WRE
	pragma unset vec; pragma unset val; pragma unused M_WRE

  Nobs = rows(*pXex)
  NClustVar = cols(*pID)
  k  = cols(*pR)
  kEx = cols(*pXex)
  if (cols(*pZexcl)==0) pZexcl = &J(Nobs,0,0)
  if (cols(*pXend )==0) pXend  = &J(Nobs,0,0)
  D = cols(*pXend) + 1
  REst = rows(*pR1)  // base model contains restrictions?
  if (pZexcl != NULL) el = cols(*pZexcl) + kEx  // l = # of exogenous variables
  if (kappa==.) kappa = cols(*pZexcl)>0  // if kappa in kappa-class estimation not specified, it's 0 or 1 for OLS or 2SLS
  GMM = pW != NULL
  IV = kappa & GMM==0
  WRE = (kappa & scoreBS==0) | ARubin
  WREnonARubin = WRE & ARubin==0

  if (weights = rows(*pwt)>1)
    sumwt = sum(*pwt)
  else
    pwt = &(sumwt = 1)
  _Nobs = weights & wttype=="fweight"? sumwt : Nobs

  if (WREnonARubin)
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

      pinfoAllData = NClustVar > NBootClustVar? &_panelsetup(*pID,            1..NClustVar) : &infoBootData  // info for grouping by intersections of all bootstrap & clustering vars wrt data; used to speed crosstab EZAR wrt bootstrapping cluster & intersection of all error clusterings
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

      if (scoreBS | WREnonARubin==0)
        ClustShare = weights? _panelsum(*pwt, *pinfoErrData)/sumwt : ((*pinfoErrData)[,2]-(*pinfoErrData)[,1]:+ 1)/Nobs // share of observations by group 

    } else {  // if no clustering, cast "robust" as clustering by observation
      pBootClust = &(Clust = structboottestClust())
      Clust.multiplier = small? _Nobs / (_Nobs - 1) : 1
      Clust.even = 1
      sumN = Clust.N = Nobs
      Clust.info = J(Nobs, 0, 0)  // signals _panelsum not to aggregate
      NErrClustCombs = 1
      if (scoreBS | WREnonARubin==0)
        ClustShare = weights? *pwt/sumwt : 1/_Nobs
    }

    purerobust = NClustVar & (scoreBS | subcluster)==0 & NBootClust==Nobs  // do we ever error-cluster *and* bootstrap-cluster by individual?
    granular   = NClustVar & scoreBS==0 & (purerobust | (Clust.N+NBootClust)*k*B + (Clust.N-NBootClust)*B + k*B < Clust.N*k*k + Nobs*k + Clust.N * NBootClust * k + Clust.N * NBootClust)

    if (robust & purerobust==0) {
      if (subcluster | granular)
        infoErrAll = _panelsetup(*pIDAll, subcluster+1..NClustVar)  // info for error clusters wrt data collapsed to intersections of all bootstrapping & error clusters; used to speed crosstab EZAR wrt bootstrapping cluster & intersection of all error clusterings
      if (scoreBS & B)
        J_ClustN_NBootClust = J(Clust.N, NBootClust, 0)
    }
  } else
    minN = rows(infoBootData)

  if (NFE) {
    sortID = (*pFEID)[o = stableorder(*pFEID, 1)]
    i_FE = 1; FEboot = B>0; j = Nobs; _FEID = wtFE = J(Nobs, 1, 1)
    FEs = structFE(NFE)
    for (i=Nobs-1;i;i--) {
      if (sortID[i] != sortID[i+1]) {
        FEs[i_FE].is = o[|i+1\j|]
        if (weights) {
          tmp  = (*pwt)[FEs[i_FE].is]
          FEs[i_FE].wt = tmp / colsum(tmp)
        } else
          FEs[i_FE].wt = J(j-i, 1, 1/(j-i))
        if (B & robust & granular < NErrClust)
          wtFE[FEs[i_FE].is] = FEs[i_FE].wt

        j = i
        
        if (B & FEboot & NClustVar) {  // are all of this FE's obs in same bootstrapping cluster? (But no need to check if B=0 for then CT_WE in 2nd term of (62) orthogonal to v = col of 1's)
          tmp = (*pID)[FEs[i_FE].is, 1..NBootClustVar]
          FEboot = all(tmp :== tmp[1,])
        }
        ++i_FE
      }
      _FEID[o[i]] = i_FE
    }
    FEs[NFE].is = FEs[NFE].is = o[|.\j|]
    if (weights) {
      tmp  = (*pwt)[FEs[NFE].is]
      FEs[NFE].wt = tmp / colsum(tmp)
    } else
      FEs[NFE].wt = J(j-i,1,1/(j-i))
    if (B & robust & granular < NErrClust)
      wtFE[FEs[NFE].is] = FEs[NFE].wt

    if (B & FEboot & NClustVar) {  // are all of this FE's obs in same bootstrapping cluster?
      tmp = (*pID)[FEs[NFE].is, 1..NBootClustVar]
      FEboot = all(tmp :== tmp[1,])
    }

    pFEID = &_FEID  // ordinal fixed effect ID

    if (robust & FEboot==0 & granular < NErrClust & B & FEboot==0 & bootstrapt)
      infoBootAll = _panelsetup(*pIDAll, 1..NBootClustVar)  // info for bootstrapping clusters wrt data collapsed to intersections of all bootstrapping & error clusters
  }

  if (granular & (WREnonARubin & purerobust)==0)
    if (NFE)
      (void) _panelsetup(*pID   , 1..NBootClustVar, IDBootData)
    else
      (void) _panelsetup(*pIDAll, 1..NBootClustVar, IDBootAll )

  if (ML)
    df = rows(*pR)
  else {
    if (REst) {  // restricted estimation, e.g., cnsreg?
      symeigensystem(*pR1 ' invsym(*pR1 * *pR1') * *pR1, vec, val)  // make "inverse" T,t of constraint matrices; formulas adapted from [P] makecns
      L = vec[|.,.\.,q1|]  // eigenvectors not in kernel of projection onto R1
      T = q1 < k? vec[|.,q1+1\.,.|] : J(k,0,0)  // eigenvectors in kernel
      t = L * luinv(*pR1 * L) * *pr1
      if (ARubin) {
        TARubin = blockdiag(T[|.,. \ k-cols(*pXend) , k-q1-cols(*pXend)|] , I(cols(*pZexcl))) // adapt T,t from XExog, XEndog to XExog, ZEXcl. Assumes no constraints link XExog and XEndog
        tARubin = t[|.\k-cols(*pXend)|] \ J(cols(*pZexcl),1,0)
      }
    }

    // Estimation with null imposed along with any model constraints; in IV, Z is unconstrained regardless of overlap with potentially constrained X
    _pR0 = null? pR : &J(0, k, 0)
    r = REst? *pR1 \ *_pR0 : *_pR0  // combine model and hypothesis constraints to prepare to "invert" them as a group too
    if (q1q = rows(r)) {
      L0 = invsym(r * r')
      if (all(diagonal(L0))==0)
        _error(111, "A null hypothesis constraint is inconsistent or redundant.")
      symeigensystem(r ' L0 * r, vec, val)
      L0  = vec[|.,. \ .,q1q|]
      T0 = q1q < cols(vec)? vec[|.,q1q+1 \ .,.|] : J(rows(vec), 0, 0)
      L0invR0L0 = L0 * luinv(r * L0)
    } else
      T0 = J(0,0,0)

    M_DGP.parent = &this
    M_DGP.InitExog()

    if (WRE) {
      pM_Repl = &(M_WRE = M_DGP)
      pM_Repl->isDGP = 0
      M_DGP.LIML = 1; M_DGP.Fuller = 0; M_DGP.kappa = 1
      if (ARubin==0) { 
        pM_Repl->LIML = this.LIML; pM_Repl->Fuller = this.Fuller; pM_Repl->kappa = this.kappa
      }
      pM_Repl->SetAR(ARubin)
      pM_Repl->SetT(ARubin? TARubin : T)
    } else {
      M_DGP.LIML = this.LIML; M_DGP.Fuller = this.Fuller; M_DGP.kappa = this.kappa
    }

    M_DGP.InitEndog(pY, pXend)

    if (ARubin) {
      if (willplot) {  // for plotting/CI purposes get original point estimate if not normally generated
        M_DGP.SetT(T)  // no-null model in DGP
        M_DGP.InitEstimate()
        M_DGP.Estimate(t)
        confpeak = *pR * M_DGP.beta  // estimated coordinate of confidence peak
      }
      pR = &(J(cols(*pZexcl),kEx,0), I(cols(*pZexcl)))  // for ARubin test, picks out coefs on excluded exogenous variables
    }
    df = rows(*pR)
    
    if (granular) {
      euZAR = smatrix(df)
      pX = ARubin? &(*pXex, *pZexcl) : pXex
    }

    M_DGP.SetT(T0)  // (potentially) constrained model in DGP; T0 imposes constraints, pR tests hypotheses on results
    M_DGP.InitEstimate()

    if (ARubin) {
      k = el
      kappa = 0
      pM = pM_Repl
    } else {
      M_DGP.InitTestDenoms(T)
      pM = &M_DGP
    }
  }

  if (bootstrapt) {
    denom = smatrix(df,df)
    if (WREnonARubin==0 & robust) {
      if (B) Kd = smatrix(df)
      Kcd = smatrix(NErrClustCombs, df)
      pJcd = B? &smatrix(NErrClustCombs, df) : &Kcd  // if B = 0, Kcd will be multiplied by v, which is all 1's, and will constitute Jcd
    }
  }

  if (WREnonARubin & NClustVar)
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
  if (enumerate = (B & weighttype==0 & NBootClust*ln(2) < ln(B)+1e-6))  // generate full Rademacher set?
    MaxMatSize = .

  Ng = MaxMatSize == .? 1 : ceil((B+1) * max((rows(IDBootData), rows(IDBootAll), NBootClust)) * 8 / MaxMatSize / 1.0X+1E) // 1.0X+1E = giga(byte)
  if (Ng == 1) {
    makeWildWeights(B, 1)  // make all wild weights, once
    if (enumerate) B = cols(v) - 1  // replications reduced to 2^G
    WeightGrpStart = 1
  } else {
    seed = rseed()
    _reps = ceil((B+1) / Ng)
     WeightGrpStart = (0::Ng-1) * _reps :+ 1
    (WeightGrpStop  = (1::Ng  ) * _reps     )[Ng] = B+1
  }

  if (bootstrapt & (WREnonARubin | df>1 | MaxMatSize<.)) // unless nonWRE or df=1 or splitting weight matrix, code will create Dist element-by-element, so pre-allocate vector now
    pDist = &J(B+1, 1, .)
  if (Ng>1 | WREnonARubin | (null==0 & df<=2))
    numer = J(df, B+1, .)

  if (interpolable = B & WREnonARubin==0 & null & scoreBS==0 & Ng==1) {    
    dnumerdr = smatrix(q)
    if (interpolate_u = (robust | ML | GMM)==0) dudr = dnumerdr
    if (robust = robust) {
      ddenomdr = dJcddr = ssmatrix(q)
      ddenomdr2 = ssmatrix(q, q)
      for (i=q;i;i--) {
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

	if (initialized==0)
    Initialize()
  else if (null==0) {  // if not imposing null and we have returned, then df=1 or 2; we're plotting and only test stat, not distribution, changes with r
    if (WREnonARubin)
      numer[,1] = *pR * pM_Repl->beta - *pr
    else if (ARubin) {
      PrepARubin(*pr)
      numer[,1] = u_sd * pM->beta[|kEx+1\.|] // coefficients on excluded instruments in ARubin OLS
    } else
      numer[,1] = u_sd * (*pR * (ML? beta : pM->beta) - *pr) // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this.

    (*pDist)[1] = df==1? numer[1] / sqrt(statDenom) : numer[,1] ' invsym(statDenom) * numer[,1]
    return
  }

	if (Ng > 1) {
		rseed(seed)
    makeWildWeights(WeightGrpStop[1] - 1, 1)
  }

	makeInterpolables()  // make stuff that depends linearly on r, possibly by interpolating, for first weight group

  for (g=1; g<=Ng; g++) {  // do group 1 first because it includes col 1, which is all that might need updating in constructing CI in WCU
		if (g > 1)
			makeWildWeights(WeightGrpStop[g] - WeightGrpStart[g] + (g>1), g==1)

		if (WREnonARubin) {
			makeWREStats(g)
			if (bootstrapt==0)
				makeBootstrapcDenom(g)
		} else {
			if (bootstrapt)
				makeNonWREStats(g)
			else
				makeBootstrapcDenom(g)
		}
	}

	BFeast = (*pDist)[1]==.? 0 : colnonmissing(*pDist) - 1
	DistCDR = J(0,0,0)
	setdirty(0)
	initialized = 1
}

// compute bootstrap-c denominator from all bootstrap numerators
void boottestModel::makeBootstrapcDenom(real scalar g) {
	real colvector tmp
  if (g == 1) {
		tmp = numer[,1]
    statDenom = numer * numer' - tmp * tmp'
		numersum = rowsum(numer) - tmp
	} else {
		statDenom = statDenom + numer * numer'
		numersum = numersum + rowsum(numer)
	}
	if (g == Ng) {  // last weight group?
		statDenom = (statDenom - numersum * numersum' / B) / B
		pDist = &((sqrt? numer:/sqrt(statDenom) : colsum(numer :* invsym(statDenom) * numer))')
	}
}

// draw wild weight matrix of width _reps. If first=1, insert column of 1s at front
void boottestModel::makeWildWeights(real scalar _reps, real scalar first) {
	if (_reps) {
		if (enumerate)
			v = J(NBootClust,1,1), count_binary(NBootClust, -1-WREnonARubin, 1-WREnonARubin)  // complete Rademacher set
		else if (weighttype==3)
			v = rnormal(NBootClust, _reps+1, -WREnonARubin, 1)  // normal weights
		else if (weighttype==4)
			v = rgamma(NBootClust, _reps+1, 4, .5) :- (2 + WREnonARubin)  // Gamma weights
		else if (weighttype==2)
			if (WREnonARubin)
				v = sqrt(2 * ceil(runiform(NBootClust, _reps+first) * 3)) :* ((runiform(NBootClust, _reps+1):>=.5):-.5) :- 1  // Webb weights, minus 1 for convenience in WRE
			else {
				v = sqrt(    ceil(runiform(NBootClust, _reps+first) * 3)) :* ((runiform(NBootClust, _reps+1):>=.5):-.5)       // Webb weights, divided by sqrt(2)
				u_sd = 1.6a09e667f3bcdX-001 /*sqrt(.5)*/
			}
		else if (weighttype)
			if (WREnonARubin)
				v = ( rdiscrete(NBootClust, _reps+first, 1.727c9716ffb76X-001\1.1b06d1d200914X-002 /*.5+sqrt(.05)\.5-sqrt(.05)*/) :- 1.5 ) * 1.1e3779b97f4a8X+001 /*sqrt(5)*/ :- .5  // Mammen weights, minus 1 for convenience in WRE
			else {
				v = ( rdiscrete(NBootClust, _reps+first, 1.727c9716ffb76X-001\1.1b06d1d200914X-002 /*.5+sqrt(.05)\.5-sqrt(.05)*/) :- 1.5 ) :+ 1.c9f25c5bfedd9X-003 /*.5/sqrt(5)*/  // Mammen weights, divided by sqrt(5)
				u_sd = 1.c9f25c5bfedd9X-002 /*sqrt(.2)*/
			}
		else if (WREnonARubin) {
      v = runiform(NBootClust, _reps+first) :<  .5; v = (-2) * v  // Rademacher weights, minus 1 for convenience in WRE
    } else {
      v = runiform(NBootClust, _reps+first) :>= .5; v = v :- .5   // Rademacher weights, divided by 2
      u_sd = .5
    }

		if (first)
			v[,1] = J(NBootClust, 1, WREnonARubin? 0 : u_sd)  // keep original residuals in first entry to compute base model stat
	} else
		v = J(0,1,0)  // in places, cols(v) indicates number of B -- 1 for classical tests
}


void boottestModel::makeWREStats(real scalar g) {
	real scalar c, j, i
	real colvector _u, _beta, betaEnd, _v, numer_j
	real matrix Subscripts, Zi, VAR, tmp, XExi
	pointer (real matrix) scalar pustarZAR, pJ

	if (g == 1) {  // first/only weight group? initialize a couple of things
		_u = M_DGP.uddot + M_DGP.uddot2 * M_DGP.beta[|kEx+1\.|]

		if (NClustVar & NFE==0)  // prep for optimized computation for bootstrapping cluster when no FE
			for (i=NBootClust; i; i--) {
				Subscripts = (infoBootData)[i,]', (.\.)
				XExi = kEx? (*pXex)[|Subscripts|] : J(Subscripts[2,1]-Subscripts[1,1]+1,0,0)
				Zi = XExi , (*pZexcl)[|Subscripts|]  // inefficient?
				if (weights) Zi = Zi :* (*pwt)[|Subscripts|]
				XZi[i].M = cross(XExi, Zi) \ cross((*pXend)[|Subscripts|], Zi) \ cross((*pY)[|Subscripts|], Zi)
				eZi[i].M =                   cross(M_DGP.uddot2[|Subscripts|], Zi) \ cross(   _u[|Subscripts|], Zi)
			}
	}

	for (j=cols(v); j; j--) {  // WRE bootstrap
		_v = v[IDBootData,j]  // S'v^*j

		pM_Repl->InitEndog(&(*M_DGP.pY:+_u:*_v) , &(*pXend:+M_DGP.uddot2:*_v))
		pM_Repl->InitEstimate()
		pM_Repl->InitTestDenoms(T)  // prepare for replication regressions, null not imposed
		pM_Repl->Estimate(t)

		numer_j = null | g==1 & j==1? *pR * pM_Repl->beta - *pr : *pR * (pM_Repl->beta - M_DGP.beta0)

		if (bootstrapt) {
			if (robust) {  // Compute denominator for this WRE test stat
				denom = smatrix()
				if (NClustVar != 1 | NFE)  // collapse meat+sandwich  to all-Cluster-var intersections. If no collapsing needed, _panelsum() will still fold in any weights
					pustarZAR = &_panelsum(pM_Repl->ZAR, weights? *pwt :* pM_Repl->uddot : pM_Repl->uddot, *pinfoErrData)
				for (c=1; c<=NErrClustCombs; c++) {
					if (NClustVar != 1 & rows(Clust[c].order))
						pustarZAR = &((*pustarZAR)[Clust[c].order,])
					if (*pBootClust==Clust[c] & NClustVar & NFE==0) {  // optimized computation for bootstrapping Cluster when no FE
						VAR = pM_Repl->V * pM_Repl->AR; _beta = -pM_Repl->beta \ 1; betaEnd = _beta[|kEx+1\.|]
						pragma unset tmp
						for (i=1; i<=NBootClust; i++) {
							pJ = &((_beta'XZi[i].M + betaEnd'eZi[i].M * v[i,j]) * VAR) // R * V * Z_i'estar_i
							tmp = i==1? cross(*pJ,*pJ) : tmp + cross(*pJ,*pJ)
						}
            _clustAccum(denom.M, c, tmp)
					} else {
						pJ = &_panelsum(*pustarZAR, Clust[c].info)
            _clustAccum(denom.M, c, cross(*pJ,*pJ))
					}
				}
			} else
				denom.M = (*pR * pM_Repl->AR) * pM_Repl->eec

 			(*pDist)[j+WeightGrpStart[g]-1] = sqrt? numer_j/sqrt(denom.M) : cross(numer_j, invsym(denom.M) * numer_j)
		}
		numer[,j+WeightGrpStart[g]-1] = numer_j  // slight inefficiency: in usual bootstrap-t case, only need to save numerators in numer if getdist("numer") is coming because of svmat(numer)
	}

	if (g==1 & bootstrapt) statDenom = denom.M  // original-sample denominator
}


void boottestModel::PrepARubin(real colvector r) {
  pM_Repl->InitEndog(&(*pY - *pXend * r), NULL, &(*M_DGP.pZexclY - M_DGP.ZexclXend * r), &(*M_DGP.pXexY - M_DGP.XexXend * r))
  pM_Repl->InitEstimate()
  pM_Repl->Estimate(tARubin)
  pM_Repl->InitTestDenoms(TARubin)
}

// Construct stuff that depends linearly or quadratically on r, possibly by interpolation
void boottestModel::makeInterpolables() {
	real scalar h1, h2, d1, d2, c; real matrix tmp; real colvector Delta, newPole

  if (interpolable) {
    if (rows(anchor)==0) {  // first call? save current r as permanent anchor for interpolation
      _makeInterpolables(anchor = *pr)
      numer0 = numer
      if (interpolate_u) uddot0 = *puddot
      if (robust) Jcd0 = *pJcd
      return
    }

    if (rows(poles))  //  been here at least twice? interpolate unless current r stretches range > 2X in some dimension(s)
      newPole = abs(*pr - anchor) :> 2 * abs(poles)
    else {  // second call: from anchor make set of orthogonal poles, which equal anchor except in one dimension
      poles = *pr - anchor
      if (robust)  // grab quadratic denominator from *previous* (1st) evaluation
        denom0 = denom
      newPole = J(q,1,1)  // all poles new
    }

    if (any(newPole)) {  // prep interpolation
      for (h1=1;h1<=q;h1++)
        if (newPole[h1]) {
        	poles[h1] = (*pr)[h1] - anchor[h1]
        	(tmp = anchor)[h1] = (*pr)[h1]  // if q>1 this creates anchor points that are not graphed, an inefficiency. But simpler to make the deviations from 1st point orthogonal
          _makeInterpolables(tmp)  // calculate linear stuff at new anchor

          dnumerdr[h1].M = (numer - numer0) / poles[h1]
          if (interpolate_u)
            dudr[h1].M = (*puddot - uddot0) / poles[h1]
          if (robust)  // df > 1 for an ARubin test with >1 instruments. 
            for (d1=1;d1<=df;d1++) {
              for (c=1;c<=NErrClustCombs;c++) {
                dJcddr[h1].M[c,d1].M = ((*pJcd)[c,d1].M - Jcd0[c,d1].M) / poles[h1]
                for (d2=1;d2<=d1;d2++) {
                                tmp =       colsum(Jcd0 [c,d1].M :* dJcddr[h1].M[c,d2].M)
                  if (d1 != d2) tmp = tmp + colsum(Jcd0 [c,d2].M :* dJcddr[h1].M[c,d1].M)  // for diagonal items, faster to just double after the c loop
                  _clustAccum(ddenomdr[h1].M[d1,d2].M, c, tmp)
                }
              }
              ddenomdr[h1].M[d1,d1].M = ddenomdr[h1].M[d1,d1].M + ddenomdr[h1].M[d1,d1].M  // double diagonal terms
            }
        }
      if (robust)  // quadratic interaction terms
        for (h1=1;h1<=q;h1++)
          for (h2=h1;h2;h2--)
            if (newPole[h1] | newPole[h2])
              for (d1=df;d1;d1--)
                for (d2=d1;d2;d2--)
                  for (c=1;c<=NErrClustCombs;c++)
                    _clustAccum(ddenomdr2[h1,h2].M[d1,d2].M, c, colsum(dJcddr[h1].M[c,d1].M :* dJcddr[h2].M[c,d2].M))

      Delta = poles
      interpolating = 1

    } else {  // routine linear interpolation if the anchors not moved

      Delta = *pr - anchor
      numer = numer0 + dnumerdr.M * Delta[1]; if (q > 1) numer = numer + dnumerdr[2].M * Delta[2]
      if (interpolate_u) {
        puddot = &(uddot0 + dudr.M * Delta[1]); if (q > 1) puddot = &(*puddot + dudr[2].M * Delta[2])
      }
    }

    if (robust)  // even if an anchor was just moved, and linear components just compued from scratch, do the quadratic interpolation now, from the updated linear factors
      if (q==1)
        for (d1=df;d1;d1--)
          for (d2=d1;d2;d2--)
              denom[d1,d2].M = denom0[d1,d2].M + ddenomdr.M[d1,d2].M * Delta + ddenomdr2.M[d1,d2].M * (Delta * Delta)
      else  // q==2
        for (d1=df;d1;d1--)
          for (d2=d1;d2;d2--)
            denom[d1,d2].M = denom0[d1,d2].M + 
                         ddenomdr[1].M[d1,d2].M * Delta[1] + 
                         ddenomdr[2].M[d1,d2].M * Delta[2] + 
                         ddenomdr2[1,1].M[d1,d2].M * (Delta[1] * Delta[1]) + 
                         ddenomdr2[2,1].M[d1,d2].M * (Delta[1] * Delta[2]) + 
                         ddenomdr2[2,2].M[d1,d2].M * (Delta[2] * Delta[2])
  } else  // non-interpolable cases
    _makeInterpolables(*pr)
}

// Construct stuff that depends linearly or quadratically on r and doesn't depend on v. Does one bit relevant for WRE prep. No interpolation.
void boottestModel::_makeInterpolables(real colvector r) {
  pointer (real matrix) scalar puwt; real scalar d, i, c; real matrix tmp; pointer (real matrix) scalar pustarZAR, ptmp; real colvector r0

  if (ML==0) {
    if (ARubin)
      PrepARubin(r)
    else {
      r0 = null? r : J(0, 1, 0); if (REst) r0 =  *pr1 \ r0  // constant terms of model + null constraints
      t0 = rows(r0) ? L0invR0L0 * r0 : J(0,1,0)
      M_DGP.Estimate(t0)
    }
    puddot = &(pM->uddot)
    if (WREnonARubin) return
  }

  if (ML)
    uZAR = *pSc * (AR = *pA * *pR')
  else if (scoreBS | (robust & granular < NErrClustCombs))
    uZAR = *puddot :* pM->ZAR

  if (scoreBS)
    SuwtXA = B? (NClustVar? _panelsum(uZAR, *pwt, infoBootData) : (weights?       uZAR :* *pwt  :        uZAR) ) :
                                                                  (weights? cross(uZAR,   *pwt) : colsum(uZAR)')
  else {  // same calc as in score BS but broken apart to grab intermediate stuff, and assuming residuals defined
    puwt = weights? &(*puddot :* *pwt) : puddot
    ptmp = &_panelsum(*pXex, *puwt, infoBootData)
    if (ARubin)
      ptmp = &(*ptmp, _panelsum(*pZexcl, *puwt, infoBootData))
    SuwtXA = *pM->pA * (*ptmp)'
  }

  if (robust & GMM==0 & granular < NErrClustCombs) {
    pustarZAR = &_panelsum(uZAR, *pwt, *pinfoAllData)  // collapse data to all-boot & error-cluster-var intersections. If no collapsing interpolate_ed, _panelsum() will still fold in any weights
    if (B) {
      if (scoreBS)
        for (d=df;d;d--)
          Kd[d].M = J_ClustN_NBootClust  // inefficient, but not optimizing for the score bootstrap
      else
        for (d=df;d;d--) {
          tmp = pM->ZAR[,d]; if (weights) tmp = tmp :* *pwt
          if (ARubin)  // final term in (64) in paper, for c=intersection of all error clusters
            Kd[d].M = (_panelsum(*pXex, tmp, *pinfoErrData), _panelsum(*pZexcl, tmp, *pinfoErrData)) * SuwtXA
          else
            Kd[d].M =  _panelsum(*pXex, tmp, *pinfoErrData                                         ) * SuwtXA
        }

      if (NFE & FEboot==0)
        CT_WE = _panelsum(crosstab(wtFE :* *puddot), infoBootAll)'

      // subtract crosstab of E:*ZAR wrt bootstrapping cluster combo and all-cluster-var intersections
      if (*pBootClust == Clust[1])  // crosstab c,c* is square
        for (d=df;d;d--) {  // if bootstrapping on all-cluster-var intersections (including one-way clustering), the base crosstab is diagonal
          for (i=Clust.N;i;i--)
            Kd[d].M[i,i] = Kd[d].M[i,i] - (*pustarZAR)[i,d]
          if (scoreBS)
            Kd[d].M = Kd[d].M :+ ClustShare * (*pustarZAR)[,d]' // for score bootstrap, recenter; "+" because we subtracted *pustarZAR

        }
      else
        if (subcluster) // crosstab c,c* is wide
          for (d=df;d;d--) {
            for (i=Clust.N;i;i--) {
              tmp = infoErrAll[i,]'
              Kd.M[|(i\i), tmp|] = Kd.M[|(i\i), tmp|] - (*pustarZAR)[|tmp, (d\d)|]'
            }
            if (scoreBS)
	            Kd[d].M = Kd[d].M - ClustShare * colsum(Kd[d].M) // for score bootstrap, recenter
          }
        else // crosstab c,c* is tall
          for (d=df;d;d--) {
            for (i=NBootClust;i;i--) {
              tmp = pBootClust->info[i,]'
              Kd[d].M[|tmp, (i\i)|] = Kd[d].M[|tmp, (i\i)|] - (*pustarZAR)[|tmp, (d\d)|]
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
            Kcd[c,d].M = Kcd[c,d].M + _panelsum(pM->CT_ZAR[c,d].M, Clust[c].info) * CT_WE
        }
      }
    } else {  // B = 0. In this case, only 1st term of (64) is non-zero after multiplying by v* (= all 1's), and it is then a one-way sum by c
      if (scoreBS)
        pustarZAR = &(*pustarZAR :- ClustShare * colsum(*pustarZAR))  // recenter if OLS

      for (c=1; c<=NErrClustCombs; c++) {
        if (rows(Clust[c].order))
          pustarZAR = &((*pustarZAR)[Clust[c].order,])
        ptmp = &_panelsum(*pustarZAR, Clust[c].info)  // when c=1 (unless subcluster bootstrap), two args have same # of rows, &_panelsum() returns 1st arg by reference. Using & then prevents unnecessary cloning.
        for (d=df;d;d--)
          Kcd[c,d].M = (*ptmp)[,d]
      }
    }
  }

  makeNumerAndJ(1, r)  // compute J = kappa * v; if Ng > 1, then this is for 1st group; if interpolating, it is only group, and may be needed now to prep interpolation

  if      ( ARubin) numer[,1] = u_sd * pM->beta[|kEx+1\.|] // coefficients on excluded instruments in ARubin OLS
  else if (null==0) numer[,1] = u_sd * (*pR * (ML? beta : pM->beta) - r) // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this.
}

// compute stuff depending linearly on v, needed to prep for interpolation
void boottestModel::makeNumerAndJ(real scalar g, real scalar r) {  // called to *prepare* interpolation, or when g>1, in which case there is no interpolation
  real scalar c, d

  if (Ng==1)
    if (scoreBS)
      numer                                         = B? cross(SuwtXA, v) : SuwtXA * u_sd
    else if (robust==0 | granular)
      numer                                         = *pR * (betadev = SuwtXA * v)
    else
      numer                                         = (*pR * SuwtXA) * v
  else
    if (scoreBS)
      numer[|WeightGrpStart[g] \ WeightGrpStop[g]|] = B? cross(SuwtXA, v) : SuwtXA * u_sd
    else if (robust==0 | granular)
      numer[|WeightGrpStart[g] \ WeightGrpStop[g]|] = *pR * (betadev = SuwtXA * v)
    else
      numer[|WeightGrpStart[g] \ WeightGrpStop[g]|] = (*pR * SuwtXA) * v

  if (B & robust & GMM==0) {
    if (granular)  // prep optimized treatment when bootstrapping by many/small groups
      if (purerobust)
        ustar = *M_DGP.partialFE(&(*puddot :* v)) - *pX * betadev
      else {  // clusters small but not all singletons
        if (NFE) {
          ustar = *M_DGP.partialFE(&(*puddot :* v[IDBootData,]))
          for (d=df;d;d--)
            (*pJcd)[1,d].M = _panelsum(ustar,             pM->WZAR[d].M, *pinfoErrData)                                   - _panelsum(*pX, pM->WZAR[d].M, *pinfoErrData) * betadev
        } else
          for (d=df;d;d--)
            (*pJcd)[1,d].M = _panelsum(_panelsum(*puddot, pM->WZAR[d].M, *pinfoAllData) :* v[IDBootAll,], infoErrAll) - _panelsum(*pX, pM->WZAR[d].M, *pinfoErrData) * betadev
      }

    for (c=NErrClustCombs; c>granular; c--)
      for (d=df;d;d--)
        (*pJcd)[c,d].M = Kcd[c,d].M * v
  }
}

void boottestModel::makeNonWREStats(real scalar g) {
	real scalar i, c, j, l; real matrix ustar2, tmp; real colvector numer_l; pointer (real matrix) scalar pAR; real rowvector t1, t2, t12

  if (g > 1) makeNumerAndJ(g, *pr)

	if (robust & GMM==0) {
    if (interpolating==0) {  // these quadratic computation are not needed to *prepare* for interpolation but are superseded by interpolation once it is going
      if (purerobust)
        ustar2 = ustar :* ustar
      for (i=df;i;i--)
        for (j=i;j;j--) {
          if (c = purerobust) _clustAccum(denom[i,j].M, c, cross(pM->WZAR[i].M, pM->WZAR[j].M, ustar2))  // c=purerobust not a bug
          for (c++; c<=NErrClustCombs; c++)
            _clustAccum(denom[i,j].M, c, colsum((*pJcd)[c,i].M :* (*pJcd)[c,j].M))  // (60)
        }
    }

		if (df == 1) {
			if (Ng > 1)
				(*pDist)[|WeightGrpStart[g] \ WeightGrpStop[g]|] =   (numer[|WeightGrpStart[g] \ WeightGrpStop[g]|] :/ sqrt(denom.M))'
			else
          pDist                                          = &((numer                                         :/ sqrt(denom.M))')
			if (g==1)
				statDenom = denom.M[1]  // original-sample denominator
		} else if (df==2) {  // hand-code 2D numer'inv(denom)*numer
    	t1 = numer[1,]; t2 = numer[2,]; t12 = t1:*t2
			if (Ng > 1)
				(*pDist)[|WeightGrpStart[g] \ WeightGrpStop[g]|] =   ( (t1:*t1:*denom[2,2].M - (t12+t12):*denom[2,1].M + t2:*t2:*denom[1,1].M) :/ (denom[1,1].M:*denom[2,2].M - denom[2,1].M:*denom[2,1].M) )'
			else                                                                                                                         
          pDist                                          = &(( (t1:*t1:*denom[2,2].M - (t12+t12):*denom[2,1].M + t2:*t2:*denom[1,1].M) :/ (denom[1,1].M:*denom[2,2].M - denom[2,1].M:*denom[2,1].M) )')
			if (g==1)
				statDenom = denom[1,1].M[1], denom[2,1].M[1] \ denom[2,1].M[1], denom[2,2].M[1]  // original-sample denominator
    } else {  // build each replication's denominator from vectors that hold values for each position in denominator, all replications
			tmp = J(df,df,.)
			for (l=cols(v); l; l--) {
				for (i=df;i;i--)
					for (j=i;j;j--)
						tmp[i,j] = denom[i,j].M[l]
				_makesymmetric(tmp)
				numer_l = numer[,l]
				(*pDist)[l+WeightGrpStart[g]-1] = numer_l ' invsym(tmp) * numer_l  // in degenerate cases, cross() would turn cross(.,.) into 0
			}
			if (g==1)
				statDenom = tmp  // original-sample denominator
		}

	} else { // non-robust or GMM

		pAR = ML? &AR : &(pM->AR)
		if (df == 1) {  // optimize for one null constraint
			denom.M = *pR * *pAR

			if ((ML | GMM)==0) {
                     ustar = B? v :* *puddot : *puddot
				if (scoreBS) ustar = ustar :- (weights? cross(ClustShare, ustar) : colsum(ustar) * ClustShare)  // Center variance if interpolated
				        else ustar = ustar  - (*pXex, *pZexcl) * betadev  // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
				denom.M = denom.M :* (weights? cross(*pwt, ustar :* ustar) : colsum(ustar :* ustar))
			}
			if (Ng > 1)
				(*pDist)[|WeightGrpStart[g] \ WeightGrpStop[g]|] =   (numer :/ sqrt(denom.M))'
			else
				  pDist                                          = &((numer :/ sqrt(denom.M))')
			if (g==1)
				statDenom = denom.M[1]  // original-sample denominator
		} else {
			denom.M = *pR * *pAR

			if (ML | GMM) {
				for (l=cols(v); l; l--) {
					numer_l = numer[,l]
					(*pDist)[l+WeightGrpStart[g]-1] = cross(numer_l, invsym(denom.M), numer_l)
				}
				if (g==1)
					statDenom = denom.M  // original-sample denominator
			} else {
				for (l=cols(v); l; l--) {
					numer_l = numer[,l]
					(*pDist)[l+WeightGrpStart[g]-1] = cross(numer_l, invsym(denom.M) * numer_l)
                       ustar = B? v[,l] :* *puddot : *puddot
					if (scoreBS) ustar = ustar :- (weights? cross(*pwt, ustar) : colsum(ustar)) * ClustShare  // Center variance if interpolated
					        else ustar = ustar  - (*pXex, *pZexcl) * betadev[,l]  // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
					(*pDist)[l+WeightGrpStart[g]-1] = (*pDist)[l+WeightGrpStart[g]-1] / (tmp = cross(ustar, *pwt, ustar))
				}
				if (g==1)
					statDenom = denom.M * tmp  // original-sample denominator
			}
		}
	}
}


// like panelsetup() but can group on multiple columns, like sort(), and faster. But doesn't take minobs, maxobs arguments.
// Does take optional third argument, a matrix in which to store standardized ID variable, starting from 1
real matrix _panelsetup(real matrix X, real rowvector cols, | real colvector ID) {
	real matrix info; real scalar i, N; real scalar p; real rowvector tmp, id
	N = rows(X)
	info = J(N, 2, N); if (args()>2) ID = J(N, 1, 1)
	info[1,1] = p = 1
	id = X[1, cols]
	for (i=2; i<=N; i++) {
		if ((tmp=X[i,cols]) != id) {
			info[  p,2] = i - 1
			info[++p,1] = i
			id = tmp
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
	real matrix tmp
	if (N<=1) return (lo , hi)
	tmp = count_binary(N-1, lo, hi)
	return (J(1, cols(tmp), lo), J(1, cols(tmp), hi) \ tmp, tmp)
}

// cross-tab sum of a column vector w.r.t. intersection-of-error & bootstrap-clustering-vars and fixed-effect var
real matrix boottestModel::crosstab(real colvector v) {
	real matrix retval; real scalar i, j, tmp; real colvector _FEID, _v
	retval = J(rows(*pinfoAllData), NFE, 0)
	for (i=rows(*pinfoAllData);i;i--) {
		_FEID = panelsubmatrix(*pFEID, i, *pinfoAllData)
		_v    = panelsubmatrix(v     , i, *pinfoAllData)
		for (j=rows(_FEID);j;j--) {
			tmp = _FEID[j] 
			retval[i,tmp] = retval[i,tmp] + _v[j]
		}
	}
	return(retval)
}

// given a pre-configured boottest linear model with one-degree null imposed, compute distance from target p value of boostrapped one associated with given value of r
// used with optimize() to construct confidence intervals
// performs no error checking
real scalar boottestModel::r_to_p(real colvector r) {
	pr = &r
	setdirty(1, 1) // set dirty = 1, but leave initialized=0, which we want when only changing r
	return (getpadj())
}

real scalar boottestModel::search(real scalar alpha, real scalar p_lo, real scalar lo, real scalar p_hi, real scalar hi) {
	real scalar mid, _p
	mid = lo + (alpha-p_lo)/(p_hi-p_lo)*(hi-lo)  // interpolate linearly
//	mid = lo + ( reldif(p_hi-p_lo, 1/BFeast)<1e-6 & BFeast? (hi - lo)/2 : (alpha-p_lo)/(p_hi-p_lo)*(hi-lo) ) // interpolate linearly until bracketing a single bump-up; then switch to binary search
	if (reldif(lo,mid)<ptol | reldif(hi,mid)<ptol | (BFeast & abs(p_hi-p_lo)<(1+(ptype==1))/BFeast*1.000001))
		return (mid)
	if ( ((_p = r_to_p(mid)) < alpha) == (p_lo < alpha) )
		return(search(alpha, _p, mid, p_hi, hi))
	return(search(alpha, p_lo, lo, _p, mid))
}


// derive wild bootstrap-based CI, for case of linear model with one-degree null imposed.
void boottestModel::plot() {
	real scalar tmp, alpha, _quietly, c, d, i, j, halfwidth, p_lo, p_hi, p_confpeak; real colvector lo, hi; pointer (real colvector) scalar _pr

	_quietly = quietly; _pr = pr
	setquietly(1)
  alpha = 1 - level*.01
  _editmissing(gridpoints, 25)

  boottest()
  if (ARubin==0) {
    halfwidth = (-1.5 * invnormal(alpha/2)) * sqrt(diagonal(getV()))
    confpeak = getb() + *pr
  } else
    halfwidth = abs(confpeak) * invnormal(getpadj(1)/2) / invnormal(alpha/2)

  if (q==2) {
    lo = hi = J(2, 1, .)
    for(d=df;d;d--) {
      lo[d] = editmissing(gridmin[d], confpeak[d] - halfwidth[d])
      hi[d] = editmissing(gridmin[d], confpeak[d] + halfwidth[d])

      stata("_natscale " + strofreal(lo[d]) + " " + strofreal(hi[d]) + " 4")  // using Stata algorithm for setting graph bounds ensures good-looking contour plot
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
    if (alpha<=0) alpha = .05  // if level=100, no CI constructed, but we need a reasonable alpha to choose graphing bounds

    if (alpha > 0 & cols(v)-1 <= 1/alpha-1e6) {
      setquietly(_quietly)
      if (quietly==0) errprintf("\nError: need at least %g replications to resolve a %g%% two-sided confidence interval.\n", ceil(1/alpha), level)
      return
    }
    
    if (gridmin[1]==. | gridmax[1]==.) {
      if (B) {
        lo = editmissing(gridmin[1], confpeak - halfwidth) // initial guess based on classical distribution
        hi = editmissing(gridmax[1], confpeak + halfwidth)
      } else {
        tmp = sqrt(statDenom) * (small? invttail(df_r, alpha/2) : -invnormal(alpha/2))
        lo = editmissing(gridmin[1], confpeak - tmp)
        hi = editmissing(gridmax[1], confpeak + tmp)
        if (scoreBS & (null | willplot)==0) {  // if doing simple Wald test with no graph, we're done
          CI = lo, hi
          return
        }
      }
      
      if (abs(lo - *pr) > abs(hi - *pr)) {  // brute force way to ensure that first trial bound tested is the farther one from *pr, for better interpolation
        if (gridmin[1]==. & ptype!=2)  // unless upper-tailed p value, try at most 10 times to bracket confidence set by symmetrically widening
          for (i=10; i & -(p_lo=r_to_p(lo)) < -alpha; i--) {
            tmp = hi - lo
            lo = lo - tmp
            if (gridmax[1]==. & twotailed) hi = hi + tmp  // maintain rough symmetry unless user specified upper bound
          }
        if (gridmax[1]==. & ptype!=3)  // ditto for high side
          for (i=10; i & -(p_hi=r_to_p(hi)) < -alpha; i--) {
            tmp = hi - lo
            if (gridmin[1]==. & twotailed) lo = lo - tmp
            hi = hi + tmp
          }
      } else {
        if (gridmax[1]==. & ptype!=3)  // ditto for high side
          for (i=10; i & -(p_hi=r_to_p(hi)) < -alpha; i--) {
            tmp = hi - lo
            if (gridmin[1]==. & twotailed) lo = lo - tmp
            hi = hi + tmp
          }
        if (gridmin[1]==. & ptype!=2)  // unless upper-tailed p value, try at most 10 times to bracket confidence set by symmetrically widening
          for (i=10; i & -(p_lo=r_to_p(lo)) < -alpha; i--) {
            tmp = hi - lo
            lo = lo - tmp
            if (gridmax[1]==. & twotailed) hi = hi + tmp  // maintain rough symmetry unless user specified upper bound
          }
      }
    } else {  // both grid bounds pre-specified
      lo = gridmin[1]
      hi = gridmax[1]
    }

    plotX = rangen(lo, hi, gridpoints[1])
    plotY = J(rows(plotX), 1, .)
    plotY[1]  = p_lo; plotY[rows(plotX)]  = p_hi
    p_confpeak = WREnonARubin? . : (twotailed? 1 : .5)
    if (confpeak < lo) { // insert original point estimate into grid
      if (gridmin[1] == .) {
        plotX =   confpeak \ plotX
        plotY = p_confpeak \ plotY
        c = 1
      }
    } else if (confpeak > hi) {
      if (gridmax[1] == .) {
        plotX = plotX \   confpeak
        plotY = plotY \ p_confpeak
        c = gridpoints[1] + 1
      }
    } else {
      c = floor((confpeak - lo)/(hi - lo)*(gridpoints[1] - 1)) + 2
      plotX = plotX[|.\c-1|] \   confpeak \ plotX[|c\.|]
      plotY = plotY[|.\c-1|] \ p_confpeak \ plotY[|c\.|]
    }
  }  // end 1D plot

  i = 1
  if (_quietly==0 & WREnonARubin) {
    printf("{txt}")
    do {  // loop so that 1st 2 points are extremes, for efficient interpolation
      if (plotY[i] == .) plotY[i] = r_to_p(plotX[i,]')
      printf(".")
      if (mod(i-rows(plotX)-2, 50)) displayflush()
        else printf("\n")
    } while (1 < (i = mod(i-2,rows(plotX))+1))
    printf("\n")
  } else
    do {
      if (plotY[i] == .) plotY[i] = r_to_p(plotX[i,]')
    } while (1 < (i = mod(i-2,rows(plotX))+1))

  if (q==1 & level<100) {  // find CI bounds
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

	setquietly(_quietly); pr = _pr; dirty = 1  // restore backups
	notplotted = 0
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
	string scalar plotname, string scalar peakname, real scalar level, real scalar ptol, real scalar ML, real scalar LIML, real scalar Fuller, 
	real scalar kappa, real scalar ARubin, real scalar null, real scalar scoreBS, string scalar weighttype, string scalar ptype, string scalar statistic, string scalar madjtype, real scalar NumH0s,
	string scalar XExnames, string scalar XEndnames, real scalar hascons, string scalar Ynames, string scalar bname, string scalar Aname, string scalar Wname, 
	string scalar ZExclnames, string scalar samplename, string scalar scnames, real scalar robust, string scalar IDnames, real scalar NBootClustVar, real scalar NErrClust, 
	string scalar FEname, real scalar NFE, string scalar wtname, string scalar wttype, string scalar C1name, string scalar Cname, real scalar B, string scalar repsname, string scalar repsFeasname, 
	real scalar small, string scalar diststat, string scalar distname, string scalar gridmin, string scalar gridmax, string scalar gridpoints, real scalar MaxMatSize, real scalar quietly,
	string scalar b0name, string scalar V0name, string scalar vname, string scalar NBootClustname) {

	real matrix C1, R1, C, R, ZExcl, ID, FEID, sc, XEnd, XEx
	real colvector r1, wt, r, Y
	class boottestModel scalar M
	pragma unset ID; pragma unset wt; pragma unset XEnd; pragma unset XEx; pragma unset Y; pragma unset ZExcl; pragma unset sc

	C = st_matrix(Cname)
	R = C[|.,.\.,cols(C)-1|]
	r = C[,cols(C)]
	C1 = st_matrix(C1name)
	if (rows(C1)) { // restricted OLS?
		R1 = C1[|.,.\.,cols(C1)-1|]
		r1 = C1[,cols(C1)]
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
	M.setR1(R1, r1)
	M.setR(R, r)
	M.setnull(null)
	M.setsmall(small)
	M.setrobust(robust)
	M.setscoreBS(scoreBS)
	M.setweighttype(weighttype)
	M.setptype(ptype)
	M.setstattype(statistic)
	M.setwttype(wttype)
	M.setB(B)
	M.setLIML(LIML)
	M.setFuller(Fuller)
	M.setkappa(kappa)
	M.setARubin(ARubin)
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
	if (Aname != "") M.setA   (st_matrix(Aname) )
	if (Wname != "") M.setW   (st_matrix(Wname) )
	M.setwillplot(plotname != "" | ciname != "")  // indicates whether to optimize for many p evaluations
	if (plotname != "" | (level<100 & ciname != "")) {
		if (plotname != "") st_matrix(plotname, M.getplot())
		if (cols(M.peak)) st_matrix(peakname, M.getpeak())
		if (level<100 & ciname != "") st_matrix(ciname, M.getCI())
	}
	st_numscalar(statname, M.getstat())
	st_numscalar(pname   , M.getp   ())
	st_numscalar(repsname, M.getreps())
	st_numscalar(repsFeasname, M.getrepsFeas())
	st_numscalar(NBootClustname, M.getNBootClust())
	st_numscalar(padjname, M.getpadj())
	st_numscalar(dfname  , M.getdf  ())
	st_numscalar(dfrname , M.getdf_r())
	st_matrix(b0name, M.getb()')
	st_matrix(V0name, M.getV())
	if (distname != "") st_matrix(distname, M.getdist(diststat))
	if (vname != "" & B) st_matrix(vname, M.getv())
	M.M_DGP.parent = NULL // break loop in data structure topology to avoid Mata garbage collection leak
}

mata mlib create lboottest, dir("`c(sysdir_plus)'l") replace
mata mlib add lboottest *(), dir("`c(sysdir_plus)'l")
mata mlib index
end
