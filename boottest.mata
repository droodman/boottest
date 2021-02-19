/*
cd "C:\Users\drood\OneDrive\Documents\Work\Econometrics\Wild cluster"
use Levitt, clear
ivregress 2sls D.lPropertypop (DL.lpris_totpop = ibnL.stage#i(1/3)L.substage) D.(lincomepop unemp lpolicepop metrop black a*pop) i.year i.state, robust
boottest DL.lpris_totpop=-1/3, cluster(state year) bootcluster(year) ptype(equaltail) reps(199) noci seed(1231)
*/
*! boottest 3.0.2 18 December 2020
*! Copyright (C) 2015-21 David Roodman

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
mata set matalnum on

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

class boottestEstimator {  // class for analyitcal OLS, 2SLS, LIML, GMM estimation--everything but iterative ML
	real scalar LIML, y1y1, uuc, Fuller, ARubin, kappa, isDGP, k, kX1, y1pary1par
	real matrix XZ, XY2, XX, H_2SLS, invH, V, AR, XAR, U2ddot, ZY2, dbetadr, X2Y2, X1Y2, R1invR1R1, R0invR0R0, R1perp, T0, Rpar, RperpX, XinvXX, PXZ, YPXY, Zperp, invZperpZperp, Zex, Z, ZperpinvZperpZperp, YY, RRpar, invXX1, Pi, gamma, ZR1ex, ZR1, ZR1ZR1, X2ZR1, ZR1Y2, X1ZR1, ZZR1, RparX, RparY, R1invR1R1X, R1invR1R1Y, Rt1
	real colvector t0, u1ddot, u1dddot, beta, beta0, PXy1
	real rowvector y1Y2, Yendog, y1ZR1
	pointer(real rowvector) py1parY2
	pointer(real colvector) scalar py1, py1par, pX2y1par, pX1y1par, pZy1par, pXy1par
	pointer(real matrix) scalar pY2, pZZ, pZy1, pX2y1, pX1y1, pA, pW, pinvXX, pX1X1, pXX1, pH, pR, pT, pX1, pX2, pinvXX2
	pointer (class boottest scalar) scalar parent
	struct smatrix matrix CT_XAR, FillingT0
	struct smatrix rowvector WXAR, SPXYZperp, SMZperpYX

	private void new(), InitExog(), InitEndog(), InitTestDenoms(), SetR(), InitEstimate(), Estimate(), SetARubin()
	pointer(real matrix) scalar partialFE()
	private real matrix _select()/*, perp()
	private struct smatrix rowvector ParPerp()*/
}

class boottest {
	real scalar scoreBS, B, small, weighttype, null, dirty, initialized, ML, GMM, Nobs, _Nobs, k, kY2, kX1, l, sumwt, NClustVar, haswt, REst, multiplier, smallsample, quietly, FEboot, NErrClustCombs, ///
		sqrt, hascons, LIML, Fuller, kappa, IV, WRE, WREnonARubin, ptype, twotailed, df, df_r, ARubin, confpeak, willplot, notplotted, NumH0s, p, NBootClustVar, NErrClust, ///
		NFE, granular, purerobust, subcluster, Ncstar, BFeas, u_sd, level, ptol, MaxMatSize, Nw, enumerate, bootstrapt, q1, q, interpolable, interpolating, interpolate_u, robust
	real matrix AR, numer, v, ustar, TARubin, T0, L0invR0L0, CI, CT_WE, infoBootData, infoBootAll, infoErrAll, JNcapNcstar, statDenom, uXAR, SuwtXA, numer0, betadev, IDErr, U2parddot, deltadenom_b, JNcapk, ZR1ex, invR1R1par
	real colvector DistCDR, t, tARubin, plotX, plotY, t0, beta, wtFE, ClustShare, IDBootData, IDBootAll, WeightGrpStart, WeightGrpStop, gridmin, gridmax, gridpoints, numersum, uddot0, anchor, poles
	real rowvector peak
	string scalar wttype, madjtype, seed
	pointer (real matrix) scalar pX2, pR1, pR, pID, pFEID, pY2, pX1, pX, pinfoAllData, pinfoErrData, pIDAll
	pointer (real colvector) scalar pr1, pr, py1, pSc, pwt, pW, pA, puddot, pDist
	class boottestEstimator scalar DGP, Repl
	pointer (class boottestEstimator scalar) scalar pM
	struct structboottestClust colvector Clust
	pointer (struct structboottestClust scalar) scalar pBootClust
	struct smatrix matrix denom, Kcd, denom0, Jcd0, CTuX, SCTcapuXinvXX, FillingR2nonCTZperp, FillingR2nonCTX, Scstaruu
	struct smatrix rowvector Kd, euXAR, dudr, dnumerdr, IDCTCapcstar, infoCTCapcstar, ScstaruX, ScstaruXinvXX, ScstaruZperpinvZperpZperp, FillingT2, deltadenom, ZZg, Zyg, ScstaruY, ScstaruPZperp, ScstaruPX

  struct ssmatrix colvector ddenomdr, dJcddr
  struct ssmatrix matrix ddenomdr2
	pointer(struct smatrix matrix) scalar pJcd
	struct structFE rowvector FEs
	
	void new(), setsqrt(), setXEx(), setptype(), setdirty(), setXEnd(), setY(), setZExcl(), setwt(), setsc(), setML(), setLIML(), setARubin(),
		setFuller(), setkappa(), setquietly(), setbeta(), setA(), setW(), setsmall(), sethascons(), setscoreBS(), setB(), setnull(), setWald(), setRao(), setwttype(), setID(), setFEID(), setlevel(), setptol(), 
		setrobust(), setR1(), setR(), setwillplot(), setgrid(), setmadjust(), setweighttype(), setMaxMatSize(), setstattype(), close()
  private void makeNumerAndJ(), _clustAccum(), makeWREStats(), PrepARubin(), makeInterpolables(), _makeInterpolables(), makeNonWREStats(), makeBootstrapcDenom(), Init(), plot(), makeWildWeights(), boottest(), crosstabCapcstarMinus()
	real matrix getplot(), getCI(), getV(), getv()
	real scalar getp(), getpadj(), getstat(), getdf(), getdf_r(), getreps(), getrepsFeas(), getNBootClust()
	real rowvector getpeak()
	real colvector getdist(), getb()
	private real scalar r_to_p(), search()
  private real matrix count_binary(), crosstabAllFE(), ProductTermkappa()
	private pointer(real matrix) scalar Filling()
  private static real matrix combs()
  private static real colvector stableorder()
	private real vector _selectindex()
  static void _st_view()
}

void boottestEstimator::new() {
	ARubin = 0
  isDGP = 1  // by default, created object is for DGP rather than replication regression
}

void boottestEstimator::SetARubin(real scalar _AR) {
	if (ARubin = _AR) LIML = Fuller = kappa = 0
}

// do select() but handle case that both cases are scalar and second is 0 by interpreting second arg as rowvector and returning J(1,0,0)
// if v = 0 (so can't tell if row or col vector), returns J(1, 0, 0) 
real matrix boottestEstimator::_select(real matrix X, real rowvector v)
	return (rows(X)==1 & cols(X)==1 & v==0? J(1,0,0) : select(X,v))

function P(X) return(X*invsym(X'X)*X')
function M(X) return(I(rows(X))-P(X))


struct smatrix rowvector ParPerp(real matrix A) {
	real matrix vec; real rowvector val; struct smatrix rowvector retval
	pragma unset vec; pragma unset val
	symeigensystem(A'invsym(A*A')*A, vec, val)
	_edittozero(val, 10)
	retval = smatrix(2)
	retval[1].M = select(vec,  val)
	retval[2].M = select(vec, !val)
	return (retval)
}

real matrix perp(real matrix A) {
	real matrix vec; real rowvector val
	pragma unset vec; pragma unset val
	symeigensystem(A'invsym(A*A')*A, vec, val)
	_edittozero(val, 10)
	return (select(vec, !val))
}
real matrix par(real matrix A) {
	real matrix vec; real rowvector val
	pragma unset vec; pragma unset val
	symeigensystem(A'invsym(A*A')*A, vec, val)
	_edittozero(val, 10)
	return (select(vec, val))
}

// for DGP regression R has 0 rows, R1 for maintained constraints + null; for WRE replication regression R is for null, R1 for maintained constraints
// for non-WRE, R should have zero rows or not be passed if null not imposed on DGP
void boottestEstimator::SetR(real matrix R1, real matrix R) {
	real matrix RR1perp, vec; pointer (real matrix) scalar _pR; real rowvector val; pragma unset vec; pragma unset val

	if (parent->WREnonARubin==0) {
		// for OLS-based tests, this is called for DGP only. FWL is not exploited to shrink attack surface for null. 
		// No boottestEstimator() is created for replication regressions [change?] so construct R1perp, T0 at in this instance, for replication, DGP regressions respectively
		if (rows(R1)) {
			R1invR1R1 = R1 ' invsym(R1 * R1')
			symeigensystem(R1invR1R1 * R1, vec, val); _edittozero(val, 10)
			R1perp = _select(vec, !val)
			if (parent->ARubin) {
				R1perp = blockdiag(R1perp[|.,. \ parent->kX1 , parent->kX1-rows(R1)|] , I(parent->kY2)) // adapt R1perp,R1invR1R1 from XExog, XEndog to XExog, ZEXcl. Assumes no constraints link XExog and XEndog
				R1invR1R1 = R1invR1R1[|.,.\parent->kX1,.|] \ J(parent->kY2,1,0)
			}
		} else
			R1perp = J(0,0,0)

		if (rows(R1))
			_pR = rows(R)? &(R1 \ R) : &R1
		else if (rows(R))
			_pR = &R
		else {
			T0 = J(0,0,0)
			return
		}
		R0invR0R0 = invsym(*_pR * *_pR')
		if (all(diagonal(R0invR0R0))==0)
			_error(111, "Null hypothesis or model constraints are inconsistent or redundant.")
		R0invR0R0 = *_pR ' R0invR0R0
		symeigensystem(R0invR0R0 * *_pR, vec, val); _edittozero(val, 10)
		T0 = _select(vec, !val)

		return
	}

	this.pR = &R
	if (rows(R1)) {
		R1invR1R1 = invsym(R1 * R1')
		if (all(diagonal(R1invR1R1))==0)
			_error(111, "Null hypothesis or model constraints are inconsistent or redundant.")
		R1invR1R1 = R1 ' R1invR1R1
		symeigensystem(R1invR1R1 * R1, vec, val); _edittozero(val, 10)
		R1perp = _select(vec, !val)  // eigenvectors orthogonal to span of R1; foundation for parameterizing subspace compatible with constraints
	} else
		R1invR1R1 = J(parent->k,0,0)  // and R1perp = I

	// prepare to reduce regression via FWL
	RR1perp = R \ J(parent->kY2, parent->kX1, 0), I(parent->kY2)  // rows to prevent partialling out of endogenous regressors

	if (rows(R1))
		RR1perp = RR1perp * R1perp 
	symeigensystem(RR1perp ' invsym(RR1perp * RR1perp') * RR1perp, vec, val); _edittozero(val, 10)

	Rpar  = _select(vec,  val)
	RperpX = _select(vec, !val)
	if (rows(R1)) {  // fold model constraint factors into Rpar, RperpX
		Rpar  = R1perp * Rpar
		RperpX = R1perp * RperpX
	}
	RRpar = R * Rpar

	if (parent->kX1) {
		RperpX     = cols(RperpX   )? RperpX   [|.,.\parent->kX1,.|] : J(parent->kX1,0,0) // Zperp=Z*RperpX; though formally a multiplier on Z, it will only extract exogenous components, in X1, since all endogenous ones will be retained
		RparX      = cols(Rpar     )? Rpar     [|.,.\parent->kX1,.|] : J(parent->kX1,0,0)  // part of Rpar that refers to X1
		R1invR1R1X = cols(R1invR1R1)? R1invR1R1[|.,.\parent->kX1,.|] : J(parent->kX1,0,0)
	}
	RparY      = cols(Rpar     )? Rpar     [|parent->kX1+1,.\.,.|] : J(parent->kY2,0,0)  // part of Rpar that refers to Y2
	R1invR1R1Y = cols(R1invR1R1)? R1invR1R1[|parent->kX1+1,.\.,.|] : J(parent->kY2,0,0)
}

// stuff that can be done before r set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
void boottestEstimator::InitExog() {
	real matrix X2X1, S; real scalar i

	if (cols(*parent->pX2)) {  // GMM, 2SLS, LIML
		pX1 = partialFE(parent->pX1)  // <--drop modification of parent
		pX2 = partialFE(parent->pX2)

		Zex = *pX1 * RparX         // exogenous component of Zpar
		ZR1ex = *pX1 * R1invR1R1X  // exogenous component of ZR1

		Zperp = *parent->pX1 * RperpX
		invZperpZperp = invsym(cross(Zperp, *parent->pwt, Zperp))
		ZperpinvZperpZperp = Zperp * invZperpZperp

		pX1 = &(*parent->pX1 * perp(RperpX'))  // FWL-process X1
		pX1 = &(*pX1 - ZperpinvZperpZperp * cross(Zperp, *parent->pwt, *pX1))
		pX2 = &(*pX2 - ZperpinvZperpZperp * cross(Zperp, *parent->pwt, *pX2))

		kX1 = cols(*pX1)
		k   = cols(Rpar)

		pX1X1 = &cross(*pX1, *parent->pwt, *pX1)
		 X2X1 =  cross(*pX2, *parent->pwt, *pX1)
		pXX1 = &(*pX1X1 \ X2X1)
		if (parent->IV) {  // non-GMM
			pinvXX = &invsym((XX = *pXX1, (X2X1' \ cross(*pX2, *parent->pwt, *pX2))))
			if (kX1) {
				invXX1 = (*pinvXX)[|.,.\kX1,.|]
				pinvXX2 = &((*pinvXX)[|kX1+1,.\.,.|])
				XinvXX = *pX1 * invXX1 + *pX2 * *pinvXX2
			} else {
				invXX1 = J(0,parent->l,0)
				pinvXX2 = pinvXX
				XinvXX = *pX2 * (*pinvXX)
			}
		}

		if (isDGP==0)
			SMZperpYX = SPXYZperp = smatrix(k+1)
	} else {
		pX1 = parent->pX1	= partialFE(parent->pX1)  // <--drop modification of parent?
		pZZ = pX1X1 = &cross(*pX1, *parent->pwt, *pX1)
	}
}

// stuff that can be done before R & r set, but depend on endogenous variables, which are bootstrapped in WRE
void boottestEstimator::InitEndog(pointer (real colvector) scalar _py1, pointer (real matrix) scalar _pY2) {
	real matrix ZperpZ, ZperpY2, ZperpZR1; real scalar ind1, ind2, i; pointer (real colvector) scalar puwt

	py1 = partialFE(_py1)
	pY2 = partialFE(_pY2)
	if (ARubin)
		pX2y1 = &cross(*pX2, *parent->pwt, *py1)

	if (kappa) {  // **** check which computations needed for DGp vs Repl
		Z   = Zex   + *pY2 * RparY // Zpar
		ZR1 = ZR1ex + *pY2 * R1invR1R1Y
		
		ZperpY2  = cross(Zperp, *parent->pwt, *pY2 )
		ZperpZ   = cross(Zperp, *parent->pwt, Zex  ) + ZperpY2 * RparY
		ZperpZR1 = cross(Zperp, *parent->pwt, ZR1ex) + ZperpY2 * R1invR1R1Y

		Z   = Z   - ZperpinvZperpZperp * ZperpZ
		ZR1 = ZR1 - ZperpinvZperpZperp * ZperpZR1

		pY2 = &(*pY2 - ZperpinvZperpZperp * ZperpY2)
		py1 = &(*py1 - ZperpinvZperpZperp * cross(Zperp, *parent->pwt, *py1))

		X1Y2 = cross(*pX1, *parent->pwt, *pY2)
		X2Y2 = cross(*pX2, *parent->pwt, *pY2)
		XY2 = X1Y2 \ X2Y2
		y1Y2 =   cross(*py1, *parent->pwt, *pY2)
		pX2y1 = &cross(*pX2, *parent->pwt, *py1)
		pX1y1 = &cross(*pX1, *parent->pwt, *py1)
		y1y1 =   cross(*py1, *parent->pwt, *py1)
		pZy1 =  &cross(Z, *parent->pwt, *py1)		
		XZ =     cross(*pX1, *parent->pwt, Z) \ cross(*pX2, *parent->pwt, Z)
		ZY2 =    cross(Z, *parent->pwt, *pY2)
		pZZ =   &cross(Z, *parent->pwt, Z)

		if (cols(R1invR1R1)) {
			X2ZR1  = cross(*pX2, *parent->pwt, ZR1 )
			X1ZR1  = cross(*pX1, *parent->pwt, ZR1 )
			ZZR1   = cross(Z   , *parent->pwt, ZR1 )
			y1ZR1  = cross(*py1, *parent->pwt, ZR1 )
			ZR1ZR1 = cross(ZR1 , *parent->pwt, ZR1 )
			ZR1Y2  = cross(ZR1 , *parent->pwt, *pY2)
		} else {
			py1parY2   = &y1Y2 
			pX2y1par   = pX2y1
			pX1y1par   = pX1y1
			pZy1par    = pZy1 
			y1pary1par = y1y1 
			pXy1par = &(*pX1y1 \ *pX2y1)
			py1par = py1
		}

		// if GMM weight matrix not provided, prepare 2SLS one
		V =  (parent->IV? *pinvXX : *parent->pW) * XZ // in 2SLS case, estimator is (V' XZ)^-1 * (V'Xy1). Also apparently used in k-class and LIML robust VCV by Stata convention

		H_2SLS = V ' XZ  // Hessian

		if (isDGP==0 & parent->robust) {  // for replication regression, prepare for CRVE
			Yendog = 1, colsum(Rpar[|parent->kX1+1,.\.,.|] :!= 0)  // columns of Y = [y1 Zpar] that are endogenous (normally all)

			PXZ  = XinvXX * XZ

			FillingT0 = smatrix(k+1, k+1)  // fixed component of groupwise term in sandwich filling
			for (ind1=k; ind1; ind1--) {
				puwt = parent->haswt? *parent->pwt :* PXZ[,ind1] : &PXZ[,ind1]
				for (ind2=k; ind2; ind2--)
					FillingT0[ind1+1,ind2+1].M = _panelsum(Z[,ind2], *puwt, *parent->pinfoErrData)
			}

			for (i=k; i; i--) {  // precompute various clusterwise sums
				puwt = &(PXZ[,i]); if (parent->haswt) puwt = &(*puwt :* *parent->pwt)
				SPXYZperp[i+1].M = _panelsum(Zperp, *puwt, *parent->pinfoErrData)  // S_cap(P_(MZperpX) * Z :* Zperp)
				
				puwt = &(Z[,i]); if (parent->haswt) puwt = &(*puwt :* *parent->pwt)
				SMZperpYX[i+1].M = _panelsum(*pX1, *puwt, *parent->pinfoErrData), _panelsum(*pX2, *puwt, *parent->pinfoErrData)  // S_cap(M_Zperp[Z or y1] :* P_(MZperpX)])
			}
		}

	} else {  // OLS / ARubin
		pZy1 = pX1y1 = &cross(*pX1, *parent->pwt, *py1)
		if (ARubin) {
			pZy1 = &(*pZy1 \ *pX2y1)
			pZZ = &XX
		}
	}
}

// do most of estimation; for LIML r1 must be passed now in order to solve eigenvalue problem involving it
// For OLS, compute beta0 (beta when r=0) and dbetadr without knowing r1, for efficiency
// For WRE, should only be called once for the replication regressions, since for them r1 is the unchanging model constraints
void boottestEstimator::InitEstimate(real colvector r1) {
	real rowvector val; real colvector tmp; real matrix vec, XMZperpYpar, YparPXYpar; pointer (real matrix) scalar pbetadenom; real scalar i; pointer (real colvector) scalar puwt
	pragma unset vec; pragma unset val

	if (kappa) {
		if (cols(R1invR1R1)) {
			y1pary1par = y1y1 - 2 * y1ZR1 * r1 + r1 ' ZR1ZR1 * r1
			py1parY2 =  &(y1Y2  - r1 ' ZR1Y2)
			pX2y1par = &(*pX2y1 - X2ZR1 * r1)
			pX1y1par = &(*pX1y1 - X1ZR1 * r1)
			pZy1par  = &(*pZy1  - ZZR1  * r1)
			pXy1par = &(*pX1y1par \ *pX2y1par)

			py1par = &(*py1 - ZR1 * r1) 
		}

		if (isDGP==0 & parent->robust) {  // for replication regression, prepare for CRVE
			Rt1 = *pR * R1invR1R1 * r1

			PXy1 = XinvXX * *pXy1par

			puwt = parent->haswt? &(PXy1    :* *parent->pwt) : &PXy1
			for (i=k; i; i--)
				FillingT0[1,i+1].M = _panelsum(Z[,i], *puwt, *parent->pinfoErrData)
			SPXYZperp.M = _panelsum(Zperp, *puwt, *parent->pinfoErrData)  // S_cap(P_(MZperpX) * y1 :* Zperp)

			puwt = parent->haswt? &(*py1par :* *parent->pwt) : py1par
			for (i=k; i; i--)
				FillingT0[i+1,1].M = _panelsum(PXZ[,i], *puwt, *parent->pinfoErrData)
			FillingT0.M          = _panelsum(PXy1   , *puwt, *parent->pinfoErrData)
			SMZperpYX.M = _panelsum(*pX1,                                *puwt       , *parent->pinfoErrData), _panelsum(*pX2, *puwt, *parent->pinfoErrData)  // S_cap(M_Zperp*y1 :* P_(MZperpX)])
		}
		
		
		YY = y1pary1par, *pZy1par' \ *pZy1par, *pZZ

		tmp = *pinvXX * *pXy1par
		YPXY = tmp ' XZ 
		YPXY = tmp ' (*pXy1par) , YPXY \ YPXY' , XZ ' (*pinvXX) * XZ 

		if (LIML) {
			XMZperpYpar = *pXy1par, XZ
			YparPXYpar = XMZperpYpar ' (*pinvXX) * XMZperpYpar

			V = *pinvXX * XZ
			H_2SLS = V ' XZ  // Hessian
			eigensystemselecti(invsym(YY) * YparPXYpar, rows(YY)\rows(YY), vec, val)  // eigensystemselecti(invsym(YY) * YPXY, rows(YY)\rows(YY), ... gives 1 - the eigenvalue, but can cause eigensystem() to return all missing
			kappa = 1/(1 - Re(val)) // sometimes a tiny imaginary component sneaks into val
			if (Fuller) kappa = kappa - 1 / (parent->_Nobs - parent->l)
		}
	}

  pH = kappa? (kappa==1? &H_2SLS : &((1-kappa)* *pZZ + kappa*H_2SLS)) : pZZ
	pbetadenom = &(invH = invsym(*pH))

  if (kappa)
		beta0 = *pbetadenom *  (kappa==1? V ' (*pXy1par) : (kappa * V ' (*pXy1par) + (1-kappa) * *pZy1par))  // concentrated estimator, for computing residuals and, via left-multiplication by R*Rpar, generating null LHS
	else {  // OLS / ARubin
		beta0 = *pbetadenom * *pZy1 * R0invR0R0
		dbetadr = *pbetadenom * *pZZ * R0invR0R0 - R0invR0R0
	}
}

// stuff that depends on r and endogenous variables: compute beta and residuals
void boottestEstimator::Estimate(real colvector r1) {
	real colvector negXuinvuu; real matrix Xu; real scalar uu

	beta = rows(r1) & parent->WREnonARubin==0? beta0 - dbetadr * r1 : beta0  // change so WRE call to Estimate() doesn't pass r1, so "rows(r1)" suffices

	if (isDGP | parent->bootstrapt | parent->WREnonARubin==0 | parent->robust==0) {  // don't need residuals in replication regressions in bootstrap-c on WRE/non-ARubin
		if (ARubin) {
			u1ddot = *py1 - *pX2 * beta[|cols(*pX1)+1\.|]
			if (cols(*pX1))
				u1ddot =  u1ddot - *pX1 * beta[|.\cols(*pX1)|]
		} else
  			u1ddot = *py1par - Z * beta
		if ((parent->robust | parent->scoreBS)==0 | (isDGP & LIML))  // useful in non-robust, residual-based bootstrap, and in computing u1ddot^2 in LIML (just below)
			uu = (1 \ -beta) ' YY * (1 \ -beta)
		if ((parent->robust | parent->scoreBS)==0 & parent->bootstrapt==0)
			uuc = parent->hascons? uu : uu - (parent->haswt? cross(u1ddot, *parent->pwt) : sum(u1ddot))^2 / parent->_Nobs  // sum of squares after centering, N * Var
	}

	if (isDGP & LIML) {  // after IV/GMM DGP regression, compute Y2 residuals by regressing Y2 on X while controlling for y1 residuals, done through FWL

		Xu = *pXy1par - XZ * beta
		negXuinvuu = Xu / -uu
		Pi = invsym(XX + negXuinvuu * Xu') * (negXuinvuu * (*py1parY2 - beta ' ZY2) + XY2)
		U2ddot = *pY2 - *pX2 * Pi[|kX1+1,.\.,.|]; if (kX1) U2ddot = U2ddot - *pX1 * Pi[|.,.\kX1,.|]

		gamma = R1invR1R1Y * r1; if (cols(Rpar)) gamma = gamma + RparY * beta
		u1dddot = u1ddot + U2ddot * gamma
	}
}


// stuff that doesn't depend on r, for test stat denominators in replication regressions
// but, confusingly, since the non-AR OLS code never creates an object for replication regresssions, in that case this is called on the DGP regression object
void boottestEstimator::InitTestDenoms() {
	real matrix VAR; real scalar d, c; struct smatrix rowvector _CT_XAR; pointer (real matrix) scalar pWXAR
	pA = &(rows(invH)? invH : invsym(*pH))
	AR = *pA * Rpar ' (*parent->pR)'

	if (parent->scoreBS | parent->robust) {
		if (kappa) {
			VAR = V * AR
			XAR = *pX2 * VAR[|kX1+1,.\.,.|]; if (kX1) XAR = XAR + *pX1 * VAR[|.,.\kX1,.|]
		} else if (ARubin) {
			XAR = *pX2 * AR[|kX1+1,.\.,.|]
			if (cols(*pX1))
				XAR = XAR + *pX1 * AR[|.,.\kX1,.|]
		} else
			XAR = *pX1 * AR

			if (parent->bootstrapt == 0) return

			if (parent->granular) {
				pWXAR = &(parent->haswt? XAR :* *parent->pwt : XAR)
				WXAR = smatrix(parent->df)
				for (d=parent->df;d;d--)
					WXAR[d].M = (*pWXAR)[,d]
			}
			
      if (parent->NFE & parent->robust & (parent->FEboot | parent->scoreBS)==0 & parent->granular < parent->NErrClustCombs) {
				if (pWXAR==NULL) pWXAR = &(parent->haswt? XAR :* *parent->pwt : XAR)
				CT_XAR = smatrix(parent->NErrClustCombs, parent->df)
				for (d=parent->df;d;d--)
					CT_XAR[1,d].M = parent->crosstabAllFE((*pWXAR)[,d])
				if (parent->NClustVar > 1) {
					_CT_XAR = CT_XAR[1,]
					for (c=2; c<=parent->NErrClustCombs; c++)
						for (d=parent->df;d;d--) {
							if (rows(parent->Clust[c].order))
								_CT_XAR[d].M = _CT_XAR[d].M[parent->Clust[c].order,]
							CT_XAR[c,d].M = _panelsum(_CT_XAR[d].M, parent->Clust[c].info)
						}
				}
			}
	}
}


// partial fixed effects out of a data matrix
pointer(real matrix) scalar boottestEstimator::partialFE(pointer(real matrix) scalar pIn) {
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

void boottest::new() {
	ARubin = LIML = Fuller = WRE = small = scoreBS = weighttype = ML = initialized = quietly = sqrt = hascons = IV = ptype = robust = NFE = FEboot = granular = NErrClustCombs = subcluster = B = BFeas = interpolating = 0
	twotailed = null = dirty = willplot = u_sd = bootstrapt = notplotted = 1
	level = 95
  ptol = 1e-6
	confpeak = MaxMatSize = .
	pY2 = pX1 = pX2 = py1 = pSc = pID = pFEID = pR1 = pR = pwt = &J(0,0,0)
	pr1 = pr = &J(0,1,0)
	IDBootData = .
}

// import to call this when done: break loops in data structure topology to enable Mata garbage collection
void boottest::close() {
	DGP.parent = NULL
	if (WRE)
		Repl.parent = NULL
}

void boottest::setdirty(real scalar _dirty, | real scalar noinitialize) {
	dirty = _dirty
	if (_dirty & noinitialize!=1)
		initialized = 0
}

void boottest::setsqrt(real scalar _sqrt) {
	if (_sqrt < sqrt) {
		if (dirty==0) {
    	pDist = &(*pDist :* *pDist)
      multiplier = multiplier * multiplier
    }
	} else
		setdirty(1)
	sqrt = _sqrt
}

void boottest::setptype(string scalar ptype) {
	real scalar p
	p = cross( (strtrim(strlower(ptype)) :== ("symmetric"\"equaltail"\"lower"\"upper")), 1::4 ) - 1
	if (p<0) 
		_error(198, `"p-value type must be "symmetric", "equaltail", "lower", or "upper"."')
	this.ptype = p
	this.twotailed = p<=1
}
void boottest::setstattype(string scalar stattype) {
	real scalar p
	p = cross( (strtrim(strlower(stattype)) :== ("c"\"t")), 1::2 ) - 1
	if (p<0) 
		_error(198, `"statistic type must be "t" or "c"."')
	this.bootstrapt = p
	setdirty(1)
}	
void boottest::setXEnd    (real matrix X ) {
	this.pY2  = &X; setdirty(1)
}
void boottest::setXEx     (real matrix X ) {
	this.pX1  = &X; setdirty(1)
}
void boottest::setY       (real matrix Y ) {
	this.py1  = &Y; setdirty(1)
}
void boottest::setZExcl   (real matrix Z ) {
	this.pX2  = &Z; setdirty(1)
}
void boottest::setwt      (real matrix wt) {
	this.pwt  = &wt; setdirty(1)
}
void boottest::setsc(real matrix Sc) {
	this.pSc  = &Sc
	setdirty(1)
}
void boottest::setML(real scalar ML) {
	this.ML  = ML; setdirty(1)
	if (ML) setscoreBS(1)
}
void boottest::setLIML(real scalar LIML) {
	this.LIML = LIML; setdirty(1)
}
void boottest::setARubin(real scalar ARubin) {
	this.ARubin = ARubin; setdirty(1)
}
void boottest::setFuller    (real scalar Fuller) {
	this.Fuller = Fuller; setdirty(1)
}
void boottest::setkappa(real scalar kappa) {  // kappa as in k-class
	this.kappa = kappa; setdirty(1)
}
void boottest::setquietly(real scalar quietly )
	this.quietly = quietly
void boottest::setbeta(real colvector beta) {
	this.beta = beta; setdirty(1)
}
void boottest::setA(real matrix V) {
	this.pA = &V; setdirty(1)
}
void boottest::setW(real matrix W) {
	this.pW = &W; setdirty(1)
}
void boottest::setsmall(real scalar small) {
	this.small = small; setdirty(1)
}
void boottest::sethascons(real scalar hascons) {
	this.hascons = hascons; setdirty(1)
}
void boottest::setscoreBS (real scalar scoreBS) {
	this.scoreBS = scoreBS; setdirty(1)
}
void boottest::setB(real scalar B) {
	this.B = B
	if (B==0)
		setscoreBS(1)
	setdirty(1)
}
void boottest::setnull    (real scalar null) {
	this.null = null; setdirty(1)
}
void boottest::setWald() { // set-up for classical Wald test
	this.scoreBS = 1; this.B = 0; this.null = 0; setdirty(1)
}
void boottest::setRao() { // set-up for classical Rao test
	this.scoreBS = 1; this.B = 0; this.null = 1; setdirty(1)
}
void boottest::setwttype  (string scalar wttype) {
	this.wttype = wttype; setdirty(1)
}
void boottest::setID      (real matrix ID, | real scalar NBootClustVar, real scalar NErrClust) {
	this.pID = &ID; this.NBootClustVar = editmissing(NBootClustVar,1); this.NErrClust=editmissing(NErrClust,1); setdirty(1)
	if (cols(ID)) this.robust = 1
}
void boottest::setFEID(real matrix ID, real scalar NFE) {
	this.pFEID = &ID; this.NFE = NFE; setdirty(1)
}
void boottest::setlevel(real scalar level)
	this.level = level
void boottest::setptol(real scalar ptol)
	this.ptol = ptol
void boottest::setrobust(real scalar robust) {
	this.robust = robust
	if (robust==0) setID(J(0,0,0), 1, 1)
	setdirty(1)
}
void boottest::setR1(real matrix R1, real colvector r1) {
	this.pR1 = &R1; 	this.pr1  = &r1; q1 = rows(R1); setdirty(1)
}
void boottest::setR(real matrix R, real colvector r) {
	this.pR = &R; this.pr = &r; q = rows(R); setdirty(1)  // q can differ from df in ARubin test
}
void boottest::setwillplot(real scalar willplot) {
	this.willplot = willplot
}
void boottest::setgrid(real rowvector gridmin, real rowvector gridmax, real rowvector gridpoints) {
	this.gridmin = gridmin; this.gridmax = gridmax; this.gridpoints = gridpoints
}
void boottest::setmadjust(string scalar madjtype, real scalar NumH0s) {
	this.madjtype = strlower(madjtype)
	this.NumH0s = NumH0s
	if (this.madjtype != "bonferroni" & this.madjtype != "sidak" & this.madjtype != "")
		_error(198, `"Multiple-hypothesis adjustment type must be "Bonferroni" or "Sidak"."')
}
void boottest::setweighttype(string scalar weighttype) {
	weighttype = strlower(weighttype)
	if (.==(this.weighttype = weighttype=="rademacher" ? 0 : (weighttype=="mammen" ? 1 : (weighttype=="webb" ? 2 : (weighttype=="normal" ? 3 : (weighttype=="gamma" ? 4 : .))))))
		_error(198, `"Wild type must be "Rademacher", "Mammen", "Webb", "Normal", or "Gamma"."')
	setdirty(1)
}
void boottest::setMaxMatSize(real scalar MaxMatSize) {
	this.MaxMatSize = MaxMatSize; setdirty(1)
}

real colvector boottest::getdist(| string scalar diststat) {
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

// get p value. Robust to missing bootstrapped values interpreted as +infinity.
real scalar boottest::getp(|real scalar classical) {
	real scalar tmp
	if (dirty) boottest()
	tmp = (*pDist)[1]
	if (tmp == .) return (.)
	if (B & classical==.)
		if (sqrt & ptype != 3) {
			if (ptype==0)
				p = colsum(-abs(tmp) :> -abs(*pDist)) / BFeas  // symmetric p value; do so as not to count missing entries in *pDist
			else if (ptype==1)  // equal-tail p value
				p = 2 * min((colsum(tmp :> *pDist) , colsum(-tmp:>- *pDist))) / BFeas
			else
				p = colsum( tmp :>   *pDist) / BFeas  // lower-tailed p value
		} else
				p = colsum(-tmp :> - *pDist) / BFeas  // upper-tailed p value or p value based on squared stats
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
real colvector boottest::getb() {
	if (dirty) boottest()
	return(u_sd == 1? numer[,1] : numer[,1] / u_sd)
}

// denominator for full-sample test stat
real matrix boottest::getV() {
	if (dirty) boottest()
	return (statDenom / ((u_sd == 1? smallsample : u_sd * u_sd * smallsample)  * (sqrt? multiplier*multiplier : multiplier) * df))
}

// wild weights
real matrix boottest::getv()
	return(u_sd==1? v[|.,2\.,.|] : v[|.,2\.,.|] / u_sd)

// Return number of bootstrap replications with feasible results
// Returns 0 if getp() not yet accessed, or doing non-bootstrapping tests
real scalar boottest::getrepsFeas()
	return (BFeas)

real scalar boottest::getNBootClust()
	return (Ncstar)

// return number of replications, possibly reduced to 2^G
real scalar boottest::getreps()
	return (B)

real scalar boottest::getpadj(|real scalar classical) {
	(void) getp(classical)
	if (madjtype=="bonferroni") return(min((1, NumH0s*p)))
	if (madjtype=="sidak"     ) return(1 - (1 - p)^NumH0s)
	return(p)
}

real scalar boottest::getstat() {
	if (dirty) boottest()
	return(multiplier * (*pDist)[1])
}
real scalar boottest::getdf() {
	if (dirty) boottest()
	return(df)
}
real scalar boottest::getdf_r() {
	if (dirty) boottest()
	return(df_r)
}
real matrix boottest::getplot() {
	if (notplotted) plot()
	return((plotX,plotY))
}
real rowvector boottest::getpeak() {  // x and y values of confidence curve peak (at least in OLS & ARubin)
	if (notplotted) plot()
	return(peak)
}
real matrix boottest::getCI() {
	if (notplotted) plot()
	return(CI)
}

void boottest::_st_view(real matrix V, real scalar i, string rowvector j, string scalar selectvar) {
	if (favorspeed() | 1)
		V = length(tokens(j))? st_data(i, j, selectvar) : st_data(i, J(1,0,0), selectvar)
	else
		st_view(V, i, j, selectvar)
}

// helper for summing over clusterings while factoring in clustering-specific parity and small-sample adjustments
// replace X with Y if c=1; otherwise add it
void boottest::_clustAccum(real matrix X, real scalar c, real matrix Y)
  X = c == 1?
          (Clust.   even?
            (Clust.   multiplier != 1?   Clust.   multiplier  * Y :  Y) :
            (Clust.   multiplier != 1? (-Clust.   multiplier) * Y : -Y)) :
      X + (Clust[c].even?                
            (Clust[c].multiplier != 1?   Clust[c].multiplier  * Y :  Y) :
            (Clust[c].multiplier != 1? (-Clust[c].multiplier) * Y : -Y))

						
void boottest::Init() {  // for efficiency when varying r repeatedly to make CI, do stuff once that doesn't depend on r
  real colvector sortID, o, _FEID, IDAllData, IDErrData
	real rowvector ClustCols
	real matrix Combs, tmp
	real scalar i, j, c, minN, sumN, _B, i_FE
	pragma unset IDAllData; pragma unset IDErrData

  Nobs = rows(*pX1)
  NClustVar = cols(*pID)
  k  = cols(*pR)
  kX1 = cols(*pX1)
  if (cols(*pX2)==0) pX2 = &J(Nobs,0,0)
  if ((kY2 = cols(*pY2))==0) pY2 = &J(Nobs,0,0)
	if (LIML & cols(*pX2)==kY2) {  // exactly identified LIML = 2SLS
		kappa = 1
		LIML = 0
	}
  REst = rows(*pR1)  // base model contains restrictions?
	if (REst == 0) pR1 = &J(0,k,0)
  if (pX2 != NULL) l = cols(*pX2) + kX1  // l = # of exogenous variables
  if (kappa==.) kappa = cols(*pX2)>0  // if kappa in kappa-class estimation not specified, it's 0 or 1 for OLS or 2SLS
  GMM = pW != NULL
  IV = kappa & GMM==0
  WRE = (kappa & scoreBS==0) | ARubin
  WREnonARubin = WRE & ARubin==0

  if (haswt = rows(*pwt)>1)
    sumwt = sum(*pwt)
  else
    pwt = &(sumwt = 1)
  _Nobs = haswt & wttype=="fweight"? sumwt : Nobs

  if (WREnonARubin)
    if (NClustVar)
      infoBootData = _panelsetup(*pID, 1..NBootClustVar, IDBootData)
    else
      pinfoErrData = &(infoBootData = J(Nobs,0,0))
  else if (NClustVar)
    if (NClustVar > NBootClustVar)  // bootstrap Cluster grouping defs rel to original data
      infoBootData = _panelsetup(*pID, 1..NBootClustVar)
    else
      infoBootData = _panelsetup(*pID, 1..NClustVar)
  else
    pinfoErrData = pinfoAllData = &(infoBootData = J(Nobs,0,0))  // causes no collapsing of data in _panelsum() calls, only multiplying by weights if any
  Ncstar = rows(infoBootData)

  if (bootstrapt) {
    if (NClustVar) {
      minN = .; sumN = 0

      Combs = combs(NErrClust)  // represent all error clustering combinations. First is intersection of all error clustering vars
      Clust = structboottestClust(rows(Combs)-1)  // leave out no-cluster combination
      NErrClustCombs = length(Clust)
      subcluster = NClustVar - NErrClust

      pinfoAllData = NClustVar > NBootClustVar? &_panelsetup(*pID,            1..NClustVar) : &infoBootData  // info for grouping by intersections of all bootstrap & clustering vars wrt data; used to speed crosstab EXAR wrt bootstrapping cluster & intersection of all error clusterings
      pinfoErrData = NClustVar > NErrClust    ? &_panelsetup(*pID, subcluster+1..NClustVar) : pinfoAllData  // info for intersections of error clustering wrt data
       IDErr = rows(*pinfoErrData)==Nobs? *pID :   (*pID)[(*pinfoErrData)[,1],]   // version of ID matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
      pIDAll = rows(*pinfoAllData)==Nobs?  pID : &((*pID)[(*pinfoAllData)[,1],])  // version of ID matrix with one row for each all-bootstrap & error cluster-var intersection instead of 1 row for each obs

      if (subcluster) {  // for subcluster bootstrap, bootstrapping cluster is not among error clustering combinations
        pBootClust = &(structboottestClust())
        pBootClust->info = _panelsetup(*pIDAll, 1..NBootClustVar)  // bootstrapping cluster info w.r.t. all-bootstrap & error-cluster intersections
        Ncstar = rows(pBootClust->info)
      } else
        pBootClust = &(Clust[2^(NClustVar - NBootClustVar)])  // *** inefficient because copies large objections? *** location of bootstrap clustering within list of cluster combinations

			for (c=1; c<=NErrClustCombs; c++) {  // for each error clustering combination
        ClustCols = subcluster :+ _selectindex(Combs[c,])
        Clust[c].even = mod(cols(ClustCols),2)

        if (c == 1)
          if (subcluster) {
            IDErr = IDErr[ Clust.order = stableorder(IDErr, ClustCols), ]
            Clust.info       = _panelsetup(IDErr, ClustCols)
          } else
            Clust.info       = J(rows(*pinfoAllData),0,0)  // causes no collapsing of data in _panelsum() calls
        else {
          if (any( Combs[|c, min(_selectindex(Combs[c,] :!= Combs[c-1,])) \ c,.|])) // if this sort ordering same as last to some point and missing thereafter, no need to re-sort
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
        ClustShare = haswt? _panelsum(*pwt, *pinfoErrData)/sumwt : ((*pinfoErrData)[,2]-(*pinfoErrData)[,1]:+ 1)/Nobs // share of observations by group 
if (WREnonARubin) {
(void) _panelsetup(*pID,            1..NClustVar, IDAllData)
(void) _panelsetup(*pID, subcluster+1..NClustVar, IDErrData)
				IDCTCapcstar = infoCTCapcstar = smatrix(Ncstar)
				for (i=Ncstar;i;i--) {
					t = IDAllData[|infoBootData[i,]'|]                        // ID numbers w.r.t. intersection of all bootstrap/error clusterings contained in bootstrap cluster i
					infoCTCapcstar[i].M = (*pinfoAllData)[t[1]::t[rows(t)],]  // for those ID's, panel info for the all-bootstrap/error-clusterings data row groupings
					IDCTCapcstar[i].M = IDErrData[infoCTCapcstar[i].M[,1]]    // ID numbers of those groupings w.r.t. the all-error-clusterings grouping
				}
}

    } else {  // if no clustering, cast "robust" as clustering by observation
      pBootClust = &(Clust = structboottestClust())
      Clust.multiplier = small? _Nobs / (_Nobs - 1) : 1
      Clust.even = 1
      sumN = Clust.N = Nobs
      Clust.info = J(Nobs, 0, 0)  // signals _panelsum not to aggregate
      NErrClustCombs = 1
      if (scoreBS | WREnonARubin==0)
        ClustShare = haswt? *pwt/sumwt : 1/_Nobs
    }

    purerobust = NClustVar & (scoreBS | subcluster)==0 & Ncstar==Nobs  // do we ever error-cluster *and* bootstrap-cluster by individual?
    granular   = NClustVar & scoreBS==0 & (purerobust | (Clust.N+Ncstar)*k*B + (Clust.N-Ncstar)*B + k*B < Clust.N*k*k + Nobs*k + Clust.N * Ncstar * k + Clust.N * Ncstar)

    if (robust & purerobust==0) {
      if (subcluster | granular)
        infoErrAll = _panelsetup(*pIDAll, subcluster+1..NClustVar)  // info for error clusters wrt data collapsed to intersections of all bootstrapping & error clusters; Fused to speed crosstab UXAR wrt bootstrapping cluster & intersection of all error clusterings
      if ((scoreBS & B) | (WREnonARubin & (NClustVar != NBootClustVar | subcluster)))
        JNcapNcstar = J(Clust.N, Ncstar, 0)
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
        if (haswt) {
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
    if (haswt) {
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
    DGP.parent = &this
		DGP.SetR(null? *pR1 \ *pR : *pR1, J(0,k,0))
		DGP.InitExog()

    if (WRE) {
			Repl.parent = &this
      Repl.isDGP = 0
      DGP.LIML = 1; DGP.Fuller = 0; DGP.kappa = 1
      if (ARubin==0) { 
        Repl.LIML = this.LIML; Repl.Fuller = this.Fuller; Repl.kappa = this.kappa
      }
      Repl.SetARubin(ARubin)
			Repl.SetR(*pR1, *pR)
			Repl.InitExog()
			Repl.InitEndog(py1, pY2)
			Repl.InitEstimate(*pr1)

			if (ARubin==0) {
			  ZZg = smatrix(Repl.k, Repl.k)
				deltadenom_b = J(Repl.k, Repl.k, 0)
				ScstaruPZperp = ScstaruPX = ScstaruX = ScstaruXinvXX = ScstaruZperpinvZperpZperp = deltadenom = ScstaruY = Zyg = smatrix(Repl.k+1)
				Scstaruu = smatrix(Repl.k+1, Repl.k+1)
				CTuX = smatrix(Repl.k+1, Clust.N); for (i=Repl.k; i>=0; i--) for (j=Clust.N;j;j--) CTuX[i+1,j].M = J(Ncstar, cols(*Repl.pX1)+cols(*Repl.pX2), 0)
				SCTcapuXinvXX = smatrix(Repl.k+1, Ncstar)
				FillingR2nonCTZperp = FillingR2nonCTX = smatrix(Repl.k+1, Clust.N)
				JNcapk = J(Clust.N, Repl.k, 0)
				FillingT2 = smatrix(Clust.N)
			}
		} else {
      DGP.LIML = this.LIML; DGP.Fuller = this.Fuller; DGP.kappa = this.kappa
    }

    DGP.InitEndog(py1, pY2)

    if (ARubin) {
      if (willplot) {  // for plotting/CI purposes get original point estimate if not normally generated
        DGP.SetR(*pR1, J(0,k,0))  // no-null model in DGP
        DGP.InitEstimate(*pr1)
        DGP.Estimate(*pr1)
        confpeak = *pR * DGP.beta  // estimated coordinate of confidence peak
      }
      pR = &(J(cols(*pX2),kX1,0), I(cols(*pX2)))  // for ARubin test, picks out coefs on excluded exogenous variables
    }
    df = rows(*pR)
    
    if (granular) {
      euXAR = smatrix(df)
      pX = ARubin? &(*pX1, *pX2) : pX1
    }

    DGP.InitEstimate(null? *pr1 \ *pr : *pr1)

    if (ARubin) {
      k = l
      kappa = 0
      pM = &Repl
    } else {
      if (WREnonARubin == 0) DGP.InitTestDenoms()
      pM = &DGP
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
  if (enumerate = (B & weighttype==0 & Ncstar*ln(2) < ln(B)+1e-6))  // generate full Rademacher set?
    MaxMatSize = .

  Nw = MaxMatSize == .? 1 : ceil((B+1) * max((rows(IDBootData), rows(IDBootAll), Ncstar)) * 8 / MaxMatSize / 1.0X+1E) // 1.0X+1E = giga(byte)
  if (Nw == 1) {
    makeWildWeights(B, 1)  // make all wild weights, once
    if (enumerate) B = cols(v) - 1  // replications reduced to 2^G
    WeightGrpStart = 1
  } else {
    seed = rseed()
    _B = ceil((B+1) / Nw)
     WeightGrpStart = (0::Nw-1) * _B :+ 1
    (WeightGrpStop  = (1::Nw  ) * _B     )[Nw] = B+1
  }

  if (bootstrapt & (WREnonARubin | df>1 | MaxMatSize<.)) // unless nonWRE or df=1 or splitting weight matrix, code will create Dist element-by-element, so pre-allocate vector now
    pDist = &J(B+1, 1, .)
  if (Nw>1 | WREnonARubin | (null==0 & df<=2))
    numer = J(df, B+1, .)

  if (interpolable = B & WREnonARubin==0 & null & scoreBS==0 & Nw==1) {    
    dnumerdr = smatrix(q)
    if (interpolate_u = (robust | ML | GMM)==0) dudr = dnumerdr
    if (robust) {
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
void boottest::boottest() {
	real scalar w

	if (initialized==0)
    Init()
  else if (null==0) {  // if not imposing null and we have returned, then df=1 or 2; we're plotting and only test stat, not distribution, changes with r
    if (WREnonARubin)
      numer[,1] = *pR * Repl.beta - *pr
    else if (ARubin) {
      PrepARubin(*pr)
      numer[,1] = u_sd * pM->beta[|kX1+1\.|] // coefficients on excluded instruments in ARubin OLS
    } else
      numer[,1] = u_sd * (*pR * (ML? beta : pM->beta) - *pr) // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this

    (*pDist)[1] = df==1? numer[1] / sqrt(statDenom) : numer[,1] ' invsym(statDenom) * numer[,1]
    return
  }

	if (Nw > 1) {
		rseed(seed)
    makeWildWeights(WeightGrpStop[1] - 1, 1)
  }
	makeInterpolables()  // make stuff that depends linearly on r, possibly by interpolating, for first weight group

  for (w=1; w<=Nw; w++) {  // do group 1 first because it includes col 1, which is all that might need updating in constructing CI in WCU
		if (w > 1)
			makeWildWeights(WeightGrpStop[w] - WeightGrpStart[w] + (w>1), w==1)

		if (WREnonARubin) {
			makeWREStats(w)
			if (bootstrapt==0)
				makeBootstrapcDenom(w)
		} else {
			if (bootstrapt)
				makeNonWREStats(w)
			else
				makeBootstrapcDenom(w)
		}
	}

	BFeas = (*pDist)[1]==.? 0 : colnonmissing(*pDist) - 1
	DistCDR = J(0,0,0)
	setdirty(0)
	initialized = 1
}

// compute bootstrap-c denominator from all bootstrap numerators
void boottest::makeBootstrapcDenom(real scalar w) {
	real colvector tmp
  if (w == 1) {
		tmp = numer[,1]
    statDenom = numer * numer' - tmp * tmp'
		numersum = rowsum(numer) - tmp
	} else {
		statDenom = statDenom + numer * numer'
		numersum = numersum + rowsum(numer)
	}
	if (w == Nw) {  // last weight group?
		statDenom = (statDenom - numersum * numersum' / B) / B
		pDist = &((sqrt? numer:/sqrt(statDenom) : colsum(numer :* invsym(statDenom) * numer))')
	}
}

// draw wild weight matrix of width _B. If first=1, insert column of 1s at front
void boottest::makeWildWeights(real scalar _B, real scalar first) {
	if (_B) {
		if (enumerate)
			v = J(Ncstar,1,1), count_binary(Ncstar, -1-WREnonARubin, 1-WREnonARubin)  // complete Rademacher set
		else if (weighttype==3)
			v = rnormal(Ncstar, _B+1, -WREnonARubin, 1)  // normal weights
		else if (weighttype==4)
			v = rgamma(Ncstar, _B+1, 4, .5) :- (2 + WREnonARubin)  // Gamma weights
		else if (weighttype==2)
			if (WREnonARubin)
				v = sqrt(2 * ceil(runiform(Ncstar, _B+first) * 3)) :* ((runiform(Ncstar, _B+1):>=.5):-.5) :- 1  // Webb weights, minus 1 for convenience in WRE
			else {
				v = sqrt(    ceil(runiform(Ncstar, _B+first) * 3)) :* ((runiform(Ncstar, _B+1):>=.5):-.5)       // Webb weights, divided by sqrt(2)
				u_sd = 1.6a09e667f3bcdX-001 /*sqrt(.5)*/
			}
		else if (weighttype)
			if (WREnonARubin)
				v = ( rdiscrete(Ncstar, _B+first, 1.727c9716ffb76X-001\1.1b06d1d200914X-002 /*.5+sqrt(.05)\.5-sqrt(.05)*/) :- 1.5 ) * 1.1e3779b97f4a8X+001 /*sqrt(5)*/ :- .5  // Mammen weights, minus 1 for convenience in WRE
			else {
				v = ( rdiscrete(Ncstar, _B+first, 1.727c9716ffb76X-001\1.1b06d1d200914X-002 /*.5+sqrt(.05)\.5-sqrt(.05)*/) :- 1.5 ) :+ 1.c9f25c5bfedd9X-003 /*.5/sqrt(5)*/  // Mammen weights, divided by sqrt(5)
				u_sd = 1.c9f25c5bfedd9X-002 /*sqrt(.2)*/
			}
		else if (WREnonARubin) {
      v = runiform(Ncstar, _B+first) :<  .5; v = (-2) * v  // Rademacher weights, minus 1 for convenience in WRE
    } else {
      v = runiform(Ncstar, _B+first) :>= .5; v = v :- .5   // Rademacher weights, divided by 2
      u_sd = .5
    }

		if (first)
			v[,1] = J(Ncstar, 1, WREnonARubin? 0 : u_sd)  // keep original residuals in first entry to compute base model stat		
	} else
		v = J(0,1,0)  // in places, cols(v) indicates number of B -- 1 for classical tests
}


// With reference to Y = [y1 Z], given 0-based columns indexes within it, ind1, ind2, return all bootstrap realizations of Y[,ind1]'(M_Zperp-kappa*M_Xpar)*Y[,ind2] for kappa CONSTANT across replications,
// ind1 can be a rowvector of indices, which then refers to Z only; in this case ind1 vars should all be in X1 or Y2 but not both
real matrix boottest::ProductTermkappa(real rowvector ind1, real scalar ind2, real scalar kappa) {
	real scalar i; real matrix retval, T1, _T1; real colvector T0; struct smatrix rowvector T2

	T0 = Repl.YPXY[ind1:+1, ind2+1]

	_T1 = ind2? Repl.XZ[,ind2] : *Repl.pXy1par
	if (cols(ind1)==1)
		T1 = ScstaruXinvXX[ind1+1].M ' _T1
	else {
	  T1 = J(Ncstar, cols(ind1), 0)
		for (i=cols(ind1);i;i--)
			T1[,i] = ScstaruXinvXX[i+1].M ' _T1
	}

	if (Repl.Yendog[ind2+1]) {
		T1 = T1 + ScstaruXinvXX[ind2+1].M ' (ind1!=0? Repl.XZ[,ind1] : *Repl.pXy1par)

		T2 = smatrix(cols(ind1))
		for (i=cols(ind1);i;i--)  // quadratic term can't be vectorized; vectorize manually
			T2[i].M = ScstaruXinvXX[i+1].M ' ScstaruX[ind2+1].M
	}

	if (kappa != 1) {
		T0 = kappa * T0
		T1 = kappa * T1
		for (i=cols(ind1);i;i--)
			T2[i].M = kappa * T2[i].M
  }

	if (kappa != 1) {  // add stuff for I term
		kappa = 1 - kappa

		T0 = T0 + kappa * (ind1!=0? (ind2? (*Repl.pZZ )[ind1,ind2] : (*Repl.pZy1par)[ind1]) : 
			                          (ind2? (*Repl.pZy1par)[ind2]'  :   Repl.y1pary1par    ) )

		_T1 = ScstaruY[ind2+1].M[,ind1:+1]
		for (i=cols(_T1);i;i--)
			_T1[,i] = _T1[,i] + ScstaruY[ind1[i]+1].M[,ind2+1]
		T1 = T1 + kappa * _T1

		for (i=cols(ind1);i;i--)
			_diag(T2[i].M, diagonal(T2[i].M) + kappa * (ind1[i] <= ind2? Scstaruu[ind2+1, ind1[i]+1].M : Scstaruu[ind1[i]+1, ind2+1].M))
	}

	retval = T0 :+ T1'v
	if (Repl.Yendog[ind2+1])
		for (i=rows(retval);i;i--)
			retval[i,] = retval[i,] + colsum(v :* T2[i].M * v) 
	return(retval)
}


// Workhorse for CRVE sandwich filling
// With reference to notional Y = [y1 Z], given 0-based columns indexes within it, ind1, ind2, return all bootstrap realizations of (P_X-P_Zperp)*Y[,ind1]_g' (M_Zperp*Y[,ind2])_g,
// for all groups in the intersection of all error clusterings
pointer(real matrix) scalar boottest::Filling(real scalar ind1, real scalar ind2) {
	real scalar i; real matrix retval, T1, S

	if (Repl.Yendog[ind1+1]==0 & Repl.Yendog[ind2+1]==0)  // product involves purely exogenous variables, doesn't depend on bootstrap
		return (&(Repl.FillingT0[ind1+1,ind2+1].M))

	if (Repl.Yendog[ind1+1])  // S_∩ (M_(Z_⊥ ) y_(1∥):*X) (X^' X)^(-1) (S_(c^* ) (U ̈_(2∥i):*X))^'
  	T1 = Repl.SMZperpYX[ind2+1].M * ScstaruXinvXX[ind1+1].M

	if (Repl.Yendog[ind2+1]) {  // add CT of u2 * PXY1
		if (NClustVar == NBootClustVar & !subcluster)  // simple case of one clustering: full crosstab is diagonal
			if (Repl.Yendog[ind1+1])
				_diag(T1, diagonal(T1) + ScstaruXinvXX[ind2+1].M ' (ind1? Repl.XZ[,ind1] : *Repl.pXy1par))
			else
				T1 =                     ScstaruXinvXX[ind2+1].M ' (ind1? Repl.XZ[,ind1] : *Repl.pXy1par)  // keep T1 as vector if it's just going to be a diagonal matrix
		else {
		  if (Repl.Yendog[ind1+1]==0)
				T1 = JNcapNcstar
			for (i=Ncstar;i;i--)
				T1[IDCTCapcstar[i].M, i] = T1[IDCTCapcstar[i].M, i] + SCTcapuXinvXX[ind2+1,i].M * (ind1? Repl.XZ[,ind1] : *Repl.pXy1par)
		}

		if (cols(Repl.Zperp))
      T1 = T1 - Repl.SPXYZperp[ind1+1].M * ScstaruZperpinvZperpZperp[ind2+1].M
	}

  retval = Repl.FillingT0[ind1+1,ind2+1].M :+ (cols(T1)==1? T1 :* v : T1 * v)

	if (Repl.Yendog[ind1+1] & Repl.Yendog[ind2+1]) {  // pu1 = NULL => ind1 points to an exogenous control with no varying error component
		if (cols(Repl.RperpX))
			for (i=Clust.N;i;i--) {  // innermost loop
				S = (*pinfoErrData)[i,]', (.\.)
  			FillingT2[i].M = CTuX[ind2+1,i].M * ScstaruXinvXX[ind1+1].M - ScstaruPZperp[ind2+1].M[|S|] ' ScstaruPX[ind1+1].M[|S|]
			}
		else
			for (i=Clust.N;i;i--)
				FillingT2[i].M = CTuX[ind2+1,i].M * ScstaruXinvXX[ind1+1].M

		for (i=rows(retval);i;i--)
			retval[i,] = retval[i,] + colsum(v :* FillingT2[i].M ' v)
	}

	return(&retval)
}


void boottest::makeWREStats(real scalar w) {
	real scalar c, b, i, j, g, kappa
	real colvector numer_b, beta_b
	real rowvector betas, As, _numer
	real matrix ScstaruX1, ScstaruX2, deltanumer, A_b, Jcap, J_b, SuX1g, SuX2g, tmp, Jcaps
	pointer (real colvector) scalar puwt

//real matrix U
//U = (*DGP.py1par, *DGP.pY2) - ((*DGP.pX1,*DGP.pY2)* Rpar * DGP.beta + *DGP.pY2 * DGP.gamma - (*DGP.pX1,*DGP.pX2) * DGP.Pi * DGP.gamma, (*DGP.pX1,*DGP.pX2) * DGP.Pi * Repl.RparY)
timer_on(2)

	U2parddot = DGP.U2ddot * Repl.RparY

	if (w==1 & NClustVar & NFE==0) {  // first/only weight group? prep for optimized computation for bootstrapping cluster when no FE
kappa=1
		for (i=Repl.k; i>=0; i--) {  // precompute various clusterwise sums
			puwt = i? &(U2parddot[,i]) : &DGP.u1dddot; if (haswt) puwt = &(*puwt :* *pwt)

			if (kappa != 1) {
				ScstaruY[i+1].M = _panelsum(*Repl.py1par, *puwt, infoBootData), _panelsum(Repl.Z, *puwt, infoBootData)
				for (j=i; j>=0; j--)
					Scstaruu[i+1,j+1].M = _panelsum(j? U2parddot[,j] : DGP.u1dddot, *puwt, infoBootData)
			}

			// S_cstar(u :* X), S_cstar(u :* Zperp) for residuals u for each endog var; store transposed
			ScstaruX1 = _panelsum(*DGP.pX1, *puwt, infoBootData)'
			ScstaruX2 = _panelsum(*DGP.pX2, *puwt, infoBootData)'

			ScstaruX                 [i+1].M = ScstaruX1 \ ScstaruX2
			ScstaruXinvXX            [i+1].M = *Repl.pinvXX * ScstaruX[i+1].M
			ScstaruZperpinvZperpZperp[i+1].M = Repl.invZperpZperp * _panelsum(Repl.Zperp, *puwt, infoBootData)'
			ScstaruPZperp            [i+1].M =                                                        Repl.Zperp * ScstaruZperpinvZperpZperp[i+1].M
			ScstaruPX                [i+1].M = Repl.kX1? *Repl.pX1*ScstaruX1 + *Repl.pX2*ScstaruX2 : *Repl.pX2   * ScstaruXinvXX            [i+1].M

			// Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u:*X and u:*Zperp (and times invXX or invZperpZperp)
			for (g=Ncstar;g;g--) {
			  SuX1g = _panelsum(*DGP.pX1, *puwt, infoCTCapcstar[g].M)
				SuX2g = _panelsum(*DGP.pX2, *puwt, infoCTCapcstar[g].M)

				SCTcapuXinvXX[i+1,g].M = Repl.kX1? SuX1g * Repl.invXX1 + SuX2g * *Repl.pinvXX2 : SuX2g * *Repl.pinvXX

				// CT_cstar,cap[u:*X]. smatrix has 1 row for each u, 1 col for each all-error-cluster intersection. Elements have 1 row for each bootstrapping cluster, 1 col for each X var
				for (j=rows(SuX1g); j; j--)
					CTuX[i+1, IDCTCapcstar[g].M[j]].M[g,] = SuX1g[j,], SuX2g[j,]
			}
		}

		deltanumer = ProductTermkappa(1..Repl.k, 0, kappa)  // *** handle LIML exception: kappa varies

		for (i=Repl.k;i;i--)
			deltadenom[i].M = ProductTermkappa(1..i, i, kappa)
		for (i=Repl.k;i;i--) {
			Zyg[i].M = *Filling(i,0)
			for(j=Repl.k;j;j--)
  			ZZg[i,j].M = *Filling(i,j)
		}
	}

	if (Repl.k == 1) {  // WRE bootstrap optimized for for Repl.k (# of coefficients in bootstrap regression) = 1
		betas = deltanumer :/ (As = deltadenom.M)
		if (null)
			if (Nw == 1)
				numer = betas :+ (Repl.Rt1 - *pr) / Repl.RRpar
			else
				numer[|WeightGrpStart[w] \ WeightGrpStop[w]|] = _numer = betas :+ (Repl.Rt1 - *pr) / Repl.RRpar
		else {
			if (Nw == 1)
				numer = betas :- DGP.beta0
			else
				numer[|WeightGrpStart[w] \ WeightGrpStop[w]|] = _numer = betas :- DGP.beta0
			if (w==1) {
				numer[1] = betas[1] + (Repl.Rt1 - *pr) / Repl.RRpar
				if (Nw == 1)
					_numer[1] = numer[1]
			}
		}
		if (bootstrapt) {
			if (robust) {
				Jcaps = (Zyg.M :- ZZg.M :* betas) :/ As  // if null tests exog vars only, then ZZ is same in all replications, and ZZg.M has just 1 col
				for (c=1; c<=NErrClustCombs; c++) {  // sum sandwich over error clusterings
					if (NClustVar != 1 & rows(Clust[c].order))
						Jcaps = Jcaps[Clust[c].order,]
					Jcap = _panelsum(Jcaps, Clust[c].info)
					_clustAccum(denom.M, c, colsum(Jcap:*Jcap))
				}
			} else  // spherical errors
				denom.M = Repl.uuc :/ As

			if (Nw==1)
				pDist                                            = &((sqrt?  numer:/sqrt(denom.M) :  numer :*  numer :/ denom.M)')
			else
				(*pDist)[|WeightGrpStart[w] \ WeightGrpStop[w]|] =   (sqrt? _numer:/sqrt(denom.M) : _numer :* _numer :/ denom.M)'
			denom.M = Repl.RRpar * Repl.RRpar * denom.M[1]
		}

	} else {  // WRE bootstrap for Repl.k > 1

		for (b=cols(v); b; b--) {  // WRE bootstrap
			for (i=Repl.k;i;i--)
				deltadenom_b[|.,i\i,i|] = deltadenom[i].M[,b]  // fill uppper triangle, which is all that invsym() looks at
			beta_b  = (A_b = invsym(deltadenom_b)) * deltanumer[,b]
			numer_b = null | w==1 & b==1? (Repl.RRpar * beta_b + Repl.Rt1) - *pr : Repl.RRpar * (beta_b - DGP.beta0)

			if (bootstrapt) {
				if (robust) {  // Compute denominator for this WRE test stat
					Jcap = JNcapk
					for(i=Repl.k;i;i--) {
						tmp = Zyg[i].M[,b]
						for(j=Repl.k-1;j;j--)
							tmp    = tmp - (cols(ZZg[i,     j].M)==1? ZZg[i,     j].M : ZZg[i,     j].M[,b]) * beta_b[     j]  // if a row of null tests exog vars only, then some elements of ZZ are same in all replications, and ZZg[i,j].M has just 1 col
						Jcap[,i] = tmp - (cols(ZZg[i,Repl.k].M)==1? ZZg[i,Repl.k].M : ZZg[i,Repl.k].M[,b]) * beta_b[Repl.k]
					}
					Jcap = Jcap * (A_b * Repl.RRpar')

					for (c=1; c<=NErrClustCombs; c++) {
						if (NClustVar != 1 & rows(Clust[c].order))
							Jcap = Jcap[Clust[c].order,]
						J_b = _panelsum(Jcap, Clust[c].info)
						_clustAccum(denom.M, c, cross(J_b,J_b))
					}
				} else  // spherical errors
					denom.M = (Repl.RRpar * A_b * Repl.RRpar') * Repl.uuc

				(*pDist)[b+WeightGrpStart[w]-1] = sqrt? numer_b/sqrt(denom.M) : cross(numer_b, invsym(denom.M) * numer_b)  // hand-code for 2-dimensional?
			}
			numer[,b+WeightGrpStart[w]-1] = numer_b  // slight inefficiency: in usual bootstrap-t case, only need to save numerators in numer if getdist("numer") is coming because of svmat(numer)
		}
	}
	if (w==1 & bootstrapt) statDenom = denom.M  // original-sample denominator
timer_off(2)
}


void boottest::PrepARubin(real colvector r) {
  Repl.InitEndog(&(*Repl.py1par - *pY2 * r), NULL/*, &(*DGP.pX2y1 - DGP.X2Y2 * r), &(*DGP.pX1y1 - DGP.X1Y2 * r)*/)
  Repl.InitEstimate(*pr1)
  Repl.Estimate(*pr1)
  Repl.InitTestDenoms()
}

// Construct stuff that depends linearly or quadratically on r, possibly by interpolation
void boottest::makeInterpolables() {
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
void boottest::_makeInterpolables(real colvector r) {
  pointer (real matrix) scalar puwt; real scalar d, c; real matrix tmp; pointer (real matrix) scalar pustarXAR, ptmp; real colvector r0

  if (ML)
		uXAR = *pSc * (AR = *pA * *pR')
	else {
    if (ARubin)
      PrepARubin(r)
    else {
      r0 = null? r : J(0, 1, 0); if (REst) r0 =  *pr1 \ r0  // constant terms of model + null constraints
			DGP.InitEstimate(r0)
      DGP.Estimate(r0)
    }

    if (WREnonARubin) return

    puddot = &(pM->u1ddot)
		if (scoreBS | (robust & granular < NErrClustCombs))
			uXAR = *puddot :* pM->XAR
	}
    

  if (scoreBS)
    SuwtXA = B? (NClustVar? _panelsum(uXAR, *pwt, infoBootData) : (haswt?       uXAR :* *pwt  :        uXAR) ) :
                                                                  (haswt? cross(uXAR,   *pwt) : colsum(uXAR)')
  else {  // same calc as in score BS but broken apart to grab intermediate stuff, and assuming residuals defined
    puwt = haswt? &(*puddot :* *pwt) : puddot
    ptmp = &_panelsum(*pX1, *puwt, infoBootData)
    if (ARubin)
      ptmp = &(*ptmp, _panelsum(*pX2, *puwt, infoBootData))
    SuwtXA = *pM->pA * (*ptmp)'
  }

  if (robust & GMM==0 & granular < NErrClustCombs) {
    pustarXAR = &_panelsum(uXAR, *pwt, *pinfoAllData)  // collapse data to all-boot & error-cluster-var intersections. If no collapsing needed, _panelsum() will still fold in any weights
    if (B) {
      if (scoreBS)
        for (d=df;d;d--)
          Kd[d].M = JNcapNcstar  // inefficient, but not optimizing for the score bootstrap
      else
        for (d=df;d;d--) {
          tmp = pM->XAR[,d]; if (haswt) tmp = tmp :* *pwt
          if (ARubin)  // final term in (64) in paper, for c=intersection of all error clusters
            Kd[d].M = (_panelsum(*pX1, tmp, *pinfoErrData), _panelsum(*pX2, tmp, *pinfoErrData)) * SuwtXA
          else
            Kd[d].M =  _panelsum(*pX1, tmp, *pinfoErrData                                      ) * SuwtXA
        }

      if (NFE & FEboot==0)
        CT_WE = _panelsum(crosstabAllFE(wtFE :* *puddot), infoBootAll)'

			for (d=df;d;d--) {  // subtract crosstab of u:*XAR wrt bootstrapping cluster combo and all-cluster-var intersections
				crosstabCapcstarMinus(Kd[d].M, (*pustarXAR)[,d])
				if (scoreBS)
					Kd[d].M = Kd[d].M - ClustShare * colsum(Kd[d].M) // recenter
			}

      for (c=1+granular; c<=NErrClustCombs; c++) {
        if (rows(Clust[c].order))
          for (d=df;d;d--)
            Kd[d].M = Kd[d].M[Clust[c].order,]
        for (d=df;d;d--) {
          Kcd[c,d].M = _panelsum(Kd[d].M, Clust[c].info)
          if (NFE & FEboot==0)
            Kcd[c,d].M = Kcd[c,d].M + _panelsum(pM->CT_XAR[c,d].M, Clust[c].info) * CT_WE
        }
      }
    } else {  // B = 0. In this case, only 1st term of (64) is non-zero after multiplying by v* (= all 1's), and it is then a one-way sum by c
      if (scoreBS)
        pustarXAR = &(*pustarXAR :- ClustShare * colsum(*pustarXAR))  // recenter if OLS

      for (c=1; c<=NErrClustCombs; c++) {
        if (rows(Clust[c].order))
          pustarXAR = &((*pustarXAR)[Clust[c].order,])
        ptmp = &_panelsum(*pustarXAR, Clust[c].info)  // when c=1 (unless subcluster bootstrap), two args have same # of rows, &_panelsum() returns 1st arg by reference. Using & then prevents unnecessary cloning.
        for (d=df;d;d--)
          Kcd[c,d].M = (*ptmp)[,d]
      }
    }
  }

  makeNumerAndJ(1)  // compute J = kappa * v; if Nw > 1, then this is for 1st group; if interpolating, it is only group, and may be needed now to prep interpolation

  if      ( ARubin) numer[,1] = u_sd * pM->beta[|kX1+1\.|] // coefficients on excluded instruments in ARubin OLS
  else if (null==0) numer[,1] = u_sd * (*pR * (ML? beta : pM->beta) - r) // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this.
}

// compute stuff depending linearly on v, needed to prep for interpolation
void boottest::makeNumerAndJ(real scalar w) {  // called to *prepare* interpolation, or when w>1, in which case there is no interpolation
  real scalar c, d

  if (Nw==1)
    if (scoreBS)
      numer                                         = B? cross(SuwtXA, v) : SuwtXA * u_sd
    else if (robust==0 | granular)
	    numer                                         = *pR * (betadev = SuwtXA * v)
    else
	    numer                                         = (*pR * SuwtXA) * v
  else
    if (scoreBS)
      numer[|WeightGrpStart[w] \ WeightGrpStop[w]|] = B? cross(SuwtXA, v) : SuwtXA * u_sd
    else if (robust==0 | granular)
      numer[|WeightGrpStart[w] \ WeightGrpStop[w]|] = *pR * (betadev = SuwtXA * v)
    else
      numer[|WeightGrpStart[w] \ WeightGrpStop[w]|] = (*pR * SuwtXA) * v

  if (B & robust & GMM==0) {
    if (granular)  // prep optimized treatment when bootstrapping by many/small groups
      if (purerobust)
        ustar = *DGP.partialFE(&(*puddot :* v)) - *pX * betadev
      else {  // clusters small but not all singletons
        if (NFE) {
          ustar = *DGP.partialFE(&(*puddot :* v[IDBootData,]))
          for (d=df;d;d--)
            (*pJcd)[1,d].M = _panelsum(ustar,             pM->WXAR[d].M, *pinfoErrData)                                   - _panelsum(*pX, pM->WXAR[d].M, *pinfoErrData) * betadev
        } else
          for (d=df;d;d--)
            (*pJcd)[1,d].M = _panelsum(_panelsum(*puddot, pM->WXAR[d].M, *pinfoAllData) :* v[IDBootAll,], infoErrAll) - _panelsum(*pX, pM->WXAR[d].M, *pinfoErrData) * betadev
      }

    for (c=NErrClustCombs; c>granular; c--)
      for (d=df;d;d--)
        (*pJcd)[c,d].M = Kcd[c,d].M * v
  }
}

void boottest::makeNonWREStats(real scalar w) {
	real scalar i, c, j, l; real matrix ustar2, tmp; real colvector numer_l; pointer (real matrix) scalar pAR; real rowvector t1, t2, t12

  if (w > 1) makeNumerAndJ(w)

	if (robust & GMM==0) {
    if (interpolating==0) {  // these quadratic computation needed to *prepare* for interpolation but are superseded by interpolation once it is going
      if (purerobust)
        ustar2 = ustar :* ustar
      for (i=df;i;i--)
        for (j=i;j;j--) {
          if (c = purerobust) _clustAccum(denom[i,j].M, c, cross(pM->WXAR[i].M, pM->WXAR[j].M, ustar2))  // c=purerobust not a bug
          for (c++; c<=NErrClustCombs; c++)
            _clustAccum(denom[i,j].M, c, colsum((*pJcd)[c,i].M :* (*pJcd)[c,j].M))  // (60)
        }
    }

		if (df == 1) {
			if (Nw > 1)
				(*pDist)[|WeightGrpStart[w] \ WeightGrpStop[w]|] =   (numer[|WeightGrpStart[w] \ WeightGrpStop[w]|] :/ sqrt(denom.M))'
			else
          pDist                                          = &((numer                                         :/ sqrt(denom.M))')
			if (w==1)
				statDenom = denom.M[1]  // original-sample denominator
		} else if (df==2) {  // hand-code 2D numer'inv(denom)*numer
    	t1 = numer[1,]; t2 = numer[2,]; t12 = t1:*t2
			if (Nw > 1)
				(*pDist)[|WeightGrpStart[w] \ WeightGrpStop[w]|] =   ( (t1:*t1:*denom[2,2].M - (t12+t12):*denom[2,1].M + t2:*t2:*denom[1,1].M) :/ (denom[1,1].M:*denom[2,2].M - denom[2,1].M:*denom[2,1].M) )'
			else                                                                                                                         
          pDist                                          = &(( (t1:*t1:*denom[2,2].M - (t12+t12):*denom[2,1].M + t2:*t2:*denom[1,1].M) :/ (denom[1,1].M:*denom[2,2].M - denom[2,1].M:*denom[2,1].M) )')
			if (w==1)
				statDenom = denom[1,1].M[1], denom[2,1].M[1] \ denom[2,1].M[1], denom[2,2].M[1]  // original-sample denominator
    } else {  // build each replication's denominator from vectors that hold values for each position in denominator, all replications
			tmp = J(df,df,0)
			for (l=cols(v); l; l--) {
				for (i=df;i;i--)
					for (j=i;j;j--)
						tmp[j,i] = denom[i,j].M[l]  // fill upper triangle, which is all invsym() looks at
				numer_l = numer[,l]
				(*pDist)[l+WeightGrpStart[w]-1] = numer_l ' invsym(tmp) * numer_l  // in degenerate cases, cross() would turn cross(.,.) into 0
			}
			if (w==1)
				statDenom = tmp  // original-sample denominator
		}

	} else { // non-robust or GMM

		pAR = ML? &AR : &(pM->AR)
		if (df == 1) {  // optimize for one null constraint
			denom.M = *pR * *pAR

			if ((ML | GMM)==0) {
                     ustar = B? v :* *puddot : *puddot
				if (scoreBS) ustar = ustar :- (haswt? cross(ClustShare, ustar) : colsum(ustar) * ClustShare)  // Center variance if interpolated
				        else ustar = ustar  - (*pX1, *pX2) * betadev  // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
				denom.M = denom.M :* (haswt? cross(*pwt, ustar :* ustar) : colsum(ustar :* ustar))
			}
			if (Nw > 1)
				(*pDist)[|WeightGrpStart[w] \ WeightGrpStop[w]|] =   (numer :/ sqrt(denom.M))'
			else
				  pDist                                          = &((numer :/ sqrt(denom.M))')
			if (w==1)
				statDenom = denom.M[1]  // original-sample denominator
		} else {
			denom.M = *pR * *pAR

			if (ML | GMM) {
				for (l=cols(v); l; l--) {
					numer_l = numer[,l]
					(*pDist)[l+WeightGrpStart[w]-1] = cross(numer_l, invsym(denom.M), numer_l)
				}
				if (w==1)
					statDenom = denom.M  // original-sample denominator
			} else {
				for (l=cols(v); l; l--) {
					numer_l = numer[,l]
					(*pDist)[l+WeightGrpStart[w]-1] = cross(numer_l, invsym(denom.M) * numer_l)
                       ustar = B? v[,l] :* *puddot : *puddot
					if (scoreBS) ustar = ustar :- (haswt? cross(*pwt, ustar) : colsum(ustar)) * ClustShare  // Center variance if interpolated
					        else ustar = ustar  - (*pX1, *pX2) * betadev[,l]  // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
					(*pDist)[l+WeightGrpStart[w]-1] = (*pDist)[l+WeightGrpStart[w]-1] / (tmp = cross(ustar, *pwt, ustar))
				}
				if (w==1)
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
// if v = 0 (so can't tell if row or col vector), returns J(1, 0, 0) 
real vector boottest::_selectindex(real vector v) {
	real scalar rows
	if (v==0) return(J(1,0,0))
	rows = rows(v)
	return(select(rows>1? 1::rows : 1..cols(v), v))
}

// Return matrix that counts from 0 to 2^N-1 in binary, one column for each number, one row for each binary digit
// except use provided lo and hi values for 0 and 1
real matrix boottest::count_binary(real scalar N, real scalar lo, real scalar hi) {
	real matrix tmp
	if (N<=1) return (lo , hi)
	tmp = count_binary(N-1, lo, hi)
	return (J(1, cols(tmp), lo), J(1, cols(tmp), hi) \ tmp, tmp)
}

// cross-tab sum of a column vector w.r.t. intersection-of-error & bootstrap-clustering-vars and fixed-effect var
real matrix boottest::crosstabAllFE(real colvector v) {
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

// subtract crosstab of v wrt bootstrapping cluster and all-cluster-var intersections from M
// M should have one row for each all-cluster-var (including bootstrap cluster) intersection and one col for each bootstrap cluster
// *** v needs to have been panelsum'd with pinfoAllData
void boottest::crosstabCapcstarMinus(real matrix M, real colvector v) {
	real colvector tmp; real scalar i

	if (subcluster)  // crosstab c,c* is wide
		for (i=Clust.N;i;i--) {
			tmp = infoErrAll[i,]'
			M[|(i\i), tmp|] = M[|(i\i), tmp|] - v[|tmp|]'
		}
	else if (NClustVar == NBootClustVar)  // crosstab c,c* is square
		_diag(M, diagonal(M) + v)
	else  // crosstab c,c* is tall
		for (i=Ncstar;i;i--) {
			tmp = pBootClust->info[i,]'
			M[|tmp, (i\i)|] = M[|tmp, (i\i)|] - v[|tmp|]
		}
}

// given a pre-configured boottest linear model with one-degree null imposed, compute distance from target p value of boostrapped one associated with given value of r
// used with optimize() to construct confidence intervals
// performs no error checking
real scalar boottest::r_to_p(real colvector r) {
	pr = &r
	setdirty(1, 1) // set dirty = 1, but leave initialized=0, which we want when only changing r
	return (getpadj())
}

real scalar boottest::search(real scalar alpha, real scalar p_lo, real scalar lo, real scalar p_hi, real scalar hi) {
	real scalar mid, _p
	mid = lo + (alpha-p_lo)/(p_hi-p_lo)*(hi-lo)  // interpolate linearly
//	mid = lo + ( reldif(p_hi-p_lo, 1/BFeas)<1e-6 & BFeas? (hi - lo)/2 : (alpha-p_lo)/(p_hi-p_lo)*(hi-lo) ) // interpolate linearly until bracketing a single bump-up; then switch to binary search
//	mid = lo + (hi - lo)/2  // binary search
	if (reldif(lo,mid)<ptol | reldif(hi,mid)<ptol | (BFeas & abs(p_hi-p_lo)<(1+(ptype==1))/BFeas*1.000001))
		return (mid)
	if ( ((_p = r_to_p(mid)) < alpha) == (p_lo < alpha) )
		return(search(alpha, _p, mid, p_hi, hi))
	return(search(alpha, p_lo, lo, _p, mid))
}


// derive wild bootstrap-based CI, for case of linear model with one-degree null imposed.
void boottest::plot() {
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
    lo = _selectindex(CI:== 1)
    hi = _selectindex(CI:==-1)
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
real matrix boottest::combs(real scalar d) {
	real matrix retval; real scalar i
	retval = J(2^d, 0, 0)
	for (i=d;i;i--)
		retval = retval, J(2^(d-i),1,1) # (1\0) # J(2^(i-1),1,1) 
	return (retval)
}

// Like Mata's order() but does a stable sort
real colvector boottest::stableorder(real matrix X, real rowvector idx)
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
timer_clear()
timer_on(1)
	real matrix C1, R1, C, R, ZExcl, ID, FEID, sc, XEnd, XEx
	real colvector r1, wt, r, Y
	class boottest scalar M
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
	M.close()
timer_off(1)
timer()
}

mata mlib create lboottest, dir("`c(sysdir_plus)'l") replace
mata mlib add lboottest *(), dir("`c(sysdir_plus)'l")
mata mlib index
end
boottest DL.lpris_totpop=-1/3, cluster(state year) bootcluster(year) ptype(equaltail) reps(199) noci seed(1231)