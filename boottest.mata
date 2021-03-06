*! boottest 3.1.0 13 March 2021
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

// return pointer to chosen columns of a matrix, but don't duplicate data if return value is whole matrix
pointer (real matrix) pcol(real matrix A, real vector p)
	return(length(p)==cols(A)? &A : &A[,p])

// return [X1 X2]B where X1 or X2 may have 0 cols
real matrix X12B(real matrix X1, real matrix X2, real matrix B)
	return(cols(X1)? (cols(X2)? X1 * B[|.,.\cols(X1),.|] + X2 * B[|cols(X1)+1,.\.,.|] : X1 * B) : X2 * B)

// return &X[|S|] with appropiate behavior if cols(X)=0 or S requests no rows (S[2,1]=0), presuming that in degenerate cases S does not specify columns
// if retval = X, doesn't duplicate data
// S should be 2x1 because the function is only for selecting rows
pointer(real matrix) scalar pXS(real matrix X, real matrix S)
	return(cols(X)? (S[2,1]? ((S[1,1]==. | S[1,1]==1) & (S[2,1]==. | S[2,1]==rows(X))? &X : &X[|S,(.\.)|]) : &J(0,cols(X),0)) : &J(editmissing(S[2,1],rows(X))-editmissing(S[1,1],1)+1,0,0))

// :* operation, handling case that either argument is just 1 without duplicating data
pointer (real colvector) scalar pvHadw(real matrix v, real matrix w)
	return(w==1? &v : (v==1? &w : &(v :* w)))

class boottestOLS {  // class for analyitcal OLS, 2SLS, LIML, GMM estimation--everything but iterative ML
	real scalar LIML, uu, Fuller, ARubin, kappa, isDGP, kZ, kX1, kX2, kX
	real colvector u1ddot, u1dddot, beta, beta0, PXy1, invXXXy1par
	real rowvector Yendog
  real matrix Zperp, invZperpZperp, ZperpinvZperpZperp, XZ, PXZ, Z, YPXY, R1invR1R1, R1perp, Rpar, RperpX, RRpar, R1invR1R1X, RparX, RparY, R1AR1, RR1invR1R1, invH, dbetadr, YY, AR, XAR, R1invR1R1Y, invXXXZ, U2ddot, invXX1, XinvXX, Rt1
	pointer(real colvector) scalar py1, py1par, pZy1par, pXy1par, pr1, pZy1, pX1y1
	pointer(real matrix) scalar pY2, pA, pinvXX, pH, pX1, pX2, pinvXX2
	pointer (class boottest scalar) scalar parent
	struct smatrix matrix FillingT0
	struct smatrix rowvector WXAR, SPXYZperp, SMZperpYX, CT_XAR, CT_FEcapPY

	private void new(), InitTestDenoms()
  private virtual void InitVars(), SetR(), InitEstimate(), Estimate(), MakeResiduals()
	real matrix _select(), perp()
}

class boottestARubin extends boottestOLS {
  real matrix XY2
	private void new()
	private virtual void InitVars(), InitEstimate(), Estimate()
}

class boottestIVGMM extends boottestOLS {
	real matrix XY2, XX, H_2SLS, V, ZY2, X2Y2, X1Y2, Zex, ZR1ex, ZR1, ZR1ZR1, X2ZR1, ZR1Y2, X1ZR1, ZZR1
  real colvector ZXinvXXXy1par
  real rowvector y1Y2, y1ZR1
  real scalar y1y1, y1pary1par
  pointer(real colvector) scalar pX2y1par, pX1y1par, pX2y1
	pointer(real rowvector) scalar py1parY2
  pointer(real matrix) scalar pZZ, _pRperp
	private void new()
	private virtual void InitVars(), InitEstimate(), Estimate(), MakeResiduals()
}

class boottest {
	real scalar scoreBS, B, small, weighttype, null, dirty, initialized, ML, Nobs, _Nobs, kZ, kY2, kX1, l, sumwt, NClustVar, haswt, REst, multiplier, smallsample, quietly, FEboot, NErrClustCombs, ///
		sqrt, hascons, LIML, Fuller, kappa, IV, WRE, WREnonARubin, ptype, twotailed, df, df_r, ARubin, confpeak, willplot, notplotted, NumH0s, p, NBootClustVar, NErrClust, BootClust, ///
		NFE, granular, purerobust, subcluster, Nstar, BFeas, u_sd, level, ptol, MaxMatSize, Nw, enumerate, bootstrapt, q, interpolable, interpolating, interpolate_u, robust, kX2, kX
	real matrix AR, v, ustar, CI, CT_WE, infoBootData, infoBootAll, infoErrAll, JNcapNstar, statDenom, uXAR, SuwtXA, numer0, betadev, IDCap, U2parddot, deltadenom_b, _Jcap, YYstar_b, YPXYstar_b, numerw
	real colvector DistCDR, plotX, plotY, beta, ClustShare, WeightGrpStart, WeightGrpStop, gridmin, gridmax, gridpoints, numersum, uddot0, anchor, poles, invFEwt
	real rowvector peak, betas, As
	string scalar wttype, madjtype, seed
	pointer (real matrix) scalar pX2, pR1, pR, pID, pFEID, pY2, pX1, pinfoAllData, pinfoCapData, pIDAll, pnumer
	pointer (real colvector) scalar pr1, pr, py1, pSc, pwt, pA, puddot, pDist, pIDBootData, pIDBootAll
	class boottestOLS scalar DGP, Repl
  pointer(class boottestOLS scalar) scalar pM
	struct structboottestClust rowvector Clust
	struct smatrix matrix denom, Kcd, denom0, Jcd0, SCTcapuXinvXX, Sstaruu, CTuX
	struct smatrix rowvector Kd, dudr, dnumerdr, IDCTCapstar, infoCTCapstar, SstaruX, SstaruXinvXX, SstaruZperpinvZperpZperp, deltadenom, Zyg, SstaruY, SstaruMZperp, SstaruPX, SstaruZperp, YYstar, YPXYstar, CTFEu
  struct ssmatrix rowvector ddenomdr, dJcddr
  struct ssmatrix matrix ddenomdr2
	pointer(struct smatrix matrix) scalar pJcd
	struct structFE rowvector FEs
	
	void new(), setsqrt(), setX1(), setptype(), setdirty(), setY2(), setY(), setX2(), setwt(), setsc(), setML(), setLIML(), setARubin(),
		setFuller(), setkappa(), setquietly(), setbeta(), setA(), setsmall(), sethascons(), setscoreBS(), setB(), setnull(), setWald(), setRao(), setwttype(), setID(), setFEID(), setlevel(), setptol(), 
		setrobust(), setR1(), setR(), setwillplot(), setgrid(), setmadjust(), setweighttype(), setMaxMatSize(), setstattype(), close()
  private void makeNumerAndJ(), _clustAccum(), makeWREStats(), makeInterpolables(), _makeInterpolables(), makeNonWREStats(), makeBootstrapcDenom(), Init(), plot(), makeWildWeights(), boottest(), crosstabCapstarMinus(), PrepWRE(), storeWtGrpResults()
	real matrix getplot(), getCI(), getV(), getv()
	real scalar getp(), getpadj(), getstat(), getdf(), getdf_r(), getreps(), getrepsFeas(), getNBootClust()
	real rowvector getpeak()
	real colvector getdist(), getb()
	private real scalar r_to_p(), search()
  private real matrix count_binary(), crosstabFE(), HessianFixedkappa()
	private pointer(real matrix) scalar Filling(), partialFE()
  private static real matrix combs()
  private static real colvector stableorder()
	private real vector _selectindex()
  static void _st_view()
}

void boottestOLS::new() {
	kappa = ARubin = 0
  isDGP = 1
}

void boottestARubin::new() {
	ARubin = isDGP = 1
  LIML = Fuller = kappa = 0    
}

void boottestIVGMM::new() {
	ARubin = 0
  kappa = isDGP = 1
}


// do select() but handle case that both cases are scalar and second is 0 by interpreting second arg as rowvector and returning J(1,0,0)
// if v = 0 (so can't tell if row or col vector), returns J(1, 0, 0) 
real matrix boottestOLS::_select(real matrix X, real rowvector v)
	return (rows(X)==1 & cols(X)==1 & v==0? J(1,0,0) : select(X,v))

real matrix boottestOLS::perp(real matrix A) {
	real matrix vec; real rowvector val
	symeigensystem(A*invsym(A'A)*A', vec, val)
	_edittozero(val, 1000)
	return (_select(vec, !val))
}


// R1 is constraints. R is attack surface for null; only needed when using FWL for WRE
// for DGP regression, R1 is maintained constraints + null if imposed while R should have 0 rows
// for replication regressions R1 is maintained constraints, R is null
void boottestOLS::SetR(real matrix R1, | real matrix R) {
	real matrix RR1perp, vec, S; real rowvector val
	pragma unset vec; pragma unset val

	if (rows(R1)) {
		R1invR1R1 = invsym(R1 * R1')
		if (all(diagonal(R1invR1R1))==0)
			_error(111, "Null hypothesis or model constraints are inconsistent or redundant.")
		R1invR1R1 = R1 ' R1invR1R1
		symeigensystem(R1invR1R1 * R1, vec, val); _edittozero(val, 10)
		R1perp = _select(vec, !val)  // eigenvectors orthogonal to span of R1; foundation for parameterizing subspace compatible with constraints
	} else
		R1invR1R1 = J(parent->kZ,0,0)  // and R1perp = I

  if (kappa) {			
		// prepare to reduce regression via FWL
		RR1perp = R \ J(parent->kY2, parent->kX1, 0), I(parent->kY2)  // rows to prevent partialling out of endogenous regressors

		if (rows(R1))
			RR1perp = RR1perp * R1perp 
		symeigensystem(RR1perp ' invsym(RR1perp * RR1perp') * RR1perp, vec, val); _edittozero(val, 10)
		Rpar   = _select(vec,  val)
		RperpX = _select(vec, !val)

		if (rows(R1)) {  // fold model constraint factors into Rpar, RperpX
			Rpar   = R1perp * Rpar
			RperpX = R1perp * RperpX
		}
		RRpar = R * Rpar

		S = .\parent->kX1
		RperpX     = *pXS(RperpX   , S)  // Zperp=Z*RperpX; though formally a multiplier on Z, it will only extract exogenous components, in X1, since all endogenous ones will be retained
		RparX      = *pXS(Rpar     , S)  // part of Rpar that refers to X1
		R1invR1R1X = *pXS(R1invR1R1, S)

		S = parent->kX1+1\.
		RparY      = *pXS(Rpar     , S)  // part of Rpar that refers to Y2
		R1invR1R1Y = *pXS(R1invR1R1, S)

    RR1invR1R1 = R * R1invR1R1
	}
}

// stuff that can be done before r set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
void boottestOLS::InitVars(pointer(real matrix) scalar pRperp) {  // Rperp is for replication regression--no null imposed
	real matrix X2X1

  pX1 = parent->pX1
  pX2 = parent->pX2
	py1par = parent->py1
  pZy1par = &cross(*pX1, *parent->pwt, *py1par)
  pH = &cross(*pX1, *parent->pwt, *pX1)

  R1AR1 = rows( R1perp )?    R1perp * invsym( R1perp ' (*pH) *  R1perp) *  R1perp'   :  invsym(*pH)
	pA    = rows(*pRperp )? &(*pRperp * invsym(*pRperp ' (*pH) * *pRperp) * *pRperp' ) : &invsym(*pH)
  AR = *pA * *parent->pR'
	if (parent->scoreBS | parent->robust)
		XAR = X12B(*pX1, *pX2, AR)

  beta0   = R1AR1 * *pZy1par
  dbetadr = R1AR1 * *pH * R1invR1R1 - R1invR1R1
}

// stuff that can be done before r set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
void boottestARubin::InitVars(| pointer(real matrix) pRperp) {
  pragma unused pRperp

	real matrix X2X1

  pX1 = parent->pX1
  pX2 = parent->pX2
	py1 = parent->py1
	pY2 = parent->pY2
  X2X1 = cross(*pX2, *parent->pwt, *pX1)
  pH = &(cross(*pX1, *parent->pwt, *pX1), X2X1' \ X2X1, cross(*pX2, *parent->pwt, *pX2))
  pA = &invsym(*pH)
  AR = *pA * *parent->pR'
	if (parent->scoreBS | parent->robust)
		XAR = X12B(*pX1, *pX2, AR)

  pZy1 = &cross(*pX1, *parent->pwt, *py1)
  XY2 =    cross(*pX1, *parent->pwt, *pY2) \ cross(*pX2, *parent->pwt, *pY2) // for use on LHS-- Y - Y2 *r  "ZR1" = Y2
  pZy1 = &(*pZy1 \ cross(*pX2, *parent->pwt, *py1))  // X1 & X2 become "Z" (RHS) in A-R replication regression

  R1AR1 = rows(R1perp)? R1perp * invsym(R1perp ' (*pH) * R1perp) * R1perp' : *pA
}

void boottestIVGMM::InitVars(|pointer(real matrix) scalar pRperp) {
	real matrix X2X1, XX1, ZperpZ, ZperpY2, ZperpZR1; real scalar i,j; pointer (real colvector) scalar puwt
  
  _pRperp = pRperp  // can't process the replication regressions' R matrix yet because A potenntially depends on kappa, which depends on r; so save it for now

  Zex   = *parent->pX1 * RparX       // exogenous component of Zpar
  ZR1ex = *parent->pX1 * R1invR1R1X  // exogenous component of ZR1

  Zperp = *parent->pX1 * RperpX
  invZperpZperp = invsym(cross(Zperp, *parent->pwt, Zperp))
  ZperpinvZperpZperp = Zperp * invZperpZperp

  pX1 = &(*parent->pX1 * perp(RperpX))  // FWL-process X1
  pX1 = &(*pX1 - ZperpinvZperpZperp * cross(Zperp, *parent->pwt, *pX1))
  pX2 = &(*parent->pX2 - ZperpinvZperpZperp * cross(Zperp, *parent->pwt, *parent->pX2))
  kX = (kX1 = cols(*pX1)) + (kX2 = cols(*pX2))
  kZ = cols(Rpar)

  X2X1 = cross(*pX2, *parent->pwt, *pX1)
  XX1 = cross(*pX1, *parent->pwt, *pX1) \ X2X1
  pinvXX = &invsym((XX = XX1, (X2X1' \ cross(*pX2, *parent->pwt, *pX2))))

  if (isDGP==0) {
    invXX1 = *pXS(*pinvXX, .\kX1  )
    pinvXX2 = pXS(*pinvXX, kX1+1\.)
  
    SMZperpYX = SPXYZperp = smatrix(kZ+1)
    XinvXX = (*pX1) * invXX1 + *pX2 * *pinvXX2
  }

	// handle endogenous variables
	py1 = parent->py1
	pY2 = parent->pY2

  Z   = Zex    + *pY2 * RparY // Zpar
  ZR1 = ZR1ex :+ *pY2 * R1invR1R1Y  // :+ needed in case cols(ZR1ex)=0
  
  ZperpY2  = cross(Zperp, *parent->pwt, *pY2 )
  ZperpZ   = cross(Zperp, *parent->pwt, Zex  ) + ZperpY2 * RparY
  ZperpZR1 = cross(Zperp, *parent->pwt, ZR1ex) + ZperpY2 * R1invR1R1Y

  Z   =    Z   - ZperpinvZperpZperp * ZperpZ
  ZR1 =    ZR1 - ZperpinvZperpZperp * ZperpZR1
  pY2 = &(*pY2 - ZperpinvZperpZperp * ZperpY2)
  py1 = &(*py1 - ZperpinvZperpZperp * cross(Zperp, *parent->pwt, *py1))

  X1Y2 = cross(*pX1, *parent->pwt, *pY2)
  X2Y2 = cross(*pX2, *parent->pwt, *pY2)
  XY2 = X1Y2 \ X2Y2
  y1Y2 =   cross(*py1, *parent->pwt, *pY2)
  pX2y1 = &cross(*pX2, *parent->pwt, *py1)
  pX1y1 = &cross(*pX1, *parent->pwt, *py1)
  y1y1 =   cross(*py1, *parent->pwt, *py1)
  pZy1 =  &cross( Z  , *parent->pwt, *py1)		
  XZ =     cross(*pX1, *parent->pwt, Z   ) \ 
           cross(*pX2, *parent->pwt, Z   )
  ZY2 =    cross(Z   , *parent->pwt, *pY2)
  pZZ =   &cross(Z   , *parent->pwt, Z   )

  invXXXZ = *pinvXX * XZ

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

  V =  *pinvXX * XZ // in 2SLS case, estimator is (V' XZ)^-1 * (V'Xy1). Also used in kZ-class and LIML robust VCV by Stata convention
  H_2SLS = V ' XZ  // Hessian

  if (isDGP==0) {
    Yendog = 1, colsum(RparY :!= 0)  // columns of Y = [y1par Zpar] that are endogenous (normally all)

    if (parent->robust) {  // for WRE replication regression, prepare for CRVE

      if (parent->bootstrapt) {
        PXZ = X12B(*pX1, *pX2, invXXXZ)

        FillingT0 = smatrix(kZ+1, kZ+1)  // fixed component of groupwise term in sandwich filling
        if (parent->NFE)
          CT_FEcapPY = smatrix(kZ+1)
        for (i=kZ; i; i--) {
          puwt = pvHadw(*pcol(PXZ,i), *parent->pwt)
          if (parent->NFE)
            CT_FEcapPY[i+1].M = parent->crosstabFE(*puwt, *parent->pinfoCapData) :* parent->invFEwt
          for (j=kZ; j; j--)
            FillingT0[i+1,j+1].M = *_panelsum(*pcol(Z,j), *puwt, *parent->pinfoCapData)
        }

        for (i=kZ; i; i--) {  // precompute various clusterwise sums
          SPXYZperp[i+1].M = *_panelsum(Zperp, *pvHadw(*pcol(PXZ,i), *parent->pwt), *parent->pinfoCapData)  // S_cap(P_(MZperpX) * Z :* Zperp)
          if (parent->granular==0)
            SMZperpYX[i+1].M = *_panelsum2(*pX1, *pX2, *pvHadw(*pcol(Z,i), *parent->pwt), *parent->pinfoCapData)  // S_cap(M_Zperp[Z or y1] :* P_(MZperpX)])
        }
      }
    }
  }
}


// do most of estimation; for LIML r1 must be passed now in order to solve eigenvalue problem involving it
// inconsistency: for replication regression of Anderson-Rubin, r1 refers to the *null*, not the maintained constraints, because that's what affects the endogenous variables
// For OLS, compute beta0 (beta when r=0) and dbetadr without knowing r1, for efficiency
// For WRE, should only be called once for the replication regressions, since for them r1 is the unchanging model constraints
void boottestOLS::InitEstimate(real colvector _r1)
  pr1 = &_r1

void boottestARubin::InitEstimate(real colvector _r1) {
  py1par   = &(*py1  - *pY2 * _r1)
  pZy1par  = &(*pZy1  - XY2 * _r1)
}

void boottestIVGMM::InitEstimate(real colvector _r1) {
	real rowvector val; real matrix vec; real scalar i; pointer (real colvector) scalar puwt
	pragma unset vec; pragma unset val

  pr1 = &_r1

  if (cols(R1invR1R1)) {
		y1pary1par = y1y1 - 2 * y1ZR1 * *pr1 + *pr1 ' ZR1ZR1 * *pr1
		py1parY2 =  &(y1Y2  - *pr1 ' ZR1Y2)
		pX2y1par = &(*pX2y1 - X2ZR1 * *pr1)
		pX1y1par = &(*pX1y1 - X1ZR1 * *pr1)
		pZy1par  = &(*pZy1  - ZZR1  * *pr1)
		pXy1par  = &(*pX1y1par \ *pX2y1par)

		py1par = &(*py1 - ZR1 * *pr1)
	} else {
		py1par  = py1
		pZy1par = pZy1
	}

  ZXinvXXXy1par = invXXXZ ' (*pXy1par)
  YY = y1pary1par, *pZy1par' \ *pZy1par, *pZZ
  invXXXy1par = *pinvXX * *pXy1par
  YPXY = invXXXy1par ' (*pXy1par) , ZXinvXXXy1par' \ ZXinvXXXy1par , XZ ' invXXXZ

  if (isDGP) {
    if (LIML) {
      eigensystemselecti(invsym(YY) * YPXY, rows(YY)\rows(YY), vec, val)
      kappa = 1/(1 - Re(val)) // sometimes a tiny imaginary component sneaks into val
      if (Fuller) kappa = kappa - 1 / (parent->_Nobs - parent->l)
    }

    pH = kappa==1? &H_2SLS : &((1-kappa) * *pZZ + kappa * H_2SLS)
    invH = invsym(*pH)

    if (_pRperp) {  // for score bootstrap
      pA = cols(*_pRperp)? &(*_pRperp * invsym(*_pRperp ' (*pH) * *_pRperp) * *_pRperp') : &invH
      AR = *pA * (parent->scoreBS? *parent->pR' : RRpar')

      if (parent->scoreBS | parent->robust)
        XAR = X12B(*pX1, *pX2, V * AR)
    }
  } else if (parent->WREnonARubin) {
    Rt1 = RR1invR1R1 * *pr1
    
    if (parent->robust & parent->bootstrapt) {  // prepare WRE replication regressions
      PXy1 = *pX2 * invXXXy1par[|kX1+1,.\.,.|]; if (kX1) PXy1 = PXy1 + *pX1 * invXXXy1par[|.,.\kX1,.|]

      puwt = pvHadw(PXy1, *parent->pwt)
      for (i=kZ; i; i--)
        FillingT0[1,i+1].M = *_panelsum(*pcol(Z,i), *puwt, *parent->pinfoCapData)
      SPXYZperp.M = *_panelsum(Zperp, *puwt, *parent->pinfoCapData)  // S_cap(P_(MZperpX) * y1 :* Zperp)

      if (parent->NFE)
        CT_FEcapPY.M = parent->crosstabFE(*puwt, *parent->pinfoCapData) :* parent->invFEwt

      puwt = pvHadw(*py1par, *parent->pwt)
      for (i=kZ; i; i--)
        FillingT0[i+1,1].M = *_panelsum(*pcol(PXZ,i), *puwt, *parent->pinfoCapData)
      FillingT0.M          = *_panelsum(PXy1        , *puwt, *parent->pinfoCapData)
      if (parent->granular==0)
        SMZperpYX.M        = *_panelsum2(*pX1, *pX2 , *puwt, *parent->pinfoCapData)  // S_cap(M_Zperp*y1 :* P_(MZperpX)])
    }
  }
}

// stuff that depends on r1 and on Hessian (thus potentially on kappa and endog vars)
void boottestOLS::Estimate()
  beta = beta0 - dbetadr * *pr1

void boottestARubin::Estimate()
  beta = R1AR1 * *pZy1par

// stuff that depends on r1 and on Hessian (thus potentially on kappa and endog vars)
void boottestIVGMM::Estimate() {
  beta = invH * (kappa==1?  ZXinvXXXy1par : kappa * ZXinvXXXy1par + (1-kappa) * *pZy1par)
}

void boottestOLS::MakeResiduals() {
  real colvector _beta
  
  u1ddot = *py1par - X12B(*pX1, *pX2, beta)

  if ((parent->robust | parent->scoreBS | parent->bootstrapt)==0) {
  	_beta = 1 \ -beta; uu = _beta ' YY * _beta
    if (parent->hascons == 0)
      uu = uu - cross(*parent->pwt, u1ddot)^2 / parent->_Nobs  // sum of squares after centering, N * Var
  }
}

void boottestIVGMM::MakeResiduals() {
	real matrix Xu; real colvector negXuinvuu; real colvector _beta
  
  u1ddot = *py1par - Z * beta

  if (parent->scoreBS==0) {
    _beta = 1 \ -beta
    uu = _beta ' YY * _beta

    Xu = *pXy1par - XZ * beta  // after IV/GMM DGP regression, compute Y2 residuals by regressing Y2 on X while controlling for y1 residuals, done through FWL
    negXuinvuu = Xu / -uu
    U2ddot = *pY2 - X12B(*pX1, *pX2, invsym(XX + negXuinvuu * Xu') * (negXuinvuu * (*py1parY2 - beta ' ZY2) + XY2))  // complex expression is Pihat

    u1dddot = u1ddot + U2ddot * (R1invR1R1Y * *pr1 + RparY * beta)
  }
}


// non-WRE stuff that only depends on r in A-R case, for test stat denominators in replication regressions
// since the non-AR OLS code never creates an object for replication regresssions, in that case this is called on the DGP regression object
// depends on results of Estimate() only when doing OLS-style bootstrap on an overidentified IV/GMM regression--score bootstrap or A-R. Then kappa from DGP LIML affects Hessian, pH.
void boottestOLS::InitTestDenoms() {
	real scalar d; pointer (real matrix) scalar pWXAR

  if (parent->scoreBS | parent->robust) {
		if (parent->bootstrapt) {
      if (parent->granular | parent->purerobust) {
        pWXAR = pvHadw(XAR, *parent->pwt)
        WXAR = smatrix(parent->df)
        for (d=parent->df;d;d--)
          WXAR[d].M = (*pWXAR)[,d]
      }

      if (parent->NFE & parent->robust & (parent->FEboot | parent->scoreBS)==0 & parent->granular < parent->NErrClustCombs) {  // make first factor of second term of (64) for c=∩ (c=1)
        if (pWXAR==NULL)
          pWXAR = pvHadw(XAR, *parent->pwt)
        CT_XAR = smatrix(parent->df)
        for (d=parent->df;d;d--)
          CT_XAR[d].M = parent->crosstabFE((*pWXAR)[,d], *parent->pinfoCapData)
      }
    }
	}
}


// partial fixed effects out of a data matrix
pointer(real matrix) scalar boottest::partialFE(pointer(real matrix) scalar pIn) {
	real matrix Out, tmp; real scalar i
	if (NFE & pIn!=NULL) {
		Out = *pIn
		for (i=NFE;i;i--) {
			tmp = Out[FEs[i].is,]
			Out[FEs[i].is,] = tmp :- cross(FEs[i].wt, tmp)
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
	pIDBootData = &.
}

// important to call this when done: break loops in data structure topology to enable garbage collection
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

void boottest::setY2    (real matrix X ) {
	this.pY2  = &X; setdirty(1)
}
void boottest::setX1     (real matrix X ) {
	this.pX1  = &X; setdirty(1)
}
void boottest::setY       (real matrix Y ) {
	this.py1  = &Y; setdirty(1)
}
void boottest::setX2   (real matrix Z ) {
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
void boottest::setR1(real matrix R1, real matrix r1) {
	this.pR1 = &R1; 	this.pr1  = &r1; setdirty(1)
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
	pointer (real rowvector) scalar _pnumer
	if (dirty) boottest()
	if (diststat == "numer") {
		_pnumer = u_sd==1? pnumer : &(*pnumer / u_sd)
		_sort( DistCDR = (*_pnumer)[|2\.|]' :+ *pr , 1)
	} else if (rows(DistCDR)==0)
		if (cols(*pDist) > 1)
			_sort( DistCDR=(*pDist)[|2\.|]' , 1)
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
				p = rowsum(-abs(tmp) :> -abs(*pDist)) / BFeas  // symmetric p value; do so as not to count missing entries in *pDist
			else if (ptype==1)  // equal-tail p value
				p = 2 * min((rowsum(tmp :> *pDist) , rowsum(-tmp:>- *pDist))) / BFeas
			else
				p = rowsum( tmp :>   *pDist) / BFeas  // lower-tailed p value
		} else
				p = rowsum(-tmp :> - *pDist) / BFeas  // upper-tailed p value or p value based on squared stats
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
	return(u_sd == 1? (*pnumer)[,1] : (*pnumer)[,1] / u_sd)
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
	return (Nstar)

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


// efficiently store bootstrap results for one wild weight group in a matrix, one col per replication, handling Nw>1 case (matsizegb() option)
void boottest::storeWtGrpResults(pointer(real matrix) scalar pdest, real scalar w, real matrix content)
	if (Nw==1)
    pdest = &content
  else
    (*pdest)[|., WeightGrpStart[w] \ ., WeightGrpStop[w]|] = content


void boottest::Init() {  // for efficiency when varying r repeatedly to make CI, do stuff once that doesn't depend on r
  real colvector sortID, o, _FEID
  pointer (real colvector) scalar pIDAllData, pIDCapData
	real rowvector ClustCols
	real matrix Combs, tmp
	real scalar i, j, c, minN, sumN, _B, i_FE, sumFEwt
	pragma unset pIDAllData; pragma unset pIDCapData

  Nobs = rows(*pX1)
  NClustVar = cols(*pID)
  kZ  = cols(*pR)
  kX = (kX1 = cols(*pX1)) + (kX2 = cols(*pX2))
  if (kX2 == 0) pX2 = &J(Nobs,0,0)
  if ((kY2 = cols(*pY2)) == 0) pY2 = &J(Nobs,0,0)
	if (LIML & kX2 == kY2) {  // exactly identified LIML = 2SLS
		kappa = 1
		LIML = 0
	}
  REst = rows(*pR1)  // base model contains restrictions?
	if (REst == 0) {
    pR1 = &J(0,kZ,0)
    pr1 = &J(0,1 ,0)
  }
  if (pX2 != NULL) l = kX2 + kX1  // l = # of exogenous variables
  if (kappa==.) kappa = kX2>0  // if kappa in kappa-class estimation not specified, it's 0 or 1 for OLS or 2SLS
  IV = kappa  // "IV" includes 2SLS, LIML, k-class, Fuller
  WRE = (kappa & scoreBS==0) | ARubin
  WREnonARubin = WRE & ARubin==0

  if (haswt = rows(*pwt)>1)
    sumwt = sum(*pwt)
  else
    pwt = &(sumwt = 1)
  _Nobs = haswt & wttype=="fweight"? sumwt : Nobs

  if (WREnonARubin)
    if (NClustVar)
      infoBootData = _panelsetup(*pID, 1..NBootClustVar, pIDBootData)
    else
      pinfoCapData = &(infoBootData = J(Nobs,0,0))  // no clustering, so no collapsing by cluster
  else if (NClustVar)
    if (NClustVar > NBootClustVar)  // bootstrap cluster grouping defs rel to original data
      infoBootData = _panelsetup(*pID, 1..NBootClustVar)
    else
      infoBootData = _panelsetup(*pID, 1..NClustVar)
  else
    pinfoCapData = pinfoAllData = &(infoBootData = J(Nobs,0,0))  // causes no collapsing of data in _panelsum() calls, only multiplying by weights if any
  Nstar = rows(infoBootData)

  if (bootstrapt) {
    if (NClustVar) {
      minN = .; sumN = 0

      Combs = combs(NErrClust)  // represent all error clustering combinations. First is intersection of all error clustering vars
      Clust = structboottestClust(rows(Combs)-1)  // leave out no-cluster combination
      NErrClustCombs = length(Clust)
      subcluster = NClustVar - NErrClust

      if (NClustVar > NBootClustVar)  // info for grouping by intersections of all bootstrap & clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster & intersection of all error clusters
        pinfoAllData = WREnonARubin & granular==0? &_panelsetup(*pID, 1..NClustVar, pIDAllData) : 
                                                   &_panelsetup(*pID, 1..NClustVar            )
      else {
        pinfoAllData = &infoBootData  // info for grouping by intersections of all bootstrap & clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster & intersection of all error clusters
        if (WREnonARubin & granular==0)
          pIDAllData = pIDBootData
      }

      if (NClustVar > NErrClust)  // info for intersections of error clustering wrt data
        pinfoCapData = WREnonARubin & granular==0? &_panelsetup(*pID, subcluster+1..NClustVar, pIDCapData) :
                                                   &_panelsetup(*pID, subcluster+1..NClustVar            )
      else {
        pinfoCapData = pinfoAllData  // info for intersections of error clustering wrt data
        if (WREnonARubin & granular==0)
          pIDCapData = pIDAllData
      }

       IDCap = rows(*pinfoCapData)==Nobs? *pID :   (*pID)[(*pinfoCapData)[,1],]   // version of ID matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
      pIDAll = rows(*pinfoAllData)==Nobs?  pID : &((*pID)[(*pinfoAllData)[,1],])  // version of ID matrix with one row for each all-bootstrap & error cluster-var intersection instead of 1 row for each obs

      BootClust = 2^(NClustVar - NBootClustVar)  // location of bootstrap clustering within list of cluster combinations

			for (c=1; c<=NErrClustCombs; c++) {  // for each error clustering combination
        ClustCols = subcluster :+ _selectindex(Combs[c,])
        Clust[c].even = mod(cols(ClustCols),2)

        if (c == 1)
          if (subcluster) {
            IDCap = IDCap[ Clust.order = stableorder(IDCap, ClustCols), ]
            Clust.info  = _panelsetup(IDCap, ClustCols)
          } else
            Clust.info  = J(rows(*pinfoAllData),0,0)  // causes no collapsing of data in _panelsum() calls
        else {
          if (any(Combs[|c, min(_selectindex(Combs[c,] :!= Combs[c-1,])) \ c,.|])) // if this sort ordering same as last to some point and missing thereafter, no need to re-sort
            IDCap = IDCap[ Clust[c].order = stableorder(IDCap, ClustCols), ]

          Clust[c].info = _panelsetup(IDCap, ClustCols)
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
        ClustShare = haswt? *_panelsum(*pwt, *pinfoCapData)/sumwt : ((*pinfoCapData)[,2]-(*pinfoCapData)[,1]:+ 1)/Nobs // share of observations by group 

    } else {  // if no clustering, cast "robust" as clustering by observation
      Clust = structboottestClust()
      Clust.multiplier = small? _Nobs / (_Nobs - 1) : 1
      Clust.even = 1
      sumN = Clust.N = Nobs
      Clust.info = J(Nobs, 0, 0)  // signals _panelsum not to aggregate
      NErrClustCombs = 1
      if (scoreBS | WREnonARubin==0)
        ClustShare = haswt? *pwt/sumwt : 1/_Nobs
    }

    purerobust = robust & (scoreBS | subcluster)==0 & Nstar==Nobs  // do we ever error-cluster *and* bootstrap-cluster by individual?
    granular   = WREnonARubin? 2*Nobs*B*(2*Nstar+1) < Nstar*(Nstar*Nobs+Clust.N*B*(Nstar+1)) :
		                           NClustVar & scoreBS==0 & (purerobust | (Clust.N+Nstar)*kZ*B + (Clust.N-Nstar)*B + kZ*B < Clust.N*kZ*kZ + Nobs*kZ + Clust.N * Nstar * kZ + Clust.N * Nstar)

    if (robust & purerobust==0) {
      if (subcluster | granular)
        infoErrAll = _panelsetup(*pIDAll, subcluster+1..NClustVar)  // info for error clusters wrt data collapsed to intersections of all bootstrapping & error clusters; used to speed crosstab UXAR wrt bootstrapping cluster & intersection of all error clusterings
      if ((scoreBS & B) | (WREnonARubin & granular==0 & bootstrapt))
        JNcapNstar = J(Clust.N, Nstar, 0)
    }

  if (WREnonARubin & robust & bootstrapt & granular==0) {
  	if (cols(*pIDAllData) == 0) {
      (void) _panelsetup(*pID,            1..NClustVar, pIDAllData)
      (void) _panelsetup(*pID, subcluster+1..NClustVar, pIDCapData)
    }
    IDCTCapstar = infoCTCapstar = smatrix(Nstar)
    for (i=Nstar;i;i--) {
      tmp = (*pIDAllData)[|infoBootData[i,]'|]                           // ID numbers w.r.t. intersection of all bootstrap/error clusterings contained in bootstrap cluster i
      infoCTCapstar[i].M = (*pinfoAllData)[tmp[1]::tmp[rows(tmp)],]  // for each of those ID's, panel info for the all-bootstrap/error-clusterings data row groupings
      IDCTCapstar[i].M = (*pIDCapData)[infoCTCapstar[i].M[,1]]           // ID numbers of those groupings w.r.t. the all-error-clusterings grouping
    }
  }

  } else
    minN = rows(infoBootData)

  if (NFE) {
    sortID = (*pFEID)[o = stableorder(*pFEID, 1)]
    i_FE = 1; FEboot = B>0; j = Nobs; _FEID = J(Nobs, 1, 1)
    invFEwt = J(NFE,1,0)
    FEs = structFE(NFE)
    for (i=Nobs-1;i;i--) {
      if (sortID[i] != sortID[i+1]) {
        FEs[i_FE].is = o[|i+1\j|]
        if (haswt) {
          tmp  = (*pwt)[FEs[i_FE].is]
          FEs[i_FE].wt = tmp / (sumFEwt = colsum(tmp))
        } else
          FEs[i_FE].wt = J(j-i, 1, 1/(sumFEwt = j-i))
        if ((B & robust & granular < NErrClust) | (WREnonARubin & robust & granular & bootstrapt))
          invFEwt[i_FE] = 1 / sumFEwt

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
      FEs[NFE].wt = tmp / (sumFEwt = colsum(tmp))
    } else
      FEs[NFE].wt = J(j-i,1,1/(sumFEwt = j-i))
    if ((B & robust & granular < NErrClust) | (WREnonARubin & robust & granular & bootstrapt))
      invFEwt[NFE] = 1 / sumFEwt
    if (B & FEboot & NClustVar) {  // are all of this FE's obs in same bootstrapping cluster?
      tmp = (*pID)[FEs[NFE].is, 1..NBootClustVar]
      FEboot = all(tmp :== tmp[1,])
    }

    pFEID = &_FEID  // ordinal fixed effect ID

    if (robust & FEboot==0 & granular < NErrClust & B & FEboot==0 & bootstrapt)
      infoBootAll = _panelsetup(*pIDAll, 1..NBootClustVar)  // info for bootstrapping clusters wrt data collapsed to intersections of all bootstrapping & error clusters

		pX1 = partialFE(pX1)
		pX2 = partialFE(pX2)
		py1 = partialFE(py1)
		pY2 = partialFE(pY2)
	}

  if (B & robust & granular & purerobust==0 & WREnonARubin==0)
    if (NFE)
      (void) _panelsetup(*pID   , 1..NBootClustVar, pIDBootData)
    else
      (void) _panelsetup(*pIDAll, 1..NBootClustVar, pIDBootAll )

  if (enumerate = (B & weighttype==0 & Nstar*ln(2) < ln(B)+1e-6))  // generate full Rademacher set?
    MaxMatSize = .

  Nw = MaxMatSize == .? 1 : ceil((B+1) * max((rows(*pIDBootData), rows(*pIDBootAll), Nstar)) * 8 / MaxMatSize / 1.0X+1E) // 1.0X+1E = giga(byte)
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

  if (ML)
    df = rows(*pR)
  else {
    if (ARubin) {
      pR  = &(J(kX2,kX1,0), I(kX2))  // attack surface is all endog vars
      pR1 = &(kX1 & rows(*pR1)? ((*pR1)[|.,.\.,kX1|], J(rows(*pR1),kX2,0)) : J(0, kX, 0))  // and convert model constraints from referring to X1, Y2 to X1, X2
    }
    df = rows(*pR)

    if (WRE==0 & kappa==0) {  // regular OLS
      DGP.parent = Repl.parent = &this
      DGP.LIML = this.LIML; DGP.Fuller = this.Fuller; DGP.kappa = this.kappa
      DGP.SetR (null? *pR1 \ *pR : *pR1)  // DGP constraints: model constraints + null if imposed
      Repl.SetR(*pR1)  // model constraints only
      DGP.InitVars(&Repl.R1perp)
      DGP.InitTestDenoms()
      pM = &DGP  // estimator object from which to get A, AR, XAR

    } else if (ARubin) {

      if (willplot) {  // for plotting/CI purposes get original point estimate since not normally generated
        DGP = boottestIVGMM(); DGP.parent = &this
        DGP.LIML = this.LIML; DGP.Fuller = this.Fuller; DGP.kappa = this.kappa
        DGP.SetR(*pR1, J(0,kZ,0))  // no-null model
        DGP.InitVars()
        DGP.InitEstimate(*pr1)
        DGP.Estimate()
        confpeak = DGP.beta  // estimated coordinate of confidence peak
      }
       
      DGP = boottestARubin(); DGP.parent = &this
      DGP.SetR(*pR1)
      DGP.InitVars()
      DGP.InitTestDenoms()
      pM = &DGP  // estimator object from which to get A, AR, XAR
      kZ = l

    } else if (WREnonARubin) {

      DGP = boottestIVGMM(); DGP.parent = &this
      DGP.LIML = DGP.LIML = kX2!=kY2; DGP.Fuller = 0; DGP.kappa = 1
      DGP.SetR(null? *pR1 \ *pR : *pR1, J(0,kZ,0))  // DGP constraints: model constraints + null if imposed
      DGP.InitVars()
      if (null==0) {  // if not imposing null, then DGP constraints, kappa, Hessian, etc. do not vary with r and can be set now
      	DGP.InitEstimate(*pr1)
        DGP.Estimate()
        DGP.MakeResiduals()
      }

      Repl = boottestIVGMM()
      Repl.parent = &this
      Repl.isDGP = 0
      Repl.LIML = this.LIML; Repl.Fuller = this.Fuller; Repl.kappa = this.kappa
      Repl.SetR(*pR1, *pR)
      Repl.InitVars()
      Repl.InitEstimate(*pr1)

      if (LIML & Repl.kZ==1 & Nw==1) As = betas = J(1, B+1, 0)
      SstaruZperpinvZperpZperp = SstaruZperp = SstaruY = SstaruXinvXX = SstaruX = smatrix(Repl.kZ+1)
      Sstaruu = smatrix(Repl.kZ+1, Repl.kZ+1)
      if (bootstrapt) {
        deltadenom_b = J(Repl.kZ, Repl.kZ, 0)
        SstaruMZperp = SstaruPX = deltadenom = Zyg = SstaruX
        _Jcap = J(Clust.N, Repl.kZ, 0)
        if (granular==0)
          SCTcapuXinvXX = smatrix(Repl.kZ+1, Nstar)
        if (LIML | robust==0) {
          YYstar = YPXYstar = SstaruX
          YYstar_b = YPXYstar_b = J(Repl.kZ+1, Repl.kZ+1, 0)
        }
        if (NFE & (bootstrapt | kappa != 1 | LIML))
          CTFEu = SstaruX
      }

    } else {  // the score bootstrap for IV/GMM uses a IV/GMM DGP but then masquerades as an OLS test because most factors are fixed during the bootstrap. To conform, need DGP and Repl objects with different R, R1, one with FWL, one not

      DGP = boottestIVGMM(); DGP.parent = &this
      DGP.LIML = this.LIML; DGP.Fuller = this.Fuller; DGP.kappa = this.kappa
      DGP.SetR(null? *pR1 \ *pR : *pR1, J(0,kZ,0))  // DGP constraints: model constraints + null if imposed
      DGP.InitVars()
      Repl = boottestIVGMM(); Repl.parent = &this
      Repl.LIML = this.LIML; Repl.Fuller = this.Fuller; Repl.kappa = this.kappa
      Repl.SetR(*pR1, I(kZ))  // process replication restraints = model constraints only
      Repl.InitVars(&Repl.R1perp)
      Repl.InitEstimate(*pr1)
      Repl.Estimate()  // bit inefficient to estimate in both objects, but maintains the conformity
      Repl.InitTestDenoms()
      pM = &Repl  // estimator object from which to get A, AR, XAR; DGP follows WRE convention of using FWL, Repl follows OLS convention of not; scoreBS for IV/GMM mixes the two
      if (null==0) {  // if not imposing null, then DGP constraints, kappa, Hessian, etc. do not vary with r and can be set now
      	DGP.InitEstimate(*pr1)
        DGP.Estimate()
        DGP.MakeResiduals()
      }
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

  if (small) df_r = NClustVar? minN - 1 : _Nobs - kZ - NFE

  if (df==1) setsqrt(1)  // work with t/z stats instead of F/chi2

  if (small)
    multiplier = (smallsample = (_Nobs - kZ - NFE) / (_Nobs - robust)) / df  // divide by # of constraints because F stat is so defined
  else
    multiplier = smallsample = 1

  if ((robust | ML)==0)
    multiplier = multiplier * _Nobs  // will turn sum of squared errors in denom of t/z into mean
  if (sqrt) multiplier = sqrt(multiplier)

  if (bootstrapt & (WREnonARubin | df>1 | MaxMatSize<.)) // unless nonWRE or df=1 or splitting weight matrix, code will create Dist element-by-element, so pre-allocate vector now
    pDist = &J(1, B+1, .)
  if (Nw>1 | WREnonARubin | (null==0 & df<=2))
    pnumer = &J(df, B+1, .)

  if (WREnonARubin==0)
    if (interpolable = B & null & Nw==1 & (kappa==0 | ARubin)) {    
      dnumerdr = smatrix(q)
      if (interpolate_u = (robust | ML)==0) dudr = dnumerdr
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
      (*pnumer)[,1] = *pR * betas[1] - *pr
    else if (ARubin) {
      DGP.InitEstimate(*pr)
      DGP.Estimate()
      DGP.MakeResiduals()
      (*pnumer)[,1] = u_sd * Repl.beta[|kX1+1\.|] // coefficients on excluded instruments in ARubin OLS
    } else
      (*pnumer)[,1] = u_sd * (*pR * (ML? beta : pM->beta) - *pr) // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this

    (*pDist)[1] = df==1? (*pnumer)[1] / sqrt(statDenom) : (*pnumer)[,1] ' invsym(statDenom) * (*pnumer)[,1]
    return
  }

	if (Nw > 1) {
		rseed(seed)
    makeWildWeights(WeightGrpStop[1] - 1, 1)
  }
	
  if (WREnonARubin)
    PrepWRE()
   else
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

  BFeas = (*pDist)[1]==.? 0 : rownonmissing(*pDist) - 1
	DistCDR = J(0,0,0)
	setdirty(0)
	initialized = 1
}

// compute bootstrap-c denominator from all bootstrap numerators
void boottest::makeBootstrapcDenom(real scalar w) {
	real colvector tmp
  if (w == 1) {
		tmp = (*pnumer)[,1]
    statDenom = *pnumer * *pnumer' - tmp * tmp'
		numersum = rowsum(*pnumer) - tmp
	} else {
		statDenom = statDenom + *pnumer * *pnumer'
		numersum = numersum + rowsum(*pnumer)
	}
	if (w == Nw) {  // last weight group?
		statDenom = (statDenom - numersum * numersum' / B) / B
		pDist = &(sqrt? *pnumer:/sqrt(statDenom) : colsum(*pnumer :* invsym(statDenom) * *pnumer))
	}
}

// draw wild weight matrix of width _B. If first=1, insert column of 1s at front. For non-Anderson-Rubin WRE, subtract 1 from all weights
void boottest::makeWildWeights(real scalar _B, real scalar first) {

	if (_B /*| WREnonARubin*/) {  // in scoretest or waldtest WRE, still make v a col of 1's
		if (enumerate)
			v = J(Nstar,1,1), count_binary(Nstar, -1-WREnonARubin, 1-WREnonARubin)  // complete Rademacher set
		else if (weighttype==3)
			v = rnormal(Nstar, _B+1, -WREnonARubin, 1)  // normal weights
		else if (weighttype==4)
			v = rgamma(Nstar, _B+1, 4, .5) :- (2 + WREnonARubin)  // Gamma weights
		else if (weighttype==2)
			if (WREnonARubin)
				v = sqrt(2 * ceil(runiform(Nstar, _B+first) * 3)) :* ((runiform(Nstar, _B+1):>=.5):-.5) :- 1  // Webb weights, minus 1 for convenience in WRE
			else {
				v = sqrt(    ceil(runiform(Nstar, _B+first) * 3)) :* ((runiform(Nstar, _B+1):>=.5):-.5)       // Webb weights, divided by sqrt(2)
				u_sd = 1.6a09e667f3bcdX-001 /*sqrt(.5)*/
			}
		else if (weighttype)
			if (WREnonARubin)
				v = ( rdiscrete(Nstar, _B+first, 1.727c9716ffb76X-001\1.1b06d1d200914X-002 /*.5+sqrt(.05)\.5-sqrt(.05)*/) :- 1.5 ) * 1.1e3779b97f4a8X+001 /*sqrt(5)*/ :- .5  // Mammen weights, minus 1 for convenience in WRE
			else {
				v = ( rdiscrete(Nstar, _B+first, 1.727c9716ffb76X-001\1.1b06d1d200914X-002 /*.5+sqrt(.05)\.5-sqrt(.05)*/) :- 1.5 ) :+ 1.c9f25c5bfedd9X-003 /*.5/sqrt(5)*/  // Mammen weights, divided by sqrt(5)
				u_sd = 1.c9f25c5bfedd9X-002 /*sqrt(.2)*/
			}
		else if (WREnonARubin) {
      v = runiform(Nstar, _B+first) :<  .5; v = (-2) * v  // Rademacher weights, minus 1 for convenience in WRE
    } else {
      v = runiform(Nstar, _B+first) :>= .5; v = v :- .5   // Rademacher weights, divided by 2
      u_sd = .5
    }

		if (first)
			v[,1] = J(Nstar, 1, WREnonARubin? 0 : u_sd)  // keep original residuals in first entry to compute base model stat		
	} else
		v = J(0,1,0)  // in places, cols(v) indicates number of B -- 1 for classical tests
}


// For WRE, and with reference to Y = [y1 Z], given 0-based columns indexes within it, ind1, ind2, return all bootstrap realizations of Y[,ind1]'((1-kappa)*M_Zperp-kappa*M_Xpar)*Y[,ind2] for kappa constant across replications
// ind1 can be a rowvector
// (only really the Hessian when we narrow Y to Z)
real matrix boottest::HessianFixedkappa(real rowvector ind1, real scalar ind2, real scalar kappa) {
	real scalar i; real matrix retval, T2, MZperpTerm; pointer (real colvector) scalar pT1L, pT1R

	if (kappa) {
		pT1R = ind2? pcol(Repl.invXXXZ,ind2) : &Repl.invXXXy1par
		if (Repl.Yendog[ind2+1])
			pT1R = &(*pT1R :+ SstaruXinvXX[ind2+1].M * v)  // right-side linear term
		if (cols(ind1) == 1) {
			pT1L = ind1? pcol(Repl.XZ,ind1) : Repl.pXy1par
			if (Repl.Yendog[ind1+1])
				pT1L = &(*pT1L :+ SstaruX[ind1+1].M * v)
			retval = colsum(*pT1L :* *pT1R)  // multiply in the left-side linear term
		} else {
			retval = J(cols(ind1), cols(v), 0)
			for (i=cols(ind1);i;i--) {
				pT1L = ind1[i]? pcol(Repl.XZ,ind1[i]) : Repl.pXy1par
				if (Repl.Yendog[ind1[i]+1])
					pT1L = &(*pT1L :+ SstaruX[ind1[i]+1].M * v)
				retval[i,] = colsum(*pT1L :* *pT1R)  // multiply in the left-side linear terms for each ind1
			}
		}
	}

	if (kappa != 1) {
		if (cols(ind1) > 1)
			MZperpTerm = J(cols(ind1),cols(v),0)
		for (i=cols(ind1);i;i--)  // 
			if (Repl.Yendog[ind1[i]+1]) {
				T2 = SstaruZperpinvZperpZperp[ind1[i]+1].M ' SstaruZperp[ind2+1].M  // quadratic term
				_diag(T2, diagonal(T2) - (ind1[i] <= ind2? Sstaruu[ind2+1, ind1[i]+1].M : Sstaruu[ind1[i]+1, ind2+1].M))  // minus diagonal crosstab

        if (NFE)
          T2 = T2 + CTFEu[ind1[i]+1].M ' (invFEwt :* CTFEu[ind2+1].M)

				if (cols(ind1) > 1)
					MZperpTerm[i,] = (*pcol(SstaruY[ind2+1].M, ind1[i]+1) + *pcol(SstaruY[ind1[i]+1].M, ind2+1)) ' v - colsum(v :* T2 * v)
				else
					MZperpTerm     = (*pcol(SstaruY[ind2+1].M, ind1[i]+1) + *pcol(SstaruY[ind1[i]+1].M, ind2+1)) ' v - colsum(v :* T2 * v)
			}
			
		retval = kappa? kappa :* retval + (1 - kappa) :* (MZperpTerm :+ Repl.YY[ind1:+1,ind2+1]) :
		                                                  MZperpTerm :+ Repl.YY[ind1:+1,ind2+1]
	}
	return(retval)
}


// Workhorse for WRE CRVE sandwich filling
// With reference to notional Y = [y1 Z], given 0-based columns index within it, ind1, and a matrix betas of all the boostrap estimates, return all bootstrap realizations of P_X * Y[,ind1]_g ' u\hat_1g^*b
// for all groups in the intersection of all error clusterings
// return value has one row per cap cluster, one col per bootstrap replication
pointer(real matrix) scalar boottest::Filling(real scalar ind1, real matrix betas) {
	real scalar i, ind2; real matrix retval, T1; pointer (real matrix) scalar pbetav; pointer (real colvector) pPXYstar; real rowvector _beta; real colvector S
	pragma unset retval

	if (granular) {
		if (Nw == 1) {  // create or avoid NxB matrix?
			pPXYstar = ind1? pcol(Repl.PXZ, ind1) : &Repl.PXy1
			if (Repl.Yendog[ind1+1])
				pPXYstar = &(*pPXYstar :+ SstaruPX[ind1+1].M * v)

			retval = *_panelsum(*pPXYstar :* (*Repl.py1 :- SstaruMZperp.M * v), *pwt, *pinfoCapData)

			for (ind2=Repl.kZ;ind2;ind2--) {
				_beta = -betas[ind2,]
				retval = retval + *_panelsum(*pPXYstar :* (Repl.Yendog[ind2+1]? *pcol(Repl.Z,ind2) * _beta :- SstaruMZperp[ind2+1].M * (v :* _beta) :
				                                                                *pcol(Repl.Z,ind2) * _beta                                           ), *pwt, *pinfoCapData)
			}
		} else { // create pieces of each N x B matrix one at a time rather than whole thing at once
			retval = J(Clust.N, cols(v), 0)
			for (ind2=0; ind2<=Repl.kZ; ind2++) {
				if (ind2)
					pbetav = &(v :* (_beta = -betas[ind2,]))

        if (purerobust) {
          for (i=Clust.N;i;i--) {
            pPXYstar = ind1? &Repl.PXZ[i,ind1] : &Repl.PXy1[i,]
            if (Repl.Yendog[ind1+1])
              pPXYstar = &(*pPXYstar :+ SstaruPX[ind1+1].M[i,] * v)

            if (ind2)
              retval[i,] = retval[i,] + cross(*pwt, *pPXYstar :* (Repl.Yendog[ind2+1]? Repl.Z[i,ind2] * _beta :- SstaruMZperp[ind2+1].M[i,] * *pbetav :
                                                                                       Repl.Z[i,ind2] * _beta                                           ))
            else
              retval[i,] =              cross(*pwt, *pPXYstar :* ((*Repl.py1)[i] :- SstaruMZperp.M[i,] * v))
          }
        } else {
          for (i=Clust.N;i;i--) {
            S = (*pinfoCapData)[i,]'
            pPXYstar = ind1? &Repl.PXZ[|S,(ind1\ind1)|] : &Repl.PXy1[|S|]
            if (Repl.Yendog[ind1+1])
              pPXYstar = &(*pPXYstar :+ SstaruPX[ind1+1].M[|S,(.\.)|] * v)

            if (ind2)
              retval[i,] = retval[i,] + cross(*pwt, *pPXYstar :* (Repl.Yendog[ind2+1]? Repl.Z[|S,(ind2\ind2)|] * _beta :- SstaruMZperp[ind2+1].M[|S,(.\.)|] * *pbetav :
                                                                                       Repl.Z[|S,(ind2\ind2)|] * _beta                                                 ))
            else
              retval[i,] =              cross(*pwt, *pPXYstar :* ((*Repl.py1)[|S|] :- SstaruMZperp.M[|S,(.\.)|] * v))
          }
        }
			}
		}
	} else {  // non-granular
		for (ind2=0; ind2<=Repl.kZ; ind2++) {
			pbetav = ind2? &(v :* (_beta = -betas[ind2,])) : &v

			if (Repl.Yendog[ind1+1])
				T1 = Repl.SMZperpYX[ind2+1].M * SstaruXinvXX[ind1+1].M

			if (Repl.Yendog[ind2+1]) {  // add CT of U2 * PXY1
				if (NClustVar == NBootClustVar & !subcluster)  // simple case of one clustering: full crosstab is diagonal
					if (cols(T1))
						_diag(T1, diagonal(T1) + SstaruXinvXX[ind2+1].M ' (ind1? *pcol(Repl.XZ,ind1) : *Repl.pXy1par))
					else
						T1 =                     SstaruXinvXX[ind2+1].M ' (ind1? *pcol(Repl.XZ,ind1) : *Repl.pXy1par)  // keep T1 as vector if it's just going to be a diagonal matrix
				else {
					if (Repl.Yendog[ind1+1]==0)
						T1 = JNcapNstar
					for (i=Nstar;i;i--)
						T1[IDCTCapstar[i].M, i] = T1[IDCTCapstar[i].M, i] + SCTcapuXinvXX[ind2+1,i].M * *(ind1? pcol(Repl.XZ,ind1) : Repl.pXy1par)
				}
				if (cols(Repl.Zperp))
					T1 = T1 - Repl.SPXYZperp[ind1+1].M * SstaruZperpinvZperpZperp[ind2+1].M
        if (NFE)
	        T1 = T1 - Repl.CT_FEcapPY[ind1+1].M ' CTFEu[ind2+1].M
			}

			retval = ind2? retval + Repl.FillingT0[ind1+1,ind2+1].M * _beta  + (cols(T1)==1? T1 :* *pbetav : T1 * *pbetav) :   // - x*beta components
															Repl.FillingT0[ind1+1,     1].M         :+ (cols(T1)==1? T1 :* *pbetav : T1 *       v)     // y component
			if (Repl.Yendog[ind1+1] & Repl.Yendog[ind2+1])
				for (i=Clust.N;i;i--) {
					S = (*pinfoCapData)[i,]', (.\.)
					retval[i,] = retval[i,] - colsum(v :* cross(SstaruPX[ind1+1].M[|S|], haswt? (*pwt)[|S|] : 1, SstaruMZperp[ind2+1].M[|S|]) * *pbetav)  // *** efficient, or move Xg over?
				}		
		}
	}
	return(&retval)
}


void boottest::PrepWRE() {
  real scalar i, j, g; real matrix SuX1g, SuX2g; pointer (real colvector) scalar puwt

	DGP.InitEstimate(null? *pr1 \ *pr : *pr1)
	DGP.Estimate()
  DGP.MakeResiduals()
	U2parddot = DGP.U2ddot * Repl.RparY

	for (i=Repl.kZ; i>=0; i--) {  // precompute various clusterwise sums
		puwt = pvHadw(i? *pcol(U2parddot,i) : DGP.u1dddot, *pwt)

    // S_star(u :* X), S_star(u :* Zperp) for residuals u for each endog var; store transposed
    SstaruX                   [i+1].M = *_panelsum2(*Repl.pX1, *Repl.pX2, *puwt, infoBootData)'
    SstaruXinvXX              [i+1].M = *Repl.pinvXX * SstaruX[i+1].M

    if (bootstrapt | kappa!=1 | LIML) {
      SstaruZperp             [i+1].M = *_panelsum(Repl.Zperp, *puwt, infoBootData)'
      SstaruZperpinvZperpZperp[i+1].M = Repl.invZperpZperp * SstaruZperp[i+1].M
    }

    if (kappa != 1 | LIML | robust==0) {
      SstaruY[i+1].M = *_panelsum2(*Repl.py1par, Repl.Z, *puwt, infoBootData)
      for (j=i; j>=0; j--)
        Sstaruu[i+1,j+1].M = *_panelsum(j? *pcol(U2parddot,j) : DGP.u1dddot, *puwt, infoBootData)
    }

    if (NFE & (bootstrapt | kappa != 1 | LIML))
      CTFEu[i+1].M = crosstabFE(*puwt, infoBootData)

    if (robust & bootstrapt) {
      if (granular==0)  // Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u:*X and u:*Zperp (and times invXX or invZperpZperp)
        for (g=Nstar;g;g--) {
          SuX1g = *_panelsum(*DGP.pX1, *puwt, infoCTCapstar[g].M)
          SuX2g = *_panelsum(*DGP.pX2, *puwt, infoCTCapstar[g].M)
          SCTcapuXinvXX[i+1,g].M = Repl.kX1? SuX1g * Repl.invXX1 + SuX2g * *Repl.pinvXX2 : 
                                                                   SuX2g * *Repl.pinvXX
        }
      
      SstaruPX                [i+1].M = Repl.XinvXX        * SstaruX                 [i+1].M
      SstaruMZperp            [i+1].M = Repl.Zperp         * SstaruZperpinvZperpZperp[i+1].M
      if (Nobs == Nstar)  // subtract "crosstab" of observation by cap-group of u
        _diag(SstaruMZperp[i+1].M, diagonal(SstaruMZperp[i+1].M) - (i? U2parddot[,i] : DGP.u1dddot))  // case: bootstrapping by observation
      else
        for (g=Nobs;g;g--)
          SstaruMZperp[i+1].M[g,(*pIDBootData)[g]] = SstaruMZperp[i+1].M[g,(*pIDBootData)[g]] - (i? U2parddot[g,i] : DGP.u1dddot[g])
      if (NFE)
        SstaruMZperp[i+1].M = SstaruMZperp[i+1].M + (invFEwt :* CTFEu[i+1].M)[*pFEID,]
    }
	}
}

void boottest::makeWREStats(real scalar w) {
	real scalar c, b, i
	real colvector numer_b
	real rowvector numerw, val
	real matrix deltanumer, Jcap, J_b, Jcaps, vec
	struct smatrix rowvector A
	pragma unset vec; pragma unset val

	if (Repl.kZ == 1) {  // optimized code for 1 coefficient in bootstrap regression
		if (LIML) {
			if (Nw>1) As = betas = J(1, cols(v), 0)
			YYstar.M      = HessianFixedkappa(0   , 0, 0)  // kappa=0 => Y*MZperp*Y
			YPXYstar.M    = HessianFixedkappa(0..1, 0, 1)  // kappa=1 => Y*PXpar*Y
			YYstar  [2].M = HessianFixedkappa(0..1, 1, 0)  // kappa=0 => Y*MZperp*Y
			YPXYstar[2].M = HessianFixedkappa(1   , 1, 1)  // kappa=1 => Y*PXpar*Y
			for (b=cols(v); b; b--) {
				YYstar_b  [1,1] = YYstar.M     [ b]  // fill uppper triangle, which is all that invsym() looks at
				YYstar_b   [,2] = YYstar  [2].M[,b]  // fill uppper triangle, which is all that invsym() looks at
				YPXYstar_b [,1] = YPXYstar.M   [,b]  // fill lower triangle to prepare for _makesymmetric()
				YPXYstar_b[2,2] = YPXYstar[2].M[ b]  // fill lower triangle to prepare for _makesymmetric()
				YPXYstar_b[1,2] = YPXYstar_b[2,1]
				eigensystemselecti(invsym(YYstar_b) * YPXYstar_b, 2\2, vec, val)
				kappa = 1/(1 - Re(val))
				if (Fuller) kappa = kappa - 1 / (_Nobs - l)
        vec = kappa * YPXYstar_b[,2] + (1 - kappa) * YYstar_b[,2]
				betas[b] = vec[1] / (As[b] = vec[2])
			}
		} else
			betas = HessianFixedkappa(1, 0, kappa) :/ (As = HessianFixedkappa(1, 1, kappa))

    if (null)
 			numerw      = betas   :+ (Repl.Rt1 - *pr) / Repl.RRpar
		else {
			numerw = betas :- DGP.beta0
			if (w==1)
				numerw[1] = betas[1] + (Repl.Rt1 - *pr) / Repl.RRpar
		}
    storeWtGrpResults(pnumer, w, numerw)

		if (bootstrapt) {
			if (robust) {
				Jcaps = *Filling(1, betas) :/ As
				for (c=1; c<=NErrClustCombs; c++) {  // sum sandwich over error clusterings
					if (NClustVar != 1 & rows(Clust[c].order))
						Jcaps = Jcaps[Clust[c].order,]
					Jcap = *_panelsum(Jcaps, Clust[c].info)
					_clustAccum(denom.M, c, colsum(Jcap:*Jcap))
				}
			} else
        denom.M = (HessianFixedkappa(0,0,0) - 2 * betas :* HessianFixedkappa(0, 1, 0) + betas:*betas :* HessianFixedkappa(1, 1, 0)) / _Nobs :/ As  // classical error variance

      storeWtGrpResults(pDist, w, sqrt? numerw:/sqrt(denom.M) : numerw :* numerw :/ denom.M)
			denom.M = Repl.RRpar * Repl.RRpar * denom.M[1]
		}
	} else {  // WRE bootstrap for more than 1 coefficeint in bootstrap regression

		betas = J(Repl.kZ, cols(v), 0)
		A = smatrix(cols(v))

		if (LIML) {
			for (i=Repl.kZ;i>=0;i--) {
				YYstar  [i+1].M = HessianFixedkappa(0..i      , i, 0)  // kappa=0 => Y*MZperp*Y
				YPXYstar[i+1].M = HessianFixedkappa(i..Repl.kZ, i, 1)  // kappa=1 => Y*PXpar*Y
			}
			for (b=cols(v); b; b--) {
				for (i=Repl.kZ;i>=0;i--) {
					YYstar_b  [|.  ,i+1\      i+1,i+1|] = YYstar  [i+1].M[,b]  // fill uppper triangle, which is all that invsym() looks at
					YPXYstar_b[|i+1,i+1\Repl.kZ+1,i+1|] = YPXYstar[i+1].M[,b]  // fill lower triangle to prepare for _makesymmetric()
				}
				_makesymmetric(YPXYstar_b)
				eigensystemselecti(invsym(YYstar_b) * YPXYstar_b, Repl.kZ+1\Repl.kZ+1, vec, val)
				kappa = 1/(1 - Re(val)) // sometimes a tiny imaginary component sneaks into val
				if (Fuller) kappa = kappa - 1 / (_Nobs - l)
				betas[,b] = (A[b].M = invsym(kappa*YPXYstar_b[|2,2\.,.|] + (1-kappa)*YYstar_b[|2,2\.,.|])) * (kappa*YPXYstar_b[|2,1\.,1|] + (1-kappa)*YYstar_b[|1,2\1,.|]')
			}
		} else {
			deltanumer = HessianFixedkappa(1..Repl.kZ, 0, kappa)

			for (i=Repl.kZ;i;i--)
				deltadenom[i].M = HessianFixedkappa(1..i, i, kappa)

			for (b=cols(v); b; b--) {
				for (i=Repl.kZ;i;i--)
					deltadenom_b[|.,i\i,i|] = deltadenom[i].M[,b]  // fill uppper triangle, which is all that invsym() looks at
				betas[,b]  = (A[b].M = invsym(deltadenom_b)) * deltanumer[,b]
			}
		}
		
		if (bootstrapt)
      if (robust)
        for(i=Repl.kZ;i;i--)
          Zyg[i].M = *Filling(i, betas)
      else
        for (i=Repl.kZ;i>=0;i--)
          YYstar[i+1].M = HessianFixedkappa(i..Repl.kZ, i, 0)  // kappa=0 => Y*MZperp*Y

		for (b=cols(v); b; b--) {
			numer_b = null | w==1 & b==1? (Repl.RRpar * betas[,b] + Repl.Rt1) - *pr : Repl.RRpar * (betas[,b] - DGP.beta0)

			if (bootstrapt) {
				if (robust) {  // Compute denominator for this WRE test stat
          for(i=Repl.kZ;i;i--)
            _Jcap[,i] = Zyg[i].M[,b]
          Jcap = _Jcap * (A[b].M * Repl.RRpar')

					for (c=1; c<=NErrClustCombs; c++) {
						if (NClustVar != 1 & rows(Clust[c].order))
							Jcap = Jcap[Clust[c].order,]
						J_b = *_panelsum(Jcap, Clust[c].info)
						_clustAccum(denom.M, c, cross(J_b,J_b))
					}
				} else {  // non-robust
          for (i=Repl.kZ;i>=0;i--)
            YYstar_b[|i+1,i+1\Repl.kZ+1,i+1|] = YYstar[i+1].M[,b]  // fill lower triangle to prepare for makesymmetric()
          denom.M = (Repl.RRpar * A[b].M * Repl.RRpar') * ((-1 \ betas[,b]) ' makesymmetric(YYstar_b) * (-1 \ betas[,b]) / _Nobs)  // 2nd half is sig2 of errors
        }
				(*pDist)[b+WeightGrpStart[w]-1] = sqrt? numer_b/sqrt(denom.M) : cross(numer_b, invsym(denom.M) * numer_b)  // hand-code for 2-dimensional?
			}
			(*pnumer)[,b+WeightGrpStart[w]-1] = numer_b  // slight inefficiency: in usual bootstrap-t case, only need to save numerators in numer if getdist("numer") is coming because of svmat(numer)
		}
	}

	if (w==1 & bootstrapt) statDenom = denom.M  // original-sample denominator
}


// Construct stuff that depends linearly or quadratically on r, possibly by interpolation
void boottest::makeInterpolables() {
	real scalar h1, h2, d1, d2, c; real matrix tmp; real colvector Delta, newPole

  if (interpolable) {
    if (rows(anchor)==0) {  // first call? save current r as permanent anchor for interpolation
      _makeInterpolables(anchor = *pr)
      numer0 = *pnumer
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

          dnumerdr[h1].M = (*pnumer - numer0) / poles[h1]
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
      numerw   =   numer0 + dnumerdr.M * Delta[1] ; if (q > 1) numerw =    numerw + dnumerdr[2].M * Delta[2]
      if (interpolate_u) {
        puddot = &(uddot0 +     dudr.M * Delta[1]); if (q > 1) puddot = &(*puddot + dudr    [2].M * Delta[2])
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
  real scalar d, c; pointer (real matrix) scalar pustarXAR, ptmp

  if (ML)
		uXAR = *pSc * (AR = *pA * *pR')
	else {
    if (ARubin)
      DGP.InitEstimate(r)
    else if (kappa) {
      if (null) { // in score bootstrap for IV/GMM, if imposing null, then DGP constraints, kappa, Hessian, etc. do vary with r and must be set now
      	DGP.InitEstimate(*pr1 \ r)
        DGP.InitTestDenoms()
      }
    } else  // regular OLS
    	DGP.InitEstimate(null? *pr1 \ r : *pr1)

    DGP.Estimate()
    DGP.MakeResiduals()
    puddot = &DGP.u1ddot

		if (scoreBS | (robust & granular < NErrClustCombs))
      uXAR = *puddot :* pM->XAR
  }

  SuwtXA = scoreBS?
	            (B? 
		             (NClustVar? *_panelsum(uXAR, *pwt, infoBootData) : 
					                   *pvHadw(uXAR, *pwt)                  ) :
				         cross(*pwt, uXAR)')                             :
              *pM->pA * *_panelsum2(*pX1, *pX2, *pvHadw(*puddot, *pwt), infoBootData)'  // same calc as in score BS but broken apart to grab intermediate stuff, and assuming residuals defined; X2 empty except in Anderson-Rubin

  if (robust & granular < NErrClustCombs) {
    pustarXAR = _panelsum(uXAR, *pwt, *pinfoAllData)  // collapse data to all-boot & error-cluster-var intersections. If no collapsing needed, _panelsum() will still fold in any weights
    if (B) {
      if (scoreBS)
        for (d=df;d;d--)
          Kd[d].M = JNcapNstar  // inefficient, but not optimizing for the score bootstrap
      else
        for (d=df;d;d--)
          Kd[d].M = *_panelsum2(*pX1, *pX2, *pvHadw(*pcol(pM->XAR,d), *pwt), *pinfoCapData) * SuwtXA  // final term in (64), for c=intersection of all error clusters

      if (NFE & FEboot==0)
        CT_WE = crosstabFE(*pwt :* *puddot, infoBootData)

			for (d=df;d;d--) {  // subtract crosstab of u:*XAR wrt bootstrapping cluster combo and all-cluster-var intersections
				crosstabCapstarMinus(Kd[d].M, *pcol(*pustarXAR,d))
        if (NFE & FEboot==0)
          Kd[d].M = Kd[d].M + pM->CT_XAR[d].M ' (invFEwt :* CT_WE)  // middle term of (64)
        if (scoreBS)
					Kd[d].M = Kd[d].M - ClustShare * colsum(Kd[d].M) // recenter
			}

      for (c=1+granular; c<=NErrClustCombs; c++) {
        if (rows(Clust[c].order))
          for (d=df;d;d--)
            Kd[d].M = Kd[d].M[Clust[c].order,]
        for (d=df;d;d--)
          Kcd[c,d].M = *_panelsum(Kd[d].M, Clust[c].info)
      }
    } else {  // B = 0. In this case, only 1st term of (64) is non-zero after multiplying by v* (= all 1's), and it is then a one-way sum by c

      if (scoreBS)
        pustarXAR = &(*pustarXAR :- ClustShare * colsum(*pustarXAR))  // recenter if OLS

      for (c=1; c<=NErrClustCombs; c++) {
        if (rows(Clust[c].order))
          pustarXAR = &((*pustarXAR)[Clust[c].order,])
        ptmp = _panelsum(*pustarXAR, Clust[c].info)
        for (d=df;d;d--)
          Kcd[c,d].M = *pcol(*ptmp,d)
      }
    }
  }

  makeNumerAndJ(1, r)  // compute J = kappa * v; if Nw > 1, then this is for 1st group; if interpolating, it is only group, and may be needed now to prep interpolation
}

// compute stuff depending linearly on v, needed to prep for interpolation
void boottest::makeNumerAndJ(real scalar w, | real colvector r) {  // called to *prepare* interpolation, or when w>1, in which case there is no interpolation
  real scalar c, d

  numerw = scoreBS?
             (B? 
               cross(SuwtXA, v) : 
               SuwtXA * u_sd    ) :
             (robust==0 | granular | purerobust?
                *pR * (betadev = SuwtXA * v) :
               (*pR * SuwtXA) * v)

  if (w==1) {
  	if      ( ARubin) numerw[,1] = u_sd * pM->beta[|kX1+1\.|]  // coefficients on excluded instruments in ARubin OLS
    else if (null==0) numerw[,1] = u_sd * (*pR * (ML? beta : pM->beta) - r)  // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this.
  }
  storeWtGrpResults(pnumer, w, numerw)

	if (B & robust) {
    if (granular | purerobust)  // optimized treatment when bootstrapping by many/small groups
      if (purerobust)
        ustar = *partialFE(&(*puddot :* v)) - X12B(*pX1, *pX2, betadev)
      else {  // clusters small but not all singletons
        if (NFE) {
          ustar = *partialFE(&(*puddot :* v[*pIDBootData,]))
          for (d=df;d;d--)
            (*pJcd)[1,d].M = *_panelsum(ustar, pM->WXAR[d].M, *pinfoCapData)                                              - *_panelsum2(*pX1, *pX2, pM->WXAR[d].M, *pinfoCapData) * betadev
        } else 
          for (d=df;d;d--)
	          (*pJcd)[1,d].M = *_panelsum(*_panelsum(*puddot, pM->WXAR[d].M, *pinfoAllData) :* v[*pIDBootAll,], infoErrAll) - *_panelsum2(*pX1, *pX2, pM->WXAR[d].M, *pinfoCapData) * betadev
      }

		for (c=NErrClustCombs; c>granular; c--)
      for (d=df;d;d--)
        (*pJcd)[c,d].M = Kcd[c,d].M * v
  }
}

void boottest::makeNonWREStats(real scalar w) {
	real scalar i, c, j, l; real matrix ustar2, tmp; real colvector numer_l; pointer (real matrix) scalar pAR; real rowvector t1, t2, t12

  if (w > 1) makeNumerAndJ(w)

	if (robust) {
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
    	storeWtGrpResults(pDist, w, numerw :/ sqrt(denom.M))
			if (w==1)
				statDenom = denom.M[1]  // original-sample denominator
		} else if (df==2) {  // hand-code 2D numer'inv(denom)*numer
    	t1 = numerw[1,]; t2 = numerw[2,]; t12 = t1:*t2
			storeWtGrpResults(pDist, w, (t1:*t1:*denom[2,2].M - (t12+t12):*denom[2,1].M + t2:*t2:*denom[1,1].M) :/ (denom[1,1].M:*denom[2,2].M - denom[2,1].M:*denom[2,1].M))
			if (w==1)
				statDenom = denom[1,1].M[1], denom[2,1].M[1] \ denom[2,1].M[1], denom[2,2].M[1]  // original-sample denominator
    } else {  // build each replication's denominator from vectors that hold values for each position in denominator, all replications
			tmp = J(df,df,0)
			for (l=cols(v); l; l--) {
				for (i=df;i;i--)
					for (j=i;j;j--)
						tmp[j,i] = denom[i,j].M[l]  // fill upper triangle, which is all invsym() looks at
				numer_l = numerw[,l]
				(*pDist)[l+WeightGrpStart[w]-1] = numer_l ' invsym(tmp) * numer_l  // in degenerate cases, cross() would turn cross(.,.) into 0
			}
			if (w==1)
				statDenom = tmp  // original-sample denominator
		}

	} else { // non-robust

		pAR = ML? &AR : &pM->AR
		if (df == 1) {  // optimize for one null constraint
			denom.M = *pR * *pAR

			if (ML==0) {
                     ustar = B? v :* *puddot : *puddot
				if (scoreBS) ustar = ustar :- (haswt? cross(ClustShare, ustar) : colsum(ustar) * ClustShare)  // Center variance if interpolated
				        else ustar = ustar  - X12B(*pX1, *pX2, betadev)  // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
				denom.M = denom.M :* (haswt? cross(*pwt, ustar :* ustar) : colsum(ustar :* ustar))
			}
			storeWtGrpResults(pDist, w,  numerw :/ sqrt(denom.M))
			if (w==1)
				statDenom = denom.M[1]  // original-sample denominator
		} else {
			denom.M = *pR * *pAR

			if (ML) {
				for (l=cols(v); l; l--) {
					numer_l = numerw[,l]
					(*pDist)[l+WeightGrpStart[w]-1] = cross(numer_l, invsym(denom.M), numer_l)
				}
				if (w==1)
					statDenom = denom.M  // original-sample denominator
			} else {
				for (l=cols(v); l; l--) {
					numer_l = numerw[,l]
					(*pDist)[l+WeightGrpStart[w]-1] = cross(numer_l, invsym(denom.M) * numer_l)
                       ustar = B? v[,l] :* *puddot : *puddot
					if (scoreBS) ustar = ustar :- (haswt? cross(*pwt, ustar) : colsum(ustar)) * ClustShare  // Center variance if interpolated
					        else ustar = ustar  - X12B(*pX1, *pX2, betadev[,l])  // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
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
real matrix _panelsetup(real matrix X, real rowvector cols, | pointer(real colvector) scalar pID) {
	real matrix info; real scalar i, N; real scalar p; real rowvector tmp, id
	N = rows(X)
	info = J(N, 2, N); if (args()>2) pID = &J(N, 1, 1)
	info[1,1] = p = 1
	id = X[1, cols]
	for (i=2; i<=N; i++) {
		if ((tmp=X[i,cols]) != id) {
			info[  p,2] = i - 1
			info[++p,1] = i
			id = tmp
		}
		if (args()>2) (*pID)[i] = p
	}
	return (info[|.,.\p,.|])
}

// Do panelsum() except that a single missing value in X doesn't make all results missing and
// efficiently handles case when all groups have one row.
pointer(real matrix) scalar _panelsum(real matrix X, real matrix arg2, | real matrix arg3) {
	if (args()==2) {
		if (rows(arg2)==0 | rows(arg2)==rows(X))
			return(&X)
	} else if (rows(arg3)==0 | rows(arg3)==rows(X))
		return(arg2==1? &X : &(X :* arg2)) // if no collapsing called for, still fold in provided weights
	return (cols(arg3)? &panelsum(X, arg2, arg3) : &panelsum(X, arg2))
}

// concatenation of two _panelsum's
pointer(real matrix) scalar _panelsum2(real matrix X1, real matrix X2, real matrix arg2, | real matrix arg3)
	return(&(args()==2? *_panelsum(X1,arg2),*_panelsum(X2,arg2) : *_panelsum(X1,arg2,arg3),*_panelsum(X2,arg2,arg3)))

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


// cross-tab sum of a column vector w.r.t. given panel info and fixed-effect var
// one row per FE, one col per other grouping
real matrix boottest::crosstabFE(real colvector v, real matrix info) {
	real matrix retval; real scalar i, j, tmp; real colvector _FEID, _v
	retval = J(NFE, rows(info), 0)
  if (cols(info))
    for (i=cols(retval);i;i--) {
      _FEID = panelsubmatrix(*pFEID, i, info)
      _v    = panelsubmatrix(v     , i, info)
      for (j=rows(_FEID);j;j--) {
        tmp = _FEID[j] 
        retval[tmp,i] = retval[tmp,i] + _v[j]
      }
    }
  else  // "robust" case, no clustering, indicated by cols(info)=0
    for (i=cols(retval);i;i--)
      retval[(*pFEID)[i],i] = v[i]
	return(retval)
}

// subtract crosstab of v wrt bootstrapping cluster and all-cluster-var intersections from M
// M should have one row for each all-cluster-var (including bootstrap cluster) intersection and one col for each bootstrap cluster
// *** v needs to have been panelsum'd with pinfoAllData
void boottest::crosstabCapstarMinus(real matrix M, real colvector v) {
	real colvector tmp; real scalar i

	if (subcluster)  // crosstab c,c* is wide
		for (i=Clust.N;i;i--) {
			tmp = infoErrAll[i,]'
			M[|(i\i), tmp|] = M[|(i\i), tmp|] - v[|tmp|]'
		}
	else if (NClustVar == NBootClustVar)  // crosstab c,c* is square
		_diag(M, diagonal(M) - v)
	else  // crosstab c,c* is tall
		for (i=Nstar;i;i--) {
			tmp = Clust[BootClust].info[i,]'
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


// Chandrupatla 1997, "A new hybrid quadratic-bisection algorithm for finding the zero of a nonlinear function without using derivatives"
// x1, x2 must bracket the true value, with f1=f(x1) and f2=f(x2)
real scalar boottest::search(real scalar alpha, real scalar f1, real scalar x1, real scalar f2, real scalar x2) {
	real scalar x, t, fx, phi1, phi1_2, xi1, x3, f3
	
	t = 0.5
	do {
		fx = r_to_p(x = x1 + t * (x2 - x1))

		if (fx>f1 == fx>f2)   // violation of monotonicity because of precision problems? That's as good as it gets.
			return(x)

    if (fx<alpha == f1<alpha) {
			x3 = x1; x1 = x; f3 = f1; f1 = fx
		} else {
			x3 = x2; x2 = x1; x1 = x; f3 = f2; f2 = f1; f1 = fx
		}

		if ((B & abs(fx - alpha) < (1+(ptype==1))/BFeas*1.000001) | reldif(x2,x1) < ptol)
			return (abs(f1 - alpha) < abs(f2 - alpha)? x1 : x2)

		phi1 = (f1 - f2) / (f3 - f2)
		phi1_2 = phi1 * phi1
		xi1 = (x1 - x2) / (x3 - x2)
		if (phi1_2 > xi1 | xi1 > phi1 + phi1 - phi1_2)
			t = 0.5
		else {
			t = ((f3 - alpha) / (f1 - f2) + (x3 - x1) / ((x2 - x1) * (f3 - f1)) * (f2 - alpha)) * (f1 - alpha) / (f3 - f2)
			if (t < 0.000001)
				t = 0.000001
			else if (t > 0.999999)
				t = 0.999999
		}
	} while (1)
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
	string scalar X1names, string scalar Y2names, real scalar hascons, string scalar Ynames, string scalar bname, string scalar Aname, 
	string scalar X2names, string scalar samplename, string scalar scnames, real scalar robust, string scalar IDnames, real scalar NBootClustVar, real scalar NErrClust, 
	string scalar FEname, real scalar NFE, string scalar wtname, string scalar wttype, string scalar R1name, string scalar r1name, string scalar Rname, string scalar rname, real scalar B, string scalar repsname, string scalar repsFeasname, 
	real scalar small, string scalar diststat, string scalar distname, string scalar gridmin, string scalar gridmax, string scalar gridpoints, real scalar MaxMatSize, real scalar quietly,
	string scalar b0name, string scalar V0name, string scalar vname, string scalar NBootClustname) {
	real matrix X2, ID, FEID, sc, Y2, X1
	real colvector wt, Y
	class boottest scalar M
	pragma unset ID; pragma unset wt; pragma unset Y2; pragma unset X1; pragma unset Y; pragma unset X2; pragma unset sc


	M._st_view(sc, ., scnames, samplename)
	M._st_view(Y , ., Ynames , samplename)
	M._st_view(X2, ., X2names , samplename)
	if (FEname != "" ) FEID = st_data(., FEname , samplename)
	if (IDnames != "") ID   = st_data(., IDnames, samplename)
	if (wtname  != "") wt   = st_data(., wtname , samplename) // panelsum() doesn't like views as weights
	M.setMaxMatSize(MaxMatSize)
	M.sethascons(hascons)
	M.setsc(sc)
	M.setML(ML)
	M.setY (Y)
	M.setX2(X2)
	M.setwt (wt)
	M.setID(ID, NBootClustVar, NErrClust)
	M.setFEID(FEID, NFE)
	M.setR1(st_matrix(R1name), st_matrix(r1name))
	M.setR(st_matrix(Rname), st_matrix(rname))
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
	M._st_view(Y2, ., Y2names, samplename)
	M.setY2(Y2)
	M._st_view(X1, ., X1names, samplename)
	M.setX1(X1)
	if (bname != "") M.setbeta(st_matrix(bname)')
	if (Aname != "") M.setA   (st_matrix(Aname) )
	M.setwillplot(plotname != "")
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
}

mata mlib create lboottest, dir("`c(sysdir_plus)'l") replace
mata mlib add lboottest *(), dir("`c(sysdir_plus)'l")
mata mlib index
end
