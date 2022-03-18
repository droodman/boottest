*! boottest 4.0.0 18 March 2022
*! Copyright (C) 2015-22 David Roodman

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

// right-multiply a data matrix by a matrix, with efficient handling of special cases like latter is identity
pointer (real matrix) scalar pXB(real matrix X, real matrix M) {
  scalar r
  r = rows(M)
  return (all(colsum(M) :== 1) & all(colsum(!M) :== r-1)?  // is M 0's but for one 1 in each col, so it's just copying/reordering cols?
             (all(diagonal(M)) & rows(M)==cols(M)?  // is M just the identity matrix?
                &X            :
                &X[,colsum(M:*(1::r))]) :  // reorder X
             &(X * M))
}

// return [X1 X2]B where X1 or X2 may have 0 cols
pointer(real matrix) scalar pX12B(real matrix X1, real matrix X2, real matrix B)
  return(cols(B)? (cols(X1)? (cols(X2)? &(*pXB(X1, B[|.,.\cols(X1),.|]) + *pXB(X2, B[|cols(X1)+1,.\.,.|])) : pXB(X1, B)) : pXB(X2, B)) : &J(rows(X1),0,0))

// return &X[|S|] with appropiate behavior if cols(X)=0 or S requests no rows (S[2,1]=0), presuming that in degenerate cases S does not specify columns
// if retval = X, doesn't duplicate data
// S should be 2x1 because the function is only for selecting rows
pointer(real matrix) scalar pXS(real matrix X, real matrix S)
  return(cols(X)? (S[2,1]? ((S[1,1]==. | S[1,1]==1) & (S[2,1]==. | S[2,1]==rows(X))? &X : &X[|S,(.\.)|]) : &J(0,cols(X),0)) : &J(editmissing(S[2,1],rows(X))-editmissing(S[1,1],1)+1,0,0))

// :* operation, handling case that either argument is just 1 without duplicating data
pointer (real colvector) scalar pvHadw(real matrix v, real matrix w)
  return(w==1? &v : (v==1? &w : &(v :* w)))

class boottestOLS {  // class for performing OLS
  real scalar LIML, Fuller, kappa, isDGP, kZ
  real colvector y1, u1ddot, u1dddot, beta, beta0, invXXXy1par
  real rowvector Yendog
  real matrix invZperpZperp, XZ, PXZ, YPXY, R1invR1R1, R1perp, Rpar, RperpX, RRpar, RparY, RR1invR1R1, dbetadr, YY, AR, XAR, R1invR1R1Y, invXXXZ, U2ddot, XinvXX, Rt1, invXX, Y2, X2, invH
  pointer(real colvector) scalar py1par, pXy1par
  pointer(real matrix) scalar pA, pZ, pZperp, pX1
  pointer (class boottest scalar) scalar parent
  struct smatrix matrix FillingT0
  struct smatrix rowvector WXAR, ScapPXYZperp, ScapYX, CT_XAR, CT_FEcapPY

  private void new(), InitTestDenoms()
  private virtual void InitVars(), SetR(), Estimate(), MakeResiduals()
  real matrix _select(), perp()
}

class boottestARubin extends boottestOLS {
  private virtual void InitVars(), Estimate()
}

class boottestIVGMM extends boottestOLS {
  real matrix ZZ, XY2, XX, H_2SLS, V, ZY2, X2Y2, X1Y2, ZR1ZR1, X2ZR1, ZR1Y2, X1ZR1, ZZR1, X2y1, X1y1, Zy1, ZXinvXXXZ, H_2SLSmZZ
  real colvector ZXinvXXXy1par, t1Y
  real rowvector y1Y2, twoy1ZR1
  real scalar y1y1, y1pary1par
  pointer(real colvector) scalar pX2y1par, pX1y1par, pZy1par
  pointer(real rowvector) scalar py1parY2
  pointer(real matrix) scalar pRperp, pZR1
  private void new()
  private virtual void InitVars(), Estimate(), MakeResiduals(), MakeH()
}

class boottest {
  real scalar scoreBS, B, small, auxwttype, null, dirty, initialized, ML, Nobs, _Nobs, kZ, kY2, kX1, sumwt, NClustVar, haswt, REst, multiplier, smallsample, quietly, FEboot, NErrClustCombs, ///
    sqrt, LIML, Fuller, kappa, WRE, WREnonARubin, ptype, twotailed, df, df_r, ARubin, willplot, notplotted, NumH0s, p, NBootClustVar, NErrClustVar, BootClust, FEdfadj, ///
    NFE, granular, purerobust, subcluster, Nstar, BFeas, v_sd, level, ptol, MaxMatSize, Nw, enumerate, bootstrapt, q, interpolable, interpolating, interpolate_u, robust, kX2, kX
  real matrix AR, v, ustar, CI, CT_WE, infoBootData, infoErrAll, JNcapNstar, statDenom, uXAR, SuwtXA, numer0, betadev, deltadenom_b, _Jcap, YYstar_b, YPXYstar_b, numerw
  real colvector DistCDR, plotX, plotY, beta, ClustShare, WeightGrpStart, WeightGrpStop, confpeak, gridmin, gridmax, gridpoints, numersum, uddot0, anchor, poles, invFEwt
  real rowvector peak, betas, As
  string scalar obswttype, madjtype, seed
  pointer (real matrix) scalar pX2, pR1, pR, pID, pFEID, pY2, pX1, pinfoAllData, pinfoCapData, pIDAll, pnumer, pU2parddot
  pointer (real colvector) scalar pr1, pr, py1, pSc, pwt, pA, puddot, pDist, pIDBootData, pIDBootAll
  class boottestOLS scalar DGP, Repl
  pointer(class boottestOLS scalar) scalar pM
  struct structboottestClust rowvector Clust
  struct smatrix matrix denom, Kcd, denom0, Jcd0, SCTcapuXinvXX, SstarUU, CTUX
  struct smatrix rowvector Kd, dudr, dnumerdr, IDCTCapstar, infoCTCapstar, SstarUX, SstarUXinvXX, SstarUZperpinvZperpZperp, deltadenom, Zyg, SstaruY, SstarUMZperp, SstarUPX, SstarUZperp, YYstar, YPXYstar, CTFEU
  struct ssmatrix rowvector ddenomdr, dJcddr
  struct ssmatrix matrix ddenomdr2
  pointer(struct smatrix matrix) scalar pJcd
  struct structFE rowvector FEs
  
  void new(), setsqrt(), setX1(), setptype(), setdirty(), setY2(), setY(), setX2(), setobswt(), setsc(), setML(), setLIML(), setARubin(), setauxwttype(),
    setFuller(), setkappa(), setquietly(), setbeta(), setA(), setsmall(), sethascons(), setscoreBS(), setB(), setnull(), setID(), setFEID(), setlevel(), setptol(), 
    setrobust(), setR1(), setR(), setwillplot(), setgrid(), setmadjust(), setMaxMatSize(), setstattype(), close()
  private void MakeNumerAndJ(), _clustAccum(), MakeWREStats(), MakeInterpolables(), _MakeInterpolables(), MakeNonWREStats(), UpdateBootstrapcDenom(), Init(), plot(), MakeWildWeights(), boottest(), crosstabCapstarMinus(), PrepWRE(), storeWtGrpResults(), NoNullUpdate()
  real matrix getplot(), getCI(), getV(), getv()
  real scalar getp(), getpadj(), getstat(), getdf(), getdf_r(), getreps(), getrepsFeas(), getNBootClust()
  real rowvector getpeak()
  real colvector getdist(), getb()
  private real scalar r_to_p(), search()
  private real matrix count_binary(), crosstabFE(), HessianFixedkappa()
  private real rowvector _HessianFixedkappa()
  private pointer(real matrix) scalar Filling(), partialFE()
  private static real matrix combs()
  private static real colvector stableorder()
  private real vector _selectindex()
  static void _st_view()
}

void boottestOLS::new() {
  LIML = Fuller = kappa = 0    
  Rpar = isDGP = 1
}

void boottestIVGMM::new() {
  Fuller = 0
  kappa = isDGP = 1
}


// do select() but handle case that both args are scalar and second is 0 by interpreting second arg as rowvector and returning J(1,0,0)
// if v = 0 (so can't tell if row or col vector), returns J(1, 0, 0) 
real matrix boottestOLS::_select(real matrix X, real rowvector v)
  return (rows(X)==1 & cols(X)==1 & v==0? J(1,0,0) : select(X,v))

real matrix boottestOLS::perp(real matrix A) {
  real matrix vec; real rowvector val; pragma unset vec; pragma unset val
  symeigensystem(A*invsym(A'A)*A', vec, val); _edittozero(val, 1000)
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
    symeigensystem(R1invR1R1 * R1, vec, val); _edittozero(val, 1000)
    R1perp = _select(vec, !val)  // eigenvectors orthogonal to span of R1; foundation for parameterizing subspace compatible with constraints
  } else
    R1invR1R1 = J(parent->kZ,0,0)  // and R1perp = I

  if (kappa) {      
    // prepare to reduce regression via FWL
    RR1perp = R \ J(parent->kY2, parent->kX1, 0), I(parent->kY2)  // rows to prevent partialling out of endogenous regressors
    if (rows(R1))
      RR1perp = RR1perp * R1perp 
    symeigensystem(RR1perp ' invsym(RR1perp * RR1perp') * RR1perp, vec, val); _edittozero(val, 1000)
    Rpar   = _select(vec,  val)  // perp and par of RR₁perp
    RperpX = _select(vec, !val)

    if (rows(R1)) {  // fold model constraint factors into Rpar, RperpX
      Rpar   = R1perp * Rpar
      RperpX = R1perp * RperpX
    }
    RRpar = R * Rpar
    if (isDGP & parent->WREnonARubin) RperpX = parent->Repl.RperpX
    RperpX = *pXS(RperpX, .\parent->kX1)  // Zperp=Z*RperpX; though formally a multiplier on Z, it will only extract exogenous components, in X1, since all endogenous ones will be retained

    S = parent->kX1+1\.
    RparY      = *pXS(Rpar     , S)  // part of Rpar that refers to Y2
    R1invR1R1Y = *pXS(R1invR1R1, S)
    RR1invR1R1 = R * R1invR1R1
  }
}

// stuff that can be done before r set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
void boottestOLS::InitVars(pointer(real matrix) scalar pRperp) {  // Rperp is for replication regression--no null imposed
  real matrix H; pointer(real matrix) scalar pR1AR1

  py1par = parent->py1
  invH = invsym(H = cross(*parent->pX1, *parent->pwt, *parent->pX1))

  pR1AR1 = rows(R1perp)? &( R1perp * invsym( R1perp ' H *  R1perp) *  R1perp') : &invH  // for DGP regression
  beta0   = *pR1AR1 * cross(*parent->pX1, *parent->pwt, *py1par)
  dbetadr = *pR1AR1 * H * R1invR1R1 - R1invR1R1

  pA = rows(*pRperp)? &(*pRperp * invsym(*pRperp ' H * *pRperp) * *pRperp') : &invH  // for replication regression
  AR = *pA * *parent->pR'
  if (parent->scoreBS | parent->robust)
    XAR = *parent->pX1 * AR
}

void boottestARubin::InitVars(| pointer(real matrix) pRperp) {
  pragma unused pRperp
  real matrix H, X2X1; pointer(real matrix) scalar pR1AR1

  X2X1 = cross(*parent->pX2, *parent->pwt, *parent->pX1)
  H = cross(*parent->pX1, *parent->pwt, *parent->pX1), X2X1' \ X2X1, cross(*parent->pX2, *parent->pwt, *parent->pX2)
  pA = &invsym(H)
  AR = *pA * *parent->pR'
  if (parent->scoreBS | parent->robust)
    XAR = *pX12B(*parent->pX1, *parent->pX2, AR)

  pR1AR1 = rows(R1perp)? &(R1perp * invsym(R1perp ' H * R1perp) * R1perp') : pA
  beta0   = *pR1AR1 * (cross(*parent->pX1, *parent->pwt, *parent->py1) \ cross(*parent->pX2, *parent->pwt, *parent->py1))
  dbetadr = *pR1AR1 * (cross(*parent->pX1, *parent->pwt, *parent->pY2) \ cross(*parent->pX2, *parent->pwt, *parent->pY2))
}

void boottestIVGMM::InitVars(|pointer(real matrix) scalar pRperp) {
  real matrix X2X1; real scalar i,j; pointer (real colvector) scalar puwt
  
  this.pRperp = pRperp

  pZperp = pXB(*parent->pX1, RperpX)
  invZperpZperp = invsym(cross(*pZperp, *parent->pwt, *pZperp))
  pX1 = pXB(*parent->pX1, perp(RperpX)); pX1 = &(*pX1 - *pZperp * (invZperpZperp * cross(*pZperp, *parent->pwt, *pX1)))  // FWL-process X1
  X2 = *parent->pX2 - *pZperp * (invZperpZperp * cross(*pZperp, *parent->pwt, *parent->pX2))                // FWL-process X2
  X2X1 = cross(X2, *parent->pwt, *pX1)
  invXX = invsym((XX = cross(*pX1, *parent->pwt, *pX1), X2X1' \ X2X1, cross(X2, *parent->pwt, X2)))
  pZ   = pX12B(*parent->pX1, *parent->pY2, Rpar     )  // Zpar
  pZR1 = pX12B(*parent->pX1, *parent->pY2, R1invR1R1)  // Z*R1

  pZ   =    &(*pZ   - *pZperp * (invZperpZperp * cross(*pZperp, *parent->pwt,        *pZ  )))  // partialling out
  pZR1 =    &(*pZR1 - *pZperp * (invZperpZperp * cross(*pZperp, *parent->pwt,        *pZR1)))
  Y2 = *parent->pY2 - *pZperp * (invZperpZperp * cross(*pZperp, *parent->pwt, *parent->pY2))
  y1 = *parent->py1 - *pZperp * (invZperpZperp * cross(*pZperp, *parent->pwt, *parent->py1))

  X1Y2 = cross(*pX1, *parent->pwt, Y2)
  X2Y2 = cross(X2  , *parent->pwt, Y2)
  XY2 = X1Y2 \ X2Y2
  y1Y2 = cross(y1  , *parent->pwt, Y2)
  X2y1 = cross(X2  , *parent->pwt, y1)
  X1y1 = cross(*pX1, *parent->pwt, y1)
  y1y1 = cross(y1  , *parent->pwt, y1)
  Zy1  = cross(*pZ , *parent->pwt, y1)    
  XZ   = cross(*pX1, *parent->pwt, *pZ) \ 
         cross(X2  , *parent->pwt, *pZ)
  ZY2 =  cross(*pZ , *parent->pwt, Y2)
  ZZ  =  cross(*pZ , *parent->pwt, *pZ )
  
  ZXinvXXXZ = XZ ' (invXXXZ = invXX * XZ)

  if (cols(R1invR1R1)) {
    X2ZR1  = cross(X2   , *parent->pwt, *pZR1)
    X1ZR1  = cross(*pX1 , *parent->pwt, *pZR1)
    ZZR1   = cross(*pZ  , *parent->pwt, *pZR1)
    twoy1ZR1  = cross( y1  , *parent->pwt, *pZR1); twoy1ZR1 = twoy1ZR1 + twoy1ZR1
    ZR1ZR1 = cross(*pZR1, *parent->pwt, *pZR1)
    ZR1Y2  = cross(*pZR1, *parent->pwt, Y2   )
  } else {
    py1parY2   = &y1Y2 
    pX2y1par   = &X2y1
    pX1y1par   = &X1y1
    pZy1par    = &Zy1
    y1pary1par = y1y1 
    pXy1par = &(X1y1 \ X2y1)
    py1par = &y1
  }

  V =  invXX * XZ  // in 2SLS case, estimator is (V' XZ)^-1 * (V'Xy1). Also used in kZ-class and LIML robust VCV by Stata convention
  H_2SLS = V ' XZ  // Hessian
  if (kappa != 1 | LIML) H_2SLSmZZ = H_2SLS - ZZ

  if (isDGP) {
    if (LIML==0)  // DGP is LIML except possibly when getting confidence peak for A-R plot; but LIML=0 when exactly id'd, for then kappa=1 always and Hessian doesn't depend on r1 and can be computed now
      MakeH()
  } else {
    kZ = cols(Rpar)
    Yendog = 1, colsum(RparY :!= 0)  // columns of Y = [y1par Zpar] that are endogenous (normally all)

    if (parent->robust & parent->bootstrapt) {  // for WRE replication regression, prepare for CRVE
      ScapYX = ScapPXYZperp = smatrix(kZ+1)
      XinvXX = *pX12B(*pX1, X2, invXX); PXZ = *pX12B(*pX1, X2, invXXXZ)
      if (parent->haswt) {
        XinvXX = XinvXX :* *parent->pwt
        PXZ = PXZ :* *parent->pwt
      }

      FillingT0 = smatrix(kZ+1, kZ+1)  // fixed component of groupwise term in sandwich filling
      if (parent->NFE)
        CT_FEcapPY = smatrix(kZ+1)
      for (i=kZ; i; i--) {
        puwt = pcol(PXZ,i)
        if (parent->NFE)
          CT_FEcapPY[i+1].M = parent->crosstabFE(*puwt, *parent->pinfoCapData) :* parent->invFEwt
        for (j=kZ; j; j--)
          FillingT0[i+1,j+1].M = *_panelsum(*pcol(*pZ,j), *puwt, *parent->pinfoCapData)
      }

      for (i=kZ; i; i--) {  // precompute various clusterwise sums
        ScapPXYZperp[i+1].M = *_panelsum(*pZperp, *pcol(PXZ,i), *parent->pinfoCapData)  // S_cap(P_(MZperpX) * Z :* Zperp)
        if (parent->granular==0)
          ScapYX[i+1].M = *_panelsum2(*pX1, X2, *pvHadw(*pcol(*pZ,i), *parent->pwt), *parent->pinfoCapData)  // S_cap(M_Zperp[Z or y1] :* P_(MZperpX)])
      }
    }
  }
}


// do most of estimation; for LIML r1 must be passed now in order to solve eigenvalue problem involving it
// inconsistency: for replication regression of Anderson-Rubin, r1 refers to the *null*, not the maintained constraints, because that's what affects the endogenous variables
// For OLS, compute beta0 (beta when r=0) and dbetadr without knowing r1, for efficiency
// For WRE, should only be called once for the replication regressions, since for them r1 is the unchanging model constraints
void boottestOLS::Estimate(real colvector r1)
  beta = beta0 - dbetadr * r1

void boottestARubin::Estimate(real colvector r1) {
  beta = beta0 - dbetadr * r1
  py1par = &(*parent->py1  - *parent->pY2 * r1)
}

void boottestIVGMM::MakeH() {
  pointer(real matrix) scalar pH

  pH = kappa==1? &H_2SLS : &(ZZ + kappa * H_2SLSmZZ)
  invH = invsym(*pH)

  if (pRperp) {  // for score bootstrap
    pA = cols(*pRperp)? &(*pRperp * invsym(*pRperp ' (*pH) * *pRperp) * *pRperp') : &invH
    AR = *pA * (Rpar ' (*parent->pR'))
    XAR = *pX12B(*pX1, X2, V * AR)
  }
}

void boottestIVGMM::Estimate(real colvector r1) {
  real rowvector val; real matrix vec; real scalar i; pointer (real colvector) scalar puwt
  pragma unset vec; pragma unset val

  if (cols(R1invR1R1)) {
    y1pary1par = y1y1 - twoy1ZR1 * r1 + r1 ' ZR1ZR1 * r1
    py1par   = &(y1 - *pZR1 * r1)
    py1parY2 = &(y1Y2  - r1 ' ZR1Y2)
    pX2y1par = &(X2y1 - X2ZR1 * r1)
    pX1y1par = &(X1y1 - X1ZR1 * r1)
    pZy1par  = &( Zy1 -  ZZR1 * r1)
    pXy1par  = &(*pX1y1par \ *pX2y1par)
  }

  ZXinvXXXy1par = XZ ' (invXXXy1par = invXX * *pXy1par)
  YY = y1pary1par, *pZy1par' \ *pZy1par, ZZ
  YPXY = invXXXy1par ' (*pXy1par) , ZXinvXXXy1par' \ ZXinvXXXy1par , ZXinvXXXZ

  if (isDGP) {
    if (LIML) {
      eigensystemselecti(invsym(YY) * YPXY, rows(YY)\rows(YY), vec, val)
      kappa = 1/(1 - Re(val)) // sometimes a tiny imaginary component sneaks into val
      if (Fuller) kappa = kappa - Fuller / (parent->_Nobs - parent->kX)
      MakeH()
    }

    beta = invH * (kappa==1?  ZXinvXXXy1par : kappa * (ZXinvXXXy1par - *pZy1par) + *pZy1par)
    t1Y = R1invR1R1Y * r1
  } else if (parent->WREnonARubin) {  // if not score bootstrap of IV/GMM...
    Rt1 = RR1invR1R1 * r1
    
    if (parent->robust & parent->bootstrapt) {  // prepare WRE replication regressions
      for (i=kZ; i; i--)
        FillingT0[i+1,1].M = *_panelsum (*pcol(PXZ,i), *py1par, *parent->pinfoCapData)
      if (parent->granular==0)
        ScapYX.M           = *_panelsum2(*pX1, X2    , *pvHadw(*py1par, *parent->pwt), *parent->pinfoCapData)  // S_cap(M_Zperp*y1 :* P_(MZperpX)])
    }
  }
}


void boottestOLS::MakeResiduals()
  u1ddot = *py1par - *pX12B(*parent->pX1, *parent->pX2, beta)

void boottestIVGMM::MakeResiduals() {
  real matrix Xu; real colvector negXuinvuu, _beta; real scalar uu
  
  u1ddot = *py1par - *pZ * beta

  if (parent->scoreBS==0) {
    _beta = 1 \ -beta
    uu = _beta ' YY * _beta

    Xu = *pXy1par - XZ * beta  // after DGP regression, compute Y2 residuals by regressing Y2 on X while controlling for y1 residuals, done through FWL
    negXuinvuu = Xu / -uu
    U2ddot = Y2 - *pX12B(*pX1, X2, invsym(XX + negXuinvuu * Xu') * (negXuinvuu * (*py1parY2 - beta ' ZY2) + XY2))  // large expression is Pihat

    u1dddot = u1ddot + U2ddot * (t1Y + RparY * beta)
  }
}


// non-WRE stuff that only depends on r in A-R case, for test stat denominators in replication regressions
// since the non-AR OLS code never creates an object for replication regresssions, in that case this is called on the DGP regression object
// depends on results of Estimate() only when doing OLS-style bootstrap on an overidentified IV/GMM regression--score bootstrap or A-R. Then kappa from DGP LIML affects Hessian, pH.
void boottestOLS::InitTestDenoms() {
  real scalar d; pointer (real matrix) scalar pWXAR

  if (parent->bootstrapt & (parent->scoreBS | parent->robust)) {
    if (parent->granular | parent->purerobust) {
      pWXAR = pvHadw(XAR, *parent->pwt)
      WXAR = smatrix(parent->df)
      for (d=parent->df;d;d--)
        WXAR[d].M = (*pWXAR)[,d]
    }

    if (parent->NFE & parent->robust & (parent->FEboot | parent->scoreBS)==0 & parent->granular < parent->NErrClustCombs) {  // make first factor of second term of (64) for c=∩ (c=1)
      if (pWXAR == NULL)
        pWXAR = pvHadw(XAR, *parent->pwt)
      CT_XAR = smatrix(parent->df)
      for (d=parent->df;d;d--)
        CT_XAR[d].M = parent->crosstabFE((*pWXAR)[,d], *parent->pinfoCapData)
    }
  }
}


// partial fixed effects out of a data matrix
pointer(real matrix) scalar boottest::partialFE(pointer(real matrix) scalar pIn) {
  real matrix Out, tmp; real scalar i
  if (NFE & pIn) {
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
  ARubin = LIML = Fuller = WRE = small = scoreBS = auxwttype = ML = initialized = quietly = sqrt = ptype = robust = NFE = FEboot = granular = NErrClustCombs = subcluster = B = BFeas = interpolating = 0
  twotailed = null = dirty = willplot = v_sd = notplotted = FEdfadj = bootstrapt = 1
  level = 95
  ptol = 1e-6
  confpeak = MaxMatSize = .
  pY2 = pX1 = pX2 = py1 = pSc = pID = pFEID = pR1 = pR = pwt = &J(0,0,0)
  pr1 = pr = &J(0,1,0)
  pIDBootData = pIDBootAll = &.
}

// important to call this when done: break loops in data structure topology to enable garbage collection
void boottest::close() {
  DGP.parent = Repl.parent = NULL
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

void boottest::setX1(real matrix X1) {
  this.pX1  = &X1; setdirty(1)
}
void boottest::setX2(real matrix X2) {
  this.pX2  = &X2; setdirty(1)
}
void boottest::setY(real matrix y1) {
  this.py1  = &y1; setdirty(1)
}
void boottest::setY2(real matrix Y2) {
  this.pY2  = &Y2; setdirty(1)
}
void boottest::setobswt(real matrix wt, string scalar obswttype) {
  this.pwt  = &wt; this.obswttype = obswttype; setdirty(1)
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
void boottest::setFuller(real scalar Fuller) {
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
void boottest::setA(real matrix A) {
  this.pA = &A; setdirty(1)
}
void boottest::setsmall(real scalar small) {
  this.small = small; setdirty(1)
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
void boottest::setID      (real matrix ID, | real scalar NBootClustVar, real scalar NErrClustVar) {
  this.pID = &ID; this.NBootClustVar = editmissing(NBootClustVar,1); this.NErrClustVar=editmissing(NErrClustVar,editmissing(NBootClustVar,1)); setdirty(1)
  if (cols(ID)) this.robust = 1
}
void boottest::setFEID(real matrix ID, real scalar NFE, | real scalar FEdfadj) {
  this.pFEID = &ID; this.NFE = NFE; this.FEdfadj = editmissing(FEdfadj,1); setdirty(1)
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
  this.pR1 = &R1;   this.pr1  = &r1; setdirty(1)
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
void boottest::setauxwttype(string scalar auxwttype) {
  auxwttype = strlower(auxwttype)
  if (.==(this.auxwttype = auxwttype=="rademacher" ? 0 : (auxwttype=="mammen" ? 1 : (auxwttype=="webb" ? 2 : (auxwttype=="normal" ? 3 : (auxwttype=="gamma" ? 4 : .))))))
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
    _pnumer = v_sd==1? pnumer : &(*pnumer / v_sd)
    _sort( DistCDR = (*_pnumer)[|2\.|]' :+ *pr , 1)
  } else if (rows(DistCDR)==0)
    if (cols(*pDist) > 1)
      _sort( DistCDR = multiplier * (*pDist)[|2\.|]' , 1)
    else
      DistCDR = J(0,1,0)
  return(DistCDR)
}

// get p value. Robust to missing bootstrapped values interpreted as +infinity.
real scalar boottest::getp(|real scalar classical) {
  real scalar tmp; real scalar _p
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
    _p = small? Ftail(df, df_r, sqrt? tmp*tmp : tmp) : chi2tail(df, sqrt? tmp*tmp : tmp)
    if (sqrt & twotailed==0) {
      _p = _p / 2
      if ((ptype==3) == (tmp<0))
        _p = 1 - _p
    }
    if (classical != .)
      return(_p)
    p = _p  // only save as the official p if this was not requested by plot() for AR case
  }
  return(p)
}

// numerator for full-sample test stat
real colvector boottest::getb() {
  if (dirty) boottest()
  return(v_sd == 1? (*pnumer)[,1] : (*pnumer)[,1] / v_sd)
}

// denominator for full-sample test stat
real matrix boottest::getV() {
  if (dirty) boottest()
  return (statDenom / ((v_sd == 1? smallsample : v_sd * v_sd * smallsample)  * (sqrt? multiplier*multiplier : multiplier) * df))
}

// wild weights
real matrix boottest::getv()
  return(v_sd==1? v[|.,2\.,.|] : v[|.,2\.,.|] / v_sd)

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
  real scalar _p
  _p = dirty | classical != . ? getp(classical) : p
  if (madjtype=="bonferroni") return(min((1, NumH0s*_p)))
  if (madjtype=="sidak"     ) return(1 - (1 - _p)^NumH0s)
  return(_p)
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
          (Clust.even?
                 (Clust.   multiplier != 1?   Clust.   multiplier  * Y :  Y) :
                 (Clust.   multiplier != 1? (-Clust.   multiplier) * Y : -Y)) :
          (Clust[c].even?
             X + (Clust[c].multiplier != 1?   Clust[c].multiplier  * Y :  Y) :
             X - (Clust[c].multiplier != 1?   Clust[c].multiplier  * Y :  Y))


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
  real matrix Combs, tmp, IDCap
  real scalar i, j, c, minN, _B, i_FE, sumFEwt
  pragma unset pIDAllData; pragma unset pIDCapData

  Nobs = rows(*pX1)
  NClustVar = cols(*pID)
  kX = (kX1 = cols(*pX1)) + (kX2 = cols(*pX2))
  if (kX2 == 0) pX2 = &J(Nobs,0,0)
  if ((kY2 = cols(*pY2)) == 0) pY2 = &J(Nobs,0,0)
  kZ  = kX1 + kY2
  if (LIML & kX2 == kY2) {  // exactly identified LIML = 2SLS
    kappa = 1
    LIML = 0
  }
  if ((REst = rows(*pR1)) == 0) {  // base model contains no restrictions?
    pR1 = &J(0,kZ,0)
    pr1 = &J(0,1 ,0)
  }
  if (kappa==.) kappa = kX2>0  // if kappa in kappa-class estimation not specified, it's 0 or 1 for OLS or 2SLS
  WRE = (kappa & scoreBS==0) | ARubin
  WREnonARubin = WRE & ARubin==0

  if (haswt = rows(*pwt)>1)
    sumwt = sum(*pwt)
  else
    pwt = &(sumwt = 1)
  _Nobs = haswt & obswttype=="fweight"? sumwt : Nobs

  if (WREnonARubin)
    if (NClustVar)
      infoBootData = _panelsetup(*pID, 1..NBootClustVar, pIDBootData)
    else
      pinfoCapData = &(infoBootData = J(Nobs,0,0))  // no clustering, so no collapsing by cluster
  else if (NClustVar)
    infoBootData = _panelsetup(*pID, 1..min((NClustVar, NBootClustVar)))  // bootstrap cluster grouping defs rel to original data
  else
    pinfoCapData = pinfoAllData = &(infoBootData = J(Nobs,0,0))  // causes no collapsing of data in _panelsum() calls, only multiplying by weights if any
  Nstar = rows(infoBootData)

  if (bootstrapt) {
    if (NClustVar) {
      minN = .

      Combs = combs(NErrClustVar)  // represent all error clustering combinations. First is intersection of all error clustering vars
      Clust = structboottestClust(rows(Combs)-1)  // leave out no-cluster combination
      NErrClustCombs = length(Clust)
      subcluster = NClustVar - NErrClustVar

      if (NClustVar > NBootClustVar) {  // info for grouping by intersections of all bootstrap & clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster & intersection of all error clusters
        pinfoAllData = WREnonARubin & granular==0? &_panelsetup(*pID, 1..NClustVar, pIDAllData) : 
                                                   &_panelsetup(*pID, 1..NClustVar            )
        if (subcluster & rows(infoBootData) != rows(*pinfoAllData)) {
          errprintf("\nboottest can only perform the subcluster bootstrap when the bootstrap clusters are nested within the (intersections of the) error clusters.\n")
          _error(499)
        }
      } else {
        pinfoAllData = &infoBootData  // info for grouping by intersections of all bootstrap & clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster & intersection of all error clusters
        if (WREnonARubin & granular==0)
          pIDAllData = pIDBootData
      }

      if (NClustVar > NErrClustVar)  // info for intersections of error clustering wrt data
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

        Clust[c].N = rows(Clust[c].info)

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
      Clust.N = Nobs
      Clust.info = J(Nobs, 0, 0)  // signals _panelsum not to aggregate
      NErrClustCombs = 1
      if (scoreBS | WREnonARubin==0)
        ClustShare = haswt? *pwt/sumwt : 1/_Nobs
    }

    purerobust = robust & (scoreBS | subcluster)==0 & Nstar==Nobs  // do we ever error- *and* bootstrap-cluster by individual?
    granular   = WREnonARubin? 2*Nobs*B*(2*Nstar+1) < Nstar*(Nstar*Nobs+Clust.N*B*(Nstar+1)) :
                               robust & scoreBS==0 & (purerobust | (Clust.N+Nstar)*kZ*B + (Clust.N-Nstar)*B + kZ*B < Clust.N*kZ*kZ + Nobs*kZ + Clust.N * Nstar*kZ + Clust.N*Nstar)

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
        tmp = (*pIDAllData)[|infoBootData[i,]'|]                       // ID numbers w.r.t. intersection of all bootstrap/error clusterings contained in bootstrap cluster i
        infoCTCapstar[i].M = (*pinfoAllData)[tmp[1]::tmp[rows(tmp)],]  // for each of those ID's, panel info for the all-bootstrap/error-clusterings data row groupings
        IDCTCapstar[i].M = (*pIDCapData)[infoCTCapstar[i].M[,1]]       // ID numbers of those groupings w.r.t. the all-error-clusterings grouping
      }
    }

  } else
    minN = rows(infoBootData)

  if (NFE) {  // identify FE groups
    sortID = (*pFEID)[o = stableorder(*pFEID, 1)]
    i_FE = 1; FEboot = B>0 & WREnonARubin==0 & NClustVar; j = Nobs; _FEID = J(Nobs, 1, 1)
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
        if ((B & robust & granular < NErrClustVar) | (WREnonARubin & robust & granular & bootstrapt))
          invFEwt[i_FE] = 1 / sumFEwt

        j = i
        
        if (FEboot) {  // are all of this FE's obs in same bootstrapping cluster? (But no need to check if B=0 for then CT_WE in 2nd term of (62) orthogonal to v = col of 1's)
          tmp = (*pID)[FEs[i_FE].is, 1..NBootClustVar]
          FEboot = all(tmp :== tmp[1,])
        }
        ++i_FE
      }
      _FEID[o[i]] = i_FE
    }
    FEs[NFE].is = o[|.\j|]
    if (haswt) {
      tmp  = (*pwt)[FEs[NFE].is]
      FEs[NFE].wt = tmp / (sumFEwt = colsum(tmp))
    } else
      FEs[NFE].wt = J(j,1,1/(sumFEwt = j))
    if (robust & ((B & granular < NErrClustVar) | (WREnonARubin & granular & bootstrapt)))
      invFEwt[NFE] = 1 / sumFEwt
    if (FEboot) {  // are all of this FE's obs in same bootstrapping cluster?
      tmp = (*pID)[FEs[NFE].is, 1..NBootClustVar]
      FEboot = all(tmp :== tmp[1,])
    }

    pFEID = &_FEID  // ordinal fixed effect ID
    pX1 = partialFE(pX1)
    pX2 = partialFE(pX2)
    py1 = partialFE(py1)
    pY2 = partialFE(pY2)
  }

  if (B & robust & granular & purerobust==0 & bootstrapt & WREnonARubin==0)
    if (NFE & FEboot==0)
      (void) _panelsetup(*pID   , 1..NBootClustVar, pIDBootData)
    else
      (void) _panelsetup(*pIDAll, 1..NBootClustVar, pIDBootAll )

  if (enumerate = (B & auxwttype==0 & Nstar*ln(2) < ln(B)+1e-6))  // generate full Rademacher set?
    MaxMatSize = .

  Nw = MaxMatSize == .? 1 : ceil((B+1) * max((rows(*pIDBootData), rows(*pIDBootAll), Nstar)) * 8 / MaxMatSize / 1.0X+1E) // 1.0X+1E = giga(byte)
  if (Nw == 1) {
    MakeWildWeights(B, 1)  // make all wild weights, once
    if (enumerate) B = cols(v) - 1  // replications reduced to 2^G
    WeightGrpStart = 1
  } else {
    seed = rseed()
    _B = ceil((B+1) / Nw)
    Nw = ceil((B+1) / _B)
     WeightGrpStart = (0::Nw-1) * _B :+ 1
    (WeightGrpStop  = (1::Nw  ) * _B     )[Nw] = B+1
  }

  if (ML)
    df = rows(*pR)
  else {
    if (ARubin) {
      pR  = &(J(kX2,kX1,0), I(kX2))  // attack surface is all endog vars
      pR1 = kX1 & rows(*pR1)? &((*pR1)[|.,.\.,kX1|], J(rows(*pR1),kX2,0)) : &J(0, kX, 0)  // and convert model constraints from referring to X1, Y2 to X1, X2
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
        DGP.Estimate(*pr1)
        confpeak = DGP.beta  // estimated coordinate of confidence peak
      }
       
      DGP = boottestARubin(); DGP.parent = &this
      DGP.SetR(*pR1)
      DGP.InitVars()
      DGP.InitTestDenoms()
      pM = &DGP  // estimator object from which to get A, AR, XAR
      kZ = kX

    } else if (WREnonARubin) {

      Repl = boottestIVGMM()
      Repl.parent = &this
      Repl.isDGP = 0
      Repl.LIML = this.LIML; Repl.Fuller = this.Fuller; Repl.kappa = this.kappa
      Repl.SetR(*pR1, *pR)
      Repl.InitVars()
      Repl.Estimate(*pr1)

      DGP = boottestIVGMM(); DGP.parent = &this
      DGP.LIML = kX2!=kY2
      DGP.SetR(null? *pR1 \ *pR : *pR1, J(0,kZ,0))  // DGP constraints: model constraints + null if imposed
      DGP.InitVars()
      if (null==0) {  // if not imposing null, then DGP constraints, kappa, Hessian, etc. do not vary with r and can be set now
        DGP.Estimate(*pr1)
        DGP.MakeResiduals()
      }

      if (LIML & Repl.kZ==1 & Nw==1) As = betas = J(1, B+1, 0)
      SstarUZperpinvZperpZperp = SstarUZperp = SstaruY = SstarUXinvXX = SstarUX = smatrix(Repl.kZ+1)
      SstarUU = smatrix(Repl.kZ+1, Repl.kZ+1)
      if (bootstrapt) {
        deltadenom_b = J(Repl.kZ, Repl.kZ, 0)
        Zyg = deltadenom = smatrix(Repl.kZ+1)
        SstarUMZperp = SstarUPX = SstarUX
        _Jcap = J(Clust.N, Repl.kZ, 0)
        if (granular==0)
          SCTcapuXinvXX = smatrix(Repl.kZ+1, Nstar)
        if (LIML | robust==0) {
          YYstar = YPXYstar = SstarUX
          YYstar_b = YPXYstar_b = J(Repl.kZ+1, Repl.kZ+1, 0)
        }
        if (NFE & (bootstrapt | kappa != 1 | LIML))
          CTFEU = SstarUX
      }

    } else {  // the score bootstrap for IV/GMM uses a IV/GMM DGP but then masquerades as an OLS test because most factors are fixed during the bootstrap. To conform, need DGP and Repl objects with different R, R1, one with FWL, one not

      Repl = boottestIVGMM(); Repl.parent = &this
      Repl.LIML = this.LIML; Repl.Fuller = this.Fuller; Repl.kappa = this.kappa
      Repl.SetR(*pR1, I(kZ))  // process replication restraints = model constraints only
      Repl.InitVars(&Repl.R1perp)
      Repl.Estimate(*pr1)  // bit inefficient to estimate in both objects, but maintains the conformity
      Repl.InitTestDenoms()
      DGP = boottestIVGMM(); DGP.parent = &this
      DGP.LIML = this.LIML; DGP.Fuller = this.Fuller; DGP.kappa = this.kappa
      DGP.SetR(null? *pR1 \ *pR : *pR1, J(0,kZ,0))  // DGP constraints: model constraints + null if imposed
      DGP.InitVars()
      pM = &Repl  // estimator object from which to get A, AR, XAR; DGP follows WRE convention of using FWL, Repl follows OLS convention of not; scoreBS for IV/GMM mixes the two
      if (null==0) {  // if not imposing null, then DGP constraints, kappa, Hessian, etc. do not vary with r and can be set now
        DGP.Estimate(*pr1)
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
    multiplier = (smallsample = (_Nobs - kZ - FEdfadj) / (_Nobs - robust)) / df  // divide by # of constraints because F stat is so defined
  else
    multiplier = smallsample = 1

  if ((robust | ML)==0)
    multiplier = multiplier * _Nobs  // will turn sum of squared errors in denom of t/z into mean
  if (sqrt) multiplier = sqrt(multiplier)

  if (bootstrapt & (WREnonARubin | df > 1+robust | MaxMatSize<.))  // unless nonWRE or df=1 or splitting weight matrix, code will create Dist element-by-element, so pre-allocate vector now
    pDist = &J(1, B+1, .)
  if (Nw>1 | WREnonARubin | (null==0 & df<=2))
    pnumer = &J(df, B+1, .)

  if (WREnonARubin==0) {
    poles = anchor = J(0,0,0)
    if (interpolable = bootstrapt & B & null & Nw==1 & (kappa==0 | ARubin)) {    
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
}

// main routine
void boottest::boottest() {
  real scalar w

  if (initialized==0)
    Init()
  else if (null==0) {
    NoNullUpdate()
    return
  }

  if (Nw > 1) {
    rseed(seed)
    MakeWildWeights(WeightGrpStop[1] - 1, 1)
  }
  
  if (WREnonARubin)
    PrepWRE()
  else
    MakeInterpolables()  // make stuff that depends linearly on r, possibly by interpolating, for first weight group

  for (w=1; w<=Nw; w++) {  // do group 1 first because it includes col 1, which is all that might need updating in constructing CI in WCU
    if (w > 1)
      MakeWildWeights(WeightGrpStop[w] - WeightGrpStart[w] + 1, 0)

    if (WREnonARubin)
      MakeWREStats(w)
    else
      MakeNonWREStats(w)

    if (bootstrapt==0)
      UpdateBootstrapcDenom(w)
  }
  BFeas = (*pDist)[1]==.? 0 : rownonmissing(*pDist) - 1
  DistCDR = J(0,0,0)
  setdirty(0)
  initialized = 1
}

// if not imposing null and we have returned to boottest(), then df=1 or 2; we're plotting or finding CI, and only test stat, not distribution, changes with r
void boottest::NoNullUpdate() {
  if (WREnonARubin)
    (*pnumer)[,1] = *pR * DGP.Rpar * betas[1] - *pr
  else if (ARubin) {
    DGP.Estimate(*pr)
    (*pnumer)[,1] = v_sd * DGP.Rpar * DGP.beta[|kX1+1\.|] // coefficients on excluded instruments in ARubin OLS
  } else
    (*pnumer)[,1] = v_sd * (*pR * (ML? beta : pM->Rpar * pM->beta) - *pr) // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this

  (*pDist)[1] = df==1? (*pnumer)[1] / sqrt(statDenom) : (*pnumer)[,1] ' invsym(statDenom) * (*pnumer)[,1]
}

// compute bootstrap-c denominator from all bootstrap numerators
void boottest::UpdateBootstrapcDenom(real scalar w) {
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
    pDist = sqrt? &(*pnumer:/sqrt(statDenom)) : &colsum(*pnumer :* invsym(statDenom) * *pnumer)
  }
}

// draw wild weight matrix of width _B. If first=1, insert column of 1s at front. For non-Anderson-Rubin WRE, subtract 1 from all weights
void boottest::MakeWildWeights(real scalar _B, real scalar first) {

  if (_B) {  // in scoretest or waldtest WRE, still make v a col of 1's
    if (enumerate)
      v = J(Nstar,1,1), count_binary(Nstar, -1-WREnonARubin, 1-WREnonARubin)  // complete Rademacher set
    else if (auxwttype==3)
      v = rnormal(Nstar, _B+first, -WREnonARubin, 1)  // normal weights
    else if (auxwttype==4)
      v = rgamma(Nstar, _B+first, 4, .5) :- (2 + WREnonARubin)  // Gamma weights
    else if (auxwttype==2)
      if (WREnonARubin)
        v = sqrt(2 * ceil(runiform(Nstar, _B+first) * 3)) :* ((runiform(Nstar, _B+first):>=.5):-.5) :- 1  // Webb weights, minus 1 for WRE
      else {
        v = sqrt(    ceil(runiform(Nstar, _B+first) * 3)) :* ((runiform(Nstar, _B+first):>=.5):-.5)       // Webb weights, divided by sqrt(2)
        v_sd = 1.6a09e667f3bcdX-001 /*sqrt(.5)*/
      }
    else if (auxwttype)
      if (WREnonARubin)
        v = ( rdiscrete(Nstar, _B+first, 1.727c9716ffb76X-001\1.1b06d1d200914X-002 /*.5+sqrt(.05)\.5-sqrt(.05)*/) :- 1.5 ) * 1.1e3779b97f4a8X+001 /*sqrt(5)*/ :- .5  // Mammen weights, minus 1 for convenience in WRE
      else {
        v = ( rdiscrete(Nstar, _B+first, 1.727c9716ffb76X-001\1.1b06d1d200914X-002 /*.5+sqrt(.05)\.5-sqrt(.05)*/) :- 1.5 ) :+ 1.c9f25c5bfedd9X-003 /*.5/sqrt(5)*/  // Mammen weights, divided by sqrt(5)
        v_sd = 1.c9f25c5bfedd9X-002 /*sqrt(.2)*/
      }
    else if (WREnonARubin) {
      v = runiform(Nstar, _B+first) :<  .5; v = (-2) * v  // Rademacher weights, minus 1 for WRE
    } else {
      v = runiform(Nstar, _B+first) :>= .5; v = v :- .5   // Rademacher weights, divided by 2
      v_sd = .5
    }

    if (first)
      v[,1] = J(Nstar, 1, WREnonARubin? 0 : v_sd)  // keep original residuals in first entry to compute base model stat    
  } else
    v = J(0,1,0)  // in places, cols(v) indicates B -- 1 for classical tests
}


// For WRE, and with reference to Y = [y1 Z], given 0-based column indexes within it, ind1, ind2, return all bootstrap realizations of Y[,ind1]'((1-kappa)*M_Zperp-kappa*M_Xpar)*Y[,ind2] for kappa constant across replications
// ind1 can be a rowvector
// (only really the Hessian when we narrow Y to Z)
real matrix boottest::HessianFixedkappa(real rowvector ind1, real scalar ind2, real scalar kappa) {
  real matrix retval; real scalar i
  if (cols(ind1) > 1) {
    retval = J(cols(ind1),cols(v),0)
    for (i=cols(ind1);i;i--)
      retval[i,] = _HessianFixedkappa(ind1[i], ind2, kappa)
    return(retval)
  }
  return(_HessianFixedkappa(ind1, ind2, kappa))
}
real rowvector boottest::_HessianFixedkappa(real scalar ind1, real scalar ind2, real scalar kappa) {
  real matrix retval, T2; pointer (real colvector) scalar pT1L, pT1R

  if (kappa) {
    pT1L = ind1? pcol(Repl.XZ,ind1) : Repl.pXy1par
    if (Repl.Yendog[ind1+1])
      pT1L = &(*pT1L :+ SstarUX[ind1+1].M * v)
    pT1R = ind2? pcol(Repl.invXXXZ,ind2) : &Repl.invXXXy1par
    if (Repl.Yendog[ind2+1])
      pT1R = &(*pT1R :+ SstarUXinvXX[ind2+1].M * v)  // right-side linear term
    retval = colsum(*pT1L :* *pT1R)  // multiply in the left-side linear term
  }

  if (kappa != 1) {
    if (Repl.Yendog[ind1+1]) {
      T2 = SstarUZperpinvZperpZperp[ind1+1].M ' SstarUZperp[ind2+1].M  // quadratic term
      _diag(T2, diagonal(T2) - (ind1 <= ind2? SstarUU[ind2+1, ind1+1].M : SstarUU[ind1+1, ind2+1].M))  // minus diagonal crosstab
      if (NFE)
        T2 = T2 + CTFEU[ind1+1].M ' (invFEwt :* CTFEU[ind2+1].M)

      retval = kappa? kappa :* retval  + (1 - kappa) :* (Repl.YY[ind1+1,ind2+1] :+ (*pcol(SstaruY[ind2+1].M, ind1+1) + *pcol(SstaruY[ind1+1].M, ind2+1)) ' v - colsum(v :* T2 * v)) :
                                                         Repl.YY[ind1+1,ind2+1] :+ (*pcol(SstaruY[ind2+1].M, ind1+1) + *pcol(SstaruY[ind1+1].M, ind2+1)) ' v - colsum(v :* T2 * v)
    } else
      retval = kappa? kappa :* retval :+ (1 - kappa) * Repl.YY[ind1+1,ind2+1] :
                                                       Repl.YY[ind1+1,ind2+1]
  }
  return(cols(retval)>1? retval : J(1,cols(v),retval))  // if both vars exogenous, term is same for all b; this duplication is a bit inefficient, but only arises when exog vars involved in null
}


// Workhorse for WRE CRVE sandwich filling
// Given column index ind1 and a matrix betas of all the bootstrap estimates, return all bootstrap realizations of P_X * Z[,ind1]_g ' u\hat_1g^*
// for all groups in the intersection of all error clusterings
// return value has one row per cap cluster, one col per bootstrap replication
pointer(real matrix) scalar boottest::Filling(real scalar ind1, real matrix betas) {
  real scalar i, ind2; real matrix retval, T1; pointer (real matrix) scalar pbetav; pointer (real colvector) pPXYstar; real rowvector _beta; real colvector S
  pragma unset retval

  if (granular) {
    if (Nw == 1) {  // create or avoid NxB matrix?
      pPXYstar = pcol(Repl.PXZ, ind1)
      if (Repl.Yendog[ind1+1])
        pPXYstar = &(*pPXYstar :+ SstarUPX[ind1+1].M * v)

      retval = *_panelsum(*pPXYstar :* (Repl.y1 :- SstarUMZperp.M * v), *pinfoCapData)

      for (ind2=Repl.kZ;ind2;ind2--) {
        _beta = betas[ind2,]
        retval = retval - *_panelsum(*pPXYstar :* (Repl.Yendog[ind2+1]? *pcol(*Repl.pZ,ind2) * _beta :- SstarUMZperp[ind2+1].M * (v :* _beta) :
                                                                        *pcol(*Repl.pZ,ind2) * _beta                                           ), *pinfoCapData)
      }
    } else { // create pieces of each N x B matrix one at a time rather than whole thing at once
      retval = J(Clust.N, cols(v), 0)
      for (ind2=0; ind2<=Repl.kZ; ind2++) {
        if (ind2)
          pbetav = &(v :* (_beta = betas[ind2,]))

        if (purerobust) {
          for (i=Clust.N;i;i--) {
            pPXYstar = pcol(Repl.PXZ, ind1)
            if (Repl.Yendog[ind1+1])
              pPXYstar = &(*pPXYstar :+ SstarUPX[ind1+1].M[i,] * v)

            if (ind2)
              retval[i,] = retval[i,] - colsum(*pPXYstar :* (Repl.Yendog[ind2+1]? (*Repl.pZ)[i,ind2] * _beta :- SstarUMZperp[ind2+1].M[i,] * *pbetav :
                                                                                       (*Repl.pZ)[i,ind2] * _beta                                           ))
            else
              retval[i,] =              colsum(*pPXYstar :* (Repl.y1[i] :- SstarUMZperp.M[i,] * v))
          }
        } else {
          for (i=Clust.N;i;i--) {
            S = (*pinfoCapData)[i,]'
            pPXYstar = &Repl.PXZ[|S,(ind1\ind1)|]
            if (Repl.Yendog[ind1+1])
              pPXYstar = &(*pPXYstar :+ SstarUPX[ind1+1].M[|S,(.\.)|] * v)

            if (ind2)
              retval[i,] = retval[i,] - colsum(*pPXYstar :* (Repl.Yendog[ind2+1]? (*Repl.pZ)[|S,(ind2\ind2)|] * _beta :- SstarUMZperp[ind2+1].M[|S,(.\.)|] * *pbetav :
                                                                                       (*Repl.pZ)[|S,(ind2\ind2)|] * _beta                                                 ))
            else
              retval[i,] =              colsum(*pPXYstar :* (Repl.y1[|S|] :- SstarUMZperp.M[|S,(.\.)|] * v))
          }
        }
      }
    }
  } else {  // coarse error clustering
    for (ind2=0; ind2<=Repl.kZ; ind2++) {
      pbetav = ind2? &(v :* (_beta = -betas[ind2,])) : &v

      // T1 * v will be 1st-order terms
      T1 = Repl.Yendog[ind1+1]? Repl.ScapYX[ind2+1].M * SstarUXinvXX[ind1+1].M : J(0,0,0) //  S_∩ (Y_(∥j):*X_∥ ) (X_∥^' X_∥ )^(-1) [S_* (U ̈_(∥i):*X_∥ )]^'

      if (Repl.Yendog[ind2+1]) {  // add CT_(∩,*) (P_(X_∥ ) Y_(∥i):*U ̈_(∥j) )
        if (NClustVar == NBootClustVar & !subcluster)  // simple case of one clustering: full crosstab is diagonal
          if (cols(T1))
            _diag(T1, diagonal(T1) + SstarUXinvXX[ind2+1].M ' (*pcol(Repl.XZ,ind1)))
          else
            T1 =                     SstarUXinvXX[ind2+1].M ' (*pcol(Repl.XZ,ind1))  // keep T1 as vector if it's just going to be a diagonal matrix
        else {
          if (Repl.Yendog[ind1+1]==0)
            T1 = JNcapNstar
          for (i=Nstar;i;i--)
            T1[IDCTCapstar[i].M, i] = T1[IDCTCapstar[i].M, i] + SCTcapuXinvXX[ind2+1,i].M * *pcol(Repl.XZ,ind1)
        }
        if (cols(*Repl.pZperp))  // subtract S_∩ (P_(X_∥ ) Y_(∥i):*Z_⊥ ) (Z_⊥^' Z_⊥ )^(-1) [S_* (U ̈_(∥j):*Z_⊥ )]^'
          T1 = T1 :- Repl.ScapPXYZperp[ind1+1].M * SstarUZperpinvZperpZperp[ind2+1].M
        if (NFE)
          T1 = T1 :- Repl.CT_FEcapPY[ind1+1].M ' CTFEU[ind2+1].M
      }

      if (ind2) {
        retval = retval + Repl.FillingT0[ind1+1,ind2+1].M * _beta
        if (cols(T1))
          retval = retval + (cols(T1)==1? T1 :* *pbetav : T1 * *pbetav)  // - x*beta components
      } else
        retval = Repl.FillingT0[ind1+1,1].M :+ (cols(T1)==1? T1 :* *pbetav : T1 * v)  // y component

      if (Repl.Yendog[ind1+1] & Repl.Yendog[ind2+1])
        for (i=Clust.N;i;i--) {
          S = (*pinfoCapData)[i,]', (.\.)
          retval[i,] = retval[i,] - colsum(v :* cross(SstarUPX[ind1+1].M[|S|], SstarUMZperp[ind2+1].M[|S|]) * *pbetav)
        }
    }
  }
  return(&retval)
}


void boottest::PrepWRE() {
  real scalar i, j, g; pointer (real colvector) scalar puwt, pu

  DGP.Estimate(null? *pr1 \ *pr : *pr1)
  DGP.MakeResiduals()
  pU2parddot = pXB(DGP.U2ddot, Repl.RparY)

  for (i=Repl.kZ; i>=0; i--) {  // precompute various clusterwise sums
    pu = i? pcol(*pU2parddot,i) : &DGP.u1dddot
    puwt = pvHadw(*pu, *pwt)

    // S_star(u :* X), S_star(u :* Zperp) for residuals u for each endog var; store transposed
    SstarUX                   [i+1].M = *_panelsum2(*Repl.pX1, Repl.X2, *puwt, infoBootData)'
    SstarUXinvXX              [i+1].M = Repl.invXX * SstarUX[i+1].M

    if (kappa!=1 | LIML | bootstrapt) {
      SstarUZperp             [i+1].M = *_panelsum(*Repl.pZperp, *puwt, infoBootData)'
      SstarUZperpinvZperpZperp[i+1].M = Repl.invZperpZperp * SstarUZperp[i+1].M
      if (NFE)
        CTFEU[i+1].M = crosstabFE(*puwt, infoBootData)
    }


    if (kappa!=1 | LIML | robust==0) {
      SstaruY[i+1].M = *_panelsum2(*Repl.py1par, *Repl.pZ, *puwt, infoBootData)
      for (j=i; j>=0; j--)
        SstarUU[i+1,j+1].M = *_panelsum(j? *pcol(*pU2parddot,j) : DGP.u1dddot, *puwt, infoBootData)
    }

    if (robust & bootstrapt) {
      if (granular==0 & Repl.Yendog[i+1] & !(NClustVar == NBootClustVar & !subcluster))  // Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u:*X and u:*Zperp (and times invXX or invZperpZperp)
        for (g=Nstar;g;g--)
          SCTcapuXinvXX[i+1,g].M = *_panelsum(Repl.XinvXX, *pu, infoCTCapstar[g].M)
      
      if (i) SstarUPX[i+1].M = Repl.XinvXX * SstarUX[i+1].M
      SstarUMZperp[i+1].M = *Repl.pZperp * SstarUZperpinvZperpZperp[i+1].M
      if (Nobs == Nstar)  // subtract "crosstab" of observation by cap-group of u
        _diag(SstarUMZperp[i+1].M, diagonal(SstarUMZperp[i+1].M) - (i? (*pU2parddot)[,i] : DGP.u1dddot))  // case: bootstrapping by observation
      else
        if (i)
          for (g=Nobs;g;g--)
            SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] = SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] - (*pU2parddot)[g,i]
        else
          for (g=Nobs;g;g--)
            SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] = SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] - DGP.u1dddot[g]

      if (NFE)
        SstarUMZperp[i+1].M = SstarUMZperp[i+1].M + (invFEwt :* CTFEU[i+1].M)[*pFEID,]  // CT_(*,FE) (U ̈_(∥j) ) (S_FE S_FE^' )^(-1) S_FE
    }
  }
}

void boottest::MakeWREStats(real scalar w) {
  real scalar c, b, i
  real colvector numer_b
  real rowvector numerw, val, YY11, YY12, YY22, YPXY11, YPXY12, YPXY22, x11, x12, x21, x22, kappas, YY12YPXY12
  real matrix deltanumer, Jcap, J_b, Jcaps, vec
  struct smatrix rowvector A
  pragma unset vec; pragma unset val

  if (Repl.kZ == 1) {  // optimized code for 1 coefficient in bootstrap regression
    if (LIML) {
      YY11   = HessianFixedkappa(0, 0, 0)  // kappa=0 => Y*MZperp*Y
      YY12   = HessianFixedkappa(0, 1, 0)
      YY22   = HessianFixedkappa(1, 1, 0)
      YPXY11 = HessianFixedkappa(0, 0, 1)  // kappa=1 => Y*PXpar*Y
      YPXY12 = HessianFixedkappa(0, 1, 1)
      YPXY22 = HessianFixedkappa(1, 1, 1)
      YY12YPXY12 = YY12 :* YPXY12
      x11 = YY22 :* YPXY11 - YY12YPXY12      // elements of YYstar^-1 * YPXYstar up to factor of det(YYstar)
      x12 = YY22 :* YPXY12 - YY12 :* YPXY22
      x21 = YY11 :* YPXY12 - YY12 :* YPXY11
      x22 = YY11 :* YPXY22 - YY12YPXY12
      kappas = .5 * (x11 + x22); kappas = 1 :/ (1 :- (kappas - sqrt(kappas:*kappas - x11:*x22 + x12:*x21)) :/ (YY11 :* YY22 - YY12 :* YY12))  // solve quadratic equation for smaller eignenvalue; last term is det(YYstar)
      if (Fuller) kappas = kappas :- Fuller / (_Nobs - kX)
      betas = (kappas :* (YPXY12 - YY12) + YY12) :/ (As = kappas :* (YPXY22 - YY22) + YY22)
    } else
      betas = HessianFixedkappa(1, 0, kappa) :/ (As = HessianFixedkappa(1, 1, kappa))

    if (null)
      numerw = betas :+ (Repl.Rt1 - *pr) / Repl.RRpar
    else {
      numerw = betas :- DGP.beta
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
        if (Fuller) kappa = kappa - Fuller / (_Nobs - kX)
        betas[,b] = (A[b].M = invsym(kappa*YPXYstar_b[|2,2\.,.|] + (1-kappa)*YYstar_b[|2,2\.,.|])) * (kappa*YPXYstar_b[|2,1\.,1|] + (1-kappa)*YYstar_b[|1,2\1,.|]')
      }
    } else {
      deltanumer = HessianFixedkappa(1..Repl.kZ, 0, kappa)

      for (i=Repl.kZ;i;i--)
        deltadenom[i].M = HessianFixedkappa(1..i, i, kappa)

      for (b=cols(v); b; b--) {
        for (i=Repl.kZ;i;i--)
          deltadenom_b[|.,i\i,i|] = deltadenom[i].M[,b] // fill uppper triangle, which is all that invsym() looks at
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
            YYstar_b[|i+1,i+1\Repl.kZ+1,i+1|] = YYstar[i+1].M[,b]  // fill lower triangle for makesymmetric()
          denom.M = (Repl.RRpar * A[b].M * Repl.RRpar') * ((-1 \ betas[,b]) ' makesymmetric(YYstar_b) * (-1 \ betas[,b]) / _Nobs)  // 2nd half is sig2 of errors
        }
        (*pDist)[b+WeightGrpStart[w]-1] = sqrt? numer_b/sqrt(denom.M) : cross(numer_b, invsym(denom.M) * numer_b)  // hand-code for 2-dimensional?
      }
      (*pnumer)[,b+WeightGrpStart[w]-1] = numer_b  // slight inefficiency: in usual bootstrap-t case, only need to save numerators in numer if getdist("numer") is coming because of svmat(numer)
    }
  }

  if (w==1 & bootstrapt) statDenom = denom.M  // original-sample denominator
}


// For non-WRE, construct stuff that depends linearly or quadratically on r, possibly by interpolation
void boottest::MakeInterpolables() {
  real scalar h1, h2, d1, d2, c; real matrix tmp; real colvector Delta, newPole

  if (interpolable) {
    if (rows(anchor)==0) {  // first call? save current r as permanent anchor for interpolation
      _MakeInterpolables(anchor = *pr)
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
          _MakeInterpolables(tmp)  // calculate linear stuff at new anchor

          dnumerdr[h1].M = (*pnumer - numer0) / poles[h1]
          if (interpolate_u)
            dudr[h1].M = (*puddot - uddot0) / poles[h1]
          if (robust & !purerobust)  // df > 1 for an ARubin test with >1 instruments. 
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
      if (robust & !purerobust)  // quadratic interaction terms
        for (h1=1;h1<=q;h1++)
          for (h2=h1;h2;h2--)
            if (newPole[h1] | newPole[h2])
              for (d1=df;d1;d1--)
                for (d2=d1;d2;d2--)
                  for (c=1;c<=NErrClustCombs;c++)
                    _clustAccum(ddenomdr2[h1,h2].M[d1,d2].M, c, colsum(dJcddr[h1].M[c,d1].M :* dJcddr[h2].M[c,d2].M))
      Delta = poles
      interpolating = 1

      if (q==2) {  // in this case we haven't yet actually computed interpolables at *pr, so interpolate them
        numerw = numer0 + dnumerdr.M * Delta[1] + dnumerdr[2].M * Delta[2]
        if (interpolate_u) {
          puddot = &(uddot0 + dudr.M * Delta[1] + dudr    [2].M * Delta[2])
        }
      }
    } else {  // routine linear interpolation if the anchors not moved
      Delta = *pr - anchor
      numerw = numer0 + dnumerdr.M * Delta[1]; if (q > 1) numerw = numerw + dnumerdr[2].M * Delta[2]
      if (interpolate_u) {
        puddot = &(uddot0 + dudr.M * Delta[1]); if (q > 1) puddot = &(*puddot + dudr    [2].M * Delta[2])
      }
    }

    if (robust & !purerobust)  // even if an anchor was just moved, and linear components just computed from scratch, do the quadratic interpolation now, from the updated linear factors
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
    _MakeInterpolables(*pr)
}

// Construct stuff that depends linearly or quadratically on r and doesn't depend on v. No interpolation.
void boottest::_MakeInterpolables(real colvector r) {
  real scalar d, c; pointer (real matrix) scalar pustarXAR, ptmp

  if (ML)
    uXAR = *pSc * (AR = *pA * *pR')
  else {
    if (ARubin)
      DGP.Estimate(r)
    else if (kappa) {
      if (null) {  // in score bootstrap for IV/GMM, if imposing null, then DGP constraints, kappa, Hessian, etc. do vary with r and must be set now
        DGP.Estimate(*pr1 \ r)
        DGP.InitTestDenoms()
      }
    } else  // regular OLS
      DGP.Estimate(null? *pr1 \ r : *pr1)

    DGP.MakeResiduals()
    puddot = &DGP.u1ddot

    if (scoreBS | (robust & granular < NErrClustCombs))
      uXAR = DGP.u1ddot :* pM->XAR
  }

  SuwtXA = scoreBS?
              (B? 
                 (NClustVar? *_panelsum(uXAR, *pwt, infoBootData) : 
                             *pvHadw(uXAR, *pwt)                  ) :
                 cross(*pwt, uXAR)')                             :
              *DGP.pA * *_panelsum2(*pX1, *pX2, *pvHadw(*puddot, *pwt), infoBootData)'  // same calc as in score BS but broken apart to grab intermediate stuff, and assuming residuals defined; X2 empty except in Anderson-Rubin

  if (robust & granular < NErrClustCombs & bootstrapt) {
    pustarXAR = _panelsum(uXAR, *pwt, *pinfoAllData)  // collapse data to all-boot & error-cluster-var intersections. If no collapsing needed, _panelsum() will still fold in any weights

    if (B) {
      if (scoreBS)
        for (d=df;d;d--)
          Kd[d].M = JNcapNstar  // inefficient, but not optimizing for the score bootstrap
      else
        for (d=df;d;d--)
          Kd[d].M = *_panelsum2(*pX1, *pX2, *pvHadw(*pcol(DGP.XAR,d), *pwt), *pinfoCapData) * SuwtXA  // final term in (64), for c=intersection of all error clusters
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
  MakeNumerAndJ(1, r)  // compute J = kappa * v; if Nw > 1, then this is for 1st group; if interpolating, it is only group, and may be needed now to prep interpolation
}

// compute stuff depending linearly on v, needed to prep for interpolation
void boottest::MakeNumerAndJ(real scalar w, | real colvector r) {  // called to *prepare* interpolation, or when w>1, in which case there is no interpolation
  real scalar c, d; real matrix _v

  numerw = scoreBS?
             (B? 
               cross(SuwtXA, v) : 
               SuwtXA * v_sd    ) :
             (robust==0 | granular | purerobust?
                *pR * (betadev = SuwtXA * v) :
               (*pR * SuwtXA) * v)

  if (w==1) {
    if      ( ARubin) numerw[,1] = v_sd * DGP.Rpar * DGP.beta[|kX1+1\.|]  // coefficients on excluded instruments in ARubin OLS
    else if (null==0) numerw[,1] = v_sd * (*pR * (ML? beta : pM->Rpar * pM->beta) - r)  // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this.
  }

  storeWtGrpResults(pnumer, w, numerw)

  if (B & robust & bootstrapt) {
    if (granular | purerobust)  // optimized treatment when bootstrapping by many/small groups
      if (purerobust)
        ustar = *partialFE(&(*puddot :* v)) - *pX12B(*pX1, *pX2, betadev)
      else {  // clusters small but not all singletons
        if (NFE & FEboot==0) {
          ustar = *partialFE(&(*puddot :* v[*pIDBootData,]))
          for (d=df;d;d--)
            (*pJcd)[1,d].M = *_panelsum(ustar, pM->WXAR[d].M, *pinfoCapData)                                 - *_panelsum2(*pX1, *pX2, pM->WXAR[d].M, *pinfoCapData) * betadev
        } else {
          _v = v[*pIDBootAll,]
          for (d=df;d;d--)
            (*pJcd)[1,d].M = *_panelsum(*_panelsum(*puddot, pM->WXAR[d].M, *pinfoAllData) :* _v, infoErrAll) - *_panelsum2(*pX1, *pX2, pM->WXAR[d].M, *pinfoCapData) * betadev
        }
      }
    for (c=NErrClustCombs; c>granular; c--)
      for (d=df;d;d--)
        (*pJcd)[c,d].M = Kcd[c,d].M * v
  }
}

void boottest::MakeNonWREStats(real scalar w) {
  real scalar i, c, j, k; real matrix ustar2, tmp, invdenom; real colvector numer_l; pointer (real matrix) scalar pAR; real rowvector t1, t2, t12

  if (w > 1) MakeNumerAndJ(w)

  if (bootstrapt == 0) return

  if (robust) {
    if (interpolating==0) {  // these quadratic computations needed to *prepare* for interpolation but are superseded by interpolation once it is going
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
      for (k=cols(v); k; k--) {
        for (i=df;i;i--)
          for (j=i;j;j--)
            tmp[j,i] = denom[i,j].M[k]  // fill upper triangle, which is all invsym() looks at
        numer_l = numerw[,k]
        (*pDist)[k+WeightGrpStart[w]-1] = numer_l ' invsym(tmp) * numer_l  // in degenerate cases, cross() would turn cross(.,.) into 0
      }
      if (w==1)
        statDenom = makesymmetric(tmp')  // original-sample denominator
    }

  } else { // non-robust

    pAR = ML? &AR : &pM->AR
    if (df == 1) {  // optimize for one null constraint
      denom.M = *pR * *pAR

      if (ML==0) {
                     ustar = B? v :* *puddot : *puddot
        if (scoreBS) ustar = ustar :- (haswt? cross(ClustShare, ustar) : colsum(ustar) * ClustShare)  // Center variance if interpolated
                else ustar = ustar  - *pX12B(*pX1, *pX2, betadev)  // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
        denom.M = denom.M :* (haswt? cross(*pwt, ustar :* ustar) : colsum(ustar :* ustar))
      }
      storeWtGrpResults(pDist, w,  numerw :/ sqrt(denom.M))
      if (w==1)
        statDenom = denom.M[1]  // original-sample denominator
    } else {
      denom.M = *pR * *pAR

      if (ML) {
        for (k=cols(v); k; k--) {
          numer_l = numerw[,k]
          (*pDist)[k+WeightGrpStart[w]-1] = cross(numer_l, invsym(denom.M), numer_l)
        }
        if (w==1)
          statDenom = denom.M  // original-sample denominator
      } else {
        invdenom = invsym(denom.M)
        for (k=cols(v); k; k--) {
          numer_l = numerw[,k]
          (*pDist)[k+WeightGrpStart[w]-1] = cross(numer_l, invdenom * numer_l)
                       ustar = B? v[,k] :* *puddot : *puddot
          if (scoreBS) ustar = ustar :- (haswt? cross(*pwt, ustar) : colsum(ustar)) * ClustShare  // Center variance if interpolated
                  else ustar = ustar  - *pX12B(*pX1, *pX2, betadev[,k])  // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
          (*pDist)[k+WeightGrpStart[w]-1] = (*pDist)[k+WeightGrpStart[w]-1] / (tmp = cross(ustar, *pwt, ustar))
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
    return(arg2==1? &X : &(X :* arg2))  // if no collapsing called for, still fold in provided weights
  return (cols(arg3)? &panelsum(X, arg2, arg3) : &panelsum(X, arg2))
}

// concatenation of two _panelsum's
pointer(real matrix) scalar _panelsum2(real matrix X1, real matrix X2, real matrix arg2, | real matrix arg3)
  return(args()==2? &(*_panelsum(X1,arg2),*_panelsum(X2,arg2)) : &(*_panelsum(X1,arg2,arg3),*_panelsum(X2,arg2,arg3)))

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

// given a pre-configured boottest linear model with one-degree null imposed, compute p value for given value of r
// used with optimize() to construct confidence intervals
// performs no error checking
real scalar boottest::r_to_p(real colvector r) {
  pr = &r
  setdirty(1, 1) // set dirty = 1, but leave initialized=0, which we want when only changing r
  return (getpadj())
}


// Chandrupatla 1997, "A new hybrid quadratic-bisection algorithm for finding the zero of a nonlinear function without using derivatives"
// x1, x2 must bracket the true value, with f1=f(x1) and f2=f(x2). alpha is target p value.
real scalar boottest::search(real scalar alpha, real scalar f1, real scalar x1, real scalar f2, real scalar x2) {
  real scalar x, t, fx, phi1, phi1_2, xi1, x3, f3
  
  t = 0.5
  do {
    fx = r_to_p(x = x1 + t * (x2 - x1))

    if (fx>f1 == fx>f2)  // violation of monotonicity because of precision problems? That's as good as it gets.
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

// derive wild bootstrap-based CI, for case of linear model with one-degree null imposed
// also prepare data to plot confidence curve or surface
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

  if (q==2) {  // 2D plot
    lo = hi = J(2, 1, .)
    for(d=df;d;d--) {
      lo[d] = editmissing(gridmin[d], confpeak[d] - halfwidth[d])
      hi[d] = editmissing(gridmax[d], confpeak[d] + halfwidth[d])

      stata("_natscale " + strofreal(lo[d]) + " " + strofreal(hi[d]) + " 4")  // using Stata logic for setting graph bounds ensures good-looking contour plot
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
        tmp = sqrt(statDenom[1,1]) * (small? invttail(df_r, alpha/2) : -invnormal(alpha/2))
        lo = editmissing(gridmin[1], confpeak - tmp)
        hi = editmissing(gridmax[1], confpeak + tmp)
        if (scoreBS & (null | willplot)==0) {  // if doing simple Wald test with no graph, we're done
          CI = lo, hi
          return
        }
      }
     
      if (abs(lo - *pr) > abs(hi - *pr)) {  // brute force way to ensure that first trial bound tested is the farther one from *pr, for better interpolation
        if (gridmin[1]==. & ptype!=2)  // unless lower-tailed p value, try at most 10 times to bracket confidence set by symmetrically widening
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
    plotY[1] = p_lo; plotY[rows(plotX)] = p_hi
    p_confpeak = WRE? . : (twotailed? 1 : .5)
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
  }  // end 1D plot prep

  i = 1
  do {
    if (plotY[i] == .) plotY[i] = r_to_p(plotX[i,]')
  } while (1 < (i = mod(i-2,rows(plotX))+1))

  if (hasmissing(plotY))
    CI = . , .
  else if (q==1 & level<100) {  // find CI bounds
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
  real scalar kappa, real scalar ARubin, real scalar null, real scalar scoreBS, string scalar auxweighttype, string scalar ptype, string scalar statistic, string scalar madjtype, real scalar NumH0s,
  string scalar X1names, string scalar Y2names, string scalar Ynames, string scalar bname, string scalar Aname, 
  string scalar X2names, string scalar samplename, string scalar scnames, real scalar robust, string scalar IDnames, real scalar NBootClustVar, real scalar NErrClustVar, 
  string scalar FEname, real scalar NFE, real scalar FEdfadj, string scalar wtname, string scalar obswttype, string scalar R1name, string scalar r1name, string scalar Rname, string scalar rname, real scalar B, string scalar repsname, string scalar repsFeasname, 
  real scalar small, string scalar diststat, string scalar distname, string scalar gridmin, string scalar gridmax, string scalar gridpoints, real scalar MaxMatSize, real scalar quietly,
  string scalar b0name, string scalar V0name, string scalar vname, string scalar NBootClustname) {
  real matrix X2, ID, FEID, sc, Y2, X1
  real colvector wt, Y
  class boottest scalar M
  pragma unset ID; pragma unset wt; pragma unset Y2; pragma unset X1; pragma unset Y; pragma unset X2; pragma unset sc

  M._st_view(sc, ., scnames, samplename)
  M._st_view(Y , ., Ynames , samplename)
  M._st_view(X2, ., X2names, samplename)
  if (FEname != "" ) FEID = st_data(., FEname , samplename)
  if (IDnames != "") ID   = st_data(., IDnames, samplename)
  if (wtname  != "") wt   = st_data(., wtname , samplename) // panelsum() doesn't like views as weights
  M.setMaxMatSize(MaxMatSize)
  M.setsc(sc)
  M.setML(ML)
  M.setY (Y)
  M.setX2(X2)
  M.setobswt(wt, obswttype)
  M.setID(ID, NBootClustVar, NErrClustVar)
  M.setFEID(FEID, NFE, FEdfadj)
  M.setR1(st_matrix(R1name), st_matrix(r1name))
  M.setR (st_matrix(Rname ), st_matrix(rname ))
  M.setnull(null)
  M.setsmall(small)
  M.setrobust(robust)
  M.setscoreBS(scoreBS)
  M.setauxwttype(auxweighttype)
  M.setptype(ptype)
  M.setstattype(statistic)
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
  M._st_view(Y2, ., Y2names, samplename); M.setY2(Y2)
  M._st_view(X1, ., X1names, samplename); M.setX1(X1)
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
