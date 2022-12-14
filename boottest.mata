*! boottest 4.3.1 10 December 2022
*! Copyright (C) 2015-22 David Roodman

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
  real colvector is, wt, sqrtwt
}

// return pointer to chosen columns of a matrix, but don't duplicate data if return value is whole matrix
pointer (real matrix) pcol(real matrix A, real vector p)
  return(length(p)==cols(A)? &A : &A[,p])

// right-multiply a data matrix by a matrix, with efficient handling of special cases like latter is identity
pointer (real matrix) scalar pXB(real matrix X, real matrix M) {
  real scalar r
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
pointer(real matrix) scalar pXS(real matrix X, real colvector S)
  return(cols(X)? (S[2,1]? ((S[1,1]==. | S[1,1]==1) & (S[2,1]==. | S[2,1]==rows(X))? &X : &X[|S,(.\.)|]) : &J(0,cols(X),0)) : &J(editmissing(S[2,1],rows(X))-editmissing(S[1,1],1)+1,0,0))

// do X[|S|] = Y while allowing X to have no cols and S to be a colvector
void setXS(real matrix X, real colvector S, real matrix Y) if (cols(X)) X[|S,(.\.)|] = Y;;

matrix fold(matrix X) return(uppertriangle(X) + lowertriangle(X,0)')  // fold matrix diagonally; returns same values as a quad form, but runs faster because of all the 0's

class boottestOLS {  // class for performing OLS
  real scalar y1y1, LIML, Fuller, kappa, isDGP, kZ, kX1, kX, y1bary1bar
  real colvector y1, invXXXy1par, X1y1, dbetadr, beta0, y1bar, Zperpy1, t1, t1Y, deltadddot, X2y1, Xy1bar, deltaX, deltaY
  real rowvector y1Y2, Yendog, y1barU2ddot
  real matrix Z, ZZ, invZperpZperp, XZ, XX, PXZ, R1invR1R1, R1perp, Rpar, RperpX, RRpar, RparX, RparY, RR1invR1R1, YY, AR, XAR, invXXXZ, XinvXX, Rt1, invXX, Y2, X2, invH, Deltadddot, Y2bar, perpRperpX, ZperpZperp, ZperpX1, ZperpX2, ZperpY2, XY2, XU2ddot, X1X, X2X, U2ddotU2ddot, Y2Y2, Piddot, ZY2
  pointer(real colvector) scalar py1par, pXy1par
  pointer(real matrix) scalar pA, pZperp, pX1, pX1perpRperpX, pX1par, pR1AR1, pZR1
  pointer (class boottest scalar) scalar parent
  struct smatrix rowvector WXAR, CT_XAR, beta, u1ddot, XinvHg, invMg, Xg, XXg, u1dddot, U2ddot

  private void new(), InitTestDenoms()
  private virtual void InitVars(), SetR(), Estimate(), MakeResiduals()
  real matrix _select(), perp()
}

class boottestARubin extends boottestOLS {
  private virtual void InitVars(), Estimate()
}

class boottestIVGMM extends boottestOLS {
  real matrix H_2SLS, V, ZR1ZR1, X2ZR1, ZR1Y2, X1ZR1, ZZR1, ZXinvXXXZ, H_2SLSmZZ, X2jk, Y2jk, X1jk, Zjk, ZR1jk
  real colvector Zy1, y1jk
  real rowvector twoy1ZR1
  real scalar y1pary1par
  pointer(real rowvector) scalar pZy1par
  pointer(real rowvector) scalar py1parY2
  pointer(real matrix) scalar pRperp

  pointer(struct smatrix rowvector) scalar py1parY2g, pZy1parg, pXy1parg  // jk stuff
  struct smatrix rowvector XY2g, ZY2g, XXg, XZg, YYg, Zy1g, X1y1g, X2y1g, y1Y2g, ZZg, invXXg, H_2SLSg, H_2SLSmZZg, ZR1Y2g, ZR1ZR1g, twoy1ZR1g, ZZR1g, X1ZR1g, X2ZR1g, ZXinvXXXZg, invHg
  real colvector _Nobsg, y1pary1parg, y1y1g
  pointer(real colvector) scalar py1pary1parg, py1parjk

  private void new()
  private virtual void InitVars(), Estimate(), MakeResiduals(), MakeH()
}

class boottest {
  real scalar scoreBS, B, small, auxwttype, null, dirty, initialized, ML, Nobs, _Nobs, kZ, kY2, kX1, sumwt, NClustVar, haswt, REst, multiplier, smallsample, quietly, FEboot, NErrClustCombs, ///
    sqrt, LIML, Fuller, kappa, WRE, WREnonARubin, ptype, twotailed, df, df_r, ARubin, willplot, notplotted, NumH0s, p, NBootClustVar, NErrClustVar, BootClust, FEdfadj, jk, ///
    NFE, granular, granularjk, purerobust, subcluster, Nstar, BFeas, v_sd, level, ptol, MaxMatSize, Nw, enumerate, bootstrapt, q, interpolable, interpolating, interpolate_u, robust, kX2, kX, Ncap, bootneqcap
  real matrix AR, v, CI, CT_WE, infoBootData, infoErrAll, JNcapNstar, statDenom, SuwtXA, numer0, deltadenom_b, _Jcap, YYstar_b, YPXYstar_b, numerw, ustar0, YbarYbar, XYbar, invXXXZbar, PXZbar, Zbar, invXXXX1par
  real colvector DistCDR, plotX, plotY, beta, ClustShare, WeightGrpStart, WeightGrpStop, confpeak, gridmin, gridmax, gridpoints, numersum, anchor, poles, invFEwt, sqrtwt
  real rowvector peak, betas, As
  string scalar obswttype, madjtype, seed
  pointer (real matrix) scalar pX2, pR1, pR, pID, pFEID, pY2, pX1, pinfoAllData, pinfoCapData, pIDAll, pnumer, pustar, pU2parddot
  pointer (real colvector) scalar pr1, pr, py1, pSc, pwt, pA, pDist, pIDBootData, pIDBootAll
  class boottestOLS scalar DGP, Repl
  pointer(class boottestOLS scalar) scalar pM
  struct structboottestClust rowvector Clust
  struct smatrix matrix denom, Kcd, denom0, Jcd0, CTUX, FillingT0, SstarUU
  struct smatrix rowvector Kd, dudr, dnumerdr, IDCTCapstar, infoCTCapstar, deltadenom, Zyi, SstarUPX, YYstar, YPXYstar, ScapYbarX, ScapPXYbarZperp, CT_FEcapYbar, SstarUX, SstarUXinvXX, SstarUZperpinvZperpZperp, SstarUZperp, CTFEU, SstarUYbar, SCTcapUXinvXX, SstarUMZperp, ScapXX, SbootY2X, SbootY2parX, ScapXZperp, ScapX1parReplX, SbootZR1X, Sbooty1X, SbootZX
  pointer(struct smatrix rowvector) pSbootXX
  struct ssmatrix rowvector ddenomdr, dJcddr
  
  struct ssmatrix matrix ddenomdr2
  pointer(struct smatrix matrix) scalar pJcd
  struct structFE rowvector FEs

  void new(), setsqrt(), setX1(), setptype(), setdirty(), setY2(), setY(), setX2(), setobswt(), setsc(), setML(), setLIML(), setARubin(), setauxwttype(),
    setFuller(), setkappa(), setquietly(), setbeta(), setA(), setsmall(), sethascons(), setjk(), setscoreBS(), setB(), setnull(), setID(), setFEID(), setlevel(), setptol(), 
    setrobust(), setR1(), setR(), setwillplot(), setgrid(), setmadjust(), setMaxMatSize(), setstattype(), close()
  private void MakeNumerAndJ(), _clustAccum(), MakeWREStats(), MakeInterpolables(), _MakeInterpolables(), MakeNonWREStats(), UpdateBootstrapcDenom(), Init(), plot(), MakeWildWeights(), boottest(), crosstabCapstarMinus(), InitWRE(), PrepWRE(), storeWtGrpResults(), NoNullUpdate()
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
  beta = smatrix()
}

void boottestIVGMM::new() {
  Fuller = 0
  kappa = isDGP = 1
  u1ddot = beta = smatrix()
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
  real matrix _R, vec; real rowvector val
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
    _R = R \ J(parent->kY2, parent->kX1, 0), I(parent->kY2)  // rows to prevent partialling out of endogenous regressors
    if (rows(R1))
      _R = _R * R1perp 
    symeigensystem(_R ' invsym(_R * _R') * _R, vec, val); _edittozero(val, 1000)
    Rpar = _select(vec, val); if (rows(R1)) Rpar = R1perp * Rpar  // par of RR₁perp

    if (!(isDGP & parent->WREnonARubin)) {  // in WRE, Rperp is same DGP and Repl; might not be in WUE, but is arranged so by call to SetR()
      RperpX = _select(vec, !val)
      if (rows(R1)) RperpX = R1perp * RperpX
      RperpX = *pXS(RperpX, .\parent->kX1)
    }

    RparX = *pXS(Rpar, .\parent->kX1)  // part of Rpar that refers to X1
    RparY = *pXS(Rpar, parent->kX1+1\.)  // part of Rpar that refers to Y2
    if (isDGP==0) {
      RRpar      = R * Rpar
      RR1invR1R1 = R * R1invR1R1
    }
  }
}

// stuff that can be done before r set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
void boottestOLS::InitVars(pointer(real matrix) scalar pRperp) {  // pRperp is for replication regression--no null imposed
  real matrix H, negR1AR1; real scalar g; real colvector S

  u1ddot = smatrix(1 + parent->jk)  // for jackknife, need jk'd residuals but also non-jk'd residuals for original test stat

  py1par = parent->py1
  H = cross(*parent->pX1, *parent->pX1)
  invH = invsym(H)
  pR1AR1 = rows(R1perp)? &( R1perp * invsym( R1perp ' H *  R1perp) *  R1perp') : &invH  // for DGP regression
  X1y1 = cross(*parent->pX1, *py1par)
  beta0  = *pR1AR1 * X1y1
  dbetadr = *pR1AR1 * H * R1invR1R1 - R1invR1R1

  if (parent->jk)
    if (parent->purerobust)  // set up for optimized "robust" jk
      (invMg = smatrix()).M = 1 :/ (1 :- rowsum(*parent->pX1 * fold(*pR1AR1) :* *parent->pX1))  // standard hc3 multipliers
    else {
      u1ddot[2].M = J(parent->Nobs, 1, 0)
      Xg = smatrix(parent->Nstar)
      if (parent->granularjk)
        invMg = Xg
      else
        XXg = XinvHg = smatrix(parent->Nstar)
      negR1AR1 = - *pR1AR1
      for (g=parent->Nstar; g; g--) {
        S = parent->NClustVar? parent->infoBootData[g,1] \ parent->infoBootData[g,2] : g\g
        Xg[g].M = *pXS(*parent->pX1,S)
        if (parent->granularjk) {  // for many small clusters, faster to compute jackknife errors vaia hc3-like formula
          invMg[g].M = Xg[g].M * negR1AR1 * Xg[g].M'
          _diag(invMg[g].M, diagonal(invMg[g].M) :+ 1)
          invMg[g].M  = invsym(invMg[g].M)
        } else {
          XXg[g].M = cross(Xg[g].M, Xg[g].M)
          XinvHg[g].M = Xg[g].M * (rows(R1perp)? R1perp * invsym(R1perp ' (H - XXg[g].M) *  R1perp) *  R1perp' : invsym(H - XXg[g].M))
        }
      }
    }

  pA = rows(*pRperp)? &(*pRperp * invsym(*pRperp ' H * *pRperp) * *pRperp') : &invH  // for replication regression
  AR = *pA * *parent->pR'
  if (parent->scoreBS | parent->robust)
    XAR = *parent->pX1 * AR
}

void boottestARubin::InitVars(| pointer(real matrix) pRperp) {
  pragma unused pRperp
  real matrix H, X2X1, _X2X1, _X1, _X2, tmp; real colvector S, X1y1, X2y1; real scalar g

  u1ddot = smatrix(1 + parent->jk)  // for jackknife, need jk'd residuals but also non-jk'd residuals for original test stat

  X2X1 = cross(*parent->pX2, *parent->pX1)
  H = cross(*parent->pX1, *parent->pX1), X2X1' \ X2X1, cross(*parent->pX2, *parent->pX2)
  pA = &invsym(H)
  pR1AR1 = rows(R1perp)? &(R1perp * invsym(R1perp ' H * R1perp) * R1perp') : pA
  X1y1 = cross(*parent->pX1, *parent->py1)
  X2y1 = cross(*parent->pX2, *parent->py1)
  beta0   = *pR1AR1 * (X1y1 \ X2y1)
  dbetadr = *pR1AR1 * (cross(*parent->pX1, *parent->pY2) \ cross(*parent->pX2, *parent->pY2))
  
  if (parent->jk)
    if (parent->purerobust) {  // set up for optimized "robust" jk
      tmp = fold(*pR1AR1)
      (invMg=smatrix()).M = parent->kX1? rowsum(*parent->pX1 * tmp[|.,.\parent->kX1,parent->kX1|] :* *parent->pX1) +
                             rowsum(*parent->pX1 * tmp[|.,parent->kX1+1\parent->kX1,.|] :* *parent->pX2) +
                             rowsum(*parent->pX2 * tmp[|parent->kX1+1,parent->kX1+1\.,.|] :* *parent->pX2) :
                             rowsum(*parent->pX2 * tmp                                    :* *parent->pX2)
      invMg.M = 1 :/ (1 :- invMg.M)  // standard hc3 multipliers
    } else {
      u1ddot[2].M = J(parent->Nobs, 1, 0)
      Xg = smatrix(parent->Nstar)
      if (parent->granularjk)
        invMg = Xg
      else
        XXg = XinvHg = smatrix(parent->Nstar)

      for (g=parent->Nstar; g; g--) {
        S = parent->NClustVar? parent->infoBootData[g,1] \ parent->infoBootData[g,2] : g\g
        _X1 = *pXS(*parent->pX1, S)
        _X2 = *pXS(*parent->pX2, S)
        Xg[g].M = _X1, _X2
        if (parent->granularjk) {  // for many small clusters, faster to compute jackknife errors vaia hc3-like formula
          invMg[g].M = -(Xg[g].M * *pR1AR1 * Xg[g].M')
          _diag(invMg[g].M, diagonal(invMg[g].M) :+ 1)
          invMg[g].M  = invsym(invMg[g].M)
        } else {
          _X2X1 = cross(_X2, _X1)
          XXg[g].M = cross(_X1, _X1), _X2X1' \ _X2X1, cross(_X2, _X2)
          XinvHg[g].M = Xg[g].M * (rows(R1perp)? R1perp * invsym(R1perp ' (H - XXg[g].M) *  R1perp) *  R1perp' : invsym(H - XXg[g].M))
        }
      }
    }

  AR = *pA * *parent->pR'  // for replication regression
  if (parent->scoreBS | parent->robust)
    XAR = *pX12B(*parent->pX1, *parent->pX2, AR)
}

void boottestIVGMM::InitVars(|pointer(real matrix) scalar pRperp) {
  real matrix ZperpZ, ZperpZR1, parentY2g, parentX2g, _ZperpZperp, _invZperpZperp, _ZperpX1, _ZperpX2, _ZperpZ, _ZperpZR1, _ZperpY2, _X1, _X2, _Z, _ZR1, _Y2, X1g, X2g, Zg, Zperpg, ZR1g, Y2g, _X2X1, _X1X1, _X2X2, _X1Y2, _X2Y2, X2X1
  real colvector S, parenty1g, _Zperpy1, _y1, y1g
  real scalar g

  this.pRperp = pRperp
  if (isDGP & parent->WREnonARubin) {
    perpRperpX = parent->Repl.perpRperpX
    pZperp = parent->Repl.pZperp
    ZperpZperp = parent->Repl.ZperpZperp  // change to pointer?
    invZperpZperp = parent->Repl.invZperpZperp  // change to pointer?
    pX1perpRperpX = parent->Repl.pX1perpRperpX
    ZperpX1  = parent->Repl.ZperpX1  // change to pointer?
    ZperpX2  = parent->Repl.ZperpX2  // change to pointer?
    ZperpY2  = parent->Repl.ZperpY2  // change to pointer?
    Zperpy1  = parent->Repl.Zperpy1  // change to pointer?
    kX = (kX1 = cols(*pX1perpRperpX)) + parent->kX2
  } else {
    perpRperpX = perp(RperpX)
    pZperp = pXB(*parent->pX1, RperpX)
    invZperpZperp = invsym(ZperpZperp = cross(*pZperp, *pZperp))
    pX1 = pX1perpRperpX = pXB(*parent->pX1, perpRperpX)

    ZperpX1  = cross(*pZperp, *pX1)
    ZperpX2  = cross(*pZperp, *parent->pX2)
    ZperpY2  = cross(*pZperp, *parent->pY2)
    Zperpy1  = cross(*pZperp, *parent->py1)
  }
  kZ = cols(Rpar)

  Z   = *(pX1par = pXB(*parent->pX1, RparX)) + *parent->pY2 * RparY  // Zpar
  pZR1 = pX12B(*parent->pX1, *parent->pY2, R1invR1R1)  // Z*R1
  ZperpZ   = cross(*pZperp, Z    )
  ZperpZR1 = cross(*pZperp, *pZR1)
  

  if (isDGP)
    u1ddot = u1dddot = U2ddot = smatrix()

  if (parent->jk & isDGP & parent->WREnonARubin) {
    XY2g = ZY2g = XXg = XZg = YYg = Zy1g = X1y1g = X2y1g = y1Y2g = invHg = ZZg = invXXg = H_2SLSg = H_2SLSmZZg = ZR1Y2g = ZR1ZR1g = twoy1ZR1g = ZZR1g = X1ZR1g = X2ZR1g = ZXinvXXXZg = smatrix(parent->Nstar)
    y1y1g = J(parent->Nstar, 1, 0)
    beta = smatrix(parent->Nstar + 1)
    u1ddot.M = u1dddot.M = J(parent->Nobs, 1, 0); U2ddot.M = J(parent->Nobs, parent->kY2, 0)

    if (cols(R1invR1R1)) {
      py1pary1parg = &J(parent->Nstar,1,0)
      py1parY2g = &smatrix(parent->Nstar)
      pZy1parg  = &smatrix(parent->Nstar)
      pXy1parg  = &smatrix(parent->Nstar)
    }

    Y2jk = J(parent->Nobs, parent->kY2, 0)
    X2jk = J(parent->Nobs, parent->kX2, 0)
    X1jk = J(parent->Nobs, kX1, 0)
    Zjk = J(parent->Nobs, kZ, 0)
    ZR1jk = J(parent->Nobs, cols(*pZR1), 0)
    y1jk = J(parent->Nobs, 1, 0)

    for (g=parent->Nstar; g; g--) {
      S = parent->NClustVar? parent->infoBootData[g,]' : g\g
      parentX2g = *pXS(*parent->pX2, S)
      parentY2g = *pXS(*parent->pY2, S)
      parenty1g = *pXS(*parent->py1, S)
      Zperpg    = *pXS(*pZperp, S)
      X1g       = *pXS(*pX1perpRperpX, S)
      Zg        = *pXS(Z  , S)
      ZR1g      = *pXS(*pZR1, S)

      _ZperpX1    = ZperpX1    - cross(Zperpg, X1g)
      _ZperpX2    = ZperpX2    - cross(Zperpg, parentX2g)
      _ZperpZ     = ZperpZ     - cross(Zperpg, Zg)
      _ZperpZR1   = ZperpZR1   - cross(Zperpg, ZR1g)
      _ZperpY2    = ZperpY2    - cross(Zperpg, parentY2g)
      _Zperpy1    = Zperpy1    - cross(Zperpg, parenty1g)
      _ZperpZperp = ZperpZperp - cross(Zperpg, Zperpg)
      _invZperpZperp = invsym(_ZperpZperp)

      _X1  = *pX1perpRperpX - *pZperp * (_invZperpZperp * _ZperpX1)  // FWL-process
      _X2  = *parent->pX2   - *pZperp * (_invZperpZperp * _ZperpX2)
      _Z   = Z              - *pZperp * (_invZperpZperp * _ZperpZ)
      _ZR1 = *pZR1          - *pZperp * (_invZperpZperp * _ZperpZR1)
      _Y2  = *parent->pY2   - *pZperp * (_invZperpZperp * _ZperpY2)
      _y1  = *parent->py1   - *pZperp * (_invZperpZperp * _Zperpy1)

      X1g  = *pXS(_X1 , S)
      X2g  = *pXS(_X2 , S)
      Zg   = *pXS(_Z  , S)
      ZR1g = *pXS(_ZR1, S)
      Y2g  = *pXS(_Y2 , S)
      y1g  = *pXS(_y1 , S)

      setXS(X1jk , S, X1g )  // save partialled-out vars from each jk iteration, to compute residuals later
      setXS(X2jk , S, X2g )
      setXS(Y2jk , S, Y2g )
      setXS(y1jk , S, y1g )
      setXS(Zjk  , S, Zg  )
      setXS(ZR1jk, S, ZR1g)

      _X2X1 = cross(_X2, _X1) - cross(X2g, X1g)
      _X1X1 = cross(_X1, _X1) - cross(X1g, X1g)
      _X2X2 = cross(_X2, _X2) - cross(X2g, X2g)
      _X1Y2 = cross(_X1, _Y2) - cross(X1g, Y2g)
      _X2Y2 = cross(_X2, _Y2) - cross(X2g, Y2g)
      y1Y2g[g].M = cross(_y1, _Y2) - cross(y1g, Y2g)
      X2y1g[g].M = cross(_X2, _y1) - cross(X2g, y1g)
      X1y1g[g].M = cross(_X1, _y1) - cross(X1g, y1g)
      y1y1g[g]   = cross(_y1, _y1) - cross(y1g, y1g)
      Zy1g [g].M = cross(_Z , _y1) - cross(Zg , y1g)    
      XZg  [g].M = cross(_X1, _Z)  - cross(X1g, Zg) \ 
                   cross(_X2, _Z)  - cross(X2g, Zg)
      ZZg[g].M   = cross(_Z , _Z)  - cross(Zg , Zg)
      if (isDGP & parent->WREnonARubin)
        ZY2g[g].M = cross(_Z, _Y2)  - cross(Zg, Y2g)

      XY2g[g].M = _X1Y2 \ _X2Y2
      invXXg[g].M = invsym((XXg[g].M = _X1X1, _X2X1' \ _X2X1, _X2X2))
      ZXinvXXXZg[g].M = XZg[g].M ' invXXg[g].M * XZg[g].M

      if (cols(R1invR1R1)) {
        X2ZR1g[g].M     = cross(_X2 , _ZR1) - cross(X2g , ZR1g)
        X1ZR1g[g].M     = cross(_X1 , _ZR1) - cross(X1g , ZR1g)
        ZZR1g[g].M      = cross(_Z  , _ZR1) - cross(Zg  , ZR1g)
        twoy1ZR1g[g].M  = cross(_y1 , _ZR1) - cross(y1g , ZR1g); twoy1ZR1g[g].M = twoy1ZR1g[g].M + twoy1ZR1g[g].M
        ZR1ZR1g[g].M    = cross(_ZR1, _ZR1) - cross(ZR1g, ZR1g)
        ZR1Y2g[g].M     = cross(_ZR1, _Y2 ) - cross(ZR1g, Y2g)
      }
      H_2SLSg[g].M = XZg[g].M ' invXXg[g].M * XZg[g].M
      if (kappa!=1 | LIML) H_2SLSmZZg[g].M = H_2SLSg[g].M - ZZg[g].M
    }

    if (cols(R1invR1R1)==0) {  // DGP not constrained
      py1parY2g    = &y1Y2g 
      pZy1parg     = &Zy1g
      py1pary1parg = &y1y1g 
      pXy1parg = &smatrix(parent->Nstar); for (g=parent->Nstar;g;g--) (*pXy1parg)[g].M = X1y1g[g].M \ X2y1g[g].M
      py1parjk = &y1jk
    }

		if (Fuller)  // number of observations outside each bootstrapping cluster
	    _Nobsg = parent->_Nobs :- (parent->obswttype=="fweight"? panelsum(*parent->pwt, parent->infoBootData) : (parent->infoBootData[,2] - parent->infoBootData[,1]) :+ 1)
  }

  if (isDGP & parent->WREnonARubin) {
    pX1 = parent->Repl.pX1
    X2 = parent->Repl.X2  // change to pointer
    X1X = parent->Repl.X1X
    X2X = parent->Repl.X2X
    XX = parent->Repl.XX
    invXX = parent->Repl.invXX
    Y2 = parent->Repl.Y2
    y1 = parent->Repl.y1
    XY2  = parent->Repl.XY2
    Y2Y2  = parent->Repl.Y2Y2
    y1Y2 = parent->Repl.y1Y2   
    X2y1 = parent->Repl.X2y1   
    X1y1 = parent->Repl.X1y1   
    y1y1 = parent->Repl.y1y1   
  } else {
    pX1 = &(*pX1perpRperpX - *pZperp * (invZperpZperp * ZperpX1))      // FWL-process X1
    X2 = *parent->pX2 - *pZperp * (invZperpZperp * ZperpX2)
    X2X1 = cross(X2, *pX1)
    invXX = invsym(XX = (X1X = cross(*pX1, *pX1), X2X1') \ (X2X = X2X1, cross(X2, X2)))
    Y2 = *parent->pY2 - *pZperp * (invZperpZperp * ZperpY2)
    if (parent->scoreBS==0) Y2Y2 = cross(Y2, Y2)
    y1 = *parent->py1 - *pZperp * (invZperpZperp * Zperpy1)
    XY2  = cross(*pX1, Y2) \ cross(X2, Y2)
    y1Y2 = cross(y1  , Y2)
    X2y1 = cross(X2  , y1)
    X1y1 = cross(*pX1, y1)
    y1y1 = cross(y1  , y1)
  }
  Z    =      Z   - *pZperp * (invZperpZperp * ZperpZ  )
  pZR1 =  &(*pZR1 - *pZperp * (invZperpZperp * ZperpZR1))

  Zy1  = cross(Z   ,  y1)    
  XZ   = cross(*pX1, Z  ) \ 
         cross(  X2, Z  )
  ZZ  =  cross( Z  , Z  )
  if (parent->WREnonARubin) ZY2 = cross(Z, Y2)
  
  ZXinvXXXZ = XZ ' (invXXXZ = invXX * XZ)

  if (cols(R1invR1R1)) {
    X2ZR1    = cross(X2   , *pZR1)
    X1ZR1    = cross(*pX1 , *pZR1)
    ZZR1     = cross(Z    , *pZR1)
    twoy1ZR1 = cross( y1  , *pZR1); twoy1ZR1 = twoy1ZR1 + twoy1ZR1
    ZR1ZR1   = cross(*pZR1, *pZR1)
    ZR1Y2    = cross(*pZR1, Y2   )
  } else {
    py1parY2   = &y1Y2
    pZy1par    = &Zy1
    y1pary1par = y1y1 
    pXy1par    = &(X1y1 \ X2y1)
    py1par     = &y1
    t1  = J(parent->kZ,1,0)
    t1Y = J(parent->kY2,1,0)
  }

  V = invXX * XZ  // in 2SLS case, estimator is (V' XZ)^-1 * (V'Xy1). Also used in kZ-class and LIML robust VCV by Stata convention

  if (isDGP) {
    if (parent->jk & LIML==0 & parent->WREnonARubin)
      for (g=parent->Nstar; g; g--)
        invHg[g].M = invsym(kappa==1? H_2SLSg[g].M : ZZg[g].M + kappa * H_2SLSmZZg[g].M)
    H_2SLS = V ' XZ  // Hessian
    if (kappa!=1 | LIML) H_2SLSmZZ = H_2SLS - ZZ
    if (LIML==0)  // DGP is LIML except possibly when getting confidence peak for A-R plot; but LIML=0 when exactly id'd, for then kappa=1 always and Hessian doesn't depend on r1 and can be computed now
      MakeH()
  } else {
    pX1par = &(*pX1par - *pZperp * (invZperpZperp * cross(*pZperp, *pX1par)))
    Yendog = 1, colsum(RparY :!= 0)  // columns of Y = [y1par Zpar] that are endogenous (normally all)
    if (parent->robust & parent->bootstrapt) {  // for WRE replication regression, prepare for CRVE
      XinvXX = *pX12B(*pX1, X2, invXX)
      if (((parent->granular & parent->Nw==1) | parent->jk) & parent->WREnonARubin)
        PXZ = *pX12B(*pX1, X2, invXXXZ)
    }
  }
}


// do most of estimation; for LIML r1 must be passed now in order to solve eigenvalue problem involving it
// inconsistency: for replication regression of Anderson-Rubin, r1 refers to the *null*, not the maintained constraints, because that's what affects the endogenous variables
// For WRE, should only be called once for the replication regressions, since for them r1 is the unchanging model constraint
void boottestOLS::Estimate(real scalar _jk, real colvector r1) {
  beta.M = beta0 - dbetadr * r1
  if (_jk & rows(R1perp))
    t1 = R1invR1R1 * r1
}

void boottestARubin::Estimate(real scalar _jk, real colvector r1) {
  beta.M = beta0 - dbetadr * r1
  py1par = &(*parent->py1 - *parent->pY2 * r1)
  if (_jk & rows(R1perp))
    t1 = R1invR1R1 * r1
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

void boottestIVGMM::Estimate(real scalar _jk, real colvector r1) {
  real rowvector val; real matrix vec; real scalar g, kappag; real colvector ZXinvXXXy1par, invXXXy1parg; real matrix YPXY
  pragma unset vec; pragma unset val

  if (cols(R1invR1R1)) {
    t1 = R1invR1R1 * r1

    if (isDGP)
      py1par = &(y1 - *pZR1 * r1)

    y1pary1par = y1y1 - twoy1ZR1 * r1 + r1 ' ZR1ZR1 * r1
    py1parY2 = &(y1Y2  - r1 ' ZR1Y2)
    pZy1par  = &( Zy1 -  ZZR1 * r1)
    pXy1par  = &(X1y1 - X1ZR1 * r1 \ X2y1 - X2ZR1 * r1)

    if (_jk) {
      t1Y = t1[|parent->kX1+1,.\.,.|]
      for (g=parent->Nstar; g; g--) {
        py1parjk = &(y1jk - ZR1jk * r1)
        (*py1pary1parg)[g]   = y1y1g[g]   - twoy1ZR1g[g].M * r1 + r1 ' ZR1ZR1g[g].M * r1
        (*py1parY2g)   [g].M = y1Y2g[g].M - r1 ' ZR1Y2g[g].M
        (*pZy1parg)    [g].M =  Zy1g[g].M -  ZZR1g[g].M * r1
        (*pXy1parg)    [g].M = X1y1g[g].M - X1ZR1g[g].M * r1 \ X2y1g[g].M - X2ZR1g[g].M * r1
      }
    }
  }
    
  invXXXy1par = invXX * *pXy1par
  YY = y1pary1par, *pZy1par' \ *pZy1par, ZZ

  if (isDGP) {
    ZXinvXXXy1par = XZ ' invXXXy1par

    if (LIML) {
      YPXY = invXXXy1par ' (*pXy1par) , ZXinvXXXy1par' \ ZXinvXXXy1par , ZXinvXXXZ
      eigensystemselecti(invsym(YY) * YPXY, rows(YY)\rows(YY), vec, val)
      kappa = 1/(1 - Re(val))  // sometimes a tiny imaginary component sneaks into val
      if (Fuller) kappa = kappa - Fuller / (parent->_Nobs - parent->kX)
      MakeH()
    }
    beta.M = invH * (kappa==1? ZXinvXXXy1par : kappa * (ZXinvXXXy1par - *pZy1par) + *pZy1par)

    if (_jk)
      for (g=parent->Nstar; g; g--) {
        ZXinvXXXy1par = XZg[g].M ' (invXXXy1parg = invXXg[g].M * (*pXy1parg)[g].M)
        YYg[g].M = (*py1pary1parg)[g], (*pZy1parg)[g].M' \ (*pZy1parg)[g].M, ZZg[g].M

        if (LIML) {
          YPXY = invXXXy1parg ' (*pXy1parg)[g].M , ZXinvXXXy1par' \ ZXinvXXXy1par , ZXinvXXXZg[g].M
          eigensystemselecti(invsym(YYg[g].M) * YPXY, rows(YYg[g].M)\rows(YYg[g].M), vec, val)
          kappag = 1/(1 - Re(val))
          if (Fuller) kappag = kappag - Fuller / (_Nobsg[g] - parent->kX)
          invHg[g].M = invsym(ZZg[g].M + kappag * H_2SLSmZZg[g].M)
        } else
          kappag = kappa

        beta[g+1].M = invHg[g].M * (kappag==1? ZXinvXXXy1par : kappag * (ZXinvXXXy1par - (*pZy1parg)[g].M) + (*pZy1parg)[g].M)
      }
  } else if (parent->WREnonARubin)  // if not score bootstrap of IV/GMM...
    Rt1 = RR1invR1R1 * r1
}

void boottestOLS::MakeResiduals(real scalar _jk) {
  real scalar g, m; real colvector S, u1ddotg, Xt1

  u1ddot.M = *py1par - *pX12B(*parent->pX1, *parent->pX2, beta.M)

  if (_jk) {
    m = parent->small? sqrt((parent->Nstar - 1) / parent->Nstar) : 1
    
    if (parent->purerobust)  // somwhat faster handling classic hc3
      if (rows(R1perp)) {
        Xt1 = *pX12B(*parent->pX1, *parent->pX2, t1)
        u1ddot[2].M = m * ((invMg.M :* (u1ddot.M + Xt1)) - Xt1)
      } else
        u1ddot[2].M = m * invMg.M :* u1ddot.M
    else if (parent->granularjk) {
      for (g=parent->Nstar; g; g--) {
        S = parent->NClustVar? parent->infoBootData[g,]' : g\g
        if (rows(R1perp)) {
          Xt1 = Xg[g].M * t1
          u1ddot[2].M[|S|] = m * (cross(invMg[g].M, u1ddot.M[|S|] + Xt1) - Xt1)
        } else
          u1ddot[2].M[|S|] = m * cross(invMg[g].M, u1ddot.M[|S|])
      }
    } else
      for (g=parent->Nstar; g; g--) {
        S = parent->infoBootData[g,]'
        u1ddotg = u1ddot.M[|S|]
        u1ddot[2].M[|S|] = m * (u1ddotg + XinvHg[g].M * (rows(R1perp)? cross(Xg[g].M, u1ddotg) + XXg[g].M * t1 : cross(Xg[g].M, u1ddotg)))
      }
  }
}

void boottestIVGMM::MakeResiduals(real scalar _jk) {
  pragma unused _jk
  real matrix tmp; real colvector Xu, negXuinvuu, delta, _beta, S, _t1Y; real scalar uu; real scalar g
  
  if (parent->jk)
    for (g=parent->Nstar; g; g--) {
      S = parent->NClustVar? parent->infoBootData[g,1] \ parent->infoBootData[g,2] : g\g
      u1ddot.M[|S|] = *pXS(*py1parjk,S) - *pXS(Zjk,S) * beta[1+g].M
    }
  else
    u1ddot.M = *py1par - Z * beta.M

  if (parent->scoreBS==0) {
    _beta = 1 \ -beta.M
    uu = _beta ' YY * _beta
    Xu = *pXy1par - XZ * beta.M  // after DGP regression, compute Y2 residuals by regressing Y2 on X while controlling for y1 residuals, done through FWL
    negXuinvuu = Xu / -uu
    Piddot = invsym(XX + negXuinvuu * Xu') * (negXuinvuu * (*py1parY2 - beta.M ' ZY2) + XY2)
    Y2bar = *pX12B(*pX1, X2, Piddot)
    U2ddotU2ddot = Y2Y2 - Piddot'XY2 - XY2'Piddot + Piddot'XX*Piddot

    delta = Rpar * beta.M + t1 - parent->Repl.t1
    deltaX = *pXS(delta, .\parent->kX1)
    deltaY = *pXS(delta, parent->kX1+1\.)
    y1bar = *pX1 * (tmp = perpRperpX'deltaX) + Y2bar * deltaY
    deltadddot = Piddot * deltaY; if (kX1) deltadddot[|.\kX1|] = deltadddot[|.\kX1|] + tmp  // (X_∥'X_∥)^-1 * X_∥'y1bar
    Xy1bar = XX * deltadddot
    y1bary1bar = deltadddot'Xy1bar
    XU2ddot = XY2 - XX * Piddot
    y1barU2ddot = deltadddot'XU2ddot

    if (parent->jk) {
      _t1Y = t1Y - parent->Repl.t1Y
      for (g=parent->Nstar; g; g--) {
        S = parent->NClustVar? parent->infoBootData[g,1] \ parent->infoBootData[g,2] : g\g
        _beta = 1 \ -beta[1+g].M
        uu = _beta ' YYg[g].M * _beta
        Xu = (*pXy1parg)[g].M - XZg[g].M * beta[1+g].M
        negXuinvuu = Xu / -uu
        Piddot = invsym(XXg[g].M + negXuinvuu * Xu') * (negXuinvuu * ((*py1parY2g)[g].M - beta[1+g].M ' ZY2g[g].M) + XY2g[g].M)
        U2ddot.M [|S,(.\.)|] = tmp = *pXS(Y2jk,S) - *pX12B(*pXS(X1jk,S), *pXS(X2jk,S), Piddot)
        u1dddot.M[|S      |] = u1ddot.M[|S|] + tmp * (RparY * beta[1+g].M + _t1Y)
      }
    } else {
      U2ddot.M = Y2 - Y2bar
      u1dddot.M = u1ddot.M + U2ddot.M * deltaY
    }
  }
}


// non-WRE stuff that only depends on r in A-R case, for test stat denominators in replication regressions
// since the non-AR OLS code never creates an object for replication regresssions, in that case this is called on the DGP regression object
// depends on results of Estimate() only when doing OLS-style bootstrap on an overidentified IV/GMM regression--score bootstrap or A-R. Then kappa from DGP LIML affects Hessian, pH.
void boottestOLS::InitTestDenoms() {
  real scalar d

  if (parent->bootstrapt & (parent->scoreBS | parent->robust)) {
    if (parent->granular | parent->purerobust) {
      WXAR = smatrix(parent->df)
      for (d=parent->df;d;d--)
        WXAR[d].M = XAR[,d]
    }

    if (parent->NFE & parent->robust & (parent->FEboot | parent->scoreBS)==0 & parent->granular < parent->NErrClustCombs) {  // make first factor of second term of (64) for c=∩ (c=1)
      CT_XAR = smatrix(parent->df)
      for (d=parent->df;d;d--)
        CT_XAR[d].M = parent->crosstabFE(XAR[,d], *parent->pinfoCapData)
    }
  }
}


// partial fixed effects out of a data matrix
pointer(real matrix) scalar boottest::partialFE(pointer(real matrix) scalar pIn) {
  real matrix Out, tmp; real scalar i
  if (NFE & pIn) {
    Out = *pIn
    if (haswt)
      for (i=NFE;i;i--) {
        tmp = Out[FEs[i].is,]
        Out[FEs[i].is,] = tmp :- FEs[i].sqrtwt * cross(FEs[i].wt, tmp)
      }
    else
      for (i=NFE;i;i--) {
        tmp = Out[FEs[i].is,]
        Out[FEs[i].is,] = tmp :- cross(FEs[i].wt, tmp)
      }
    return(&Out)
  }
  return (pIn)
}

void boottest::new() {
  ARubin = LIML = Fuller = WRE = small = scoreBS = auxwttype = ML = initialized = quietly = sqrt = ptype = robust = NFE = FEboot = granular = NErrClustCombs = subcluster = B = BFeas = interpolating = jk = 0
  twotailed = null = dirty = willplot = v_sd = notplotted = FEdfadj = bootstrapt = 1
  level = 95
  ptol = 1e-6
  confpeak = MaxMatSize = .
  pY2 = pX1 = pX2 = py1 = pSc = pID = pFEID = pR1 = pR = pwt = &J(0,0,0)
  pr1 = pr = &J(0,1,0)
  pIDBootData = pIDBootAll = &.
}

// break loops in data structure topology to enable garbage collection
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
void boottest::setjk (real scalar jk) {
  this.jk = jk; setdirty(1)
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
  real colvector sortID, o, _FEID, _sqrtwt
  pointer (real colvector) scalar pIDAllData, pIDCapData
  real rowvector ClustCols
  real matrix Combs, tmp, IDCap
  real scalar i, j, c, minN, _B, i_FE, sumFEwt, _df
  pragma unset pIDAllData; pragma unset pIDCapData

  Nobs = rows(*pX1)
  NClustVar = cols(*pID)
  kX = (kX1 = cols(*pX1)) + (kX2 = cols(*pX2))
  if (kX2 == 0) pX2 = &J(Nobs,0,0)
  if ((kY2 = cols(*pY2)) == 0) pY2 = &J(Nobs,0,0)
  kZ = kX1 + kY2

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
  
  if (haswt = rows(*pwt)>1) {
    sumwt = sum(*pwt)
    sqrtwt = sqrt(*pwt)
    py1 = &(*py1 :* sqrtwt)
    if (kY2) pY2 = &(*pY2 :* sqrtwt)
    if (kX1) pX1 = &(*pX1 :* sqrtwt)
    if (kX2) pX2 = &(*pX2 :* sqrtwt)
  } else
    pwt = &(sumwt = 1)
  _Nobs = haswt & obswttype=="fweight"? sumwt : Nobs

  bootneqcap = NBootClustVar!=NClustVar | subcluster  // bootstrapping and intersection of error clusterings different?

  if (WREnonARubin)
    if (NClustVar)
      infoBootData = _panelsetup(*pID, 1..NBootClustVar, pIDBootData)
    else {
      pIDBootData = &(1::Nobs)
      pinfoCapData = &(infoBootData = *pIDBootData,*pIDBootData)  // no clustering, so no collapsing by cluster
    }
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
        pinfoAllData = &_panelsetup(*pID, 1..NClustVar, pIDAllData)
        if (subcluster & rows(infoBootData) != rows(*pinfoAllData)) {
          errprintf("\nboottest can only perform the subcluster bootstrap when the bootstrap clusters are nested within the (intersections of the) error clusters.\n")
          _error(499)
        }
      } else {
        pinfoAllData = &infoBootData  // info for grouping by intersections of all bootstrap & clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster & intersection of all error clusters
        if (WREnonARubin)
          pIDAllData = pIDBootData
      }

      if (NClustVar > NErrClustVar)  // info for intersections of error clustering wrt data
        pinfoCapData = &_panelsetup(*pID, subcluster+1..NClustVar, pIDCapData)
      else {
        pinfoCapData = pinfoAllData  // info for intersections of error clustering wrt data
        if (WREnonARubin)
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

    Ncap = rows(*pinfoCapData)

    purerobust = robust & (scoreBS | subcluster)==0 & Nstar==Nobs  // do we ever error- *and* bootstrap-cluster by individual?
    granular   = WREnonARubin? 2*Nobs*B*(2*Nstar+1) < Nstar*(Nstar*Nobs+Clust.N*B*(Nstar+1)) :
                               !jk & robust & scoreBS==0 & (purerobust | (Clust.N+Nstar)*kZ*B + (Clust.N-Nstar)*B + kZ*B < Clust.N*kZ*kZ + Nobs*kZ + Clust.N * Nstar*kZ + Clust.N*Nstar)
    if (jk & !WREnonARubin)
      granularjk = kZ^3 + Nstar * (Nobs/Nstar*kZ^2 + (Nobs/Nstar)^2*kZ + (Nobs/Nstar)^2 + (Nobs/Nstar)^3) < Nstar * (kZ^2*Nobs/Nstar + kZ^3 + 2*kZ*(kZ + Nobs/Nstar))

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
          _sqrtwt  = sqrtwt[FEs[i_FE].is]
          FEs[i_FE].wt = _sqrtwt / (sumFEwt = colsum((*pwt)[FEs[i_FE].is]))
          FEs[i_FE].sqrtwt = _sqrtwt
        } else
          FEs[i_FE].wt = J(j-i, 1, 1/(sumFEwt = j-i))
        if ((B & robust & granular < NErrClustVar) | (WREnonARubin & robust & granular & bootstrapt))
          invFEwt[i_FE] = 1 / sumFEwt

        j = i
        
        if (!WREnonARubin & FEboot) {  // are all of this FE's obs in same bootstrapping cluster? (But no need to check if B=0 for then CT_WE in 2nd term of (62) orthogonal to v = col of 1's)
          tmp = (*pID)[FEs[i_FE].is, 1..NBootClustVar]
          FEboot = all(tmp :== tmp[1,])
        }
        ++i_FE
      }
      _FEID[o[i]] = i_FE
    }
    FEs[NFE].is = o[|.\j|]
    if (haswt) {
      _sqrtwt  = sqrtwt[FEs[NFE].is]
      FEs[NFE].wt = _sqrtwt / (sumFEwt = colsum((*pwt)[FEs[NFE].is]))
      FEs[NFE].sqrtwt = _sqrtwt
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
      DGP.SetR(null? *pR1 \ *pR : *pR1)  // DGP constraints: model constraints + null if imposed
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
        DGP.Estimate(0, *pr1)
        confpeak = DGP.beta.M  // estimated coordinate of confidence peak
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
      Repl.Estimate(0, *pr1)

      DGP = boottestIVGMM(); DGP.parent = &this
      DGP.LIML = 1
      if (null)
        DGP.SetR(*pR1 \ *pR, J(0,kZ,0))  // DGP constraints: model constraints + imposed null
      else
        DGP.SetR(*pR1 , *pR)  // when null not imposed, keep it in the attack surface, though not used there, so Zperp is same in DGP and Repl
      
      DGP.InitVars()
      if (null==0) {  // if not imposing null, then DGP constraints, kappa, Hessian, etc. do not vary with r and can be set now
        DGP.Estimate(jk, *pr1)
        DGP.MakeResiduals(0)
      }

      if (LIML & Repl.kZ==1 & Nw==1) As = betas = J(1, B+1, 0)

      SstarUZperpinvZperpZperp = SstarUZperp = SstarUXinvXX = SstarUX = smatrix(Repl.kZ+1)
      SstarUU = smatrix(Repl.kZ+1, Repl.kZ+1)

      if (bootstrapt) {
        deltadenom_b = J(Repl.kZ, Repl.kZ, 0)
        Zyi = deltadenom = smatrix(Repl.kZ+1)
        SstarUMZperp = SstarUPX = SstarUX
        _Jcap = J(Clust.N, Repl.kZ, 0)
        if (granular==0)
          SCTcapUXinvXX = smatrix(Repl.kZ+1, Nstar)
        if (LIML | robust==0) {
          YYstar = YPXYstar = SstarUX
          YYstar_b = YPXYstar_b = J(Repl.kZ+1, Repl.kZ+1, 0)
        }
        if (NFE & (bootstrapt | kappa != 1 | LIML))
          CTFEU = SstarUX
      }
      InitWRE()

    } else {  // the score bootstrap for IV/GMM uses a IV/GMM DGP but then masquerades as an OLS test because most factors are fixed during the bootstrap. To conform, need DGP and Repl objects with different R, R1, one with FWL, one not

      Repl = boottestIVGMM(); Repl.parent = &this
      Repl.LIML = this.LIML; Repl.Fuller = this.Fuller; Repl.kappa = this.kappa
      Repl.SetR(*pR1, I(kZ))  // process replication restraints = model constraints only
      Repl.InitVars(&Repl.R1perp)
      Repl.Estimate(0, *pr1)  // bit inefficient to estimate in both objects, but maintains the conformity
      Repl.InitTestDenoms()
      DGP = boottestIVGMM(); DGP.parent = &this
      DGP.LIML = this.LIML; DGP.Fuller = this.Fuller; DGP.kappa = this.kappa
      DGP.SetR(null? *pR1 \ *pR : *pR1, J(0,kZ,0))  // DGP constraints: model constraints + null if imposed
      DGP.InitVars()
      pM = &Repl  // estimator object from which to get A, AR, XAR; DGP follows WRE convention of using FWL, Repl follows OLS convention of not; scoreBS for IV/GMM mixes the two
      if (null==0) {  // if not imposing null, then DGP constraints, kappa, Hessian, etc. do not vary with r and can be set now
        DGP.Estimate(jk, *pr1)
        DGP.MakeResiduals(0)
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


  if (df==1) setsqrt(1)  // work with t/z stats instead of F/chi2

  if (small) {
    _df = _Nobs - kZ + cols(Repl.R1invR1R1) - FEdfadj
    df_r = NClustVar? minN - 1 : _df
    multiplier = (smallsample = _df / (_Nobs - robust)) / df  // divide by # of constraints because F stat is so defined
  } else
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
    interpolate_u = (robust | ML)==0
    if (interpolable = bootstrapt & B & null & Nw==1 & (kappa==0 | ARubin)) {    
      dnumerdr = smatrix(q)
      if (interpolate_u) dudr = dnumerdr
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
  }
  if (bootstrapt==0)
    UpdateBootstrapcDenom()

  BFeas = (*pDist)[1]==.? 0 : rownonmissing(*pDist) - 1
  DistCDR = J(0,0,0)
  setdirty(0)
  initialized = 1
}

// if not imposing null and we have returned to boottest(), then df=1 or 2; we're plotting or finding CI, and only test stat, not distribution, changes with r
void boottest::NoNullUpdate() {
  if (WREnonARubin)
    (*pnumer)[,1] = Repl.RRpar * betas[,1] - *pr
  else if (ARubin) {
    DGP.Estimate(0, *pr)
    (*pnumer)[,1] = v_sd * DGP.Rpar * DGP.beta.M[|kX1+1\.|]  // coefficients on excluded instruments in ARubin OLS
  } else
    (*pnumer)[,1] = v_sd * (*pR * (ML? beta : pM->Rpar * pM->beta.M) - *pr)  // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this

  (*pDist)[1] = df==1? (*pnumer)[1] / sqrt(statDenom) : (*pnumer)[,1] ' invsym(statDenom) * (*pnumer)[,1]
}

// compute bootstrap-c denominator from all bootstrap numerators
// for 1-dimensional hypothesis, can ignore denominators (Young 2022, eq. 16), except that user might request actual distribution
void boottest::UpdateBootstrapcDenom() {
  real colvector numer1
  numer1 = (*pnumer)[,1]
  numersum = rowsum(*pnumer) - numer1
  statDenom = (*pnumer * *pnumer' - numer1 * numer1' - numersum * numersum' / B) / B
  pDist = sqrt? &(*pnumer:/sqrt(statDenom)) : &colsum(*pnumer :* invsym(statDenom) * *pnumer)
}

// draw wild weight matrix of width _B. If first=1, insert column of 1s at front. For non-Anderson-Rubin WRE, subtract 1 from all weights
void boottest::MakeWildWeights(real scalar _B, real scalar first) {
  real scalar m
  m = WREnonARubin & jk & small? sqrt(1 - 1 / Nstar) : 1  // finite-sample adjuster for JK regressions

  if (_B) {  // in scoretest or waldtest WRE, still make v a col of 1's
    if (enumerate)
      v = J(Nstar,1,1), count_binary(Nstar, -m, m)  // complete Rademacher set
    else if (auxwttype==3)
      v = rnormal(Nstar, _B+first, 0, m)  // normal weights
    else if (auxwttype==4)
      v = m * (rgamma(Nstar, _B+first, 4, .5) :- 2)  // Gamma weights
    else if (auxwttype==2)
      if (WREnonARubin)
        v = sqrt((2*m) * ceil(runiform(Nstar, _B+first) * 3)) :* ((runiform(Nstar, _B+first):>=.5):-.5)  // Webb weights
      else {
        v = sqrt(        ceil(runiform(Nstar, _B+first) * 3)) :* ((runiform(Nstar, _B+first):>=.5):-.5)  // Webb weights, divided by sqrt(2)
        v_sd = 1.6a09e667f3bcdX-001 /*sqrt(.5)*/
      }
    else if (auxwttype)
      if (WREnonARubin)
        v = ( rdiscrete(Nstar, _B+first, 1.727c9716ffb76X-001 \ 1.1b06d1d200914X-002 /*.5+sqrt(.05) \ .5-sqrt(.05)*/) :- 1.5 ) * (m * 1.1e3779b97f4a8X+001) /*sqrt(5)*/ :+ .5  // Mammen weights, minus 1 for convenience in WRE
      else {
        v = ( rdiscrete(Nstar, _B+first, 1.727c9716ffb76X-001 \ 1.1b06d1d200914X-002 /*.5+sqrt(.05) \ .5-sqrt(.05)*/) :- 1.5 ) :+ 1.c9f25c5bfedd9X-003 /*.5/sqrt(5)*/  // Mammen weights, divided by sqrt(5)
        v_sd = 1.c9f25c5bfedd9X-002 /*sqrt(.2)*/
      }
    else if (WREnonARubin)
      v = jk? m :+ (runiform(Nstar, _B+first) :< .5) * (-2 * m) : (runiform(Nstar, _B+first) :<  .5) * -2 :+ 1  // Rademacher weights, minus 1 for WRE; 1st exp simplifies to 2nd when jk=0 but is slower
    else {
      v = (runiform(Nstar, _B+first) :>= .5) :- .5   // Rademacher weights, divided by 2
      v_sd = .5
    }

    if (first)
      v[,1] = J(Nstar, 1, WREnonARubin? m : v_sd)  // keep original residuals in first entry to compute base model stat    
  } else
    v = J(0,1,0)  // in places, cols(v) indicates B -- 1 for classical tests
}


// For WRE, and with reference to Y = [y1 Z], given 0-based column indexes within it, ind1, ind2, return all bootstrap realizations of Y[,ind1]'((1-kappa)*M_Zperp-kappa*P_Xpar)*Y[,ind2] for kappa constant across replications
// ind1 can be a rowvector
// (only really the Hessian when we narrow Y to Z)
real matrix boottest::HessianFixedkappa(real rowvector ind1, real scalar ind2, real scalar kappa, real scalar _jk) {
  real matrix retval; real scalar i
  if (cols(ind1) > 1) {
    retval = J(cols(ind1),cols(v),0)
    for (i=cols(ind1);i;i--)
      retval[i,] = _HessianFixedkappa(ind1[i], ind2, kappa, _jk)
    return(retval)
  }
  return(_HessianFixedkappa(ind1, ind2, kappa, _jk))
}
real rowvector boottest::_HessianFixedkappa(real scalar ind1, real scalar ind2, real scalar kappa, real scalar _jk) {
  real matrix retval, _retval; pointer (real colvector) scalar pT1L, pT1R

  if (kappa) {
    pT1L = pcol(XYbar,ind1+1)  // X_∥^' Y_(∥i)
    if (Repl.Yendog[ind1+1])
      pT1L = &(*pT1L :+ SstarUX[ind1+1].M * v)

    pT1R = ind2? pcol(invXXXZbar,ind2) : &DGP.deltadddot
    if (Repl.Yendog[ind2+1])
      pT1R = &(*pT1R :+ SstarUXinvXX[ind2+1].M * v)  // right-side linear term
    retval = colsum(*pT1L :* *pT1R)  // multiply in the left-side linear term
  }

  if (kappa != 1) {
    _retval = YbarYbar[ind1+1,ind2+1] :+ 
           (*pcol(SstarUYbar[ind2+1].M, ind1+1) + *pcol(SstarUYbar[ind1+1].M, ind2+1)) ' v + 
           colsum(v :* (ind1 <= ind2? SstarUU[ind2+1, ind1+1].M : SstarUU[ind1+1, ind2+1].M) :* v)

    if (Repl.Yendog[ind1+1] & Repl.Yendog[ind2+1])
      _retval = _retval - colsum((SstarUZperpinvZperpZperp[ind1+1].M * v) :* (SstarUZperp[ind2+1].M * v))

    if (NFE)
      _retval = _retval - colsum(CTFEU[ind1+1].M * v :* (invFEwt :* CTFEU[ind2+1].M) * v)

    retval = kappa? kappa*retval:+(1-kappa)*_retval : _retval
  }

  if (_jk) {
    if (kappa)
      retval[,1] = cross(ind1? *pcol(Repl.XZ     ,ind1) : *Repl.pXy1par, 
                         ind2? *pcol(Repl.invXXXZ,ind2) :  Repl.invXXXy1par)
    if (kappa != 1)
      retval[,1] = kappa? kappa :* retval[,1] :+ (1 - kappa) * Repl.YY[ind1+1,ind2+1] :
                                                               Repl.YY[ind1+1,ind2+1]
  }
  return(cols(retval)>1? retval : J(1,cols(v),retval))  // if both vars exogenous, term is same for all b; this duplication is a bit inefficient, but only arises when exog vars involved in null
}
 

// Workhorse for WRE CRVE sandwich filling
// Given column index ind1 and a matrix betas of all the bootstrap estimates, return all bootstrap realizations of P_X * Z[,ind1]_g ' u\hat_1g^*
// for all groups in the intersection of all error clusterings
// return value has one row per cap cluster, one col per bootstrap replication
pointer(real matrix) scalar boottest::Filling(real scalar ind1, real matrix betas, real scalar _jk) {
  real scalar i, ind2; real matrix retval, CTFEUv, T1, SstarUXv, SstarUZperpinvZperpZperp_v; pointer (real matrix) scalar pbetav; pointer (real colvector) pPXYstar; real rowvector _beta, SstarUMZperp_ind2_i; real colvector S
  pragma unset retval

  if (granular) {  // create pieces of each N x B matrix one at a time rather than whole thing at once
    if (Nw == 1) {  // create or avoid NxB matrix?
      pPXYstar = pcol(PXZbar, ind1)
      if (Repl.Yendog[ind1+1])
        pPXYstar = &(*pPXYstar :+ SstarUPX[ind1+1].M * v)

      retval = *_panelsum(*pPXYstar :* (DGP.y1bar :- SstarUMZperp.M * v), *pinfoCapData)

      for (ind2=Repl.kZ;ind2;ind2--) {
        _beta = betas[ind2,]
        retval = retval - *_panelsum(*pPXYstar :* (Repl.Yendog[ind2+1]? *pcol(Zbar,ind2) * _beta :- SstarUMZperp[ind2+1].M * (v :* _beta) :
                                                                        *pcol(Zbar,ind2) * _beta                                           ), *pinfoCapData)
      }
    } else { // create pieces of each N x B matrix one at a time rather than whole thing at once--slower, but less memory-intensive
      retval = J(Clust.N, cols(v), 0)
      SstarUXv = SstarUX[ind1+1].M * v
      for (ind2=0; ind2<=Repl.kZ; ind2++) {
        if (ind2) {
          _beta = betas[ind2,]
          pbetav = &(v :* _beta)
        } else
          pbetav = &v

        if (Repl.Yendog[ind2+1])
          SstarUZperpinvZperpZperp_v = SstarUZperpinvZperpZperp[ind2+1].M * *pbetav

        if (NFE)
          CTFEUv = invFEwt :* (CTFEU[ind2+1].M * *pbetav)

        for (i=Clust.N;i;i--) {
          S = (*pinfoCapData)[i,]'
          pPXYstar = &PXZbar[|S,(ind1\ind1)|]
          if (Repl.Yendog[ind1+1])
            pPXYstar = &(*pPXYstar :+ *pXS(Repl.XinvXX,S) * SstarUXv)
          if (Repl.Yendog[ind2+1]) {
            SstarUMZperp_ind2_i = *pXS(*Repl.pZperp,S) * SstarUZperpinvZperpZperp_v - (ind2? (*pU2parddot)[|S, (ind2\ind2)|] :* (*pbetav)[*pXS(*pIDBootData,S),] :
                                                                                             DGP.u1dddot.M[|S             |] :*        v [*pXS(*pIDBootData,S),])

            if (NFE)
              SstarUMZperp_ind2_i = SstarUMZperp_ind2_i + CTFEUv[(*pFEID)[|S|],]  // CT_(*,FE) (U ̈_(∥j) ) (S_FE S_FE^' )^(-1) S_FE

            if (ind2)
              retval[i,] = retval[i,] - colsum(*pPXYstar :* (Zbar[|S,(ind2\ind2)|] * _beta :- SstarUMZperp_ind2_i))
            else
              retval[i,] =              colsum(*pPXYstar :* (DGP.y1bar[|S|]                :- SstarUMZperp_ind2_i))
          } else
              retval[i,] = retval[i,] - colsum(*pPXYstar :* (Zbar[|S,(ind2\ind2)|] * _beta))
        }
      }
    }
  } else {  // coarse error clustering
    for (ind2=0; ind2<=Repl.kZ; ind2++) {
      if (ind2) {
        _beta = -betas[ind2,]
        pbetav = &(v :* _beta)
      } else
        pbetav = &v

      // T1 * v will be 1st-order terms
      T1 = Repl.Yendog[ind1+1]? ScapYbarX[ind2+1].M * SstarUXinvXX[ind1+1].M : J(0,0,0) //  S_∩ (Ybar_(∥j)) (X_∥^' X_∥ )^(-1) [S_* (U ̈_(∥i):*X_∥ )]^'

      if (Repl.Yendog[ind2+1]) {  // add CT_(∩,*) (P_(X_∥ ) Y_(∥i):*U ̈_(∥j) )
        if (NClustVar == NBootClustVar & !subcluster)  // simple case of one clustering: full crosstab is diagonal
          if (cols(T1))
            _diag(T1, diagonal(T1) + SstarUXinvXX[ind2+1].M ' (*pcol(XYbar,ind1+1)))
          else
            T1 =                diag(SstarUXinvXX[ind2+1].M ' (*pcol(XYbar,ind1+1)))

        else {
          if (Repl.Yendog[ind1+1]==0)
            T1 = JNcapNstar
          for (i=Nstar;i;i--)
            T1[IDCTCapstar[i].M, i] = T1[IDCTCapstar[i].M, i] + SCTcapUXinvXX[ind2+1,i].M * *pcol(XYbar,ind1+1)
        }

        if (cols(*Repl.pZperp))  // subtract S_∩ (P_(X_∥ ) Y_(∥i):*Z_⊥ ) (Z_⊥^' Z_⊥ )^(-1) [S_* (U ̈_(∥j):*Z_⊥ )]^'
          T1 = T1 :- ScapPXYbarZperp[ind1+1].M * SstarUZperpinvZperpZperp[ind2+1].M
        if (NFE)
          T1 = T1 :- CT_FEcapYbar[ind1+1].M ' CTFEU[ind2+1].M
      }
      if (ind2) {
        retval = retval + FillingT0[max((ind1,ind2))+1,min((ind1,ind2))+1].M * _beta
        if (cols(T1))
          retval = retval + T1 * *pbetav  // - x*beta components
      } else
        retval = FillingT0[ind1+1,1].M :+ T1 * v  // y component

      if (Repl.Yendog[ind1+1] & Repl.Yendog[ind2+1])
        for (i=Clust.N;i;i--) {
          S = (*pinfoCapData)[i,]', (.\.)
          retval[i,] = retval[i,] - colsum(v :* cross(SstarUPX[ind1+1].M[|S|], SstarUMZperp[ind2+1].M[|S|]) * *pbetav)
        }

    }
  }
  if (_jk) retval[,1] = *_panelsum(*pcol(Repl.PXZ, ind1), Repl.y1 - Repl.Z * betas[,1], *pinfoCapData)
  return(&retval)
}

void boottest::InitWRE() {  // stuff done only once that knits together results from DGP and Repl regression prep
  real scalar i, g, j; real matrix X1g, X2g, Zperpg, X1parg, Zg, ZR1g, Y2g, Y2parg, S; real colvector y1g
  if (granular==0) {
    if (jk==0) {
      SbootZR1X = Sbooty1X = SbootZX = SbootY2X = SbootY2parX = smatrix(Nstar)
      if (!bootneqcap) pSbootXX = &ScapXX
      for (i=Repl.kZ; i>=0; i--)
        SstarUX[i+1].M = J(DGP.kX, Nstar, 0)
    }

    if (robust & bootstrapt) {
      CT_FEcapYbar = ScapYbarX = ScapPXYbarZperp = smatrix(Repl.kZ+1)
      FillingT0 = smatrix(Repl.kZ+1, Repl.kZ+1)  // fixed component of groupwise term in sandwich filling

      ScapYbarX.M = J(Ncap, DGP.kX, 0)
      invXXXX1par = DGP.invXX * (cross(*DGP.pX1, *Repl.pX1par) \ cross(DGP.X2, *Repl.pX1par))
      for (i=Repl.kZ; i; i--) {
        ScapYbarX[i+1].M = J(Ncap, DGP.kX, 0)
        ScapPXYbarZperp[i+1].M = J(Ncap, cols(*DGP.pZperp), 0)
        for (j=i; j>=0; j--)
          FillingT0[i+1,j+1].M = J(Ncap,1,0)
      }

      ScapXX = ScapXZperp = ScapX1parReplX = smatrix(Ncap)
      
      for (g=Ncap;g;g--) {
        S = (*pinfoCapData)[g,]'
        X1g = *pXS(*DGP.pX1,S)
        X2g = *pXS( DGP.X2 ,S)
        X1parg = *pXS(*Repl.pX1par,S)
        Zperpg = *pXS(*DGP.pZperp,S)
        ScapXX[g].M = cross(X1g, X2g); ScapXX[g].M = cross(X1g, X1g), ScapXX[g].M \ ScapXX[g].M', cross(X2g, X2g)
        ScapXZperp[g].M = cross(X1g,Zperpg) \ cross(X2g,Zperpg)
        ScapX1parReplX[g].M = cross(X1parg, X1g) , cross(X1parg, X2g)
        
        if (jk==0 & !bootneqcap) {  // if bootstrapping and intersection-of-error clusterings same, and using optimized non-jk code, then integrate by-boot and by-cap constructions for speed
          y1g = *pXS(DGP.y1,S)
          Zg = *pXS(DGP.Z,S)
          ZR1g = *pXS(*DGP.pZR1,S)
          Y2g = *pXS( DGP.Y2 ,S)
          SbootZR1X[g].M = cross(ZR1g, X1g) , cross(ZR1g, X2g)
          Sbooty1X[g].M = cross(y1g, X1g) , cross(y1g, X2g)
          SbootZX[g].M = cross(Zg, X1g) , cross(Zg, X2g)
          SbootY2parX[g].M = Repl.RparY ' (SbootY2X[g].M = cross(Y2g, X1g) , cross(Y2g, X2g))
        }
      }
    }

    if (jk==0 & (bootneqcap | !(robust & bootstrapt))) {
      pSbootXX = &smatrix(Nstar)
      for (g=Nstar;g;g--) {
        S = infoBootData[g,]'

        X1g = *pXS(*DGP.pX1,S)
        X2g = *pXS( DGP.X2 ,S)
        (*pSbootXX)[g].M = cross(X1g, X2g); (*pSbootXX)[g].M = cross(X1g, X1g), (*pSbootXX)[g].M \ (*pSbootXX)[g].M', cross(X2g, X2g)

        y1g = *pXS(DGP.y1,S)
        Zg = *pXS(DGP.Z,S)
        ZR1g = *pXS(*DGP.pZR1,S)
        Y2g = *pXS( DGP.Y2 ,S)
        SbootZR1X[g].M = cross(ZR1g, X1g) , cross(ZR1g, X2g)
        Sbooty1X[g].M = cross(y1g, X1g) , cross(y1g, X2g)
        SbootZX[g].M = cross(Zg, X1g) , cross(Zg, X2g)
        SbootY2parX[g].M = Repl.RparY ' (SbootY2X[g].M = cross(Y2g, X1g) , cross(Y2g, X2g))
      }
    }
  }
}

void boottest::PrepWRE() {
  real scalar i, j, g; pointer (real colvector) scalar pu; real matrix ZU2ddotpar, PiddotRparY, Pidddot; real colvector PiddotRparYi, tmp, ScapXXi, _r

  if (null) {
    DGP.Estimate(jk, _r = *pr1 \ *pr)
    DGP.MakeResiduals(0)
  } else
    _r = *pr1  // slightly inefficient to reassign every time

  XYbar = DGP.Xy1bar, (invXXXZbar = Repl.XZ - DGP.XU2ddot * Repl.RparY)
  invXXXZbar = Repl.invXX * invXXXZbar
  ZU2ddotpar = (Repl.ZY2 - Repl.XZ'DGP.Piddot) * Repl.RparY
  YbarYbar = DGP.deltadddot'Repl.XZ - DGP.y1barU2ddot * Repl.RparY
  YbarYbar = DGP.y1bary1bar, YbarYbar \ 
             YbarYbar'     , (Repl.ZZ - ZU2ddotpar' - ZU2ddotpar + Repl.RparY ' DGP.U2ddotU2ddot * Repl.RparY)

  if (granular==0)
    PiddotRparY = DGP.Piddot * Repl.RparY

  if (robust & bootstrapt) {  // for WRE replication regression, prepare for CRVE
    if (granular | NFE) PXZbar = *pX12B(*Repl.pX1, Repl.X2, invXXXZbar)

    if (granular==0) {
      Pidddot = (DGP.perpRperpX'DGP.deltaX \ J(kX2, 1, 0)) + DGP.Piddot * DGP.deltaY
      for (g=Ncap;g;g--)
        ScapYbarX.M[g,] = Pidddot'ScapXX[g].M  // S_cap(M_Zperp*y1 :* P_(MZperpX)])

      for (i=Repl.kZ; i; i--) {
        PiddotRparYi = PiddotRparY[,i]
        tmp = invXXXX1par[,i] + PiddotRparYi

        for (g=Ncap;g;g--) {
          ScapPXYbarZperp[i+1].M[g,] = tmp ' ScapXZperp[g].M  // S_cap(M_Zperp*y1 :* P_(MZperpX)])
          ScapYbarX      [i+1].M[g,] = ScapX1parReplX[g].M[i,] + PiddotRparYi ' ScapXX    [g].M  // S_cap(M_Zperp[Z or y1] :* P_(MZperpX)])

          ScapXXi = ScapXX[g].M  * tmp
          for (j=i; j; j--)
            FillingT0[i+1,j+1].M[g] = ScapX1parReplX[g].M[j,] * tmp + PiddotRparY[,j] ' ScapXXi
          FillingT0  [i+1,  1].M[g] =                                  DGP.deltadddot ' ScapXXi
        }

        if (NFE)
          CT_FEcapYbar[i+1].M = crosstabFE(*pcol(PXZbar,i), *pinfoCapData) :* invFEwt
      }
    }
  }

  if (kappa!=1 | LIML | robust==0)
    SstarUYbar = smatrix(Repl.kZ+1)
  pU2parddot = pXB(DGP.U2ddot.M, Repl.RparY)

  if (granular | jk) {  // construct objects in a way that involves explicitly constructing residuals
    if (kappa!=1 | LIML | robust==0 | granular)
      Zbar = *Repl.pX1par + DGP.Y2bar * Repl.RparY

    for (i=Repl.kZ; i>=0; i--) {  // precompute various clusterwise sums
      pu = i? pcol(*pU2parddot,i) : &DGP.u1dddot.M

      // S_star(u :* X), S_star(u :* Zperp) for residuals u for each endog var; store transposed
      SstarUX                     [i+1].M = *_panelsum2(*Repl.pX1, Repl.X2, *pu, infoBootData)'
      SstarUXinvXX                [i+1].M = Repl.invXX * SstarUX[i+1].M
      if (kappa!=1 | LIML | bootstrapt) {
        if (Repl.Yendog[i+1])
          SstarUZperpinvZperpZperp[i+1].M = Repl.invZperpZperp * (SstarUZperp[i+1].M = *_panelsum(*Repl.pZperp, *pu, infoBootData)')

        if (NFE)
          CTFEU[i+1].M = crosstabFE(*pu, infoBootData)
      }

      if (kappa!=1 | LIML | robust==0) {
        SstarUYbar[i+1].M = *_panelsum2(DGP.y1bar, Zbar, *pu, infoBootData)
        for (j=i; j>=0; j--)
          SstarUU[i+1,j+1].M = *_panelsum(j? *pcol(*pU2parddot,j) : DGP.u1dddot.M, *pu, infoBootData)
      }

      if (robust & bootstrapt & (granular==0 | Nw==1) & Repl.Yendog[i+1]) {
        SstarUPX[i+1].M = Repl.XinvXX * SstarUX[i+1].M
        SstarUMZperp[i+1].M = *Repl.pZperp * SstarUZperpinvZperpZperp[i+1].M
        if (Nobs == Nstar)  // subtract "crosstab" of observation by cap-group of u
          _diag(SstarUMZperp[i+1].M, diagonal(SstarUMZperp[i+1].M) - (i? *pcol(*pU2parddot,i) : DGP.u1dddot.M))  // case: bootstrapping by observation
        else
          if (i)
            for (g=Nobs;g;g--)
              SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] = SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] - (*pU2parddot)[g,i]
          else
            for (g=Nobs;g;g--)
              SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] = SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] - DGP.u1dddot.M[g]

        if (NFE)
          SstarUMZperp[i+1].M = SstarUMZperp[i+1].M + (invFEwt :* CTFEU[i+1].M)[*pFEID,]  // CT_(*,FE) (U ̈_(∥j) ) (S_FE S_FE^' )^(-1) S_FE
      }

      if (robust & bootstrapt & granular==0 & Repl.Yendog[i+1] & !(NClustVar == NBootClustVar & !subcluster))  // Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u:*X and u:*Zperp (and times invXX or invZperpZperp)
        for (g=Nstar;g;g--)
          SCTcapUXinvXX[i+1,g].M = *_panelsum(Repl.XinvXX, *pu, infoCTCapstar[g].M)
    }
  } else {  // for coarse clustering with no jackknife, construct objects while minimizing O(N) operations
    if (kappa!=1 | LIML | robust==0)
      Zbar = *Repl.pX1par + DGP.Y2bar * Repl.RparY

    for (i=Repl.kZ; i>=0; i--) {  // precompute various clusterwise sums
      pu = i? pcol(*pU2parddot,i) : &DGP.u1dddot.M

      // S_star(u :* X), S_star(u :* Zperp) for residuals u for each endog var; store transposed
      if (i) {
        PiddotRparYi = PiddotRparY[,i]
        for (g=Nstar;g;g--)
          SstarUX[i+1].M[,g] = SbootY2parX[g].M[i,]' - (*pSbootXX)[g].M * PiddotRparYi
      } else
        for (g=Nstar;g;g--)
          SstarUX.M[,g] = (Sbooty1X[g].M - _r ' SbootZR1X[g].M - DGP.beta.M ' SbootZX[g].M + DGP.deltaY ' (SbootY2X[g].M - DGP.Piddot ' (*pSbootXX)[g].M))'

      SstarUXinvXX[i+1].M = Repl.invXX * SstarUX[i+1].M
      if (kappa!=1 | LIML | bootstrapt) {
        if (Repl.Yendog[i+1])
          SstarUZperpinvZperpZperp[i+1].M = Repl.invZperpZperp * (SstarUZperp[i+1].M = *_panelsum(*Repl.pZperp, *pu, infoBootData)')

        if (NFE)
          CTFEU[i+1].M = crosstabFE(*pu, infoBootData)
      }

      if (kappa!=1 | LIML | robust==0) {
        SstarUYbar[i+1].M = *_panelsum2(DGP.y1bar, Zbar, *pu, infoBootData)
        for (j=i; j>=0; j--)
          SstarUU[i+1,j+1].M = *_panelsum(j? *pcol(*pU2parddot,j) : DGP.u1dddot.M, *pu, infoBootData)
      }

      if (robust & bootstrapt & (granular==0 | Nw==1) & Repl.Yendog[i+1]) {
        SstarUPX[i+1].M = Repl.XinvXX * SstarUX[i+1].M
        SstarUMZperp[i+1].M = *Repl.pZperp * SstarUZperpinvZperpZperp[i+1].M
        if (Nobs == Nstar)  // subtract "crosstab" of observation by cap-group of u
          _diag(SstarUMZperp[i+1].M, diagonal(SstarUMZperp[i+1].M) - (i? *pcol(*pU2parddot,i) : DGP.u1dddot.M))  // case: bootstrapping by observation
        else
          if (i)
            for (g=Nobs;g;g--)
              SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] = SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] - (*pU2parddot)[g,i]
          else
            for (g=Nobs;g;g--)
              SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] = SstarUMZperp[i+1].M[g,(*pIDBootData)[g]] - DGP.u1dddot.M[g]

        if (NFE)
          SstarUMZperp[i+1].M = SstarUMZperp[i+1].M + (invFEwt :* CTFEU[i+1].M)[*pFEID,]  // CT_(*,FE) (U ̈_(∥j) ) (S_FE S_FE^' )^(-1) S_FE
      }

      if (robust & bootstrapt & granular==0 & Repl.Yendog[i+1] & !(NClustVar == NBootClustVar & !subcluster))  // Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u:*X and u:*Zperp (and times invXX or invZperpZperp)
        for (g=Nstar;g;g--)
          SCTcapUXinvXX[i+1,g].M = *_panelsum(Repl.XinvXX, *pu, infoCTCapstar[g].M)
    }
  }
}

void boottest::MakeWREStats(real scalar w) {
  real scalar c, b, i, _jk
  real colvector numer_b
  real rowvector numerw, val, YY11, YY12, YY22, YPXY11, YPXY12, YPXY22, x11, x12, x21, x22, kappas, YY12YPXY12
  real matrix deltanumer, Jcap, J_b, Jcaps, vec
  struct smatrix rowvector A
  pragma unset vec; pragma unset val

  _jk = jk & w==1

  if (Repl.kZ == 1) {  // optimized code for 1 coefficient in bootstrap regression
    if (LIML) {
      YY11   = HessianFixedkappa(0, 0, 0, _jk)  // kappa=0 => Y*MZperp*Y
      YY12   = HessianFixedkappa(0, 1, 0, _jk)
      YY22   = HessianFixedkappa(1, 1, 0, _jk)
      YPXY11 = HessianFixedkappa(0, 0, 1, _jk)  // kappa=1 => Y*PXpar*Y
      YPXY12 = HessianFixedkappa(0, 1, 1, _jk)
      YPXY22 = HessianFixedkappa(1, 1, 1, _jk)
      YY12YPXY12 = YY12 :* YPXY12
      x11 = YY22 :* YPXY11 - YY12YPXY12      // elements of YYstar^-1 * YPXYstar up to factor of det(YYstar)
      x12 = YY22 :* YPXY12 - YY12 :* YPXY22
      x21 = YY11 :* YPXY12 - YY12 :* YPXY11
      x22 = YY11 :* YPXY22 - YY12YPXY12
      kappas = .5 * (x11 + x22); kappas = 1 :/ (1 :- (kappas - sqrt(kappas:*kappas - x11:*x22 + x12:*x21)) :/ (YY11 :* YY22 - YY12 :* YY12))  // solve quadratic equation for smaller eignenvalue; last term is det(YYstar)
      if (Fuller) kappas = kappas :- Fuller / (_Nobs - kX)
      betas = (kappas :* (YPXY12 - YY12) + YY12) :/ (As = kappas :* (YPXY22 - YY22) + YY22)
    } else
      betas = HessianFixedkappa(1, 0, kappa, _jk) :/ (As = HessianFixedkappa(1, 1, kappa, _jk))

    if (null)
      numerw = betas :+ (Repl.Rt1 - *pr) / Repl.RRpar
    else {
      numerw = betas :- betas[1]
      if (w==1)
        numerw[1] = betas[1] + (Repl.Rt1 - *pr) / Repl.RRpar
    }
    storeWtGrpResults(pnumer, w, numerw)

    if (bootstrapt) {
      if (robust) {
        Jcaps = *Filling(1, betas, _jk) :/ As
        for (c=1; c<=NErrClustCombs; c++) {  // sum sandwich over error clusterings
          if (NClustVar != 1 & rows(Clust[c].order))
            Jcaps = Jcaps[Clust[c].order,]
          Jcap = *_panelsum(Jcaps, Clust[c].info)
          _clustAccum(denom.M, c, colsum(Jcap:*Jcap))
        }
      } else
        denom.M = (HessianFixedkappa(0,0,0, _jk) - 2 * betas :* HessianFixedkappa(0, 1, 0, _jk) + betas:*betas :* HessianFixedkappa(1, 1, 0, _jk)) / _Nobs :/ As  // classical error variance

      storeWtGrpResults(pDist, w, sqrt? numerw:/sqrt(denom.M) : numerw :* numerw :/ denom.M)
      denom.M = Repl.RRpar * Repl.RRpar * denom.M[1]  // not a bug
    }

  } else {  // WRE bootstrap for more than 1 coefficeint in bootstrap regression

    betas = J(Repl.kZ, cols(v), 0)
    A = smatrix(cols(v))

    if (LIML) {
      for (i=Repl.kZ;i>=0;i--) {
        YYstar  [i+1].M = HessianFixedkappa(0..i      , i, 0, _jk)  // kappa=0 => Y*MZperp*Y
        YPXYstar[i+1].M = HessianFixedkappa(i..Repl.kZ, i, 1, _jk)  // kappa=1 => Y*PXpar*Y
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
      deltanumer = HessianFixedkappa(1..Repl.kZ, 0, kappa, _jk)
      for (i=Repl.kZ;i;i--)
        deltadenom[i].M = HessianFixedkappa(1..i, i, kappa, _jk)
      for (b=cols(v); b; b--) {
        for (i=Repl.kZ;i;i--)
          deltadenom_b[|.,i\i,i|] = deltadenom[i].M[,b] // fill uppper triangle, which is all that invsym() looks at
        betas[,b]  = (A[b].M = invsym(deltadenom_b)) * deltanumer[,b]
      }
    }

    if (robust)
      for(i=Repl.kZ;i;i--)
        Zyi[i].M = *Filling(i, betas, _jk)
    else
      for (i=Repl.kZ;i>=0;i--)
        YYstar[i+1].M = HessianFixedkappa(i..Repl.kZ, i, 0, _jk)  // kappa=0 => Y*MZperp*Y

    for (b=cols(v); b; b--) {
      numer_b = null | w==1 & b==1? (Repl.RRpar * betas[,b] + Repl.Rt1) - *pr : Repl.RRpar * (betas[,b] - betas[,1])
      if (bootstrapt) {
        if (robust) {  // Compute denominator for this WRE test stat
          for(i=Repl.kZ;i;i--)
            _Jcap[,i] = Zyi[i].M[,b]
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


//
// WCR/WCU
//

// For non-WRE, construct stuff that depends linearly or quadratically on r, possibly by interpolation
void boottest::MakeInterpolables() {
  real scalar h1, h2, d1, d2, c; real matrix tmp; real colvector Delta, newPole

  if (interpolable) {
    if (rows(anchor)==0) {  // first call? save current r as permanent anchor for interpolation
      _MakeInterpolables(anchor = *pr)
      numer0 = *pnumer
      if (interpolate_u) ustar0 = *pustar
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
            dudr[h1].M = (*pustar - ustar0) / poles[h1]
          if (robust & !purerobust)  // df > 1 for an ARubin test with >1 instruments. 
            for (d1=1;d1<=df;d1++) {
              for (c=1;c<=NErrClustCombs;c++) {
                dJcddr[h1].M[c,d1].M = ((*pJcd)[c,d1].M - Jcd0[c,d1].M) / poles[h1]
                for (d2=1;d2<=d1;d2++) {
                                tmp =       colsum(Jcd0[c,d1].M :* dJcddr[h1].M[c,d2].M)
                  if (d1 != d2) tmp = tmp + colsum(Jcd0[c,d2].M :* dJcddr[h1].M[c,d1].M)  // for diagonal items, faster to just double after the c loop
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
          pustar = &(ustar0 + dudr.M * Delta[1] + dudr    [2].M * Delta[2])
        }
      }
    } else {  // routine linear interpolation if the anchors not moved
      Delta = *pr - anchor
      numerw = numer0 + dnumerdr.M * Delta[1]; if (q > 1) numerw = numerw + dnumerdr[2].M * Delta[2]
      if (interpolate_u) {
        pustar = &(ustar0 + dudr.M * Delta[1]); if (q > 1) pustar = &(*pustar + dudr[2].M * Delta[2])
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
  real scalar d, c, _jk; real colvector uXAR; pointer (real matrix) scalar pustarXAR, ptmp

  if (ML)
    uXAR = *pSc * (AR = *pA * *pR')
  else {
    if (ARubin)
      DGP.Estimate(jk, r)
    else if (kappa) {
      if (null) {  // in score bootstrap for IV/GMM, if imposing null, then DGP constraints, kappa, Hessian, etc. do vary with r and must be set now
        DGP.Estimate(jk, *pr1 \ r)
        DGP.InitTestDenoms()
      }
    } else  // regular OLS
      DGP.Estimate(jk, null? *pr1 \ r : *pr1)

    DGP.MakeResiduals(jk)
  }

  for (_jk=jk; _jk>=0; _jk--) { // for jackknife, run 2nd time on full original sample to get test stat

    if (!ML & (scoreBS | (robust & granular * !(jk & !_jk) < NErrClustCombs)))
      uXAR = DGP.u1ddot[1+_jk].M :* pM->XAR

    SuwtXA = scoreBS?
                (B? 
                   (NClustVar? *_panelsum(uXAR, infoBootData) : 
                               uXAR                  ) :
                   colsum(uXAR)')                             :
                *DGP.pA * *_panelsum2(*pX1, *pX2, DGP.u1ddot[1+_jk].M, infoBootData)'  // same calc as in score BS but broken apart to grab intermediate stuff, and assuming residuals defined; X2 empty except in Anderson-Rubin

    if (robust & granular < NErrClustCombs & bootstrapt) {
      pustarXAR = _panelsum(uXAR, *pinfoAllData)  // collapse data to all-boot & error-cluster-var intersections. If no collapsing needed, _panelsum() will still fold in any weights

      if (B) {
        if (scoreBS)
          for (d=df;d;d--)
            Kd[d].M = JNcapNstar  // inefficient, but not optimizing for the score bootstrap
        else
          for (d=df;d;d--)
            Kd[d].M = *_panelsum2(*pX1, *pX2, *pcol(DGP.XAR,d), *pinfoCapData) * SuwtXA  // final term in (64), for c=intersection of all error clusters
        if (NFE & FEboot==0)
          CT_WE = crosstabFE(DGP.u1ddot[1+_jk].M, infoBootData)

        for (d=df;d;d--) {  // subtract crosstab of u:*XAR wrt bootstrapping cluster combo and all-cluster-var intersections
          crosstabCapstarMinus(Kd[d].M, *pcol(*pustarXAR,d))
          if (NFE & FEboot==0)
            Kd[d].M = Kd[d].M + pM->CT_XAR[d].M ' (invFEwt :* CT_WE)  // middle term of (64)
          if (scoreBS)
            Kd[d].M = Kd[d].M - ClustShare * colsum(Kd[d].M)  // recenter
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

    MakeNumerAndJ(1, _jk, r)  // compute J = K * v; if Nw > 1, then this is for 1st group; if interpolating, it is only group, and may be needed now to prep interpolation
  }
}

// compute stuff depending linearly on v, needed to prep for interpolation
// _jk = 0 for non-jackknife, or for non-jk, full original-sample estimation of test stat for jk
void boottest::MakeNumerAndJ(real scalar w, real scalar _jk, | real colvector r) {  // called to *prepare* interpolation, or when w>1, in which case there is no interpolation
  real scalar c, d; real matrix _v; real matrix betadev

  if (jk & !_jk) {
    if (!ARubin & null) 
      (*pnumer)[,1] = v_sd * (  // full original-sample test stat numerator
                          scoreBS?
                             (B? 
                               colsum(SuwtXA)' : 
                               SuwtXA         ) :
                             (robust==0 | granular | purerobust?
                               *pR * (betadev = rowsum(SuwtXA)) :
                               rowsum(*pR * SuwtXA)))
  } else {
    numerw = scoreBS?
               (B? 
                 cross(SuwtXA, v) : 
                 SuwtXA * v_sd    ) :
               (robust==0 | granular | purerobust?
                  *pR * (betadev = SuwtXA * v) :
                 (*pR * SuwtXA) * v)

    if (w==1) {
      if      ( ARubin) numerw[,1] = v_sd * DGP.Rpar * DGP.beta.M[|kX1+1\.|]  // coefficients on excluded instruments in ARubin OLS
      else if (null==0) numerw[,1] = v_sd * (*pR * (ML? beta : pM->Rpar * pM->beta.M) - r)  // Analytical Wald numerator; if imposing null then numer[,1] already equals this. If not, then it's 0 before this.
    }

    storeWtGrpResults(pnumer, w, numerw)
  }

  if (interpolate_u) {
    pustar = B? &(v :* DGP.u1ddot[1+_jk].M) : &DGP.u1ddot[1+_jk].M
    if (scoreBS) pustar = &(*pustar :- colsum(*pustar) * ClustShare)  // Center variance if interpolated
            else pustar = &(*pustar - *pX12B(*pX1, *pX2, betadev))  // residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline & Santos eq (11))
  }

  if (B & robust & bootstrapt)
    if (jk & !_jk)  // fill first col of K's with values for original, full-sample stat
      for (c=NErrClustCombs; c; c--)
        for (d=df;d;d--)
          ((*pJcd)[c,d].M)[,1] = rowsum(Kcd[c,d].M) * v_sd
    else {
      if (granular | purerobust)  // optimized treatment when bootstrapping by many/small groups
        if (purerobust)
          pustar = &(*partialFE(&(DGP.u1ddot[1+_jk].M :* v)) - *pX12B(*pX1, *pX2, betadev))
        else {  // clusters small but not all singletons
          if (NFE & FEboot==0) {
            pustar = partialFE(&(DGP.u1ddot.M :* v[*pIDBootData,]))
            for (d=df;d;d--)
              (*pJcd)[1,d].M = *_panelsum(*pustar, pM->WXAR[d].M, *pinfoCapData)                                           - *_panelsum2(*pX1, *pX2, pM->WXAR[d].M, *pinfoCapData) * betadev
          } else {
            _v = v[*pIDBootAll,]
            for (d=df;d;d--)
              (*pJcd)[1,d].M = *_panelsum(*_panelsum(DGP.u1ddot[1+_jk].M, pM->WXAR[d].M, *pinfoAllData) :* _v, infoErrAll) - *_panelsum2(*pX1, *pX2, pM->WXAR[d].M, *pinfoCapData) * betadev
          }
        }
      for (c=NErrClustCombs; c>granular; c--)
        for (d=df;d;d--)
          (*pJcd)[c,d].M = Kcd[c,d].M * v
    }
}

void boottest::MakeNonWREStats(real scalar w) {
  real scalar i, c, j, k; real matrix ustar2, tmp, invdenom; real colvector numer_k, ustar_k; pointer (real matrix) scalar pAR; real rowvector t1, t2, t12

  if (w > 1) MakeNumerAndJ(w, jk)

  if (bootstrapt == 0) return

  if (robust) {
    if (interpolating==0) {  // these quadratic computations needed to *prepare* for interpolation but are superseded by interpolation once it is going
      if (purerobust)
        ustar2 = *pustar :* *pustar
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
    } else if (df == 2) {  // hand-code 2D numer'inv(denom)*numer
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
        numer_k = numerw[,k]
        (*pDist)[k+WeightGrpStart[w]-1] = numer_k ' invsym(tmp) * numer_k  // in degenerate cases, cross() would turn cross(.,.) into 0
      }
      if (w==1)
        statDenom = makesymmetric(tmp')  // original-sample denominator
    }

  } else { // non-robust

    pAR = ML? &AR : &pM->AR
    if (df == 1) {  // optimize for one null constraint
      denom.M = *pR * *pAR

      if (ML==0)
        denom.M = denom.M :* colsum(*pustar :* *pustar)
      storeWtGrpResults(pDist, w,  numerw :/ sqrt(denom.M))
      if (w==1)
        statDenom = denom.M[1]  // original-sample denominator
    } else {
      denom.M = *pR * *pAR

      if (ML) {
        for (k=cols(v); k; k--) {
          numer_k = numerw[,k]
          (*pDist)[k+WeightGrpStart[w]-1] = cross(numer_k, invsym(denom.M), numer_k)
        }
        if (w==1)
          statDenom = denom.M  // original-sample denominator
      } else {
        invdenom = invsym(denom.M)
        for (k=cols(v); k; k--) {
          numer_k = numerw[,k]
          ustar_k  = (*pustar)[,k]
          (*pDist)[k+WeightGrpStart[w]-1] = cross(numer_k, invdenom * numer_k) / (tmp = cross(ustar_k, ustar_k))
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
  real matrix retval; real scalar i, j, tmp; real colvector _FEID, _v, vw
  retval = J(NFE, rows(info), 0)
  if (haswt) {
    vw = v :* sqrtwt
    if (cols(info) | rows(info)==rows(v))
      for (i=cols(retval);i;i--) {
        _FEID = panelsubmatrix(*pFEID, i, info)
        _v    = panelsubmatrix(vw    , i, info)
        for (j=rows(_FEID);j;j--) {
          tmp = _FEID[j] 
          retval[tmp,i] = retval[tmp,i] + _v[j]
        }
      }
    else  // "robust" case, no clustering
      for (i=cols(retval);i;i--)
        retval[(*pFEID)[i],i] = vw[i]
  } else
    if (cols(info) | rows(info)==rows(v))
      for (i=cols(retval);i;i--) {
        _FEID = panelsubmatrix(*pFEID, i, info)
        _v    = panelsubmatrix(v     , i, info)
        for (j=rows(_FEID);j;j--) {
          tmp = _FEID[j] 
          retval[tmp,i] = retval[tmp,i] + _v[j]
        }
      }
    else  // "robust" case, no clustering
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
  real scalar kappa, real scalar ARubin, real scalar null, real scalar scoreBS, real scalar jk, string scalar auxweighttype, string scalar ptype, string scalar statistic, string scalar madjtype, real scalar NumH0s,
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
  M.setjk(jk)
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
