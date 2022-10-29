// clear
// set obs 6
// drawnorm x z e2 e1
// replace e1 = e1 + e2
// gen y2 = z + e2
// gen y1 = y2 + x + e1
// gen id = floor(_n-1/2)

scalar γ = 0

use D:\OneDrive\Desktop\t, clear

gen e2hat = .
gen e1hat = .
constraint 1 [y1]y2 = γ
forvalues g=1/3 {
  cmp (y1 = x y2) (y2 = x z) if id!=`g', ind(1 1) qui constr(1) nolr
  predict _e1, resid eq(#1)
  predict _e2, resid eq(#2)
  replace e1hat = _e1 if id==`g'
  replace e2hat = _e2 if id==`g'
  drop _e?
}
// replace e2hat = e2hat * sqrt(1 - 1 / 3)
// replace e1hat = e1hat * sqrt(1 - 1 / 3)

cmp (y1 = x y2) (y2 = x z), ind(1 1) qui constr(1) nolr
predict xb1, eq(#1)
predict xb2, eq(#2)
replace xb1 = xb1 - γ * y2

cap mat drop res
mata v = st_data(.,"v*") * sqrt(1 - 1 / 3)
mata U2 = st_data(.,"e2hat")
mata u1 = st_data(.,"e1hat") + `=γ' * U2
mata xb1 = st_data(.,"xb1")
mata xb2 = st_data(.,"xb2")
forvalues b=1/8 {
  gen y2star = xb2              + e2hat * v`b' * sqrt(1 - 1 / 3)
  gen y1star = xb1 + γ * y2star + e1hat * v`b' * sqrt(1 - 1 / 3)
  ivregress 2sls y1star x (y2star = z), cluster(id)
  mata {
    y1 = st_data(.,"y1star")
    N = rows(y1)
    X1 = st_data(.,"x"), J(N,1,1)
    Y2 = st_data(.,"y2star")
    X2 = st_data(.,"z")
    
    _y1 = y1 - X1*invsym(X1'X1)*X1'y1
    _Y2 = Y2 - X1*invsym(X1'X1)*X1'Y2
    _X2 = X2 - X1*invsym(X1'X1)*X1'X2
    Z = _Y2
    X = _X2
    invsym(Y2'_X2*invsym(_X2'_X2)*_X2'Y2)*Y2'_X2*invsym(_X2'_X2)*_X2'y1
    invsym((xb2+U2+U2:*(v[,`b']:-1))'_X2*invsym(_X2'_X2)*_X2'Y2)*Y2'_X2*invsym(_X2'_X2)*_X2'y1
//     Y2'_X2*invsym(_X2'_X2)*_X2'Y2
  }
  mat res = nullmat(res) \ (_b[y2star] - γ) / _se[y2star]
  drop y?star
}
matlist res

ivregress 2sls y1 x (y2 = z), cluster(id)
boottest y2=γ, noci jk svmat
matlist r(dist)

// mata {  RLIML per RMNW
//   T = 0,0 \ 1,0 \ 0,1; t = 0 \ 0 \ 0
//   y1 = st_data(.,"y1")
//   N = rows(y2)
//   Y2 = st_data(.,"y2")
//   X1 = st_data(.,"x"), J(N,1,1)
//   X2 = st_data(.,"z")
//   X = X1, X2
//   Z = Y2, X1
//   Z1 = y1 - Z*t, Z*T
//   eigensystem(invsym(Z1'Z1) * (Z1'Z1 - Z1'X*invsym(X'X)*X'Z1), V=., L=.)
//   κ = 1/max(Re(L))
//   δ = (T * invsym(T'Z'(I(N)-κ*(I(N)-X*invsym(X'X)*X'))*Z*T) * ///
//                   T'Z'(I(N)-κ*(I(N)-X*invsym(X'X)*X'))*(y1-Z*t) + t)
//   δ
// }

   

