cap mata mata drop CT1()
cap mata mata drop CT2()

mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

real matrix CT1(real colvector X, real colvector ID1, real colvector ID2) {
	transmorphic CT, inds1, inds2, loc1, loc2; real scalar i, i1, i2, ind1, ind2; real colvector key; real matrix retval
		
	CT = asarray_create("real", 2); asarray_notfound(CT, 0)
	inds1 = asarray_create("real", 1)
	inds2 = asarray_create("real", 1)
	
	for (i=rows(X);i;i--) {
		key = (ind1 = ID1[i]), (ind2  = ID2[i])
		asarray(CT, key, asarray(CT, key) + X[i])
		asarray(inds1, ID1[i], 1) // add to catalogs of observed instances
		asarray(inds2, ID2[i], 1)
	}

	retval = J(asarray_elements(inds1), asarray_elements(inds2), 0)
	i1 = 0
	for (loc1=asarray_first(inds1); loc1!=NULL; loc1=asarray_next(inds1, loc1)) {
		ind1 = asarray_key(inds1, loc1)
		i1++
		i2 = 0
		for (loc2=asarray_first(inds2); loc2!=NULL; loc2=asarray_next(inds2, loc2)) {
			ind2 = asarray_key(inds2, loc2)
			retval[i1,++i2] = asarray(CT, (ind1,ind2))
		}
	}
	return(retval)
}

real matrix CT2(real colvector X, real colvector ID1, real colvector ID2) {
	class Factor scalar F
	F = _factor((ID1,ID2))
	return(rowshape(aggregate_sum(F, F.sort(X), ., ""), rows(uniqrows(F.keys[,1]))))
}

real matrix CT3(real colvector X, real colvector ID1, real colvector ID2) {
	real scalar i, j, t, N, N1, N2; real matrix retval; real colvector sortID1, sortID2, sortX, _ID2, _ID2i, Xi, o1, o2; real matrix info
	
	N = rows(X)
	
	sortID1 = ID1[o1 = order(ID1, 1)]
	sortID2 = ID2[o1]
	sortX   = X[o1]

	N1 = rows(info = _panelsetup(sortID1, 1))
	
	sortID2 = sortID2[o2 = order(sortID2, 1)]
	
	N2 = 1; _ID2 = J(N, 1, 1)
	for (i=N-1;i;i--) {
		if (sortID2[i] != sortID2[i+1])
			N2++
		_ID2[o2[i]] = N2 // ordinal version of ID2 expressed with respect to ID1-sorted data
	}

	retval = J(N1, N2, 0)
	for (i=N1;i;i--) {
		_ID2i = panelsubmatrix(_ID2, i, info)
		Xi    = panelsubmatrix(sortX    , i, info)
		for (j=rows(_ID2i);j;j--) {
			t = _ID2i[j] 
			retval[i,t] = retval[i,t] + Xi[j]
		}
	}
	return(retval)
}

real matrix CT4(real colvector X, real colvector ID1, real colvector ID2) {
	real colvector ID, o; real matrix info
	ID = ID1, ID2
	ID = ID[o = order(ID, (1,2)),]
	info = _panelsetup(ID, (1,2))
	return(1)
//	return(_panelsum(X[o], _panelsetup(ID, (1,2))))
}


rseed(9821374)
N = 1000000
N1 = 5
N2 = 5
X = runiform(N,1)
ID1 = rdiscrete(N,1,J(N1,1,1/N1)) :+ 5
ID2 = rdiscrete(N,1,J(N2,1,1/N2)) :+ 10
//CT1(X,ID1,ID2)
CT2(X,ID1,ID2)
CT3(X,ID1,ID2)
CT4(X,ID1,ID2)

timer_clear()
//timer_on(1);(void) CT1(X,ID1,ID2);timer_off(1)
timer_on(2);for(i=1;i;i--) (void) CT2(X,ID1,ID2);timer_off(2)
timer_on(3);for(i=1;i;i--) (void) CT3(X,ID1,ID2);timer_off(3)
timer_on(4);for(i=1;i;i--) (void) CT4(X,ID1,ID2);timer_off(4)
timer()
end
