mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

struct structFE {
	pointer (struct structFE scalar) scalar next
	real colvector is, wt
}

real colvector zoot1(real scalar NFE, real colvector FEID) {
	real colvector o, sortID, _FEID, wtFE; real scalar Nobs, i_FE, i, j; struct structFE rowvector FEs
	Nobs = rows(FEID)
	sortID = (FEID)[o = order(FEID, 1)]
	i_FE = 1; j = Nobs; _FEID = wtFE = J(Nobs, 1, 1)
	FEs = structFE(NFE)
	for (i=Nobs-1;i;i--) {
		if (sortID[i] != sortID[i+1]) {
			FEs[i_FE].is = o[|i+1\j|]
			FEs[i_FE].wt = J(j-i,1,1/(j-i))
			wtFE[FEs[i_FE].is] = FEs[i_FE].wt
			j = i
			++i_FE
		}
		_FEID[o[i]] = i_FE
	}
	FEs[NFE].is = FEs[NFE].is = o[|.\j|]
	FEs[NFE].wt = J(j-i,1,1/(j-i))
	wtFE[FEs[NFE].is] = FEs[NFE].wt
	return (_FEID)
}

real colvector zoot2(real scalar NFE, real colvector FEID) {
	real colvector _FEID, wtFE; real scalar Nobs, id, i, N; struct structFE rowvector FEs; transmorphic scalar A, loc
	Nobs = rows(FEID)
	_FEID = wtFE = J(Nobs, 1, 1)
	FEs = structFE(NFE)
	A = asarray_create("real", 1, NFE)
	asarray_notfound(A, J(0,1,0))
	for (i=Nobs;i;i--)
		asarray(A, FEID[i], asarray(A, FEID[i]) \ i)
	id = 0
	for (loc=asarray_first(A); loc!=NULL; loc=asarray_next(A, loc)) {
		N = rows(FEs[++id].is = asarray_contents(A, loc))
		_FEID[FEs[id].is] = J(N, 1, id)
		wtFE[FEs[id].is] = FEs[id].wt = J(N,1,1/N)
	}
	return (_FEID)
}

real colvector zoot3(real scalar NFE, real colvector FEID) {
	real colvector _FEID, wtFE; real scalar Nobs, id, i, N; struct structFE rowvector FEs; transmorphic scalar A, loc
	Nobs = rows(FEID)
	_FEID = wtFE = J(Nobs, 1, 1)
	FEs = structFE(NFE)
	A = asarray_create("real", 1, NFE)
	asarray_notfound(A, J(0,1,0))
	for (i=Nobs;i;i--)
		if (!rows(asarray(A, FEID[i])))
			asarray(A, FEID[i], selectindex(FEID:==FEID[i]))
	id = 0
	for (loc=asarray_first(A); loc!=NULL; loc=asarray_next(A, loc)) {
		N = rows(FEs[++id].is = asarray_contents(A, loc))
		_FEID[FEs[id].is] = J(N, 1, id)
		wtFE[FEs[id].is] = FEs[id].wt = J(N,1,1/N)
	}
	return (_FEID)
}

timer_clear()
FEID = 2 * floor(runiform(1000000,1) * 1000000/10)
timer_on(1)
t1 = zoot1(1000000/10, FEID)
timer_off(1)

timer_on(2)
t2 = zoot2(1000000/10, FEID)
timer_off(2)

timer_on(3)
t3 = zoot3(1000000/10, FEID)
timer_off(3)

timer()
end
