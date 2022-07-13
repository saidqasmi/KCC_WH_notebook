constrain_array = function(S_mean,Sigma_mod,Xo,Sigma_obs,Nres,centering_CX=T,ref_CX=NULL) {

	distribs = constrain(S_mean,Sigma_mod,Xo,Sigma_obs, Nres, centering_CX=centering_CX,ref_CX=ref_CX)

	full_names = names(S_mean)
	ns = length(full_names)
	year_names = substr(full_names,1,4)
	forc_names = sub(".*_","",full_names)
	locs_names = sub("_[a-z]{,3}$","",sub("^[0-9]{4}_","",full_names))
	year = unique(year_names)
	forc = unique(forc_names)
	locs = unique(locs_names)
	ny = length(year)
	nf = length(forc)
	nl = length(locs)

	sample_str = c("be",paste0("nres",1:Nres))

	XCons = array(NA,dim=c(ny,Nres+1,nl,nf,2),dimnames=list(#
						year=year,sample=sample_str,spatial=locs,forc=forc,#
						constrain=c("uncons","cons")))

	# Unconstrained
	Sigma_sqrt_uncons = matrix_sqrt(distribs$uncons$var)
	epsil_uncons = array(0,dim = c(ns,Nres+1))
	epsil_uncons[,-1] = rnorm(ns*Nres)
	X_uncons = distribs$uncons$mean%o%rep(1,Nres+1) + Sigma_sqrt_uncons %*% epsil_uncons
	dimnames(X_uncons) = list(year=full_names,sample=sample_str)

	# Constrained
	Sigma_sqrt_cons = matrix_sqrt(distribs$cons$var)
	epsil_cons = array(0,dim = c(ns,Nres+1))
	epsil_cons[,-1] = rnorm(ns*Nres)
	X_cons = distribs$cons$mean%o%rep(1,Nres+1) + Sigma_sqrt_cons %*% epsil_cons
	dimnames(X_cons) = list(year=full_names,sample=sample_str)


	for (iloc in locs) {
		for (iforc in forc) {
			for (iy in year) {
				tmp_name = paste(iy,iloc,iforc,sep="_")
				if (sum(tmp_name %in% full_names)) {
					XCons[iy,,iloc,iforc,"uncons"] = X_uncons[full_names==tmp_name,]
					XCons[iy,,iloc,iforc,"cons"]   = X_cons[full_names==tmp_name,]
				}
			}
		}
	}

	return(XCons)

}

