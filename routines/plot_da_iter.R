## Plotting function specific to DA paper.


# Plot DA results: 90% CI of GMT change
plot_da_iter = function(Xo,X0,X1,ref_period=1870:1900,pres_period=1870:2020,method="quantile",temporal="delta",ylim=NULL) {
	par(mar=c(2.5,6,1,1),cex=1.2,font=2)
	if (method!="variance" & method!="quantile") { message("error in plot_da(): unknown method"); return }

	# Computation of required quantities
	forcings	= c("Obs","ALL","NAT","ANT","GHG","OA")
	colorf	= c("black","green","blue","red","brown","orange")
	colorf_rgb = col2rgb(colorf)/255
	colorf_half = rgb(colorf_rgb[1,],colorf_rgb[2,],colorf_rgb[3,],alpha=.3)
	colorf_rgb_half = .3*colorf_rgb + .7
	colorf_half = rgb(colorf_rgb_half[1,],colorf_rgb_half[2,],colorf_rgb_half[3,])
	nfc = length(forcings)
	Qtable = array(NA,dim=c(nfc,3,2),dimnames=list(forcing=forcings,#
																  quantile=c("Q05","mean","Q95"),#
																  constrain=c("uncons","cons")))
	avg_period = function(X,ref_period,pres_period,temporal="delta") {
		year = as.numeric(dimnames(X)$year)
		ndim = length(dim(X))
		if (ndim!=2 & ndim!=3) {
			message("Error in avg_period(): ndim > 3")
			return
		}
		if (temporal=="delta") {
			if (ndim==2) {
				avg_ref = apply(X[year%in%ref_period,],2:ndim,mean)
				avg_pres = apply(X[year%in%pres_period,],2:ndim,mean)
			} else { # ndim==3
				avg_ref = apply(X[year%in%ref_period,,],2:ndim,mean)
				avg_pres = apply(X[year%in%pres_period,,],2:ndim,mean)
			}
			output = avg_pres-avg_ref
		} else if (temporal=="trend") {
			ctrend = function(x,period) {
				year = names(x)
				tt = 1:length(year)
				return(lm(x~tt)$coef[2])
			}
			if (ndim==2) {
				output = apply(X[year%in%pres_period,],2,ctrend) * 100
			} else { # ndim==3
				output = apply(X[year%in%pres_period,,],2:3,ctrend) * 100
			}
		} else {
			print("Error in plot_da.R: unknonw temporal treatment")
			return()
		}
		return(output)
	}
	# Best estimates
	# Observed trend
	obs_avg = avg_period(Xo,ref_period,pres_period,temporal=temporal)
	Qtable["Obs","mean","uncons"] = obs_avg["median"]
	# Forced trend
	all_avg_uncons = avg_period(X0[,,"all",],ref_period,pres_period,temporal=temporal)
	Qtable["ALL","mean","uncons"] = all_avg_uncons["be","uncons"]
	all_avg_cons = avg_period(X1[,,"all",],ref_period,pres_period,temporal=temporal)
	Qtable["ALL","mean","cons"] = all_avg_cons["be","cons"]
	# Nat trend
	nat_avg_uncons = avg_period(X0[,,"nat",],ref_period,pres_period,temporal=temporal)
	Qtable["NAT","mean","uncons"] = nat_avg_uncons["be","uncons"]
	nat_avg_cons = avg_period(X1[,,"nat",],ref_period,pres_period,temporal=temporal)
	Qtable["NAT","mean","cons"] = nat_avg_cons["be","cons"]
	# Ant trend
	ant_avg_uncons = avg_period(X0[,,"all",]-X0[,,"nat",],ref_period,pres_period,temporal=temporal)
	Qtable["ANT","mean","uncons"] = ant_avg_uncons["be","uncons"]
	ant_avg_cons = avg_period(X1[,,"all",]-X1[,,"nat",],ref_period,pres_period,temporal=temporal)
	Qtable["ANT","mean","cons"] = ant_avg_cons["be","cons"]
	# GHG trend
	ghg_avg_uncons = avg_period(X0[,,"ghg",],ref_period,pres_period,temporal=temporal)
	Qtable["GHG","mean","uncons"] = ghg_avg_uncons["be","uncons"]
	ghg_avg_cons = avg_period(X1[,,"ghg",],ref_period,pres_period,temporal=temporal)
	Qtable["GHG","mean","cons"] = ghg_avg_cons["be","cons"]
	# OA trend
	oa_avg_uncons = avg_period(X0[,,"all",]-X0[,,"nat",]-X0[,,"ghg",],ref_period,pres_period,temporal=temporal)
	Qtable["OA","mean","uncons"] = oa_avg_uncons["be","uncons"]
	oa_avg_cons = avg_period(X1[,,"all",]-X0[,,"nat",]-X1[,,"ghg",],ref_period,pres_period,temporal=temporal)
	Qtable["OA","mean","cons"] = oa_avg_cons["be","cons"]
	# IC bounds
	if (method=="quantile") {
		Qtable["Obs",c("Q05","Q95"),"uncons"] = quantile(obs_avg[-1],c(.05,.95))
		Qtable["ALL",c("Q05","Q95"),"uncons"] = quantile(all_avg_uncons[-1,"uncons"],c(.05,.95))
		Qtable["NAT",c("Q05","Q95"),"uncons"] = quantile(nat_avg_uncons[-1,"uncons"],c(.05,.95))
		Qtable["ANT",c("Q05","Q95"),"uncons"] = quantile(ant_avg_uncons[-1,"uncons"],c(.05,.95))
		Qtable["GHG",c("Q05","Q95"),"uncons"] = quantile(ghg_avg_uncons[-1,"uncons"],c(.05,.95))
		Qtable[ "OA",c("Q05","Q95"),"uncons"] = quantile( oa_avg_uncons[-1,"uncons"],c(.05,.95))
		Qtable["ALL",c("Q05","Q95"),"cons"] = quantile(all_avg_cons[-1,"cons"],c(.05,.95))
		Qtable["NAT",c("Q05","Q95"),"cons"] = quantile(nat_avg_cons[-1,"cons"],c(.05,.95))
		Qtable["ANT",c("Q05","Q95"),"cons"] = quantile(ant_avg_cons[-1,"cons"],c(.05,.95))
		Qtable["GHG",c("Q05","Q95"),"cons"] = quantile(ghg_avg_cons[-1,"cons"],c(.05,.95))
		Qtable[ "OA",c("Q05","Q95"),"cons"] = quantile( oa_avg_cons[-1,"cons"],c(.05,.95))
	} else if (method=="variance") {
		Qtable["Obs",c("Q05","Q95"),"uncons"] = Qtable["Obs","mean","uncons"] + qnorm(.95)*c(-1,1)*sd(obs_avg[-1])
		Qtable["ALL",c("Q05","Q95"),"uncons"] = c(1,1)%o%Qtable["ALL","mean","uncons"] + qnorm(.95)*c(-1,1)%o%apply(all_avg_uncons[-1,],2,sd)
		Qtable["NAT",c("Q05","Q95"),"uncons"] = c(1,1)%o%Qtable["NAT","mean","uncons"] + qnorm(.95)*c(-1,1)%o%apply(nat_avg_uncons[-1,],2,sd)
		Qtable["ANT",c("Q05","Q95"),"uncons"] = c(1,1)%o%Qtable["ANT","mean","uncons"] + qnorm(.95)*c(-1,1)%o%apply(ant_avg_uncons[-1,],2,sd)
		Qtable["GHG",c("Q05","Q95"),"uncons"] = c(1,1)%o%Qtable["GHG","mean","uncons"] + qnorm(.95)*c(-1,1)%o%apply(ghg_avg_uncons[-1,],2,sd)
		Qtable[ "OA",c("Q05","Q95"),"uncons"] = c(1,1)%o%Qtable[ "OA","mean","uncons"] + qnorm(.95)*c(-1,1)%o%apply( oa_avg_uncons[-1,],2,sd)
		Qtable["Obs",c("Q05","Q95"),"cons"] = Qtable["Obs","mean","uncons"] + qnorm(.95)*c(-1,1)*sd(obs_avg[-1])
		Qtable["ALL",c("Q05","Q95"),"cons"] = c(1,1)%o%Qtable["ALL","mean","cons"] + qnorm(.95)*c(-1,1)%o%apply(all_avg_cons[-1,],2,sd)
		Qtable["NAT",c("Q05","Q95"),"cons"] = c(1,1)%o%Qtable["NAT","mean","cons"] + qnorm(.95)*c(-1,1)%o%apply(nat_avg_cons[-1,],2,sd)
		Qtable["ANT",c("Q05","Q95"),"cons"] = c(1,1)%o%Qtable["ANT","mean","cons"] + qnorm(.95)*c(-1,1)%o%apply(ant_avg_cons[-1,],2,sd)
		Qtable["GHG",c("Q05","Q95"),"cons"] = c(1,1)%o%Qtable["GHG","mean","cons"] + qnorm(.95)*c(-1,1)%o%apply(ghg_avg_cons[-1,],2,sd)
		Qtable[ "OA",c("Q05","Q95"),"cons"] = c(1,1)%o%Qtable[ "OA","mean","cons"] + qnorm(.95)*c(-1,1)%o%apply( oa_avg_cons[-1,],2,sd)
	}
	Qtable["Obs",,"cons"] = Qtable["Obs",,"uncons"]

	require(plotrix)
	dx =.3
	ddx = .03
	if (is.null(ylim)) {
		ylim = range(Qtable)
	}
	plot(1,0,xlim=c(.5,nfc+.5),ylim=ylim,type="p",pch=NA,cex=1.2,cex.lab=1.1,ylab=NA,xlab=NA,xaxt="n",mgp=c(2.5,.7,0),font=2)
	axis(1,at=1:nfc,labels=forcings,font=2,font.lab=2,cex.lab=1.1,mgp=c(2.5,.7,0))
	yaxp = par("yaxp")
	yticks = seq( yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3] )
	abline(h=yticks,lty=3)
	if (temporal=="delta") {
		title(ylab="Temperature change (°C)",font.lab=2,cex.lab=1.2,mgp=c(4.5,1,0))
		periods_label=paste0(as.character(min(pres_period)),"-",#
									as.character(max(pres_period)),"  wrt  ",#
									as.character(min(ref_period)),"-",#
									as.character(max(ref_period)))
		title(ylab=periods_label,font.lab=2,cex.lab=.9,mgp=c(2.8,1,0))
	} else if (temporal=="trend") {
		title(ylab=paste0("Warming rate ",#
								as.character(min(pres_period)),"-",#
								as.character(max(pres_period))," (°C/cent.)"),#
				font.lab=2,cex.lab=1.2,mgp=c(3,1,0))
	}
	polygon(c(-dx/2,-dx/2,dx/2,dx/2)+1,#
				  c(0,Qtable[1,"mean","uncons"],Qtable[1,"mean","uncons"],0),#
				  col=colorf_half[1],border=NA)
	for (i in 2:nfc) {
		polygon(c(-dx-ddx,-dx-ddx,-ddx,-ddx)+i,#
				  c(0,Qtable[i,"mean","uncons"],Qtable[i,"mean","uncons"],0),#
				  col=colorf_half[i],border=NA)
		polygon(c(dx+ddx,dx+ddx,ddx,ddx)+i,#
				  c(0,Qtable[i,"mean","cons"],Qtable[i,"mean","cons"],0),#
				  col=colorf[i],border=NA)
	}
	abline(h=0)
	plotCI(x=1,y=Qtable[1,"mean","uncons"],li=Qtable[1,"Q05","uncons"],ui=Qtable[1,"Q95","uncons"],col="black",pch=NA,lwd=3,add=T)
	plotCI(x=2:nfc-dx/2-ddx,y=Qtable[2:nfc,"mean","uncons"],li=Qtable[2:nfc,"Q05","uncons"],ui=Qtable[2:nfc,"Q95","uncons"],col="grey50",pch=NA,lwd=3,add=T)
	plotCI(x=2:nfc+dx/2+ddx,y=Qtable[2:nfc,"mean","cons"],li=Qtable[2:nfc,"Q05","cons"],ui=Qtable[2:nfc,"Q95","cons"],col="black",pch=NA,lwd=3,add=T)

	return(Qtable)
}


# plot_x_dec_da():
#---------------
# Illustrate the decomposition of covariate X between nat/ant.
# Inputs: 
#		X_fit	: decomposition of X
#		ofile	: output pdf file
plot_x_dec_da = function(X,X_fit,ofile) {
	pdf(ofile)
	layout(matrix(c(1,2,3), 3, 1))
	year = as.numeric(dimnames(X_fit)$year)
	year_histssp = as.numeric(dimnames(X)$year)
	models = dimnames(X_fit)$model
	Nmod = length(models)
	# Add the "ant" component
	X_fit_ant = X_fit[,,"all",] - X_fit[,,"nat",]
	X_fit = abind(X_fit,X_fit_ant,along=3,use.dnns=T)
		# Set forc="ant"
		dnns_tmp = dimnames(X_fit)
		dnns_tmp$forc[4] = "ant"
		dimnames(X_fit) = dnns_tmp
	X_q95 = apply(X_fit[,-1,c("all","nat","ant"),],c(1,3,4),quantile,.95)
	X_q05 = apply(X_fit[,-1,c("all","nat","ant"),],c(1,3,4),quantile,.05)
	for (mod in models){
		# Panel1: X_all and X_full (all data))
		par(fig=c(0,1,.65,1),mgp=c(2.5,.7,0),cex=1,font=2,mar=c(1,4,2,1),font.axis=2,font.lab=2,tcl=-.4,las=1)
		plot( year_histssp, X[,mod], type="p",pch=16,xlab="",ylab="T (K)",cex=.4,main=toupper(mod),xaxt="n")
		axis(1,labels=F)
		polygon( c(year,rev(year)), c(X_q95[,"all",mod],rev(X_q05[,"all",mod])), col=rgb(1,0,0,.5), border=NA )
		lines(year,X_fit[,"be","all",mod],col="brown",lwd=1)
		mtext("ALL",adj=.02,line=-1.2)
		
		# Panel 2: X_ant
		par(fig=c(0,1,.38,.65),new=T,mar=c(1,4,0,1))
		plot( year, X_fit[,"be","ant",mod], type="l",col="forestgreen",xlab="",ylab="T (K)",font=2,font.lab=2, ylim=range(c(X_q05[,"ant",mod],X_q95[,"ant",mod])),xaxt="n")  #xaxs="i",
		axis(1,labels=F)
		polygon( c(year,rev(year)), c(X_q95[,"ant",mod],rev(X_q05[,"ant",mod])), col=rgb(0,1,0,.5), border=NA )
		abline(h=0)
		mtext("ANT",adj=.02,line=-1.2)

		# Panel 3: X_nat
		par(fig=c(0,1,0,.38),new=T,mar=c(4,4,0,1))
		plot( year_histssp, X[,mod], type="p",pch=16,xlab="",ylab="T (K)",font=2,cex=.4,font.lab=2)
		title(xlab="Years",mgp=c(2,.7,0))
		polygon( c(year,rev(year)), c(X_q95[,"nat",mod],rev(X_q05[,"nat",mod])), col=rgb(0,0,1,.5), border=NA )
		lines( year, X_fit[,"be","nat",mod], col="blue")
		mtext("NAT",adj=.02,line=-1.2)
	}
	dev.off()
}


# plot_histghg_fit()
#--------------------
# Illustrate the fit of smoothing splines to the histghg data
plot_histghg_fit = function(X_full,X_fit,ofile) {
	pdf(ofile,height=5)
	year = as.numeric(dimnames(X_full)$year)
	ny = length(year)
	year_abs = as.numeric(dimnames(X_fit)$year)
	Models_full = dimnames(X_full)$model
	Models_abs = dimnames(X_fit)$model
	if (sum(Models_full!=Models_abs)) {
		message("Error in plot_x_ghg(): inconsistent sets of models")
	}
	Models = Models_abs
	Nmod = length(Models)
	dim_X_fit = length(nnames(X_fit))
	if (dim_X_fit>=3) {
		X_q95 = apply(X_fit[as.character(year),-1,],c(1,3),quantile,.95)
		X_q05 = apply(X_fit[as.character(year),-1,],c(1,3),quantile,.05)
	}
	
	for (imod in 1:Nmod){
		par(mgp=c(3.5,.7,0),cex=1,font=2,mar=c(4,5,2,1),font.axis=2,font.lab=2,tcl=-.4,las=1)
		if (dim_X_fit>=3) {
			plot( year, X_full[,imod], type="p", pch=16, xlab="", ylab="T (K)", cex=.4, main=toupper(Models[imod]),#
				   ylim=range(X_q95[,imod],X_q05[,imod],X_full[,imod]))
			title(xlab="Year",font.lab=2,mgp=c(2.2,.7,0))
			polygon(c(year,rev(year)),#
					  c(X_q95[,imod],#
						 rev(X_q05[,imod])),#
					  col=rgb(1,0,0,.5), border=NA )
			lines(year,X_fit[as.character(year),"be",imod],col="brown",lwd=1)
		} else {
			plot( year, X_full[,imod], type="p", pch=16, xlab="", ylab="T (K)", cex=.4, main=toupper(Models[imod]),#
				   ylim=range(X_full[,imod],na.rm=T))
			title(xlab="Year",font.lab=2,mgp=c(2.2,.7,0))
			lines(year,X_fit[as.character(year),imod],col="brown",lwd=1)
		}
	}
	dev.off()
}


# scatter_mod()
#---------------
scatter_mod = function(X,ofile=NULL,col=NULL,legend=NULL,fi=1.2,ab=NULL,xlim=NULL,ylim=NULL,...) {
	if (!is.null(ofile)) { 
		pdf(ofile)
	}
	par(mar=c(4,4,1,1),mgp=c(2.7,.7,0))
	if (is.null(xlim) | is.null(ylim)) {
		dlim = apply(X,1,range)
	} else {
		dlim = cbind(xlim,ylim)
	}
	fi = 1.2	# Factor increase for xlim/ylim
	dlim = matrix(c(fi,1-fi,1-fi,fi),2,2) %*% dlim
	plot( X[1,], X[2,],#
		  col=col,lwd=2,pch=3,cex.axis=1.2,font.axis=2,font.lab=2,#
		  #xlab=paste0(toupper(dimnames(dT2010)[[1]][2])," warming (2010-2020)"),#
		  #ylab=paste0(toupper(dimnames(dT2010)[[1]][2])," warming (2010-2020)"),#
		  xlim=dlim[,1],ylim=dlim[,2],...)
	if (!is.null(ab)) {
		abline(a=ab[1],b=ab[2],lwd=1,lty=2)
	}
	text( X[1,], X[2,], labels=dimnames(X)$model, cex= 0.7, font=2, pos=3,col=col)
	if (!is.null(legend)) {
		if (length(legend)==length(unique(col))) {
			xylim = par("usr")
			legend(x = .97*xylim[1] + .03*xylim[2], y = .03*xylim[3] + .97*xylim[4], legend=legend,col=unique(col),lwd=2,lty=1)
		}
	}
	
	if (!is.null(ofile)) {
		dev.off()
	}
}


# plot_fit_mar1_cmip()
#----------------------
plot_fit_mar1_cmip = function(file_Theta_MAR1, ofile, ny_max_plot=1000, lag_max_acf=50) {
# ny_max_plot: Plot of time-series shows 'ny_max_plot' years of pictl
# lag_max_acf: Number of lags shown in ACF

	load(file_Theta_MAR1)	# Theta_MAR1, X_pictl_treated
	Models = dimnames(Theta_MAR1)$model
	nmod = length(Models)
	
	# Define lag_max, ny_max_plot
	lags = 0:lag_max_acf
	if (is.na(ny_max_plot)) {
		my_max_plot = max(apply(!is.na(X_pictl_treated),2,sum))
	}

	pdf(ofile,width=10,height=5)
	par(cex=1.2,font=2,font.axis=2,font.lab=2,cex.axis=1,mgp=c(1.8,.4,0),lwd=2,tcl=-.4)
	
	for (mod in Models) {
		#message(mod)
		y_ctl_raw = X_pictl_treated[,mod]
		y_ctl = y_ctl_raw[!is.na(y_ctl_raw)]

		# time-series
		par(fig=c(.3,1,0,1),mar=c(3.5,3,4,1))
		y_ctl_plot = c(NA,y_ctl[1:ny_max_plot])
		plot(0:ny_max_plot,	y_ctl_plot,type="l",col="red",xlab="Year",ylab="GMST (K)",main=mod,xaxs="i",cex.lab=.9)
		
		# acf
		y_ctl_acf = y_ctl[1:min(length(y_ctl),ny_max_plot)]
		acf_ctl = acf(y_ctl_acf,lag.max=lag_max_acf,plot=F)
		par(fig=c(0,.3,0,1),mar=c(3.5,3,4,1),new=T)
		plot(acf_ctl,xlab="Lag",ylab="ACF",main="",cex.lab=.9,ci.col=NA)

		# acf fitted by MAR
		theta_mod = Theta_MAR1[mod,]
		acf_mar1 = theta_mod["var_ar1"]*theta_mod["alpha_ar1"]^lags + #
						theta_mod["var_wn"]*c(1,rep(0,max(lags)))
		lines(lags,acf_mar1/acf_mar1[1],col="red")
	}
	dev.off()
}


# plot_fit_mar2_cmip()
#----------------------
plot_fit_mar2_cmip = function(file_Theta_MAR2, ofile, ny_max_plot=1000, lag_max_acf=50) {
# ny_max_plot: Plot of time-series shows 'ny_max_plot' years of pictl
# lag_max_acf: Number of lags shown in ACF

	load(file_Theta_MAR2)	# Theta_MAR1, X_pictl_treated
	Models = dimnames(Theta_MAR2)$model
	nmod = length(Models)
	
	# Define lag_max, ny_max_plot
	lags = 0:lag_max_acf
	if (is.na(ny_max_plot)) {
		my_max_plot = max(apply(!is.na(X_pictl_treated),2,sum))
	}

	pdf(ofile,width=10,height=5)
	par(cex=1.2,font=2,font.axis=2,font.lab=2,cex.axis=1,mgp=c(1.8,.4,0),lwd=2,tcl=-.4)

	for (mod in Models) {
		#message(mod)
		y_ctl_raw = X_pictl_treated[,mod]
		y_ctl = y_ctl_raw[!is.na(y_ctl_raw)]

		# time-series
		par(fig=c(.3,1,0,1),mar=c(3.5,3,4,1))
		y_ctl_plot = c(NA,y_ctl[1:ny_max_plot])
		plot(0:ny_max_plot,y_ctl_plot,type="l",col="red",xlab="Year",ylab="GMST (K)",main=mod,xaxs="i",cex.lab=.9)
		
		# acf
		y_ctl_acf = y_ctl[1:min(length(y_ctl),ny_max_plot)]
		acf_ctl = acf(y_ctl_acf,lag.max=lag_max_acf,plot=F)
		par(fig=c(0,.3,0,1),mar=c(3.5,3,4,1),new=T)
		plot(acf_ctl,xlab="Lag",ylab="ACF",main="",cex.lab=.9,ci.col=NA)

		# acf fitted by MAR
		theta_mod = Theta_MAR2[mod,]
		acf_mar2 = theta_mod["var_ar1_fast"]*theta_mod["alpha_ar1_fast"]^lags + #
					  theta_mod["var_ar1_slow"]*theta_mod["alpha_ar1_slow"]^lags
		lines(lags,acf_mar2/acf_mar2[1],col="red")
	}
	dev.off()
}



# plot_completion()
#-------------------
plot_completion = function(X,X_fit,ofile) {
	Models_compl = dimnames(X)$model
	Models_compl2 = dimnames(X_fit)$model
	if (!identical(Models_compl,Models_compl2)) {
		message("Warning in plot_completion.R: inconsistent inputs")
	}
	year = as.numeric(dimnames(X)$year)
	
	pdf(ofile,height=5)
	par(mar=c(4,4,3,1),mgp=c(2.7,.7,0))
	X_q95_compl = apply(X_fit[,-1,],c(1,3),quantile,.95)
	X_q05_compl = apply(X_fit[,-1,],c(1,3),quantile,.05)
	for (mod in Models_compl) {
		v = X[,mod]
		v_restrict = v[!is.na(v)]
		year_restrict = as.numeric(names(v_restrict))
		ny_restrict = length(year_restrict)
	
		plot(year_restrict,v_restrict,xlim=range(year),pch=16,font=2,font.axis=2,font.lab=2,xlab="Year",ylab="Temperature   (°C)",main=toupper(mod),ylim=range(v_restrict,X_q05_compl[,mod],X_q95_compl[,mod]))
		lines(year,X_fit[,"be",mod],col="red",lwd=2)
		polygon(c(year,rev(year)),#
				  c(	  X_q95_compl[,mod],#
					 rev(X_q05_compl[,mod])),#
				  col=rgb(1,0,0,.5), border=NA )
	}
	dev.off()
}


plot_Qtable_proj = function(SQtable, ofile, col=NULL, ref_period=1870:1900) {
	# Input SQtables : an array of dimensions:
	#		quantile, tables, dates
	tables = dimnames(SQtable)$table
	dates = dimnames(SQtable)$date
	nt = length(tables)
	nd = length(dates)
	col = c("darkorange1","red2","violetred","coral2")

	pdf(ofile)
	xlim = c(.5,nd+.5)
	dx = .4/(nt-1)

	ylim = range(SQtable)
	par(mar=c(4,5,1,1))
	plot(1:nd-.2,SQtable["mean",1,],pch=18,cex=.1,xaxt="n",xlim=xlim,ylim=ylim,ylab="Temperature change",xlab="",font=2,font.axis=2,font.lab=2,cex.axis=1.2,cex.lab=1.2,mgp=c(3.7,.7,0))
	title(ylab=paste0("(°C wrt ",as.character(min(ref_period)),"-",as.character(max(ref_period)),")"),font.lab=2,cex.lab=.9,mgp=c(2.2,1,0))
	axis(1,at=1:nd,labels=dates,font.lab=2,font.axis=2,cex.axis=1.2,cex.lab=1.2,mgp=c(2.7,.7,0))
	abline(h=1:10,lty=3)
	for (itab in 1:nt) {
		tab=tables[itab]
		for (idate in 1:nd) {
			lines(idate-.2+(itab-1)*dx*rep(1,2), SQtable[c("Q05","Q95"),tab,idate], lwd=4, col=col[itab])
			points(idate-.2+(itab-1)*dx, SQtable["mean",tab,idate], pch=18, cex=2,col=col[itab])
		}
	}
	lim = par("usr")
	legend("topleft",legend=tables,col=col,lwd=4,lty=1,bg="white")

	dev.off()
}





plot_ic_proj = function(XC, dates, ofile, ref_period=1870:1900, method="variance") {
# XC : an array of dimnames: (year, sample, table)

	nd = length(dates)
	tables = dimnames(XC)$table
	nt = length(tables)
	sample = dimnames(XC)$sample
	ns = length(sample)
	XC_tmp = array(NA,dim=c(ns,nt,nd),dimnames=list(sample=sample,tables=tables, dates=dates))
	XC_ic = array(NA,dim=c(4,nt,nd),dimnames=list(quantile=c("Q05","mean","Q95","ic_half_width"),tables=tables, dates=dates))

	avg_period = function(X,ref_period,futur_period) {
		year = as.numeric(dimnames(X)$year)
		ndim = length(dim(X))
		if (ndim==2) {
			avg_ref = apply(X[year%in%ref_period,,drop=F]    ,2:ndim,mean)
			avg_fut = apply(X[year%in%futur_period,,drop=F]  ,2:ndim,mean)
		} else if (ndim==3) {
			avg_ref = apply(X[year%in%ref_period,,,drop=F]   ,2:ndim,mean)
			avg_fut = apply(X[year%in%futur_period,,,drop=F] ,2:ndim,mean)
		} else if (ndim==4) {
			avg_ref = apply(X[year%in%ref_period,,,,drop=F]  ,2:ndim,mean)
		avg_fut = apply(X[year%in%futur_period,,,,drop=F],2:ndim,mean)
		} else {
			message("Error in avg_period(): ndim > 4")
			return
		}
		return(avg_fut-avg_ref)
	}


	for (date_str in dates) {
		yr_date = strsplit(date_str,"-")[[1]]
		nyr_date = length(yr_date)
		if (nyr_date==1) {
			fut_per = as.numeric(yr_date)
		} else if (nyr_date==2) {
			fut_per = as.numeric(yr_date[1]):as.numeric(yr_date[2])
		} else {
			message("Error in plot_ic_proj.R: invalid date"); stop()
		}

		XC_tmp[,,date_str] = avg_period(XC,ref_period,futur_period=fut_per)
	}
	
	XC_ic["mean",,] = XC_tmp["be",,]
	if (method=="quantile") {
		XC_ic[c("Q05","Q95"),,] = apply(XC_tmp[-1,,],c(2,3),quantile,c(.05,.95))
	} else if (method=="variance") {
		sd_tmp = apply(XC_tmp[-1,,,drop=F],c(2,3),sd)
		XC_ic[c("Q05","Q95"),,] = c(1)%o%XC_ic["mean",,] + qnorm(.95)* c(-1,1)%o%sd_tmp
	}
	XC_ic["ic_half_width",,] = (XC_ic["Q95",,] - XC_ic["Q05",,]) / 2

	plot_Qtable_proj(XC_ic,ofile)
	return(XC_ic)
}


#0970170554
#CN2308
