plot_cons = function(X_krig, Xo, ref_plot=NULL, ny=251, ylim_da=NULL, color="red", title=NULL) {
    par(font.lab=2,font.axis=2,cex.lab=1.2,mar=c(4,4,1,1),mgp=c(2.5,.7,0))
    x = X_krig[1:ny,,"all","uncons"]
    x_cons = X_krig[1:ny,,"all","cons"]
    obs_x = apply(Xo,1,median)
		year = 1850:2100
		year_obs = 1870:2021
		
    # set the ref at the preindustrial level
    ref_plot = 1870:1900
    x = x - ones(x[,1]) %o% apply(x[year %in% ref_plot,], 2, mean)
    x_cons = x_cons - ones(x_cons[,1]) %o% apply(x_cons[year %in% ref_plot,], 2, mean)
    obs_x = obs_x - mean(obs_x[year_obs %in% ref_plot])

    # compute spread
    x_q95 = apply(x[,-1],1,quantile,.95)
    x_q05 = apply(x[,-1],1,quantile,.05)
    xc_q95 = apply(x_cons[,-1],1,quantile,.95)
    xc_q05 = apply(x_cons[,-1],1,quantile,.05)

    plot(year_obs, obs_x, xlim=range(year), ylim=ylim_da, type="p", pch=16, cex=.8, xlab="Year", ylab="Temperature (K)",panel.first=abline(v=NA,col="gray"))
    title(title, line = -2, cex.main = 2)
    yaxp = par("yaxp")
    yticks = seq( yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3] )
    abline(h=yticks,lty=3)
    
    if (color == "red") {
        polygon(c(year,rev(year)), c(x_q95,rev(x_q05)),border=NA,col=rgb(1,0,0,alpha=.2))
        polygon(c(year,rev(year)), c(xc_q95,rev(xc_q05)),border=NA,col=rgb(1,0,0,alpha=.5))
        colb = col2rgb("brown")/255
    } else if (color == "blue") {
        polygon(c(year,rev(year)), c(x_q95,rev(x_q05)),border=NA,col=rgb(0,0,1,alpha=.2))
        polygon(c(year,rev(year)), c(xc_q95,rev(xc_q05)),border=NA,col=rgb(0,0,1,alpha=.5))
        colb = col2rgb("blue")/255
    }
    
    lines(year,x[,1],lwd=1.5,col=rgb(colb[1],colb[2],colb[3],alpha=.5))
    lines(year,x_cons[,1],lwd=2,col=rgb(colb[1],colb[2],colb[3],alpha=1))
}

plot_cons_iter = function(X_krig_cons, X_krig_uncons, Xo, ref_plot=NULL, ny=251, ylim_da=NULL, color="red", title=NULL) {
    par(font.lab=2,font.axis=2,cex.lab=1.2,mar=c(4,4,1,1),mgp=c(2.5,.7,0))
    x = X_krig_uncons[1:ny,,"all","uncons"]
    x_cons = X_krig_cons[1:ny,,"all","cons"]
    obs_x = apply(Xo,1,median)
		year = 1850:2100
		year_obs = 1870:2021
		
    # set the ref at the preindustrial level
    ref_plot = 1870:1900
    x = x - ones(x[,1]) %o% apply(x[year %in% ref_plot,], 2, mean)
    x_cons = x_cons - ones(x_cons[,1]) %o% apply(x_cons[year %in% ref_plot,], 2, mean)
    obs_x = obs_x - mean(obs_x[year_obs %in% ref_plot])

    # compute spread
    x_q95 = apply(x[,-1],1,quantile,.95)
    x_q05 = apply(x[,-1],1,quantile,.05)
    xc_q95 = apply(x_cons[,-1],1,quantile,.95)
    xc_q05 = apply(x_cons[,-1],1,quantile,.05)

    plot(year_obs, obs_x, xlim=range(year), ylim=ylim_da, type="p", pch=16, cex=.8, xlab="Year", ylab="Temperature (K)",panel.first=abline(v=NA,col="gray"))
    title(title, line = -2, cex.main = 2)
    yaxp = par("yaxp")
    yticks = seq( yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3] )
    abline(h=yticks,lty=3)
    
    if (color == "red") {
        polygon(c(year,rev(year)), c(x_q95,rev(x_q05)),border=NA,col=rgb(1,0,0,alpha=.2))
        polygon(c(year,rev(year)), c(xc_q95,rev(xc_q05)),border=NA,col=rgb(1,0,0,alpha=.5))
        colb = col2rgb("brown")/255
    } else if (color == "blue") {
        polygon(c(year,rev(year)), c(x_q95,rev(x_q05)),border=NA,col=rgb(0,0,1,alpha=.2))
        polygon(c(year,rev(year)), c(xc_q95,rev(xc_q05)),border=NA,col=rgb(0,0,1,alpha=.5))
        colb = col2rgb("blue")/255
    }
    
    lines(year,x[,1],lwd=1.5,col=rgb(colb[1],colb[2],colb[3],alpha=.5))
    lines(year,x_cons[,1],lwd=2,col=rgb(colb[1],colb[2],colb[3],alpha=1))
}

