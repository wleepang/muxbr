plot.od = function(od, ...) {
    plot.muxbr(od, ylim=c(0.01,2), ...)
    title(ylab='OD600')
}
plot.pd = function(pd, ...) {
    plot.muxbr(pd, ...)
    title(ylab='Amps')
}
plot.tc = function(tc, ...) {
    plot.muxbr(tc, ylim=c(35,40), ...)
    title(ylab='Temperature (degC)')
}

plot.ev = function(ev, t.scale = 1, ...) {
	xlyl = par('usr')
	pad = 0.03
	xl = xlyl[1:2]
	if (par('ylog')) {
		yl = 10^xlyl[3:4]
        yt = 10^(xlyl[4]-diff(xlyl)[3:4]*pad)
	} else {
		yl = xlyl[3:4]
        yt = yl[2]-diff(yl)*pad
	}
	
	for (i in 1:length(ev$t)) {
		lines(rep(ev$t[i], 2)/t.scale, yl, ...)
		text((ev$t[i])/t.scale+diff(xl)*pad, yt, i, ...)
	}
}

plot.muxbr = function(mxd, type='l', ltys='solid', dev=NULL, clrs = NULL, t.scale = 1, t.ix=NULL, y.ix=NULL, ...) {
    if (is.null(dev)) {
        dev.new()
    } else {
        dev.set(dev)
    }
    
    if (!is.null(t.ix)) {
        mxd$t = mxd$t[t.ix]
        mxd$y = mxd$y[t.ix,]
    }
    
    if (!is.null(y.ix)) {
        mxd$y = mxd$y[,y.ix]
    }
    
    if (is.null(clrs)) {
        clrs = rainbow(dim(mxd$y)[2])
    }
    matplot(mxd$t / t.scale, mxd$y, type=type, lty=ltys, col=clrs, ann=F, ...)
    title(xlab='Time')
}

smoother = function(x=NULL, y, f=0.08) {
    if (!is.null(x)) {
        return(lowess(x, y=y, f=f)[['y']])
    } else {
        return(lowess(y, f=f)[['y']])
    }
}

mbdSmooth = function(mbd, f=0.08) {
    mbd.s = mbd
    mbd.s$y = apply(mbd$y, 2, smoother, x=mbd$t, f=f)
    return(mbd.s)
}

plot.flask = function(id, logax='', dev=dev.cur(), evt=NULL) {
    if (!is.null(dev)) {
        dev.set(dev)
    } else {
        dev.new()
    }
	od.yl = c(0.01,3)
    par(mfcol=c(2,1))
    par(mfg=c(1,1))
    plot(od$t / 60, od$y[,id], pch=20, col='#cccccc', ylim=od.yl, log=logax, ann=F)
	lines(od$t / 60, od$s[,id], type='l')
	#points(vod[,1] / 60, vod[,id+1], pch=1)
    if (!is.null(evt)) {
        for (i in 1:length(evt)) {
            lines(c(evt[i],evt[i])/60, od.yl, lty='dotdash', col='#cc0000')
        }
    }
	title(main=sprintf('Flask %d', id), ylab='OD600')
    
    par(mfg=c(2,1))
    mu.yl = range(0,0.0025)
    plot(od$t[-1] / 60, mu[,id], pch=20, col='#cccccc', ylim=mu.yl, ann=F)
	lines(od$t[-1] / 60, mu$s[,id], type='l')
    lines(c(1380,1380)/60, mu.yl, lty='dotdash', col='#cc0000')
	title(xlab='Time (hr)')
}

getMuxBRData <- function(fname) {
    dfrm = read.delim(fname, header=F, as.is=T)
    tvec = strptime(dfrm[[1]], '%m/%d/%Y %l:%M:%S %p')
    tabs = as.numeric(difftime(tvec, tvec[1]), units='mins')
    dmat = as.matrix(dfrm[,-1])
    
    return(list(t=tabs, y=dmat))
}

getMuxBRODCal <- function(fname) {
	dfrm = read.delim(fname, comment.char='#', row.names='param')
	return(dfrm)
}

getMuxBREvents <- function(fname) {
	dfrm = read.delim(fname, header=F, sep=',', quote='"', as.is=T)
	tvec = strptime(dfrm[[1]], '%m/%d/%Y %H:%M')
	tabs = as.numeric(difftime(tvec, tvec[1]), units='mins')
	
	# trim the strings of prepended spaces
	dfrm[[2]] = sub(' ', '', dfrm[[2]])
	dfrm[[3]] = sub(' ', '', dfrm[[3]])
	
	eventType = unique(dfrm[[2]])
	
	events = vector('list', length(eventType))
	names(events) = eventType
	for (i in 1:length(eventType)) {
		ix = which(dfrm[[2]] == eventType[i])
		events[[eventType[i]]] = list(t=tabs[ix], desc=dfrm[[3]][ix])
	}
	
	return(list(eventTypes=eventType, events=events))
}

ampToOD <- function(amp, cc, rep.bad.fit=FALSE) {
	# calculate default parameter values for flow cells that do not have good
	# parameter fit statistics
	param.rows = c(-dim(cc)[1], -(dim(cc)[1]-1))
	good.fits = which(cc['r2',] >= 0.99)
	dflts = rowMeans(cc[param.rows,good.fits])
	
	# od = (amp - p0) / ((p2 / (p1 * p3)) - p3*(amp - p0))
	od = amp
	for (i in 1:dim(amp)[2]) {
		if (rep.bad.fit && is.na(match(i, good.fits))) {
			p = c(dflts['p0'], dflts['p1'], dflts['p2'], dflts['p3'])
		} else {
			p = c(cc['p0',i], cc['p1',i], cc['p2',i], cc['p3',i])
		}
		od[,i] = (amp[,i] - p[1]) / ((p[3] / (p[2] * p[4])) - p[4]*(amp[,i] - p[1]))
	}
		
	return(od)
}

despike <- function(mbd, bl.pts = 50, sdn = 2.5, bGraphic = FALSE) {
# removes large shifts in the data based on the assumption that they
# are additive artifacts

    # subtract the data baseline based on the mean of a set of
    # initial points
    mu = apply(head(mbd$y, bl.pts), 2, mean)
    mbd.b = mbd
    mbd.b$y = t(t(mbd$y) - mu)
    
    # calculate deviation statistics
    dmu = mean(diff(mbd.b$y))
    dsd = sd(as.vector(diff(mbd.b$y)))

    # determine the deviations
    delta = diff(mbd.b$y)-dmu
    
    if (bGraphic) {
        # diagnostic plot to show the deviations and the calculated
        # exclusion threshold
        matplot(mbd.b$t[-1], delta, type='l')

        xlim = par('usr')[1:2]
        lines(xlim,rep(+dsd*sdn,2), lty='dashed', col='red')
        lines(xlim,rep(-dsd*sdn,2), lty='dashed', col='red')
    }
    
    # set all deviation values of value under the threshold to zero
    delta[which(abs(delta) <= dsd*sdn)] = 0
    
    # calculate cummulative deviation over the trends
    delta = apply(delta, 2, cumsum)

    # subtract the cummulative deviations from the original trend
    # to correct the shifts
    mbd.ds = mbd.b
    mbd.ds$y[-1,] = mbd.b$y[-1,]-delta
    
    # rebaseline the data
    mu.ds = apply(head(mbd.ds$y, bl.pts), 2, mean)
    mbd.ds$y = t(t(mbd.ds$y) - mu.ds)
    
    return(mbd.ds)
}

growthRate = function(od) {
    # specific exponential growth rate:
    # rate of change of log transformed data
    mu = t(t(diff(apply(od$y, 2, log))) / diff(od$t))
    
    # there might be invalid data points (e.g. NaN values)
    # convert those to NA
    mu = matrix(sapply(mu, function(el) {if (is.nan(el)) {NA} else {el}}), nrow=dim(mu)[1])
    
    return(list(t=od$t[-1], y=mu))
}

fitAmpToOD = function(p, AMP) {
    return(AMP / ((p[1]/p[2]) - p[2] * AMP))
}

fitSSE = function(p, OD, AMP) {
    return(sum((OD - fitAmpToOD(p, AMP))^2))
}