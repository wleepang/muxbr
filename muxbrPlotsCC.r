source('muxbrFunc.r')
cc000 = getMuxBRODCal('./20110209/cal/calibrationCoeffs_Cu-0p00mM.tsv')
cc085 = getMuxBRODCal('./20110209/cal/calibrationCoeffs_Cu-0p85mM.tsv')
cc200 = getMuxBRODCal('./20110209/cal/calibrationCoeffs_Cu-2p00mM.tsv')
pd = getMuxBRData('./20110209/photodiode.txt')
ev = getMuxBREvents('./20110209/notes.txt')

od = pd

od.v = t(as.matrix(read.delim('./20110209/od600_verified.txt', header=F, row.names=1)))
od.v = apply(od.v[,-1], 2, as.numeric)
od.v = list(t=od.v[,1], y=od.v[,-1])

# use appropriate calibration constants for segments of the curve
od$y[,1:4] = ampToOD(pd$y[,1:4], cc000)

# time indices for 1st and 2nd doses
tx1st = which(od$t >= ev$events$dosing$t[1])
tx2nd = which(od$t >= ev$events$dosing$t[2])

# flask ids for 1st and 2nd dose sets and concentration sets
id1st085 = 11:13
id1st200 =  5:7
id2nd085 = 14:16
id2nd200 =  8:10

od$y[-tx1st,c(id1st085, id1st200)] = ampToOD(pd$y[-tx1st,c(id1st085, id1st200)], cc000[,c(id1st085, id1st200)])
od$y[ tx1st, id1st085] = ampToOD(pd$y[tx1st,id1st085], cc085[,id1st085])
od$y[ tx1st, id1st200] = ampToOD(pd$y[tx1st,id1st200], cc200[,id1st200])

od$y[-tx2nd,c(id2nd085, id2nd200)] = ampToOD(pd$y[-tx2nd,c(id2nd085, id2nd200)], cc000[,c(id2nd085, id2nd200)])
od$y[ tx2nd, id2nd085] = ampToOD(pd$y[tx2nd,id2nd085], cc085[,id2nd085])
od$y[ tx2nd, id2nd200] = ampToOD(pd$y[tx2nd,id2nd200], cc200[,id2nd200])

od.ds = od #despike(od, bl.pts=1)
od.ds.s = mbdSmooth(od.ds)

od.ds.s$y[which(od.ds.s$y < 0.01)] = NA

mu = growthRate(od.ds.s)

clrs = c(	rep('#000000', 4),
			rep('#0000FF', 6),
			rep('#8888FF', 6))
lwds = c(	rep(1, 2),
            rep(2, 2),
			rep(1, 3),
			rep(2, 3), 
			rep(1, 3),
			rep(2, 3))

windows(16,8)
#png(filename='muxbrMonitor.png', width=1024, height=1024)
    tx = which(od.s$t >= 2000)
    t.sc = 60
    yl.od = c(0.01, 2)
    yl.mu = c(-0.001, 0.005)
    gcf = dev.cur
    
    par(mar=c(2,4,2,0.1))
	par(mfcol=c(2,3), mfg=c(1,1))
        # OD600         Ctrls   0.85mM  2mM
        # Growth Rate   
        ix = 1:4
		plot.muxbr(od.ds.s, 
            t.ix=tx, y.ix=ix, dev=gcf(), log='y', t.scale = t.sc, ylim=yl.od, clrs=clrs[ix], lwd=lwds[ix])
		matpoints(od.v$t/t.sc, od.v$y[,ix], log='y', ylim=yl.od, col=clrs[ix], pch=1)
		plot.ev(ev$events$dosing, t.scale = t.sc, lty='solid', col='#cc0000', lwd=1)
		title(main='Controls', ylab='OD600')
	par(mfg=c(1,2))
        ix = 11:16
		plot.muxbr(od.ds.s, 
            t.ix=tx, y.ix=ix, dev=gcf(), log='y', t.scale = t.sc, ylim=yl.od, clrs=clrs[ix], lwd=lwds[ix])
		matpoints(od.v$t/t.sc, od.v$y[,ix], log='y', ylim=yl.od, col=clrs[ix], pch=1)
		plot.ev(ev$events$dosing, t.scale = t.sc, lty='solid', col='#cc0000', lwd=1)
		title(main='Dosed 0.85mM')
	par(mfg=c(1,3))
        ix = 5:10
		plot.muxbr(od.ds.s, 
            t.ix=tx, y.ix=ix, dev=gcf(), log='y', t.scale = t.sc, ylim=yl.od, clrs=clrs[ix], lwd=lwds[ix])
		matpoints(od.v$t/t.sc, od.v$y[,ix], log='y', ylim=yl.od, col=clrs[ix], pch=1)
		plot.ev(ev$events$dosing, t.scale = t.sc, lty='solid', col='#cc0000', lwd=1)
		title(main='Dosed 2.00mM')
	par(mfg=c(2,1))
		ix = 1:4
		plot.muxbr(mu, t.ix=tx-1, y.ix=ix, dev=gcf(), t.scale = t.sc, ylim=yl.mu, clrs=clrs[ix], lwd=lwds[ix])
		plot.ev(ev$events$dosing, t.scale = t.sc, lty='solid', col='#cc0000', lwd=1)
		lines(par('usr')[1:2], c(0,0), lty='dotted', col='#888888')
		title(ylab='growth rate (per min)')
	par(mfg=c(2,2))
		ix = 11:16
		plot.muxbr(mu, t.ix=tx-1, y.ix=ix, dev=gcf(), t.scale = t.sc, ylim=yl.mu, clrs=clrs[ix], lwd=lwds[ix])
		plot.ev(ev$events$dosing, t.scale = t.sc, lty='solid', col='#cc0000', lwd=1)
		lines(par('usr')[1:2], c(0,0), lty='dotted', col='#888888')
    par(mfg=c(2,3))
        ix = 5:10
		plot.muxbr(mu, t.ix=tx-1, y.ix=ix, dev=gcf(), t.scale = t.sc, ylim=yl.mu, clrs=clrs[ix], lwd=lwds[ix])
		plot.ev(ev$events$dosing, t.scale = t.sc, lty='solid', col='#cc0000', lwd=1)
		lines(par('usr')[1:2], c(0,0), lty='dotted', col='#888888')
#dev.off()

#vod = t(apply(as.matrix(read.delim('od600_verified.txt', header=F, as.is=T))[2:17,2:6], 2, as.numeric))
#colnames(vod) = c('AbsTime', as.character(seq(1,15,1)))

#mu = list(t=od$t[-1], y=t(t(diff(apply(od$ns, 2, log))) / diff(od$t)))
#mu$y = matrix(sapply(mu, function(el) {if (is.nan(el)) {0} else {el}}), nrow=dim(mu)[1])
#mu$s = apply(mu, 2, smoother, x=od$t[-1])
