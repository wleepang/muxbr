source('muxbrFunc.r')

# initial parameter guess
p0 = c(1, 0.1)

# allocate a matrix to store the parameters in
P = matrix(0, nrow=2, ncol=16)
colnames(P) = sprintf('FC%02d', 1:16)
rownames(P) = sprintf('p%d', 2:3)

# this matrix is for static parameters - e.g. baseline value and scale factor
Q = matrix(0, nrow=2, ncol=16)
colnames(Q) = colnames(P)
rownames(Q) = c('p0', 'p1')

# allocate a matrix to store the fit statistics
FS= P
rownames(FS)= c('sse', 'r2')

ODv = t(as.matrix(read.delim('./20110209/od600_verified.txt', row.names=1, header=F)))
ODv = apply(ODv[,-1], 2, as.numeric)

# only the first 5 verified readings are good to calibrate over since there is no
# copper addition yet.
ODv = list(t=ODv[,'at'], y=ODv[,-1])

cpts = matrix(0, nrow=2, ncol=16)
cpts[1,] = 2
cpts[2,1:4] = dim(ODv$y)[1]
cpts[2, c(5:7, 11:13)] = 5
cpts[2, c(8:10,14:16)] = dim(ODv$y)[1] - 1

# match the ODv times with the closest PD measurement times
PD = getMuxBRData('./20110209/photodiode.txt')
PD = despike(PD, bl.pts = 100)
PD = mbdSmooth(PD)

sf = 1e5
os = 0.005
Q[1,] = 0
Q[2,] = sf

PD$y = PD$y*sf+os

# optimization time
for (i in 1:16) {
    
    jx = cpts[1,i]:cpts[2,i]
    ix = findInterval(ODv$t[jx], PD$t)
    
    od = ODv$y[jx,i]
    pd = PD$y[ix,i]
    
    res = nlminb(p0, fitSSE, OD=od, AMP=pd, lower=c(1e-5,1e-5))
    P[,i]       = res$par
    FS['sse',i] = res$objective
    FS['r2' ,i] = cor(od, pd)
}

# make verification plots
dev.new()
split.screen(c(4,4))
for (i in 1:16) {
    od = ODv$y[,i]
    pd = PD$y[findInterval(ODv$t, PD$t),i]
    p  = P[,i]

    screen(i)
    par(mar=c(1,1,0,0)*2+.5)
    plot(pd, od, ann=F, xlim=c(0,1), ylim=c(0,1))
    
    pdrng = seq(0,1,0.05)
    lines(pdrng, fitAmpToOD(p, pdrng), col='red')
    text(par('xaxp')[1]+.05, par('yaxp')[2]-.05, i)
}

# send results to file
res = rbind(Q, P, FS)
print(res)

write.table(res, file='calibrationCoeffs.tsv', sep='\t', quote=FALSE, col.names=NA)