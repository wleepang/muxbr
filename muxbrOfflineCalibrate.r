source('muxbrFunc.r')

# list of calibration files
calPath = './20110209/cal/'
calFiles = dir(calPath, pattern='.?Cu-2p00mM.?')
calOutput= sprintf('%s%s', calPath, 'calibrationCoeffs_Cu-2p00mM.tsv')
calODs = c(0, 0.05, 0.1, 0.2, 0.5, 0.8, 1.0)

od = NULL
pd = NULL
for (i in 1:length(calFiles)) {
    cf = read.delim(sprintf('%s%s',calPath,calFiles[i]), header=F, as.is=T)
    
    # ignore the timestamp column
    cf[,1] = NULL
    od = c(od, rep(calODs[i], each=dim(cf)[1]))
    
    pd = rbind(pd, as.matrix(cf))
    
    if (i == 1) {
        nfc = dim(cf)[2]
        mu = matrix(0, nrow=length(calFiles), ncol=nfc)
        st = mu
    }
    mu[i,] = colMeans(cf)
    st[i,] = sd(cf)
}
colnames(pd) = sprintf('FC%02d', 1:dim(pd)[2])
colnames(mu) = colnames(pd)
colnames(st) = colnames(pd)

# initial parameter guess
p0 = c(1, 0.1)

# allocate a matrix to store the parameters in
P = matrix(0, nrow=2, ncol=nfc)
colnames(P) = sprintf('FC%02d', 1:nfc)
rownames(P) = sprintf('p%d', 2:3)

# this matrix is for static parameters - e.g. baseline value and scale factor
Q = P
rownames(Q) = c('p0', 'p1')

Q[1,] = mu[1,]
Q[2,] = 1e5

# allocate a matrix to store the fit statistics
FS= P
rownames(FS)= c('sse', 'r2')

# baseline and scale
pd = t(t(pd) - Q[1,])
pd = t(t(pd) * Q[2,])

# use non-linear model minimization to determine parameter fits
# loop over flow cells
for (i in 1:nfc) {
    res = nlminb(p0, fitSSE, OD=od, AMP=pd[,i], lower=c(1e-5,1e-5))
    P[,i]       = res$par
    FS['sse',i] = res$objective
    FS['r2' ,i] = cor(od, pd[,i])
}

# make verification plots
dev.new()
split.screen(c(4,4))
for (i in 1:nfc) {
    p  = P[,i]

    screen(i)
    par(mar=c(1,1,0,0)*2+.5)
    plot(pd[,i], od, ann=F, xlim=c(0,1.5), ylim=c(0,1.5))
    
    pdrng = seq(0,1.5,0.05)
    lines(pdrng, fitAmpToOD(p, pdrng), col='red')
    text(par('xaxp')[1]+.05, par('yaxp')[2]-.05, i)
}

# send results to file
res = rbind(Q, P, FS)
print(res)

write.table(res, file=calOutput, sep='\t', quote=FALSE, col.names=NA)