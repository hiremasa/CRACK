#!/usr/bin/Rscript

source("utilities.R")
source("test_synthetic_data.R")

## make crack applicable through R
## IMPORTANT: If you want to execute the tests for Crack_{\Delta} instead of Crack_{\delta} or NCI based Crack
## you have to set the '-s' parameter to '1'
crackWrapper = function(t, name="test", dir="temp", type="-t 'i'", dimX="1"){
    write.table(t, file=paste(dir, "current.txt", sep="/"), sep=" ", quote=F, row.names=F, col.names=F)
    system(paste(c("./crack.run -s 3 -o ", dir, "/crack_ -x ", dimX, " -c -d ' ' ", type, " -i ", dir, "/current.txt -a ", name), collapse=""))
    result_t = read.table(paste(dir, "crack_results.tab", sep="/"), header=F, stringsAsFactors=F, sep="\t")
    row = dim(result_t)[1]
    if(result_t[row, 1] != name){
        print("------> Output id not found -- some error occured <----------------")
        return(list(cd="--", eps=0, time=0))
    }else{
        return(list(cd=result_t[row,3], eps=result_t[row,2], time=result_t[row,4]))
    }
}

## Evaluate results on tuuebingen benchmark pairs
t = read.table(file="results/tuebingen_results.tab", sep="\t", header=F, stringsAsFactors=F)
corr = rep(0,106)
sum(t$V3 == "--")
corr[t$V3 == ref$V6] = 1
sum(corr)
sum(corr[uv])
sum(meta$V6[corr == 1]) / sum(meta$V6)
sum(meta$V6[uv[corr[uv] == 1]]) / sum(meta$V6[uv])
resCrack = data.frame(Correct=corr, Eps=t$V2, Cds=t$V3)
resCrack.uv = resCrack[uv,]

resCrack = resCrack.uv
resIGCI = read.table("results/igci_uv.tab", sep="\t", header=T)
resDC = read.table("results/dc_uv.tab", sep="\t", header=T)
resOrigo = read.table("results/origo_uv.tab", sep="\t", header=T)
resSlope = read.table("results/slope_uv.tab", sep="\t", header=T)
resCure = read.table("results/cure_uv.tab", sep="\t", header=T)
resErgo = read.table("results/ergo_uv.tab", sep="\t", header=T)

drC = decision_rate_meta(resCrack, uv=uv, oneHalf = F)
drI = decision_rate_meta(resIGCI, uv=uv, oneHalf = F)
drD = decision_rate_meta(resDC, uv=uv, oneHalf = F)
drO = decision_rate_meta(resOrigo, uv=uv, oneHalf = F)
drS = decision_rate_meta(resSlope, uv=uv, oneHalf = F)
drCu = decision_rate_meta(resCure, uv=uv, oneHalf = F)
drE = decision_rate_meta(resErgo, uv=uv, oneHalf = F)

pdf(width = 12, height=6, file="results/decision_rate_with_slope.pdf")
colors = c("black","red","darkgreen","darkblue", "gold")
plot(c(0,1), c(0,1), type="n", xlab="decision rate", ylab="% of correct decisions", bty="n", cex.lab=1.3, cex.axis=1.3)
lines(drC$S, drC$D, col=colors[1], pch=20, lty=1, lwd=2)
lines(drI$S, drI$D, col=colors[2], pch=20, lty=1, lwd=2)
lines(drD$S, drD$D, col=colors[3], pch=20, lty=1, lwd=2)
lines(drO$S, drO$D, col=colors[4], pch=20, lty=1, lwd=2)
lines(drS$S, drS$D, col=colors[5], pch=20, lty=1, lwd=2)
abline(h=0.95,col="darkgrey",lty=2)
abline(h=0.9,col="darkgrey",lty=2)
abline(h=0.85,col="darkgrey",lty=2)
abline(h=0.8,col="darkgrey",lty=2)

legend(0.8, 0.425, c("Crack", "IGCI", "DC", "Origo", "Slope"), cex=1.3, col=colors, lty=1, lwd=c(3,2,2), bty="n", y.intersp=1.2)
dev.off()

df = data.frame(CrackPos=drC$S, Crack=drC$D, IGCIPos=drI$S, IGCI=drI$D, DCPos=drD$S, DC=drD$D, OrigoPos=drO$S, Origo=drO$D, SlopePos=drS$S, Slope=drS$D, CurePos=drCu$S, Cure=drCu$D, ErgoPos=drE$S, Ergo=drE$D)
write.table(df, file="results/decision_rate_weighted.dat", row.names=F, quote=F, col.names = T)

######################
##Dependency Results##
######################

## test mixed data
dir.create("tempMixed")
corr = rep(0,21)
wrong = rep(0,21)
set.seed(1)
for(j in 0:20){
    cc = 0
    ww = 0
    dependency = j * 0.05
    for(i in 1:200){
        ttt = rep("b", 6)
        for(k in 1:6){
            tt = runif(1, min=0, max=4)
            if(tt <= 1.0){
                ttt[k] = "b"
            }else if(tt <= 2.0){
                ttt[k] = "c"
            }else{
                ttt[k] = "i"
            }
        }
        dd = generateSyntheticData(3,3,5000,dependency,0.0,ttt)
        dd = rbind(ttt, dd)
        id = round(runif(1, min=1, max=10000000))
        res = crackWrapper(dd, dir="tempMixed", name=id, type="", dimX = 3)
        if(res$cd == "->"){
            cc = cc + 1
        }else if(res$cd == "<-"){
            ww = ww + 1
        }
    }
    corr[(j+1)] = cc / 200
    wrong[(j+1)] = ww / 200
}
results = data.frame(Dependency=c(0:20*0.05), Correct=corr, Wrong=wrong, Indecisive=1-(corr+wrong))
write.table(results, file="results/crack_mixed.tab", quote=F, row.names=F)
print(results)

dir.create("tempNom")
corr = rep(0,21)
wrong = rep(0,21)
set.seed(1)
for(j in 0:20){
    cc = 0
    ww = 0
    dependency = j * 0.05
    for(i in 1:200){
        ttt = rep("b", 6)
        for(k in 1:6){
            tt = runif(1, min=0, max=1)
            if(tt <= 0.66){
                ttt[k] = "b"
            }else{
                ttt[k] = "c"
            }
        }
        dd = generateSyntheticData(3,3,5000,dependency,0.0,ttt)
        dd = rbind(ttt, dd)
        id = round(runif(1, min=1, max=10000000))
        res = crackWrapper(dd, dir="tempNom", name=id, type="", dimX = 3)
        if(res$cd == "->"){
            cc = cc + 1
        }else if(res$cd == "<-"){
            ww = ww + 1
        }
    }
    corr[(j+1)] = cc / 200
    wrong[(j+1)] = ww / 200
}
results = data.frame(Dependency=c(0:20*0.05), Correct=corr, Wrong=wrong, Indecisive=1-(corr+wrong))
write.table(results, file="results/crack_nominal.tab", quote=F, row.names=F)
print(results)

## test real valued only
dir.create("tempR")
corr = rep(0,21)
wrong = rep(0,21)
set.seed(1)
for(j in 0:20){
    cc = 0
    ww = 0
    dependency = j * 0.05
    for(i in 1:200){
        dd = generateSyntheticData(3,3,5000,dependency,0.0,c("r", "r", "r", "r", "r", "r"))
        id = round(runif(1, min=1, max=10000000))
        res = crackWrapper(dd, dir="tempR", name=id, type="-t 'i'", dimX = 3)
        if(res$cd == "->"){
            cc = cc + 1
        }else if(res$cd == "<-"){
            ww = ww + 1
        }
    }
    corr[(j+1)] = cc / 200
    wrong[(j+1)] = ww / 200
}
results = data.frame(Dependency=c(0:20*0.05), Correct=corr, Wrong=wrong, Indecisive=1-(corr+wrong))
write.table(results, file="results/crack_real.tab", quote=F, row.names=F)
print(results)

##########################
##Dimensionality Results##
##########################

## Check dependence on dimensions (do for numeric and nominal and mixed type)
# symmetric
dir.create("tempDim")
num = rep(0,5)
numun = rep(0,5)
nom = rep(0,5)
nomun = rep(0,5)
mixed = rep(0,5)
mixedun = rep(0,5)
pos = 1
set.seed(1)
for(j in c(2,3,5,7,11)){
    uu = 0
    un = 0
    oo = 0
    on = 0
    mm = 0
    mn = 0
    dims = 2*j
    dependency = 1.0
    for(i in 1:200){
        ttt = rep("b", dims)
        for(k in 1:dims){
            tt = runif(1, min=0, max=3)
            if(tt <= 1.0){
                ttt[k] = "b"
            }else if(tt <= 2.0){
                ttt[k] = "c"
            }else{
                ttt[k] = "i"
            }
        }
        dd = generateSyntheticData(j,j,5000,dependency,0,ttt)
        dd = rbind(ttt, dd)
        id = round(runif(1, min=1, max=10000000))
        res = crackWrapper(dd, dir="tempDim", name=id, type="", dimX = j)
        if(res$cd == "->"){
            mm = mm + 1
        }else if(res$cd == "--"){
            mn = mn + 1
        }
        ttt = rep("b", dims)
        for(k in 1:dims){
            tt = runif(1, min=0, max=2.0)
            if(tt <= 1.0){
                ttt[k] = "b"
            }else{
                ttt[k] = "c"
            }
        }
        dd = generateSyntheticData(j,j,5000,dependency,0,ttt)
        dd = rbind(ttt, dd)
        id = round(runif(1, min=1, max=10000000))
        res = crackWrapper(dd, dir="tempDim", name=id, type="", dimX = j)
        if(res$cd == "->"){
            oo = oo + 1
        }else if(res$cd == "--"){
            on = on + 1
        }
        ttt = rep("r", dims)
        dd = generateSyntheticData(j,j,5000,dependency,0,ttt)
        dd = rbind(ttt, dd)
        id = round(runif(1, min=1, max=10000000))
        res = crackWrapper(dd, dir="tempDim", name=id, type="", dimX = j)
        if(res$cd == "->"){
            uu = uu + 1
        }else if(res$cd == "--"){
            un = un + 1
        }
    }
    num[pos] = uu / 200
    nom[pos] = oo / 200
    mixed[pos] = mm / 200
    numun[pos] = un / 200
    nomun[pos] = on / 200
    mixedun[pos] = mn / 200
    pos = pos + 1
}
results = data.frame(Dimensions=c(2,3,5,7,11), Mixed=mixed, MixedU=mixedun, Nominal=nom, NominalU=nomun, Numeric=num, NumericU=numun)
write.table(results, file="results/crack_symmetric_dims.tab", quote=F, row.names=F)
print(results)

# asymmetric
dir.create("tempDim")
num = rep(0,5)
numun = rep(0,5)
nom = rep(0,5)
nomun = rep(0,5)
mixed = rep(0,5)
mixedun = rep(0,5)
pos = 1
set.seed(1)
for(j in c(1,3,5,7,11)){
    uu = 0
    un = 0
    oo = 0
    on = 0
    mm = 0
    mn = 0
    dims = 3+j
    d1 = 3
    d2 = j
    dependency = 1.0
    for(i in 1:200){
        if(i %% 2 == 0){
            d1 = j
            d2 = 3
        }else{
            d1 = 3
            d2 = j
        }
        ttt = rep("b", dims)
        for(k in 1:dims){
            tt = runif(1, min=0, max=4)
            if(tt <= 1.0){
                ttt[k] = "b"
            }else if(tt <= 2.0){
                ttt[k] = "c"
            }else{
                ttt[k] = "i"
            }
        }
        dd = generateSyntheticData(d1,d2,5000,dependency,0,ttt)
        dd = rbind(ttt, dd)
        id = round(runif(1, min=1, max=10000000))
        res = crackWrapper(dd, dir="tempDim", name=id, type="", dimX = d1)
        if(res$cd == "->"){
            mm = mm + 1
        }else if(res$cd == "--"){
            mn = mn + 1
        }
        ttt = rep("b", dims)
        for(k in 1:dims){
            tt = runif(1, min=0, max=2.0)
            if(tt <= 1.0){
                ttt[k] = "b"
            }else{
                ttt[k] = "c"
            }
        }
        dd = generateSyntheticData(d1,d2,5000,dependency,0,ttt)
        dd = rbind(ttt, dd)
        id = round(runif(1, min=1, max=10000000))
        res = crackWrapper(dd, dir="tempDim", name=id, type="", dimX = d1)
        if(res$cd == "->"){
            oo = oo + 1
        }else if(res$cd == "--"){
            on = on + 1
        }
        ttt = rep("r", dims)
        dd = generateSyntheticData(d1,d2,5000,dependency,0,ttt)
        dd = rbind(ttt, dd)
        id = round(runif(1, min=1, max=10000000))
        res = crackWrapper(dd, dir="tempDim", name=id, type="", dimX = d1)
        if(res$cd == "->"){
            uu = uu + 1
        }else if(res$cd == "--"){
            un = un + 1
        }
    }
    num[pos] = uu / 200
    nom[pos] = oo / 200
    mixed[pos] = mm / 200
    numun[pos] = un / 200
    nomun[pos] = on / 200
    mixedun[pos] = mn / 200
    pos = pos + 1
}
results = data.frame(Dimensions=c(1,3,5,7,11), Mixed=mixed, MixedU=mixedun, Nominal=nom, NominalU=nomun, Numeric=num, NumericU=numun)
write.table(results, file="results/crack_asymmetric_dims.tab", quote=F, row.names=F)
print(results)
