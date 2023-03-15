#!/usr/bin/Rscript

initialSampleSize = 5000

logg = function(x){
    if(x == 0){
        return(0)
    }else{
        return(log2(x))
    }
}

getYWithGauss = function(x, Sd=1, q=1, b=0, c=0,d=1,e=0,f=0, roundIt=F){
    xe = x
    xe[x == 0] = 0.01
    noise = rnorm(length(x), mean=0, sd=Sd)
    if(q != 1){
        sq0 = noise < 0
        noise = abs(noise)^q
        noise[sq0] = -noise[sq0]
    }
    if(roundIt){
        noise = round(noise)
    }
    y = b * (x^3) + c * (x^2) + e * exp(x) + d * x + f * xe^(-1) + noise
    df = data.frame(X=x, Y=y)
    return(df)
}

getYWithUnif = function(x, t=0.5, b=0, d=1,e=0,f=0, roundIt=F){
    xe = x
    xe[x == 0] = 0.01
    noise = runif(length(x), min=-t, max=t)
    if(roundIt){
        noise = round(noise)
    }
    y = b * (x^3) + e * exp(x) + d * x + f * xe^(-1) + noise
    df = data.frame(X=x, Y=y)
    return(df)
}

getYWithLocal = function(x, b=0, d=1,e=0,f=0, roundIt=F){
    xd = round(x * 10000)
    xs = unique(sort(xd))
    y = x
    for(s in xs){
        pos = (xd == s)
        val = s / 10000
        l = sum(pos)
        xhead = 1:sum(pos) - sum(pos) / 2
        par = runif(1, min=0.1, max=0.5)
        yn = xhead * d + val * d + val^2 * b
        y[pos] = y[pos] + yn
    }
    df = data.frame(X=x, Y=y)
    return(df)
}

getYWithNonAdditive = function(x, v=0.25, sigma=1, b=0, c=0, d=1,e=0,f=0){
    xe = x
    xe[x == 0] = 0.01
    y = b * (x^3) + e * exp(x) + c * (x^2) + d * x + f * xe^(-1)
    for(i in 1:length(x)){
        X = x[i]
        E = sigma * rnorm(1) * abs(sin(2*pi*v*X)) + 0.25 * sigma * rnorm(1) * abs(sin(20 * pi * v * X))
        y[i] = y[i] + E
    }
    df = data.frame(X=x, Y=y)
    return(df)
}
generateMVNominalData = function(k,l,m,connections){
    if(connections > k){
        connections = k
    }
    initialSampleSize <<- m
    # create initial df
    n = k + l
    df = data.frame(rep(0,m))
    for(i in 2:n){
        df = cbind(df, rep(0,m))
    }
    # create variables
    for(i in 1:n){
        xi = generateCategoricalData(m)
        df[i] = xi
    }
    for(j in (k+1):n){
        dependences = sample(1:k, connections)
        minDom = 10
        for(dep in dependences){
            currDom = length(unique(sort(df[, dep])))
            if(currDom < minDom){
                minDom = currDom
            }
        }
        for(i in 1:minDom){
            pos = min(i, connections)
            ind = df[, dependences[pos]] == i
            for(v in 1:connections){
                if(v != pos){
                    ind = ind & df[, dependences[v]] != i
                }
            }
            df[ind, j] = i
        }
        domJ = length(unique(sort(df[,j])))
        noise = round(rnorm(m,mean=0,sd=0.5))
        # add noise
        df[,j] = df[,j] + noise
        # make sure noise does not change the domain
        df[df[,j] > domJ, j] = domJ
        df[df[,j] < 0, j] = 0
    }
    return(df)
}

generateBinaryDataP = function(n, p1){
    p = 1000 * p1
    xr = round(runif(n, min=0, max=1000))
    xb = rep(0,n)
    xb[xr <= p] = 1
    return(xb)
}
generateBinaryData = function(n, border=0.3){
    p = round(runif(1, min=(1000*border), max=1000-(1000*border)))
    return(generateBinaryDataP(n=n, p1=(p/1000)))
}

generateCategoricalData = function(n){
    x1 = c()
    p0 = runif(1, min=0, max=6)
    categories = 3
    if(p0 > 3 & p0 <=5){
        categories = 4
    }else if(p0 > 5){
        categories = 5
    }
    px = rep(0, categories)
    for(i in 1:categories){
        px[i] = runif(1, min=0.1, max=1.0)
    }
    px = px / sum(px)
    for(i in 1:(categories-1)){
        x1 = c(x1, rep(i, floor(px[i]*n)))
    }
    x1 = c(x1, rep(categories, n-length(x1)))
    x1 = sample(x1)
    return(x1)
}

## Splits
#########
generateBinarySplit = function(p, l, r){
    oldDomain = length(unique(sort(p)))
    if(oldDomain <= 1 | length(l) <= 1 | length(r) <= 1){
        return(list(a=p, b=c()))
    }else{
        x1 = rep(0, length(l))
        x2 = rep(0, length(r))
        counter = 0
        repeat{
            p1 = runif(1, min=0, max=0.49)
            x1 = generateBinaryDataP(length(l), p1)
            x2 = generateBinaryDataP(length(r), 1-p1)
            gain = information_gain(p, x1, x2)
            counter = counter + 1
            if(gain > 0.0 | counter > 100){
                break
            }
        }
        return(list(a=x1, b=x2))
    }
}
generateCategoricalSplit = function(p,l,r){
    numsL = sample(unique(sort(l)))
    numsR = sample(unique(sort(r)))
    if(length(numsR) <= 1 | length(numsL) <= 1){
        return(list(a=p, b=c()))
    }else{
        right=T
        dummy = r
        if(length(l) > length(r) & length(numsL) > 1){
            dummy = l
            right=F
        }
        indices = 1:length(dummy)
        num = 0.5 + runif(1, min=0.1, max=0.4)
        new_ones = sample(indices)[1:round(num*length(dummy))]
        if(right){
            dummy[new_ones] = numsR[1]
        }else{
            dummy[new_ones] = numsL[1]
        }
        dummy = sample(dummy)
        if(right){
            return(list(a=l, b=dummy))
        }else{
            return(list(a=dummy, b=r))
        }
    }
}
generateCategoricalSplit_old = function(p,l,r){
    oldDomain = length(unique(sort(p)))
    domR = length(unique(sort(r)))
    nums = sample(unique(sort(p)))
    if(oldDomain <= 1 | length(l) <= 1 | length(r) <= 1){
        return(list(a=p, b=c()))
    }else{
        counter = 0
        repeat{
            categories = round(runif(1, min=(oldDomain-domR), max=oldDomain))
            px = rep(0, categories)
            for(i in 1:categories){
                px[i] = runif(1, min=0.0, max=1.0)
            }
            px = px / sum(px)
            x1 = c()
            if(categories > 1){
                for(i in 1:(categories-1)){
                    count1 = floor(px[i] * length(l))
                    x1 = c(x1, rep(nums[i], count1))
                }
            }
            x1 = c(x1, rep(nums[categories], length(l)-length(x1)))
            x1 = sample(x1)
            gain = information_gain(p, x1, r)
            counter = counter + 1
            if(gain > 0.0 | counter > 100){
                return(list(a=x1, b=r))
            }
        }
    }
}
generateRVSplit = function(p,l,r){
    mean1 = runif(1, min=1, max=5)
    shift = runif(1, min=1, max=10)
    sd1 = runif(1, min=1, max=5)
    sd2 = runif(1, min=1, max=5)
    if(length(l) <= 1 | length(r) >= 1){
        return(list(a=p, b=c()))
    }else{
        x1 = rnorm(length(l), mean=mean1, sd=sd1)
        x2 = rnorm(length(r), mean=(mean1+shift), sd=sd2)
        return(list(a=x1, b=x2))
    }
}
generateSplit = function(p,l,r,type){
    if(type == "b"){
        return(generateBinarySplit(p,l,r))
    }else if(type== "c"){
        return(generateCategoricalSplit(p,l,r))
    }else{
        return(generateRVSplit(p,l,r))
    }
}

generateSyntheticData = function(k,l,m,dependence,inbetween,types){
    if(dependence < 0.0 | dependence > 1.0){
        return(NULL)
    }
    if(inbetween < 0.0 | inbetween > 1.0){
        return(NULL)
    }
    initialSampleSize <<- m
    # create initial df
    n = k + l
    df = data.frame(rep(0,m))
    subsets = data.frame(rep(1,m))
    for(i in 2:n){
        df = cbind(df, rep(0,m))
        subsets = cbind(subsets, rep(1,m))
    }
    # create variables
    for(i in 1:n){
        if(types[i] == "b"){
            xi = generateBinaryData(m)
            df[i] = xi
        }else if(types[i] == "c"){
            xi = generateCategoricalData(m)
            df[i] = xi
        }else{
            tX = runif(1, min=1, max=10)
            gaussianity = runif(1, min=0.7, max=1.3)
            xi = rnorm(m, sd=tX)
            sq0 = xi < 0
            xi = abs(xi)^(gaussianity)
            xi[sq0] = -xi[sq0]
            df[i] = xi
        }
    }
    # create redundancies
    for(i in 1:n){
        sta = 1
        sto = i - 1
        if(i > k){
            sta = k + 1
        }
        if(sta > sto){
            next
        }
        r = runif(1, min=0, max=0.999)
        if(r < inbetween){
            j = 1
            if(sta == sto){
                j = sta
            }else{
                j = sample(sta:sto, 1)
            }
            subseti = round(runif(1, min=1, max=max(subsets[,i])))
            subset = subsets[,i] ==  subseti
            old = df[subset, i]
            if(length(old) < 5){
                next
            }
            if(types[j] == "b" | types[j] == "c"){ ## get both groups and then split also for categ and rv
                candidate = df[subset, j]
                if(types[i] != "b" & types[i] != "c"){
                    ## split
                    catInd = round(runif(1, min=1, max=(length(candidate))))
                    cat = candidate[catInd]
                    indA = subset & df[j] <= cat
                    indB = subset & df[j] > cat
                    dfA = df[indA, i]
                    dfB = df[indB, i]
                    r = generateSplit(old, dfA, dfB, types[i])
                    if(length(r$a) != length(old)){
                        df[indA, i] = r$a
                        df[indB, i] = r$b
                    }
                }else{
                    domC = unique(sort(candidate))
                    cN = length(candidate)
                    noise = runif(1, min=0.05, max=0.15)
                    cN10 = round(noise * cN)
                    transfInd = sample(1:cN, cN10)
                    for(q in transfInd){
                        candidate[q] = sample(domC, 1)
                    }
                    df[subset, i] = candidate
                }
            }else{
                candidate = df[subset, j]
                p0 = runif(1, min=0.0, max=2)
                if(types[i] == "b" | types[i] == "c"){
                    ## split
                    catInd = round(runif(1, min=1, max=(length(candidate))))
                    cat = candidate[catInd]
                    indA = subset & df[j] <= cat
                    indB = subset & df[j] > cat
                    dfA = df[indA, i]
                    dfB = df[indB, i]
                    r = generateSplit(old, dfA, dfB, types[i])
                    if(length(r$a) != length(old)){
                        df[indA, i] = r$a
                        df[indB, i] = r$b
                    }
                }else{
                    ## linear
                    sdN = runif(1, min=1, max=3)
                    dr = getYWithGauss(candidate, Sd=sdN, d=1)
                    df[subset, i] = dr[,2]
                }
            }
        }
    }
    # create Y variables
    for(j in (k+1):n){
        for(i in 1:k){
            rd = runif(1, min=0, max=0.999)
            if(rd < dependence){
                if(max(subsets[,j]) >= 2){
                    next
                }
                subseti = round(runif(1, min=1, max=max(subsets[,j])))
                subset = subsets[,j] ==  subseti
                old = df[subset, j]
                if(length(old) < 5){
                    next
                }
                if(types[i] == "b" | types[i] == "c"){ ## get both groups and then split also for categ and rv
                    candidate = df[subset, i]
                    if(types[j] != "b" & types[j] != "c"){
                        ## split
                        catInd = round(runif(1, min=1, max=(length(candidate))))
                        cat = candidate[catInd]
                        indA = subset & df[i] <= cat
                        indB = subset & df[i] > cat
                        dfA = df[indA, j]
                        dfB = df[indB, j]
                        r = generateSplit(old, dfA, dfB, types[j])
                        if(length(r$a) != length(old)){
                            df[indA, j] = r$a
                            df[indB, j] = r$b
                        }
                    }else{
                        domC = unique(sort(candidate))
                        cN = length(candidate)
                        noise = runif(1, min=0.05, max=0.15)
                        cN10 = round(noise * cN)
                        transfInd = sample(1:cN, cN10)
                        for(q in transfInd){
                            candidate[q] = sample(domC, 1)
                        }
                        df[subset, j] = candidate
                        subsets[0, j] = length(domC)
                    }
                }else{
                    candidate = df[subset, i]
                    p0 = runif(1, min=0.0, max=4)
                    if(types[j] == "b" | types[j] == "c"){
                        p0 = 0
                    }else if(k == 1){
                        p0 = 4
                    }
                    if(p0 <= 2){
                        ## split
                        catInd = round(runif(1, min=1, max=(length(candidate))))
                        cat = candidate[catInd]
                        indA = subset & df[i] <= cat
                        indB = subset & df[i] > cat
                        dfA = df[indA, j]
                        dfB = df[indB, j]
                        r = generateSplit(old, dfA, dfB, types[j])
                        if(length(r$a) != length(old)){
                            df[indA, j] = r$a
                            df[indB, j] = r$b
                            subsets[indB, j] = rep(max(subsets[,j]) + 1, length(dfB))
                        }
                    }else if(p0 > 2 & p0 <= 3){
                        ## linear
                        sdN = runif(1, min=1, max=3)
                        dr = getYWithGauss(candidate, Sd=sdN, d=1)
                        if(max(subsets[,j]) == 1){
                            df[subset, j] = dr[,2]
                        }else{
                            df[subset, j] = df[subset, j] + dr[,2]
                        }
                    }else{
                        ## quadratic
                        sdN = runif(1, min=1, max=3)
                        dr = getYWithGauss(candidate, Sd=sdN, c=1)
                        if(max(subsets[,j]) == 1){
                            df[subset, j] = dr[,2]
                        }else{
                            df[subset, j] = df[subset, j] + dr[,2]
                        }
                    }
                }
            }
        }
    }
    return(df)
}

