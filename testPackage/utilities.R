#!/usr/bin/Rscript

sysouts = T

uv = c(1:51,56:70,72:104, 106)
uvnd = c(1:46, 48:51,56:69,72:104, 106)
ref = read.table("../data/causal_tuebingen/README_polished.tab", sep="\t", header=F, stringsAsFactors = F)
meta = read.table("../data/causal_tuebingen/pairmeta.txt")
ref.uv = ref[uv, ]
ref.uvnd = ref[uvnd, ]

readI = function(i, r=ref){
    f = paste(c("../data/causal_tuebingen/", r$V1[i], ".txt"), collapse="")
    t = read.table(f, sep=" ", header=F, stringsAsFactors = F)
    return(t)
}

### na to zero
naTo0 = function(c){
    c[is.na(c)] = 0
    return(c)
}

###### Normalization
normX = function(x, n){
    if(min(x) == max(x)){
        return( rep(n, length(x)) )
    }else{
        return( ((x-min(x)) / (max(x) - min(x))) * n )
    }
}

logg = function(x){
    if(x == 0){
        return(0)
    }else{
        return(log2(x))
    }
}
log2fac = function(n){
    sum = 0
    for(i in 2:n){
        sum = sum + logg(i)
    }
    return(sum)
}
log2nChoosek = function(n, k){
    if(k > n | k == 0){
        return(0)
    }else{
        return(log2fac(n) - log2fac(k) - log2fac(n-k))
    }
}
logN = function(z){
    z = ceiling(z)
    if(z < 1){
        return(0)
    }else{
        logstar = logg(z)
        sum = logstar
        while(logstar > 0){
            logstar = logg(logstar)
            sum = sum + logstar
        }
        return(sum + logg(2.865064))
    }
}

###### Apply Span
setS = function(l){
    s = 1
    if(l > 1000){
        s = 5
    }
    if(l > 2000){
        s = 10
    }
    if(l > 5000){
        s = 20
    }
    if(l > 8000){
        s = 30
    }
    return(s)
}
getDuplicatePositions = function(x){
    last = x[1]
    pos = 0
    s = ""
    for(i in 2:length(x)){
        if(x[i] != last){
            last = x[i]
            if(i - (pos + 1) > 1){
                k = paste(pos, (i-1), sep=":")
                if(nchar(s) > 0){
                    s = paste(s, k, sep=";")
                }else{
                    s = k
                }
            }
            pos = i - 1
        }
    }
    # last elem involved
    if((pos + 1) != length(x)){
        k = paste(pos, length(x), sep=":")
        if(nchar(s) > 0){
            s = paste(s, k, sep=";")
        }else{
            s = k
        }
    }
    return(s)
}

###### Calculate the stochastic complexity
binaryComplexity = function(M){
    if(M < 1){
        return(0.0)
    }
    sum = 1.0
    b = 1.0
    p = 10
    bound = ceiling(2.0 + sqrt(2.0 * M * p * log(10)))
    for(i in 1:bound){
        b = (M - i + 1) * (b / M)
        sum = sum + b
    }
    return(sum)
}
complexityPrecal = function(M, K){
    if(K < 1){
        return(0.0)
    }else if(K == 1){
        return(1.0)
    }else{
        sum = binaryComplexity(M)
        old_sum = 1.0
        if(K > 2){
            for(j in 3:K){
                new_sum = sum + (M * old_sum) / (j-2)
                old_sum = sum
                sum = new_sum
            }
        }
        return(sum)
    }
}
complexity = function(M, K){
    costs = complexityPrecal(M, K);
    if(costs <= 0.0){
        return(0.0)
    }else{
        return(log2(costs))
    }
}
### SC
stochasticComplexity = function(x){
    n = length(x)
    loglikelihood = 0
    tx = table(x)
    for(f in tx){
        freq = as.numeric(f)
        loglikelihood = loglikelihood + freq * (logg(n) - logg(freq))
    }
    sc = loglikelihood + complexity(n, length(tx))
    return(sc)
}
### Conditional SC
conditional_sc = function(x,y){
    y.unique = unique(sort(y))
    sc = 0
    for(i in y.unique){
        xiy = x[y == i]
        sc = sc + stochasticComplexity(xiy)
    }
    return(sc)
}

###### Cumulative resitual entropy + conditional CRE
CRE_abs = function(x.uns, s=TRUE){
    x.uns = normX(x.uns,1)
    sum = 0.0
    x = x.uns
    if(s){
        x = sort(x.uns)
    }
    m = length(x)
    if(m > 1){
        for(i in 1:(m-1)){
            sum = sum - ((abs(x[i+1] - x[i]) * (i/m)) * logg(i/m))
        }
    }
    return(sum)
}
conditionalCE = function(x,y){
    y.unique = unique(sort(y))
    N = length(y)
    H = 0
    for(i in 1:length(y.unique)){
        xiy = x[y == y.unique[i]]
        H = H + ( (length(xiy) / N) * CRE_abs(xiy) )
    }
    return(H)
}

###### Entropy and conditional entropy
entropy = function(x){
    x.unique = unique(sort(x))
    N = length(x)
    H = 0
    if(N != 0){
        for(i in 1:length(x.unique)){
            frac = sum(x == x.unique[i]) / N
            H = H - (frac * logg(frac))
        }
    }
    return(H)
}
conditionalEntropy = function(x,y){
    y.unique = unique(sort(y))
    N = length(y)
    H = 0
    for(i in 1:length(y.unique)){
        xiy = x[y == y.unique[i]]
        H = H + ( (length(xiy) / N) * entropy(xiy) )
    }
    return(H)
}

information_gain = function(x, x1, x2){
    H = entropy(x)
    H1 = entropy(x1)
    H2 = entropy(x2)
    l = length(x)
    l1 = length(x1)
    l2 = length(x2)
    if(l1 == 0 | l2 == 0 | l == 0){
        return(0)
    }else{
        gain = H - (l1 / l) * H1 - (l2 / l) * H2
        return(gain)
    }
}

###### Discretization
equiWidthBinning = function(dataset, bins){
    newd = cut(dataset, bins)
    return(as.numeric(newd))
}

###### Helper
getMeanOccurance = function(){
    df = data.frame(I=uv, X1=rep(0,length(uv)), SD1=rep(0,length(uv)), X2=rep(0,length(uv)), SD2=rep(0,length(uv)), R=ref.uv$V6)
    for(i in 1:length(uv)){
        t = readI(uv[i])
        df$X1[i] = mean(table(t[,1]))
        df$SD1[i] = sd(table(t[,1]))
        df$X2[i] = mean(table(t[,2]))
        df$SD2[i] = sd(table(t[,2]))
    }
    return(df)
}
getMinOccurance = function(){
    df = data.frame(I=uv, X1=rep(0,length(uv)), SD1=rep(0,length(uv)), X2=rep(0,length(uv)), SD2=rep(0,length(uv)), R=ref.uv$V6)
    for(i in 1:length(uv)){
        t = readI(uv[i])
        df$X1[i] = min(table(t[,1]))
        df$SD1[i] = sd(table(t[,1]))
        df$X2[i] = min(table(t[,2]))
        df$SD2[i] = sd(table(t[,2]))
    }
    return(df)
}
getDims = function(){
    df = data.frame(I=uv, L=rep(0,length(uv)))
    for(i in 1:length(uv)){
        t = readI(uv[i])
        df$L[i] = dim(t)[1]
    }
    return(df)
}

###### Apply IPD
applyIPD = function(t){
    x = normX(t[,1], 1)
    y = normX(t[,2], 1)
    df= data.frame(x,y,rep(0,length(x)))
    write.table(df, file="temp/test_ipd.csv", row.names=F, col.names=F, sep=";", quote=F)
    system(paste(c("java -jar ipd.jar -FILE_INPUT temp/test_ipd.csv -FILE_CP_OUTPUT temp/cuts.txt -FILE_RUNTIME_OUTPUT temp/runtime.txt -FILE_DATA_OUTPUT temp/out.txt -NUM_ROWS ", length(x), " -NUM_MEASURE_COLS 2 -NUM_CAT_CONTEXT_COLS 0 -MAX_VAL 1.0 -METHOD 0"), collapse=""), ignore.stdout=sysouts)
    disc = read.table("temp/out.txt", sep=",", comment.char = "@")
    Xd = disc$V1
    Yd = disc$V2
    ret = list(xd=Xd, yd=Yd)
}

###### Calculate decision rate (-- = 1/2)
decision_rate_w = function(res){
    return(decision_rate(res$Correct, res$Eps, res$Cds))
}
decision_rate = function(corr, eps, cds){
    df = data.frame(A=corr, B=abs(eps), C=cds)
    df = df[with(df, order(-B)),]
    sum = 0
    dr = rep(0, dim(df)[1])
    for(i in 1:dim(df)[1]){
        if(df$C[i] == "--"){
            sum = sum + 0.5
        }else{
            sum = sum + df$A[i]
        }
        dr[i] = sum / i
    }
    return(dr)
}
decision_rate_meta = function(res, uv=uv, oneHalf=T){
    corr = res$Correct
    eps = res$Eps
    cds = res$Cds
    df = data.frame(Uv=uv, A=corr, B=abs(eps), C=cds)
    df = df[with(df, order(-B)),]
    sum = 0
    total = 0
    step = rep(0, dim(df)[1])
    dr = rep(0, dim(df)[1])
    for(i in 1:dim(df)[1]){
        val = meta$V6[df$Uv[i]]
        total = total + val
        if(df$C[i] == "--"){
            sum = sum + 0.5 * val
            if(oneHalf){
                dr[i] = sum / total
            }else{
                if(i > 1){
                    dr[i] = dr[i-1]
                }
            }
        }else{
            sum = sum + df$A[i] * val
            dr[i] = sum / total
        }
        step[i] = total
    }
    step = step / total
    return(list(D=c(1,dr), S=c(0,step)))
}
decision_rate_meta_pos = function(res, uv=uv){
    corr = res$Correct
    eps = res$Eps
    cds = res$Cds
    df = data.frame(Uv=uv, A=corr, B=abs(eps), C=cds)
    df = df[with(df, order(B)),]
    sum = 0
    total = 0
    step = rep(0, dim(df)[1])
    dr = rep(0, dim(df)[1])
    for(i in 1:dim(df)[1]){
        val = meta$V6[df$Uv[i]]
        total = total + val
        if(df$C[i] == "--"){
            sum = sum + 0.5 * val
        }else{
            sum = sum + df$A[i] * val
        }
        dr[i] = sum / total
        step[i] = total
    }
    step = step / total
    return(list(D=dr, S=step))
}
decision_rate_w_pos = function(res){
    corr = res$Correct
    eps = res$Eps
    cds = res$Cds
    df = data.frame(A=corr, B=abs(eps), C=cds)
    df = df[with(df, order(B)),]
    sum = 0
    dr = rep(0, dim(df)[1])
    for(i in 1:dim(df)[1]){
        if(df$C[i] == "--"){
            sum = sum + 0.5
        }else{
            sum = sum + df$A[i]
        }
        dr[i] = sum / i
    }
    return(dr)
}

###### SSE
SSE = function(x, yhead){
    e = sum((x-yhead)^2)
    return(e)
}

###### PlotI
plotI = function(i){plot(readI(i))}

###### MDL score
resolution = 0.01
setResolution = function(val){
    resolution <<- val
}
uniform = function(data){
    N = length(data)
    data_point = logg((max(data) - min(data))/resolution)
    return(N * data_point)
}
gaussian_score_emp = function(x){
    sse = SSE(x, mean(x))
    var = sse / length(x)
    sigma = sqrt(var)
    return(gaussian_score(sigma, x))
}
gaussian_score_old = function(sigma, x){
    sse = SSE(x, mean(x))
    n = length(x)
    sigmasq = sigma^2
    if(sse == 0.0 | sigmasq == 0.0){
        return(0.0)
    }else{
        err = (sse / (2 * sigmasq * log(2))) + ((n/2) * logg(2 * pi * sigmasq)) - n * logg(resolution)
        return(err)
    }
}
gaussian_score_emp_sse = function(sse, n){
    var = sse / n
    sigma = sqrt(var)
    return(gaussian_score_sse(sigma, sse, n))
}

gaussian_score_sse = function(sigma, sse, n){
    sigmasq = sigma^2
    if(sse == 0.0 | sigmasq == 0.0){
        return(0.0)
    }else{
        err = (sse / (2 * sigmasq * log(2))) + ((n/2) * logg(2 * pi * sigmasq)) - n * logg(resolution)
        return(max(err,0))
    }
}
best_encoding = function(resid){
    sse = SSE(resid, 0)
    gauss = gaussian_score_emp_sse(sse, length(resid))
    unif = uniform(resid)
    if(unif < gauss){
        return(unif + 1)
    }else{
        return(gauss + 1)
    }
}
gaussian_score = function(resid){
    n = length(resid)
    mresid = mean(resid)
    sdresid = sd(resid)
    xt = 1:n / (n + 1)
    xt = qnorm(xt, sd=sdresid, mean=mresid)
    sse = sum((sort(resid) - xt)^2)
    return(gaussian_score_sse(sdresid, sse, n))
}

parameterScore = function(model, n=0){
    sum = 0
    #coeff = naTo0(model$coefficients) / resolution
    coeff = naTo0(model$coefficients)
    #count = 0
    for(c in coeff){
        if(c != 0){
            #count = count + 1
            ca = abs(c)
            cdummy = ca
            prec = 1
            while(cdummy < 1000){
                cdummy = cdummy * 10
                prec = prec + 1
            }
            #sum = sum + logN(ca) + 1 ## pos or neg
            sum = sum + logN(cdummy) + logN(prec) + 1
        }
    }
    #sum = sum + count * logg(NOFC) ## which param
    return(sum)
}

### igci
S = function(x.uns){
    sum = 0.0
    x = sort(x.uns)
    m = length(x)
    for(i in 1:(m-1)){
        sum = sum + logg(abs(x[i+1] - x[i]))
    }
    sum = sum * (1/(m-1))
    sum = sum + logg(m)
    return(sum)
}
normalize = function(tt){
    t = tt
    # normalize to mean 0 and sd 1
    for(j in 1:(dim(t)[2])){
        t[,j] = as.numeric(t[,j])
        t[,j] = (t[,j] - min(t[,j])) / (max(t[,j]) - min(t[,j]))
    }
    return(t)
}

## mindiff
mindiff = function(x){
    xs = sort(x)
    diff = 0.01
    for(i in 1:(length(x)-1)){
        new_diff = xs[i+1] - xs[i]
        if(new_diff != 0 & new_diff < diff){
            diff = new_diff
        }
    }
    return(diff)
}
data_precision = function(x){
    precision = 1
    set = x != round(x)
    x = x[set]
    while(length(x) > 0){
        precision = precision / 10
        x = 10 * x
        set = x != round(x)
        x = x[set]
    }
    return(precision)
}

##### IGCI
IGCI = function(t){
    t = normalize(t)
    diffs = S(t[,2]) - S(t[,1])
    causd = "--"
    if(diffs < 0){
        causd = "->"
    }else if(diffs > 0){
        causd = "<-"
    }
    return(list(cd=causd, epsilon=diffs))
}
