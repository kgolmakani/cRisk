expit = function(x){
  return(exp(x)/(1+exp(x)))
}

mc.sim = function(P,x1) {
  sim = as.numeric(1)
  m = ncol(P)
  sim = sample(0:(m-1),1, prob = P[(x1+1),])
  sim 
}

addCol = function(df){
  n = nrow(df)
  df = df[order(df$StudyID_c,df$round),]
  df$s1 = c(NA,df$result[-n])
  df$s1[match(unique(df$StudyID_c),df$StudyID_c )] = NA 
  return (df)
}

round_1 = function(m1, m2){
  function(i, h){
    if (h==1){
      if (i==0)
        return(1/(1+exp(m1[1])+exp(m1[2])))
      if(i==1)
        return(exp(m1[1])/(1+exp(m1[1])+exp(m1[2])))
      if (i==2)
        return(exp(m1[2])/(1+exp(m1[1])+exp(m1[2])))
    } else{
      if(i==1)
        return(expit(m2[1]+m2[2]*h))
      if(i==0)
        return(1 -expit(m2[1]+m2[2]*h))
      if(i==2)
        return(0)
    }
  }
}


Pmat = function(m1,m2){
  function(h,t){
    matrix(c(1-expit(m1[1]+m1[2]*h+m1[3]*t+m1[4]*h*t), 1-expit(m2[1]+m2[2]*h+m2[3]*t+m2[4]*h*t),0,
             expit(m1[1]+m1[2]*h+m1[3]*t+m1[4]*h*t), expit(m2[1]+m2[2]*h+m2[3]*t+m2[4]*h*t), 0,
             0, 0, 1), ncol=3)}
}


Ppmat = function(m1, m2){
  function(h, t){
    matrix(c(1/(1+exp(m1[1]+m1[3]*h+m1[5]*t+m1[7]*h*t)+exp(m1[2]+m1[4]*h+m1[6]*t+m1[8]*h*t)), 1/(1+exp(m2[1]+m2[3]*h+m2[5]*t+m2[7]*h*t)+exp(m2[2]+m2[4]*h+m2[6]*t+m2[8]*h*t)), 0,
             exp(m1[1]+m1[3]*h+m1[5]*t+m1[7]*h*t)/(1+exp(m1[1]+m1[3]*h+m1[5]*t+m1[7]*h*t)+exp(m1[2]+m1[4]*h+m1[6]*t+m1[8]*h*t)),exp(m2[1]+m2[3]*h+m2[5]*t+m2[7]*h*t)/(1+exp(m2[1]+m2[3]*h+m2[5]*t+m2[7]*h*t)+exp(m2[2]+m2[4]*h+m2[6]*t+m2[8]*h*t)), 0,
             exp(m1[2]+m1[4]*h+m1[6]*t+m1[8]*h*t)/(1+exp(m1[1]+m1[3]*h+m1[5]*t+m1[7]*h*t)+exp(m1[2]+m1[4]*h+m1[6]*t+m1[8]*h*t)), exp(m2[2]+m2[4]*h+m2[6]*t+m2[8]*h*t)/(1+exp(m2[1]+m2[3]*h+m2[5]*t+m2[7]*h*t)+exp(m2[2]+m2[4]*h+m2[6]*t+m2[8]*h*t)), 1), ncol=3)
  }
}


round_1_nh = function(m1_nh, m2_nh){
  function(i, h, g, x){
    if (h==1){
      if (i==0)
        return(1/(1+exp(m1_nh[2*(g-1)+1,1]+m1_nh[2*(g-1)+1,2]*x)+exp(m1_nh[2*g,1]+m1_nh[2*g,2]*x)))
      if(i==1)
        return(exp(m1_nh[2*(g-1)+1,1]+m1_nh[2*(g-1)+1,2]*x)/(1+exp(m1_nh[2*(g-1)+1,1]+m1_nh[2*(g-1)+1,2]*x)+exp(m1_nh[2*g,1]+m1_nh[2*g,2]*x)))
      if (i==2)
        return(exp(m1_nh[2*g,1]+m1_nh[2*g,2]*x)/(1+exp(m1_nh[2*(g-1)+1,1]+m1_nh[2*(g-1)+1,2]*x)+exp(m1_nh[2*g,1]+m1_nh[2*g,2]*x)))
    } else{
      if(i==1)
        return(expit(m2_nh[g,1]+m2_nh[g,2]*h+m2_nh[g,3]*x))
      if(i==0)
        return(1-expit(m2_nh[g,1]+m2_nh[g,2]*h+m2_nh[g,3]*x))
      if(i==2)
        return(0)
    }
  }
}


Pmat_nh = function(m1, m2){
  function(h,t,g,x){
    matrix(c(1-expit(m1[g,1]+m1[g,2]*x+m1[g,3]*h+m1[g,4]*t+m1[g,5]*h*t), 1-expit(m2[g,1]+m2[g,2]*x+m2[g,3]*h+m2[g,4]*t+m2[g,5]*h*t),0,
             expit(m1[g,1]+m1[g,2]*x+m1[g,3]*h+m1[g,4]*t+m1[g,5]*h*t), expit(m2[g,1]+m2[g,2]*x+m2[g,3]*h+m2[g,4]*t+m2[g,5]*h*t), 0,
             0, 0, 1), ncol=3)}
}


Ppmat_nh = function(m1, m2){
  function(h, t, g, x){
    matrix(c(1/(1+exp(m1[2*(g-1)+1,1]+m1[2*(g-1)+1,2]*x+m1[2*(g-1)+1,3]*h+m1[2*(g-1)+1,4]*t+m1[2*(g-1)+1,5]*h*t)+exp(m1[2*g,1]+m1[2*g,2]*x+m1[2*g,3]*h+m1[2*g,4]*t+m1[2*g,5]*h*t)), 1/(1+exp(m2[2*(g-1)+1,1]+m2[2*(g-1)+1,2]*x+m2[2*(g-1)+1,3]*h+m2[2*(g-1)+1,4]*t+ m2[2*(g-1)+1,5]*h*t)+exp(m2[2*g,1]+m2[2*g,2]*x+m2[2*g,3]*h+m2[2*g,4]*t +m2[2*g,5]*h*t)), 0,
             exp(m1[2*(g-1)+1,1]+m1[2*(g-1)+1,2]*x+m1[2*(g-1)+1,3]*h+m1[2*(g-1)+1,4]*t+m1[2*(g-1)+1,5]*h*t)/(1+exp(m1[2*(g-1)+1,1]+m1[2*(g-1)+1,2]*x+m1[2*(g-1)+1,3]*h+m1[2*(g-1)+1,4]*t+m1[2*(g-1)+1,5]*h*t)+exp(m1[2*g,1]+m1[2*g,2]*x+m1[2*g,3]*h+m1[2*g,4]*t+m1[2*g,5]*h*t)),exp(m2[2*(g-1)+1,1]+m2[2*(g-1)+1,2]*x+m2[2*(g-1)+1,3]*h+m2[2*(g-1)+1,4]*t+ m2[2*(g-1)+1,5]*h*t)/(1+exp(m2[2*(g-1)+1,1]+m2[2*(g-1)+1,2]*x+m2[2*(g-1)+1,3]*h+m2[2*(g-1)+1,4]*t+ m2[2*(g-1)+1,5]*h*t)+exp(m2[2*g,1]+m2[2*g,2]*x+m2[2*g,3]*h+m2[2*g,4]*t +m2[2*g,5]*h*t)), 0,
             exp(m1[2*g,1]+m1[2*g,2]*x+m1[2*g,3]*h+m1[2*g,4]*t+m1[2*g,5]*h*t)/(1+exp(m1[2*(g-1)+1,1]+m1[2*(g-1)+1,2]*x+m1[2*(g-1)+1,3]*h+m1[2*(g-1)+1,4]*t+m1[2*(g-1)+1,5]*h*t)+exp(m1[2*g,1]+m1[2*g,2]*x+m1[2*g,3]*h+m1[2*g,4]*t+m1[2*g,5]*h*t)), exp(m2[2*g,1]+m2[2*g,2]*x+m2[2*g,3]*h+m2[2*g,4]*t +m2[2*g,5]*h*t)/(1+exp(m2[2*(g-1)+1,1]+m2[2*(g-1)+1,2]*x+m2[2*(g-1)+1,3]*h+m2[2*(g-1)+1,4]*t+ m2[2*(g-1)+1,5]*h*t)+exp(m2[2*g,1]+m2[2*g,2]*x+m2[2*g,3]*h+m2[2*g,4]*t +m2[2*g,5]*h*t)), 1), ncol=3)
  }}


cum.var = function(theta, M, G, N, N2, p){
  G = rep(G,M)
  for (j in (M-1):1){
    for (l in (j+1):M){
      for (s in 1:M){ 
        N[j,l] = N[j,l] + exp(G[j]*(s-(j+1)))*N2[l,s]
      }
    }		
  }
  var.theta = theta*(1-theta)
  cov.theta = vector(length = M)
  temp = NULL
  for (i in 1:M){
    temp = -matrix(theta[i,],ncol =1)%*%matrix(theta[i,],nrow =1)
    diag(temp) = 0
    cov.theta[i] = sum(c(temp))
  }	
  var.theta = apply(sweep(var.theta,1,p[1:M]^2/apply(N,1,sum),"*"),2,sum)-cov.theta/apply(N,1,sum)
  var.theta = sum(var.theta)
}


cens.bias = function(data,fpl, M, alpha, l, functions){
  xround_1 = functions$round_1
  xPpmat = functions$Ppmat
  xPmat = functions$Pmat
  xround_1_nh = functions$round_1_nh
  xPmat_nh = functions$Pmat_nh
  xPpmat_nh = functions$Ppmat_nh
  #data = subset(data, round <= cens2)
  theta = matrix(data = 0, nrow=M, ncol=M)
  # compute theta[h, j] for lower triangle
  cartesianProd = function(a){
    if(a==1){
      if(l==1){
        return(data.frame('i1'=1))
      } else{ stop('not possible')}
    } else{
      df = expand.grid(rep(list(c(0,1)), a-1))
      names(df) = paste0('i', 1:(a-1))
      df = subset(df, rowSums(df)==l-1)
      df[, paste0('i',a )] = 1
      return(df)}}
  for (j in l:M){
    df = cartesianProd(j)
    for (h in j:M){
      df$term = apply(df,1 , FUN = function(x){
        x = as.numeric(x)
        C = xround_1(x[1], h) #pi1
        if (j>1){
          C = C*prod(sapply(1:(j-1), function(k){
            if (k+1 == h){
              return (xPpmat(h, k+1)[x[k]+1, x[k+1]+1])
            } else{
              return (xPmat(h, k+1)[x[k]+1, x[k+1]+1])
            }
          }))}
        return (C)
      })
      theta[h, j] = sum(df$term, na.rm = T)
    }}
  
  B = list() # B[[h]][b,c] = P(S=b, nh = c-1)  
  A = list() # A[[h]][a,b,c] = P(w_l=a | S=b, nh = c-1)
  for(h in 1:(M-1)){
    B[[h]] = matrix(data=0, nrow=M, ncol=min(l-1, h)+1)  
    for (b in h:M)
      for (c in 0:min(l-1, h))
        B[[h]][b, c+1] = nrow(subset(data, S==b & get(paste0('n',h))==c))/nrow(data)
      
      A[[h]] = array(data=0, dim=c(M+1, M, 1+min(l-1, h)))
      # case I: a<=b
      for (b in max(l, h+1):M)
        for (a in l:b)
          for (c in 0:min(l-1, h)){
            df = cartesianProd(a)
            df$term = apply(df,1 , FUN = function(x){
              x = as.numeric(x)
              C = xround_1_nh(x[1], b, h, c)
              if (a>1){
                C = C*prod(sapply(1:(a-1), function(k){
                  if (k+1 == b){
                    return (xPpmat_nh(b, k+1, h, c)[x[k]+1, x[k+1]+1])
                  } else{
                    return (xPmat_nh(b, k+1, h, c)[x[k]+1, x[k+1]+1])
                  }  
                }))}
              return (ifelse(is.nan(C), 0, C))
            })
            A[[h]][a, b, c+1] = sum(df$term, na.rm = T)
          }
      # case II: a > b = M
      for (c in 0:min(l-1, h))
        A[[h]][M+1, M, c+1] = 1 - sum(sapply(l:M, function(jj) {A[[h]][jj, M, c+1]}))
      
      # case III: a >b & b < M. Also in the final round we calculate A[[h]][a,h,c+1] 
      for (b in (M-1):h)
        for (a in max(l,(b+1)):ifelse(b>h, M+1, M))
          for (c in 0:min(l-1, h)){
            if(a< M+1){
              E = diag(exp(alpha*(l-c)*seq(M+1)))
              A[[h]][a, b, c+1] = sum(A[[h]][a, , c+1]*B[[h]][, c+1])*exp(alpha*(l-c)*a)/sum(E%*%as.matrix(A[[h]][, , c+1]%*%B[[h]][, c+1, drop=F]))
            } else{
              A[[h]][a, b, c+1] = 1 - sum(sapply(l:M, function(jj) {A[[h]][jj, b, c+1]}))
            }
            
          }
  }
  
  # compute theta[h, j] for upper triangle
  for (h in 1:(M-1))
    for (j in max(h+1, l):M)
      theta[h, j] = sum(A[[h]][j,h, ]*B[[h]][h, ])*nrow(data)/nrow(subset(data, S==h))
  
  p1 = summary(factor(data$S[!duplicated(data$StudyID_c)],levels=c(1:M)))/length(data$S[!duplicated(data$StudyID_c)]) 
  D = apply(theta,1,sum)
  D = as.matrix(D, drop=F)
  risk = p1%*%D
  risk = as.numeric(risk)
  data$RESINIT_C = ifelse(fpl==data$round, 1,0)
  N1 = table(data$S[data$RESINIT_C==1], data$round[data$RESINIT_C==1])
  N2 = matrix(0, ncol=M, nrow=M)
  N2[1:nrow(N1), 1:ncol(N1)]=N1
  N = table(data$S,data$round)
  varp = cum.var(theta, M, alpha, N, N2, p1)
  se = sqrt(varp)
  return (list(risk=risk, se = se))
}


disc.surv = function(data,fpl.cens,M,l){
  theta_hat = vector(length = M)
  r = vector (length = M)
  s = vector (length =  M)
  varp = vector (length = M)
  for (k in 1:M){
    theta_hat[k] = mean(fpl.cens[data$S>=k & (fpl.cens>=k | fpl.cens==0)] == k)
    r[k] =  length(fpl.cens[data$S>=k & (fpl.cens>=k | fpl.cens==0) ]) 
    s[k]= sum(fpl.cens[data$S>=k & (fpl.cens>=k | fpl.cens==0) ]>k|fpl.cens[data$S>=k & (fpl.cens>=k | fpl.cens==0) ]==0)
  }
  theta_hat = ifelse(is.na(theta_hat),0, theta_hat)
  s = ifelse(is.na(s),0, s)
  r = ifelse(is.na(r),0, r)
  cum = sum(theta_hat*cumprod(1-c(0,theta_hat[1:M-1])))
  inner = (r-s)/r; inner = inner/s
  varp = ((1-cum)^2) * sum(inner)
  se = sqrt(varp)
  return (list(risk_ds=cum, se_ds = se))
}  


round_1_pa = function(m1_pa, m2_pa){ 
  function(i, h){
    if (h==1){
      if (i==0)
        return(1/(1+exp(m1_pa[1])+exp(m1_pa[2])))
      if(i==1)
        return(exp(m1_pa[1])/(1+exp(m1_pa[1])+exp(m1_pa[2])))
      if (i==2)
        return(exp(m1_pa[2])/(1+exp(m1_pa[1])+exp(m1_pa[2])))
    } else{
      if(i==1)
        return(expit(m2_pa[1]+m2_pa[2]*h))
      if(i==0)
        return(1-expit(m2_pa[1]+m2_pa[2]*h))
      if(i==2)
        return(0)
    }
  }
}


Pmat_pa = function(m1_pa , m2_pa){ 
  function(h){
    matrix(c(1-expit(m1_pa[1]+m1_pa[2]*h), 1-expit(m2_pa[1]+m2_pa[2]*h),0,
             expit(m1_pa[1]+m1_pa[2]*h), expit(m2_pa[1]+m2_pa[2]*h), 0,
             0, 0, 1), ncol=3)}
}


Ppmat_pa = function(m1_pa, m2_pa){
  function(h){
    matrix(c(1/(1+exp(m1_pa[1]+m1_pa[3]*h)+exp(m1_pa[2]+m1_pa[4]*h)), 1/(1+exp(m2_pa[1]+m2_pa[3]*h)+exp(m2_pa[2]+m2_pa[4]*h)), 0,
             exp(m1_pa[1]+m1_pa[3]*h)/(1+exp(m1_pa[1]+m1_pa[3]*h)+exp(m1_pa[2]+m1_pa[4]*h)),exp(m2_pa[1]+m2_pa[3]*h)/(1+exp(m2_pa[1]+m2_pa[3]*h)+exp(m2_pa[2]+m2_pa[4]*h)), 0,
             exp(m1_pa[2]+m1_pa[4]*h)/(1+exp(m1_pa[1]+m1_pa[3]*h)+exp(m1_pa[2]+m1_pa[4]*h)), exp(m2_pa[2]+m2_pa[4]*h)/(1+exp(m2_pa[1]+m2_pa[3]*h)+exp(m2_pa[2]+m2_pa[4]*h)), 1), ncol=3)
  }
}


pop.ave.var = function(S, fpl.cens, comp, delta, M){
  sig = vector(length = M)
  N = max(S)*(max(S)+3)/2-1
  # sigma_1
  eta = vector(length = N)
  k = 1
  for (i in 1:max(S)){
    for (j in i:max(S)){
      eta[k] = mean(S == j & fpl.cens == i & delta == 1 & comp>=i)
      k = k+1
    }	
  }
  for (i in 1:(max(S)-1)){
    eta[k+i-1] = mean(S == i & delta == 0) 
  }
  
  sigma = matrix(ncol = N, nrow = N)
  sigma = -outer(eta,eta)
  diag(sigma) = eta*(1-eta)
  
  sig[1] = t(c(rep(1,max(S)),rep(0,N-max(S))))%*%sigma%*%c(rep(1,max(S)),rep(0,N-max(S)))/length(S)
  
  # sigma_2 to sigma_K
  mu0 = vector(length = M)
  for (i in 1:M){
    mu0[i] = mean(S == i & delta == 0)
  }
  mu = vector("list", max(S))
  for (i in 1:max(S)){
    for (j in i:max(S)){
      mu[[i]][j-i+1] = mean(S == j & fpl.cens == i & delta == 1 & comp>=i)
    }	
  }
  
  for (i in 2:max(S)){
    k = 1
    a = rep(0,N)
    for (j in 1:i){
      a[k:(k+max(S)-j)]  = (1 + (1-i/j)*mu0[j]^(i/j)*(mu0[j]+rep(1,(max(S)-j+1))%*%mu[[j]])^(-i/j))%*%rep(1,(max(S)-j+1))
      k = k+max(S)-j+1
    }
    a[(length(a)-(max(S)-2)):length(a)] = c(rep(1,(i-1)),rep(0,(max(S)-i)))
    for (j in 1:(i-1)){
      a[(length(a)-(max(S)-2)+j-1)] = a[(length(a)-(max(S)-2)+j-1)] +
        (mu0[j]^(i/j-1)*(mu0[j]+sum(mu[[j]]))^(-i/j)*(mu0[j]+i/j*sum(mu[[j]])))
    }
    sig[i] = a%*%sigma%*%a/length(S)				
  }
  
  return(varp = sig[M])
}


createCartesian1 = function(l, j){
  L = list()
  for (i in 1:(j-1))
    L[[i]] = c(0,1)
  if (j ==1){
    df = data.frame('i1' = 1)
  } else{
    df = expand.grid(L)
    for (i in 1:(j-1))
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var==l-1)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==j-ind)})
    df = subset(df, valid==1)
    df$valid = NULL
    df[, paste0('i',j)] = 1}
  return (df)
}


f2 = function(x, h, functions2){
  round_1 = functions2$round_1
  Pmat = functions2$Pmat
  Ppmat = functions2$Ppmat
  x = as.numeric(x)
  c = round_1(x[1], h)*x[length(x)]
  for (k in 1:(length(x)-2)){
    if (k+1 == h){
      c = c*Ppmat(h,k+1)[x[k]+1, x[k+1]+1]
    } else{
      c = c*Pmat(h,k+1)[x[k]+1, x[k+1]+1]
    }
  }
  return (c)
}


f3 = function(x, h, functions1){
  round_1_pa = functions1$round_1_pa
  Pmat_pa = functions1$Pmat_pa
  Ppmat_pa = functions1$Ppmat_pa
  x = as.numeric(x)
  c = round_1_pa(x[1], h)*x[length(x)]
  for (k in 1:(length(x)-2)){
    if (k+1 >= h){
      c = c*Ppmat_pa(h)[x[k]+1, x[k+1]+1]
    } else{
      c = c*Pmat_pa(h)[x[k]+1, x[k+1]+1]
    }
  }
  return (c)
}


pop.ave = function(data,fpl.cens, functions1, functions2,l,M){
  p1 = summary(factor(data$S[!duplicated(data$StudyID_c)],levels=c(1:M)))/length(data$S[!duplicated(data$StudyID_c)]) 
  P = 0 # initial probability
  for (j in l:M)
    for (h in j:M){
      df = createCartesian1(l,j)
      df[, 'p1'] = p1[h]
      df$term = apply(df,1, f2, h=h, functions2=functions2)
      P = P + sum(df$term)
    }
  for (j in max(2,l):M)
    for (h in 1:(j-1)){
      df = createCartesian1(l,j)
      df[, 'p1'] = p1[h]
      df$term = apply(df,1, f3, h=h, functions1=functions1)
      P = P + sum(df$term)
    }
  varp = pop.ave.var(data$S, fpl.cens, data$comp, data$delta, M)
  se = sqrt(varp)
  return(list(risk_pa = P, se_pa = se))
}








