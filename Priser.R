#J'ai mis les fonctions les plus importante tout ce que j'ai utilisé n'est pas ici mais sur le fichier .pdf je vous conseille de lire le pdf seulement (tout y esr sauf les commandes de tracé de courbes)





set.seed(999) #initialisation de l'aléa
n=50
x=rnorm(n)
a=5
b=1
sigma2=2

epsilon=sqrt(2)*rnorm(n)

Z = a + b*x + epsilon
f<-function(x,lambda=0.3){
    if(x>-1/lambda){
        y=(lambda*x+1)^(1/lambda)
    }else{
        y=-(-lambda*x-1)^(1/lambda)
    }
    return(y)}

Yobs=sapply(Z,f) 
X=cbind(1,x)
h<-function(lambda){
    return(function(y){lambda
        if(lambda==0){
            return(log(y)) 
        }else{
            if(y>0){
                return((y^lambda - 1)/lambda)
            }else{
                return((-1-(-y)^lambda)/lambda)
            }
        }
        if(y==0 && lambda <= 0){
            print("error")
        }
    })
}

Q=solve(t(X)%*%X)%*%t(X)
P=diag(1,n) - X%*%Q 



theta.est<-function(lambda,Y=Yobs){
    Z=sapply(Y,h(lambda))
    return(Q%*%Z)
}
sigma2.est<-function(lambda,Y=Yobs){
    Z=sapply(Y,h(lambda))
    return(t(Z)%*%P%*%Z/n)
}

lmin<-function(lambda,Y=Yobs){
    return(n/2*log(sigma2.est(lambda,Y))-(lambda-1)*sum(log(abs(Y)))+ n/2*(log(2*pi)+1))
}
lminbis<-function(x,Y=Yobs) {  #pour pouvoir être pris en argument dans la conction curve
    return(sapply(x,lmin,Y=Y))
}

zetaObs<-function(lambda0,lambda.est,Y){
    return(2*(lmin(lambda0,Y)-lmin(lambda.est,Y)))
}


#curve(lminbis(x,Yobs),0.01,2)




resopt=nlm(lmin,p=0.3,hessian=T)
lambda.est=resopt$estimate





nivTRV<-function(k){
    niv=0
    for (i in 1:k){
        epsilon=sqrt(2)*rnorm(n)
        Z = a + b*x + epsilon
        Y=sapply(Z,f)
        resopt=nlm(lmin,p=0.3,hessian=T,Y=Y)
        lambda.est=resopt$estimate
        if(1-pchisq(zetaObs(0.3,lambda.est,Y),1)<=0.05){
            niv=niv+1
        }     
    } 

    return(niv/k)
}


pw<-function(lambda1,k=1000){
    niv=0
    for (i in 1:k){
        epsilon=sqrt(2)*rnorm(n)
        Z = a + b*x + epsilon
        Y=sapply(Z,f,lambda=lambda1)
        resopt=nlm(lmin,p=lambda1,hessian=T,Y=Y)
        lambda.est=resopt$estimate
        if(1-pchisq(zetaObs(0,lambda.est,Y),1)<0.05){
            niv=niv+1
        }     
    } 
    return(niv/k)
}

pwbis<-function(lambda1){ #pour être mis en argument de curve
    return(sapply(lambda1,pw))
}

#curve(pwbis,0.1,0.3,n=50) pour afficher la courbe de puissance

df<- read.table("NbCycleRupture.csv",header=TRUE,sep=";")


res1=lm(y~.,data=df)
res2=lm(y~x1+x2+x3+I(x1^2)+I(x2^2)+I(x3^2)+I(x1*x2)+I(x1*x3)+I(x2*x3),data=df)

res1bis=lm(log(y)~.,data=df)
res2bis=lm(log(y)~x1+x2+x3+I(x1^2)+I(x2^2)+I(x3^2)+I(x1*x2)+I(x1*x3)+I(x2*x3),data=df)