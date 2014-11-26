

fcast=function(data, method="GMDH", input=4, layer=3, f.number=10, tf="all", plotit=TRUE){


if (tf=="all"){tf_options=c(101:104)
}else if (tf=="polynomial"){tf_options=c(101)
}else if (tf=="sigmoid"){tf_options=c(102)
}else if (tf=="RBF"){tf_options=c(103)
}else if (tf=="tangent"){tf_options=c(104)

}else {stop("Transfer function you entered is not available")}


transf=function(h,dataaa){

if (h==101) dat=dataaa
if (h==102) dat=log(dataaa/(1-dataaa))
if (h==103) dat=sqrt(-log(dataaa))
if (h==104) dat=atan(dataaa)/pi*180

dat

}


back_transf=function(h,dataaa){

if (h==101) dat=dataaa
if (h==102) dat=1/(1+exp(-dataaa))
if (h==103) dat=exp(-dataaa^2)
if (h==104) dat=tan(dataaa*pi/180)

dat

}


if (min(data)<=0){
stt1=abs(min(data))+1  
}else{stt1=0}

stt2=max(data+stt1)+1  


y=(data+stt1)/stt2 


if (method=="GMDH"){


store_Astore<- list()
store_z=list()

ss=length(y)

threshold=c(rep(input,layer-1),1)
nnode=input*(input-1)/2
idn=c(1:input)


yt=y[-input:-1]
x=NULL
for (i in 1:(input-1)){
x=cbind(x,matrix(y[c(-1:-(input-i),-ss:-(ss+1-i))]))
}
x=cbind(x,matrix(y[c(-ss:-(ss-input+1))]))


for (k in 1:layer){

w=t(combn(order(idn), 2))


Astore=NULL
z=NULL

for (j in 1:nnode){

qq=cbind(1,x[,w[j,]],x[,w[j,]][,1]*x[,w[j,]][,2],x[,w[j,]]^2)


tfunc=NULL
tfunc_z=NULL
for (g in tf_options){

est_coef=t(ginv(t(qq)%*%qq)%*%t(qq)%*%transf(g,yt))

ee=as.numeric(est_coef)

est_zt=rowSums(t(ee*t(qq)))

tfunc=rbind(tfunc,c(est_coef,mean((back_transf(g,est_zt)-yt)^2),g))

tfunc_z=cbind(tfunc_z,matrix(back_transf(g,est_zt)))
}


z=cbind(z,tfunc_z[,which.min(tfunc[,7])])
Astore=rbind(Astore,tfunc[which.min(tfunc[,7]),])

}

Astore=cbind(Astore,c(1:nnode))

store_Astore[[k]]=Astore[which(Astore[,7]<=sort(Astore[,7])[threshold[k]]),]
store_z[[k]]=z[,which(Astore[,7]<=sort(Astore[,7])[threshold[k]])]

x=store_z[[k]]

if (k==layer){
store_Astore[[k]]=matrix(store_Astore[[k]],nrow=1)
store_z[[k]]=matrix(store_z[[k]],ncol=1)

}


}


for (h in 1:f.number){

yt_input=matrix(rev(tail(y,input)),nrow=1)
idn2=c(1:input)
w2=t(combn(order(idn2), 2))



for (k2 in 1:layer){


selected_coef=selected_qq2=NULL
store_qq2=NULL

for (j2 in 1:nnode){

qq2=c(1,yt_input[,w2[j2,]],yt_input[,w2[j2,]][1]*yt_input[,w2[j2,]][2],yt_input[,w2[j2,]]^2)

store_qq2=rbind(store_qq2,qq2)
}

selected_qq2=store_qq2[store_Astore[[k2]][,9],]
selected_coef=store_Astore[[k2]][,1:6]


if (k2==layer){
selected_qq2=matrix(selected_qq2,nrow=1)
selected_coef=matrix(selected_coef,nrow=1)
}

yt_input=matrix(rowSums(selected_qq2*selected_coef),nrow=1)


for (k5 in 1:threshold[k2]){

yt_input[1,k5]=back_transf(store_Astore[[k2]][k5,8],yt_input[1,k5])
}



}

y=c(y,yt_input)
}

fitted=store_z[[layer]][,1]*stt2-stt1


} 

if (method=="RGMDH"){



store_Astore<- list()
store_z=list()
store_Astore2<- list()
store_z2=list()

ss=length(y)

threshold=c(rep(input,layer-1),1)
nnode=input*(input-1)/2+input
p=input*(input-1)/2
idn=c(1:input)


yt=y[-input:-1]
x=NULL
for (i in 1:(input-1)){
x=cbind(x,matrix(y[c(-1:-(input-i),-ss:-(ss+1-i))]))
}
x=cbind(x,matrix(y[c(-ss:-(ss-input+1))]))


for (k in 1:layer){

w=t(combn(order(idn), 2))

m2 <- matrix(rep(1:input,input),input,input,byrow=T)
m2[upper.tri(m2)] <-0


Astore=NULL
z=NULL
z2=NULL
Astore2=NULL


for (j in 1:nnode){

if (j<=p){
qq=cbind(1,x[,w[j,]],x[,w[j,]][,1]*x[,w[j,]][,2],x[,w[j,]]^2)
}else{
qq=cbind(x[,m2[j-p,]])
}


tfunc=NULL
tfunc_z=NULL

tfunc2=NULL
tfunc_z2=NULL

if (j<=p){

for (g in tf_options){

est_coef=t(ginv(t(qq)%*%qq)%*%t(qq)%*%transf(g,yt))

ee=as.numeric(est_coef)

est_zt=rowSums(t(ee*t(qq)))

tfunc=rbind(tfunc,c(est_coef,mean((back_transf(g,est_zt)-yt)^2),g))

tfunc_z=cbind(tfunc_z,matrix(back_transf(g,est_zt)))
}
}else{

for (g in tf_options){

est_coef=t(ginv(t(qq)%*%qq)%*%t(qq)%*%transf(g,yt))

coef=c(est_coef,rep(0,input-length(est_coef)))

ee=as.numeric(est_coef)

est_zt=rowSums(t(ee*t(qq)))

tfunc2=rbind(tfunc2,c(coef,mean((back_transf(g,est_zt)-yt)^2),g))

tfunc_z2=cbind(tfunc_z2,matrix(back_transf(g,est_zt)))
}

}



if (j<=p){
z=cbind(z,tfunc_z[,which.min(tfunc[,7])])
Astore=rbind(Astore,tfunc[which.min(tfunc[,7]),])
}else{
z2=cbind(z2,tfunc_z2[,which.min(tfunc2[,(j-p+1)])])
Astore2=rbind(Astore2,tfunc2[which.min(tfunc2[,(j-p+1)]),])
}



}


Astore=cbind(Astore,c(1:p))
Astore2=cbind(Astore2,c((p+1):(p+input)))

checkk=rbind(Astore[,c(7,9)],Astore2[,c((input+1),(input+3))])

ord=which(checkk[,1]<=sort(checkk[,1])[threshold[k]])
ord1=ord[which(ord<=p)]
ord2=ord[which(ord>p)]


store_Astore[[k]]=Astore[ord1,]
store_Astore2[[k]]=Astore2[ord2-p,]

store_z[[k]]=z[,ord1]
store_z2[[k]]=z2[,ord2-p]

x=cbind(store_z[[k]],store_z2[[k]])


if (class(store_Astore[[k]])!="matrix"){
store_Astore[[k]]=matrix(store_Astore[[k]],nrow=1)
store_z[[k]]=matrix(store_z[[k]],ncol=1)
}

if (class(store_Astore2[[k]])!="matrix"){
store_Astore2[[k]]=matrix(store_Astore2[[k]],nrow=1)
store_z2[[k]]=matrix(store_z2[[k]],ncol=1)
}


}



for (h in 1:f.number){

yt_input=matrix(rev(tail(y,input)),nrow=1)
idn2=c(1:input)
w2=t(combn(order(idn2), 2))



for (k2 in 1:layer){


selected_coef=selected_qq2=NULL
store_qq2=NULL

selected_coef5=selected_qq5=NULL
store_qq5=NULL


for (j2 in 1:nnode){

if (j2<=p){
qq2=c(1,yt_input[,w2[j2,]],yt_input[,w2[j2,]][1]*yt_input[,w2[j2,]][2],yt_input[,w2[j2,]]^2)
store_qq2=rbind(store_qq2,qq2)
}else{
qq2=c(yt_input[,m2[j2-p,]],rep(0,input-(j2-p)))
store_qq5=rbind(store_qq5,qq2)
}

}




selected_qq2=store_qq2[store_Astore[[k2]][,9],]
selected_coef=store_Astore[[k2]][,1:6]

selected_qq5=store_qq5[store_Astore2[[k2]][,(input+3)]-p,]
selected_coef5=store_Astore2[[k2]][,1:input]



if (class(selected_qq2)!="matrix"){
selected_qq2=matrix(selected_qq2,nrow=1)
selected_coef=matrix(selected_coef,nrow=1)
}

if (class(selected_qq5)!="matrix"){
selected_qq5=matrix(selected_qq5,nrow=1)
selected_coef5=matrix(selected_coef5,nrow=1)
}


uu1=matrix(rowSums(selected_qq2*selected_coef),nrow=1)
uu2=matrix(rowSums(selected_qq5*selected_coef5),nrow=1)

d1=dim(matrix(rowSums(selected_qq2*selected_coef),nrow=1))[2]
d2=dim(matrix(rowSums(selected_qq5*selected_coef5),nrow=1))[2]

if(d1!=0){
for (k5 in 1:d1){
uu1[1,k5]=back_transf(store_Astore[[k2]][k5,8],uu1[1,k5])
}
}

if(d2!=0){
for (k5 in 1:d2){
uu2[1,k5]=back_transf(store_Astore2[[k2]][k5,(input+2)],uu2[1,k5])
}
}

yt_input=cbind(uu1,uu2)


}

y=c(y,yt_input)
}

fitted=cbind(store_z[[layer]],store_z2[[layer]])[,1]*stt2-stt1

}  

forecast_values=tail(y*stt2-stt1,f.number)


MSE=(mean((fitted-data[c(-1:-input)])^2))


if(plotit==TRUE){
plot(ts(c(data,forecast_values)),col="black",ylab="Time Series")  
abline(v=ss,lty=2)
}


out=list()
out$fitted=ts(fitted, start=input+1,end=ss)
out$MSE=MSE
out$forecasts=ts(forecast_values, start=ss+1,end=ss+f.number)

invisible(out)


}







