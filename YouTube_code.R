#Data prep
library("evd")
library("condmixt")
library("ismev")
library("mev")
library("e1071")

x1=views$Canada
x2=views$Germany
x3=views$France
x4=views$UK
x5=views$Japan
x6=views$USA
x1=na.omit(x1)
x2=na.omit(x2)
x3=na.omit(x3)
x4=na.omit(x4)
x5=na.omit(x5)
x6=na.omit(x6)

kur=matrix(NA,6,1)
kur[1,1]=kurtosis(x1)
kur[2,1]=kurtosis(x2)
kur[3,1]=kurtosis(x3)
kur[4,1]=kurtosis(x4)
kur[5,1]=kurtosis(x5)
kur[6,1]=kurtosis(x6)

x1_d=rank(x1)/length(x1)
x2_d=rank(x2)/length(x2)
x3_d=rank(x3)/length(x3)
x4_d=rank(x4)/length(x4)
x5_d=rank(x5)/length(x5)
x6_d=rank(x6)/length(x6)

#plots
plot(x1,x1_d,xlab="Canada",ylab="Density")
plot(x2,x2_d,xlab="Germany",ylab="Density")
plot(x3,x3_d,xlab="France",ylab="Density")
plot(x4,x4_d,xlab="UK",ylab="Density")
plot(x5,x5_d,xlab="Japan",ylab="Density")
plot(x6,x6_d,xlab="USA",ylab="Density")

#EVI estimation
x1_sort=sort(x1)
x2_sort=sort(x2)
x3_sort=sort(x3)
x4_sort=sort(x4)
x5_sort=sort(x5)
x6_sort=sort(x6)

k1=length(x1)*0.1
ml1=gp.fit(x1,x1_sort[length(x1)-k1],method="zhang")
gamma_1=ml1$approx.mean[2]


k2=length(x2)*0.1
ml2=gp.fit(x1,x2_sort[length(x2)-k2],method="zhang")
gamma_2=ml2$approx.mean[2]


k3=length(x3)*0.1
ml3=gp.fit(x3,x3_sort[length(x3)-k3],method="zhang")
gamma_3=ml3$approx.mean[2]

k4=length(x4)*0.1
ml4=gp.fit(x4,x4_sort[length(x4)-k4],method="zhang")
gamma_4=ml4$approx.mean[2]

k5=length(x5)*0.1
ml5=gp.fit(x5,x5_sort[length(x5)-k5],method="zhang")
gamma_5=ml5$approx.mean[2]

k6=length(x6)*0.1
ml6=gp.fit(x6,x6_sort[length(x6)-k6],method="zhang")
gamma_6=ml6$approx.mean[2]

#expected shortfall 
ind1=matrix(0,length(x1),1)
for(i in 1:(length(x1))){if(x1[i]>x1_sort[length(x1)-k1]){ind1[i]=1}}
ind_x1=ind1*x1
q1=x1_sort[length(x1)-k1]*((k1/(length(x1)*(1-0.99)))^gamma_1)
ES1=(q1/x1_sort[length(x1)-k1])*(1/k1)*sum(ind_x1)


ind2=matrix(0,length(x2),1)
for(i in 1:(length(x2))){if(x2[i]>x2_sort[length(x2)-k2]){ind2[i]=1}}
ind_x2=ind2*x2
q2=x2_sort[length(x2)-k2]*((k2/(length(x2)*(1-0.99)))^gamma_2)
ES2=(q2/x2_sort[length(x2)-k2])*(1/k2)*sum(ind_x2)


ind3=matrix(0,length(x3),1)
for(i in 1:(length(x3))){if(x3[i]>x3_sort[length(x3)-k3]){ind3[i]=1}}
ind_x3=ind3*x3
q3=x3_sort[length(x3)-k3]*((k3/(length(x3)*(1-0.99)))^gamma_3)
ES3=(q3/x3_sort[length(x3)-k3])*(1/k3)*sum(ind_x3)


ind4=matrix(0,length(x4),1)
for(i in 1:(length(x4))){if(x4[i]>x4_sort[length(x4)-k4]){ind4[i]=1}}
ind_x4=ind4*x4
q4=x4_sort[length(x4)-k4]*((k4/(length(x4)*(1-0.99)))^gamma_4)
ES4=(q4/x4_sort[length(x4)-k4])*(1/k4)*sum(ind_x4)

ind5=matrix(0,length(x5),1)
for(i in 1:(length(x5))){if(x5[i]>x5_sort[length(x5)-k5]){ind5[i]=1}}
ind_x5=ind5*x5
q5=x5_sort[length(x5)-k5]*((k5/(length(x5)*(1-0.999)))^gamma_5)
ES5=(q5/x5_sort[length(x5)-k5])*(1/k5)*sum(ind_x5)

ind6=matrix(0,length(x6),1)
for(i in 1:(length(x6))){if(x6[i]>x6_sort[length(x6)-k6]){ind6[i]=1}}
ind_x6=ind6*x6
q6=x6_sort[length(x6)-k6]*((k6/(length(x6)*(1-0.99)))^gamma_6)
ES6=(q6/x6_sort[length(x6)-k6])*(1/k6)*sum(ind_x6)




