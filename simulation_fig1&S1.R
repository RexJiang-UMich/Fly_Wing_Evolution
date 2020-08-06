library(clusterGeneration)
library(mvtnorm)
library(MBESS)
library(resample)

Ntest=100
weights=(0:10)/10
ndim=24

cor_old=matrix(0,nrow=Ntest,ncol=length(weights))
slope_old=matrix(0,nrow=Ntest,ncol=length(weights))
cor_new=matrix(0,nrow=Ntest,ncol=length(weights))
slope_new=matrix(0,nrow=Ntest,ncol=length(weights))
cor_old_error=matrix(0,nrow=Ntest,ncol=length(weights))
slope_old_error=matrix(0,nrow=Ntest,ncol=length(weights))
cor_new_error=matrix(0,nrow=Ntest,ncol=length(weights))
slope_new_error=matrix(0,nrow=Ntest,ncol=length(weights))
for(n in 1:Ntest){
	mcor=rcorrmatrix(ndim)
	mvar=rgamma(ndim,shape=0.5,scale=1)
	mmat=cor2cov(mcor,sd=sqrt(mvar))
	fcor=rcorrmatrix(ndim)
	fvar=rgamma(ndim,shape=0.5,scale=1)
	fmat=cor2cov(fcor,sd=sqrt(fvar))

	for(i in 1:length(weights)){
		rmat=(1-weights[i])*mmat+weights[i]*fmat
		
		mobs=cov(rmvnorm(100,sigma=mmat))
		robs=cov(rmvnorm(100,sigma=rmat))
		robs=robs*(sum(diag(mobs))/sum(diag(robs)))
		
		#without error
		h=(mmat+rmat)/2
		kold=eigen(h)$vectors
		xold=log10(diag((t(kold)%*%mmat)%*%kold))
		yold=log10(diag((t(kold)%*%rmat)%*%kold))
		cor_old[n,i]=cor(xold,yold)
		slope_old[n,i]=lm(yold~xold)[[1]][2]

		knew=eigen(mmat)$vectors
		xtrue=log10(diag((t(knew)%*%mmat)%*%knew))
		ytrue=log10(diag((t(knew)%*%rmat)%*%knew))
		cor_new[n,i]=cor(xtrue,ytrue)
		slope_new[n,i]=lm(ytrue~xtrue)[[1]][2]
		
		#with error
		h=(mobs+robs)/2
		kold=eigen(h)$vectors
		xold=log10(diag((t(kold)%*%mobs)%*%kold))
		yold=log10(diag((t(kold)%*%robs)%*%kold))
		cor_old_error[n,i]=cor(xold,yold)
		slope_old_error[n,i]=lm(yold~xold)[[1]][2]
		
		knew=eigen(mobs)$vectors
		xnew=log10(diag((t(knew)%*%mobs)%*%knew))
		ynew=log10(diag((t(knew)%*%robs)%*%knew))
		cor_new_error[n,i]=cor(xnew,ynew)
		slope_new_error[n,i]=lm(ynew~xnew)[[1]][2]
	}
}
out=data.frame(colMeans(slope_new),sqrt(colVars(slope_new)),colMeans(slope_old),sqrt(colVars(slope_old)),colMeans(cor_new),sqrt(colVars(cor_new)),colMeans(cor_old),sqrt(colVars(cor_old)))
out_error=data.frame(colMeans(slope_new_error),sqrt(colVars(slope_new_error)),colMeans(slope_old_error),sqrt(colVars(slope_old_error)),colMeans(cor_new_error),sqrt(colVars(cor_new_error)),colMeans(cor_old_error),sqrt(colVars(cor_old_error)))

write.table(out,file="matrix_sim_out.txt",sep="\t")
write.table(out_error,file="matrix_sim_out_error.txt",sep="\t")

write.table(cor_old,file="cor_old.txt",sep="\t")
write.table(slope_old,file="slope_old.txt",sep="\t")
write.table(cor_new,file="cor_new.txt",sep="\t")
write.table(slope_new,file="slope_new.txt",sep="\t")
write.table(cor_old_error,file="cor_old_error.txt",sep="\t")
write.table(slope_old_error,file="slope_old_error.txt",sep="\t")
write.table(cor_new_error,file="cor_new_error.txt",sep="\t")
write.table(slope_new_error,file="slope_new_error.txt",sep="\t")
