library(ggplot2);library(ggExtra)

labx=expression(paste("Estimate ","of ","the ","new ","method"))
laby=expression(paste("Estimate ","of ","Houle ",italic(et)," ",italic(al.)," 's ","method"))

#plot all data points
d1<-read.table("slope_new.txt",sep="\t")
d2<-read.table("slope_old.txt",sep="\t")
d3<-read.table("cor_new.txt",sep="\t")
d4<-read.table("cor_old.txt",sep="\t")
d=matrix(0,nrow=1,ncol=5)
for(i in 1:100){
	for(j in 1:11){
		d=rbind(d,c((j-1)/10,d1[i,j],d2[i,j],d3[i,j],d4[i,j]))
	}
}
d=d[2:nrow(d),]
colnames(d)=c("w","slope_new","slope_old","cor_new","cor_old")
d=data.frame(d)
d$w=factor(d$w,levels=(0:10)/10)

#only plot when M and R are independent
d0=d[which(d$w==1),]
g<-ggplot(d0,aes(x=slope_new,y=slope_old))+
geom_point(size=2.5)+
geom_abline(slope=0,intercept=0,linetype="dashed")+
geom_vline(xintercept=0,linetype="dashed")+
ylim(-1,1.3)+
xlab(labx)+ylab(laby)+
theme_classic()+
ggtitle("Regression slope")+
theme(plot.title=element_text(hjust=0.5),text=element_text(size=20),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),panel.border=element_rect(colour="black",fill=NA))
g<-ggMarginal(g,type="histogram")

g<-ggplot(d0,aes(x=cor_new,y=cor_old))+
geom_point(size=2.5)+
geom_abline(slope=0,intercept=0,linetype="dashed")+
geom_vline(xintercept=0,linetype="dashed")+
ylim(-1,1)+
xlab(labx)+ylab(laby)+
theme_classic()+
ggtitle("Correlation coefficient")+
theme(plot.title=element_text(hjust=0.5),text=element_text(size=20),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),panel.border=element_rect(colour="black",fill=NA))
g<-ggMarginal(g,type="histogram")

#only plot mean and sd
#plot against "true" values
d1<-read.table("matrix_sim_out.txt",sep="\t")
d2<-read.table("matrix_sim_out_error.txt",sep="\t")

d=data.frame(d1[,1:2],d2[,1:4])
colnames(d)=c("slope_true","slope_true_sd","slope_new","slope_new_sd","slope_old","slope_old_sd")

g<-ggplot(d)+
geom_point(aes(x=slope_true,y=slope_old),shape=23,fill="red",size=2)+
geom_errorbar(aes(x=slope_true,ymin=slope_old-slope_old_sd,ymax=slope_old+slope_old_sd),width=0.01)+
geom_errorbarh(aes(y=slope_old,xmin=slope_true-slope_true_sd,xmax=slope_true+slope_true_sd),height=0.01)+
geom_smooth(aes(x=slope_true,y=slope_old),method=lm,color="red",se=FALSE)+

geom_point(aes(x=slope_true,y=slope_new),shape=23,fill="blue",size=2)+
geom_errorbar(aes(x=slope_true,ymin=slope_new-slope_new_sd,ymax=slope_new+slope_new_sd),width=0.01)+
geom_errorbarh(aes(y=slope_new,xmin=slope_true-slope_true_sd,xmax=slope_true+slope_true_sd),height=0.01)+
geom_smooth(aes(x=slope_true,y=slope_new),method=lm,color="blue",se=FALSE)+

geom_abline(slope=1,intercept=0,linetype="dashed")+
theme_classic()+
xlab("True slope")+ylab("Estimated slope")+
theme(text=element_text(size=20),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),panel.border=element_rect(colour="black",fill=NA))

d=data.frame(d1[,5:6],d2[,5:8])
colnames(d)=c("cor_true","cor_true_sd","cor_new","cor_new_sd","cor_old","cor_old_sd")
g<-ggplot(d)+
geom_point(aes(x=cor_true,y=cor_old),shape=23,fill="red",size=2)+
geom_errorbar(aes(x=cor_true,ymin=cor_old-cor_old_sd,ymax=cor_old+cor_old_sd),width=0.01)+
geom_errorbarh(aes(y=cor_old,xmin=cor_true-cor_true_sd,xmax=cor_true+cor_true_sd),height=0.01)+
geom_smooth(aes(x=cor_true,y=cor_old),method=lm,color="red",se=FALSE)+

geom_point(aes(x=cor_true,y=cor_new),shape=23,fill="blue",size=2)+
geom_errorbar(aes(x=cor_true,ymin=cor_new-cor_new_sd,ymax=cor_new+cor_new_sd),width=0.01)+
geom_errorbarh(aes(y=cor_new,xmin=cor_true-cor_true_sd,xmax=cor_true+cor_true_sd),height=0.01)+
geom_smooth(aes(x=cor_new,y=cor_true),method=lm,color="blue",se=FALSE)+

geom_abline(slope=1,intercept=0,linetype="dashed")+
theme_classic()+

xlab("True correlation coefficient")+ylab("Estimated correlation coefficient")+
theme(text=element_text(size=20),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),panel.border=element_rect(colour="black",fill=NA))



