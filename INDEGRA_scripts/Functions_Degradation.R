########## Max a posteriori


EPosterior<-function(m,ninf,a,b,alpha)
{
	k=0:ninf
	Ak=lchoose(ninf,k) + k*log(alpha) + (ninf-k)*log(1-alpha) + lbeta(a+ninf-k+1, m+b)-lbeta(a,b)
	L=max(Ak)+  log(sum(exp(Ak-max(Ak))))
	return(L)

}


lT1H1<-function(m,ninf,a,b,alpha)
{
	k=0:ninf
	Ak=lchoose(ninf,k) + k*log(alpha) + (ninf-k)*log(1-alpha) + lbeta(a+ninf-k, m+b)-lbeta(a,b)
	L=max(Ak)+  log(sum(exp(Ak-max(Ak))))
	return(L)
}



lT1H0<-function(m1,ninf1,m2,ninf2,a,b,alpha1,alpha2)
{
	k1=0:ninf1
	k2=0:ninf2
	A=unlist(lapply(k1,term_likelihood, k2=k2, ninf1=ninf1,ninf2=ninf2,m1=m1,m2=m2,a=a,b=b,alpha1=alpha1,alpha2=alpha2))
	L=max(A)+  log(sum(exp(A-max(A))))
	return(L)
}

term_likelihood<-function(k1,k2,ninf1,ninf2,m1,m2,a,b,alpha1,alpha2)
{
	A=lchoose(ninf2,k2) + lchoose(ninf1,k1) + k2*log(alpha2)  + k1*log(alpha1) + (ninf1-k1)*log(1-alpha1) +(ninf2-k2)*log(1-alpha2)  + lbeta(a+ninf1-k1 +ninf2-k2, m1+m2+b)-lbeta(a,b)
	return(A)
}


TranscriptTest_Tech<-function(Df, C1,C2, a,b,thresh)
{
	L0=NA; L11=NA; L12=NA; EP1=NA; EP2=NA; Subsample1=NA; Subsample2=NA;
	if (Df$K1)
	{
		RL1=as.numeric(Df$TotalRLsample1)
		FL1=as.numeric(Df$FLsample1)
		Subsample1=FALSE
		EP1=EPosterior(RL1,FL1,a,b,C1)
		L11= lT1H1(RL1,FL1,a,b,C1)
	}	else if(!is.na(Df$ReadCount_sample1))
	{
		Norm=thresh/as.numeric(Df$ReadCount_sample1)
		RL1=round(as.numeric(Df$TotalRLsample1)*Norm)
		FL1=round(as.numeric(Df$FLsample1)*Norm)
		Subsample1=TRUE
		EP1=EPosterior(RL1,FL1,a,b,C1)
		L11= lT1H1(RL1,FL1,a,b,C1)
	}
	if (Df$K2)
	{
		RL2=as.numeric(Df$TotalRLsample2)
		FL2=as.numeric(Df$FLsample2)
		Subsample2=FALSE
		EP2=EPosterior(RL2,FL2,a,b,C2)
		L12= lT1H1(RL2,FL2,a,b,C2)
	} else if (!is.na(Df$ReadCount_sample2))
	{
		Norm=thresh/as.numeric(Df$ReadCount_sample2)
		RL2=round(as.numeric(Df$TotalRLsample2)*Norm)
		FL2=round(as.numeric(Df$FLsample2)*Norm)
		Subsample2=TRUE
		EP2=EPosterior(RL2,FL2,a,b,C2)
		L12= lT1H1(RL2,FL2,a,b,C2)		
	}
	if (!is.na(Subsample1) && !is.na(Subsample2))
		L0 = lT1H0(RL1,FL1,RL2,FL2,a,b,C1,C2)
	R=data.frame(L1_samp1=L11,L1_samp21=L12, L0=L0, Lratio=L11+L12-L0, EP1=exp(EP1-L11),EP2=exp(EP2-L12),Subsample1=Subsample1,Subsample2=Subsample2)
	return(R)
}


Test_Degradation<-function(process1,process2,p=0.1,thresh=2000)
{

	a=2
	b=1000


	samp1=read.csv(process1,sep="\t")
	samp2=read.csv(process2,sep="\t")
	EObs=c(mean(samp1$fragmentation_rate),mean(samp2$fragmentation_rate))
	Et=min(EObs)-0.0001
	alpha=1-(1-EObs)/(1-Et)
	print(alpha)

	Frag1<-cbind(samp1$transcript,samp1$read_count, samp1$fragmentation_rate, samp1$Non_full_length_reads, samp1$total_read_length)
	colnames(Frag1)=c("transcript",paste("ReadCount_","sample1",sep=""),paste("Frag","sample1",sep=""),paste("FLsample1",sep=""),paste("TotalRLsample1",sep=""))
	Frag1=as.data.frame(Frag1)
	Frag2<-cbind(samp2$transcript,samp2$read_count, samp2$fragmentation_rate, samp2$Non_full_length_reads, samp2$total_read_length)
	colnames(Frag2)=c("transcript",paste("ReadCount_","sample2",sep=""),paste("Frag","sample2",sep=""),paste("FL","sample2",sep=""),paste("TotalRL","sample2",sep=""))
	Frag2=as.data.frame(Frag2)

	mergedFrag <- merge(x=Frag1, y=Frag2, by="transcript",all=T)
	mergedFrag2=mergedFrag
	mergedFrag2$K1=is.element(mergedFrag2$transcript, Frag1$transcript[as.numeric(Frag1[,2])<thresh])
	mergedFrag2$K2=is.element(mergedFrag2$transcript, Frag2$transcript[as.numeric(Frag2[,2])<thresh])

	B= mergedFrag2 %>% group_by(transcript)  %>%
	do(data.frame(val=TranscriptTest_Tech(.,C1=alpha[1],C2=alpha[2],  a=a,b=b,thresh=thresh)))
	colnames(B)=c("transcript","logLik_samp1","logLik_samp2","logLik_E0","loglikratio", "Tau_EAP1","Tau_EAP2","Subsample1","Subsample2")


	DataResult<- merge(x=mergedFrag, y=B, by="transcript")

	DataResult$LFC=log(DataResult$Tau_EAP1/DataResult$Tau_EAP2)
	DataResult$LFC[is.na(DataResult$LFC)]<-0
	DataResult$posterior=(1-1/(1+p/(1-p)*exp(DataResult$loglikratio)))
	Degraded=as.numeric(as.numeric(DataResult$posterior)>0.5)
	Degraded[Degraded==1 & (as.numeric(DataResult$Tau_EAP2) > as.numeric(DataResult$Tau_EAP1))]<-(-1)
	DataResult$Degraded=factor(Degraded)


  return(DataResult)
}


Change_Prior<-function(ResultObject,p)
{
	NewResult=ResultObject
	NewResult$posterior=(1-1/(1+p/(1-p)*exp(ResultObject$loglikratio)))
	Degraded=as.numeric(as.numeric(NewResult$posterior)>0.5)
	Degraded[Degraded==1 & (as.numeric(NewResult$Tau_EAP2) > as.numeric(NewResult$Tau_EAP1))]<-(-1)
	NewResult$Degraded=factor(Degraded)
	return(NewResult)
}


PlotResults<-function(ResultObject,samplenames=c("sample 1","sample 2"),labels=FALSE,GeneLabel=c())
{

	G2b= ggplot(ResultObject,aes(x=LFC,y=(posterior),color=Degraded))+geom_point() + theme_minimal()  + ylim(0,1) + scale_color_manual(values=c("-1"="blue","0"="black","1"="red")) + labs(x="Biological Fragmentation\nlog-fold change",y="Posterior probability") +	geom_point(data=subset(ResultObject, Degraded==1),alpha = 0.9,size=1,color="red")+	geom_point(data=subset(ResultObject, Degraded==-1),alpha = 0.9,size=1,color="blue")+ theme(text = element_text(size = 40),legend.position = "none") 
	G3b= ggplot(ResultObject,aes(y=as.numeric(Tau_EAP1),x=as.numeric(Tau_EAP2),color=Degraded)) + theme_minimal()  + geom_point() +	geom_point(data=subset(ResultObject, Degraded==1),alpha = 0.9,size=1,color="red")+	geom_point(data=subset(ResultObject, Degraded==-1),alpha = 0.9,size=1,color="blue")+ scale_color_manual(values=c("-1"="blue","0"="black","1"="red")) + labs(x=paste("Biological Fragmentation rate\n",samplenames[2],sep=""),y=paste("Biological Fragmentation rate\n",samplenames[1],sep=""))+ theme(text = element_text(size = 40),legend.position = "none") + geom_abline(slope=1,intercept=0,linewidth=2,color="purple") 
	if (labels)
	{
		G2b = G2b +
  		geom_label_repel(aes(label = ifelse((Degraded!=0 & abs(LFC)>1),as.character(transcript),'')),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') 
		G3b = G3b +
  		geom_text(aes(label=ifelse((Degraded!=0 & abs(LFC)>1),as.character(transcript),'')),hjust=0,vjust=0)
	}
	if (length(GeneLabel)>0)
	{
		G2b = G2b +
  		geom_label_repel(aes(label = ifelse(is.element(as.character(transcript),GeneLabel),as.character(transcript),'')),
                  #box.padding   = 0.35, 
                  #point.padding = 0.5,
                  colour="deeppink",
                  segment.color = 'deeppink') +	
      geom_point(data=subset(ResultObject, is.element(as.character(transcript),GeneLabel)),alpha = 0.9,size=3,color="deeppink")
		G3b = G3b +
  		geom_text(aes(label = ifelse(is.element(as.character(transcript),GeneLabel),as.character(transcript),'')),
    		#position = position_dodge(width = 1),
    		hjust=0,vjust=0,
    		colour="deeppink") +	
      geom_point(data=subset(ResultObject, is.element(as.character(transcript),GeneLabel)),alpha = 0.9,size=3,color="deeppink")
	}

	Gt3b=ggarrange(G3b,G2b,nrow=1,ncol=2)
	return(Gt3b)
}

Get_Significant<-function(ResultObject)
{
	Result_Df=ResultObject[ResultObject$Degraded!=0,]
	Result_Df=Result_Df[order(Result_Df$posterior,decreasing=TRUE),]
	Result_Df=Result_Df[,c(1,2,6,13:17)]
	colnames(Result_Df)=c("Transcript","ReadCount_sample1","ReadCount_sample2","loglikratio","Biol_Degrad_sample1","Biol_Degrad_sample2","LogFoldChange","Posterior_Prob_Difference")
	Result_Df=Result_Df[!is.na(Result_Df$Posterior_Prob_Difference),]
	return(Result_Df)
}
