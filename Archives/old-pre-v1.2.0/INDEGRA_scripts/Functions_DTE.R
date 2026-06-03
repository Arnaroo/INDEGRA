CreateData<-function(samples,samplefiles)
{
	n=length(samples)
	if (length(samplefiles)!=n) return(print("sample names and sample files do not have the same length"))
	if (n<2) return(print("Several samples are required for differential expression analysis"))
	
	M1=read.table(samplefiles[1],header=T)	
	S1 <- M1 %>% 
  	dplyr::select(transcript, fragmentation_rate, DTI,read_count)
	colnames(S1)<-c("Transcript",paste("Frag_",samples[1],sep=""),paste("DTI_",samples[1],sep=""),paste("RC_",samples[1],sep=""))  

	for (i in 2:n)
	{
		M2=read.table(samplefiles[i],header=T)
		S2 <- M2 %>% 
  		dplyr::select(transcript, fragmentation_rate, DTI,read_count)
			colnames(S2)<-c("Transcript",paste("Frag_",samples[i],sep=""),paste("DTI_",samples[i],sep=""),paste("RC_",samples[i],sep=""))  
		S1 <- merge(x=S1, y=S2, all=F,by="Transcript")
	}

	SCounts=S1[,seq(1,n*3,by=3)+3]
	SDTI=S1[,seq(1,n*3,by=3)+2]
	SFrag=S1[,seq(1,n*3,by=3)+1]

	rownames(SCounts)=S1$Transcript
	rownames(SDTI)=S1$Transcript
	rownames(SFrag)=S1$Transcript
	
	Output=list(Df_Counts=SCounts, Df_DTI=SDTI, Df_Frag=SFrag)
	return(Output)

}

Normalize_DTI<-function(samples,samplefiles)
{

	CD<-CreateData(samples,samplefiles)
	n=length(samples)
	
	Counts= as.numeric(as.vector(as.matrix(CD$Df_Counts)))
	DTIs= as.numeric(as.vector(as.matrix(CD$Df_DTI)))
	DataGlobal=data.frame(Count=Counts, DRTIN=DTIs, Sample=factor(rep(samples,each=nrow(CD$Df_Counts))))

	CorrectedCounts<-c()
	Offset=c()
	for (sample in samples)
	{
		Data=subset(DataGlobal, Sample==sample)
		drtin <- Data$DRTIN[(Data$Count>0) &(Data$Count<=quantile(Data$Count,probs=0.99))]
		logcount <- log(Data$Count[(Data$Count>0) &(Data$Count<=quantile(Data$Count,probs=0.99))])
		fit <- loess(logcount ~ drtin)
		y.fit <- predict(fit, newdata = Data$DRTIN)
		correct=Data$Count/exp(y.fit-median(logcount))
		CorrectedCounts <- c(CorrectedCounts, correct)
		Offset=c(Offset, log(correct)-log(Data$Count))
	}

	DataGlobal$CorrectedCounts=CorrectedCounts
	DataGlobal$Offset=Offset
	NewCorrectedCounts<-c()
	NewOffset=c()
	A=matrix(DataGlobal$CorrectedCounts,ncol=n)
	retval <- normalizeQuantileRank(as.matrix(A), robust=TRUE)
	ret <- log(retval) - log(A)

	NewCorrectedCounts=as.vector(retval)
	NewOffset=DataGlobal$Offset+as.vector(ret)

	DataGlobal$NewCorrectedCounts=NewCorrectedCounts
	DataGlobal$NewOffset=NewOffset

	Matrixoffset=matrix(DataGlobal$NewOffset,ncol=n)
	colnames(Matrixoffset)=samples
	rownames(Matrixoffset)=rownames(CD$Df_Counts)
	
	return(list(MatrixOffset=Matrixoffset, DataGlobal=DataGlobal, RawCounts=CD$Df_Counts))
}

Plot_Normalisation<-function(Normalized)
{
	DataPlot=Normalized$DataGlobal
	G1=ggplot(DataPlot,aes(x=DRTIN,y=log(Count+0.1),color=Sample))+geom_smooth(se=FALSE,linewidth=2.5) +labs(x="transcript DTI",y="log(transcript count)", title="Lowess regression of counts versus DTI") +theme_minimal() + theme(text = element_text(size = 30),legend.text = element_text(size=25)) #+ annotate("text", x = 7, y = 5.5,size=12, label = "Raw counts")
	G2=ggplot(DataPlot,aes(x=DRTIN,y=log(CorrectedCounts),color=Sample))+geom_smooth(se=FALSE,linewidth=3) +labs(x="transcript DTI",y="log(transcript count)") +theme_minimal() + theme(text = element_text(size = 30),legend.text = element_text(size=30)) + annotate("text", x = 7, y = 5,size=12, label = "DTI normalized counts")
	#G3=ggplot(DataPlot,aes(x=DRTIN,y=log(NewCorrectedCounts),color=Sample))+geom_smooth(se=FALSE,linewidth=3) +labs(x="transcript DTI",y="log(transcript count)") +theme_minimal() + theme(text = element_text(size = 30),legend.text = element_text(size=30)) + annotate("text", x = 7, y = 5,size=12, label = "Global normalized counts")

	Gt3b=ggarrange(G1,G2,nrow=1,ncol=2,common.legend=TRUE)
	return(Gt3b)

}



