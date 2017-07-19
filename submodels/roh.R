options(echo=F)
### INPUT dos: three-columned matrix, 1st col for chrNum, 2nd for chrSite, 3nd col for dosage (0 for 0/0, 1 for 0/1, 2 for 1/1) #,3rd col for HW or AH statuses (0's or 1's)
### INPUT wsz: the size (in Miilions) of sliding window in which a homozygous density is calculated.
### INPUT sep: the distance (in Millions) between discrete points (ftpts) at each of which a window is formed and one density value is calculated.
hDens <- function(dos,wsz=1,sep=0.5) { #,col.shade='orange') {
	wsz = wsz*1000000
	sep=sep*1000000
	hDens.bottom <- function(site2dose,ftpts,wsz) {
		window.rb <- ftpts+wsz
		binDens <- function(ele2,site2dose) {
			a <- ele2[1]
			b <- ele2[2]
			ix <- ( site2dose[,'site']>=a & site2dose[,'site']<b )
			binDens.r <- sum(site2dose[ix,'dose']!=1)/sum(ix)
		}
		bounds <- cbind(ftpts,window.rb)
		hDens.bottom.r <- apply(bounds,1,binDens,site2dose)
		substitute <- median(hDens.bottom.r,na.rm=T)
		hDens.bottom.r[is.na(hDens.bottom.r)] <- substitute
		hDens.bottom.r
	}
	colnames(dos) <- c('chr','site','dose')
	chr.fct <- factor(dos$chr,levels=unique(dos$chr),ordered=T)
	#dos.lst <- by(dos,chr.fct,function(x) x,simplify=F)

	hDens.each <- function(dos,wsz,sep) {
		dos <- dos[order(dos$site),]
		sites <- dos[,'site']
		ftpts <- seq(from=min(sites),to=max(sites),by=sep)
		ftpts <- ftpts[1:(length(ftpts)-1)]
		hDens.each.r <- hDens.bottom(dos[,c('site','dose')],ftpts,wsz)
		hDens.each.r <- cbind(disSite=ftpts,hDens=hDens.each.r)
	}
	hDens.r <- by(dos,chr.fct,hDens.each,wsz,sep,simplify=F)
	#plot(ftpts,hDens.dis,col='blue',type='l',ylim=c(0,1),lwd=2,xlab='Genomic coordinate',ylab='Homozygosity',main=chr)
	#forPolygon <- rbind(c(min(sites),0),dos_roh[,c('site','roh')],c(max(sites),0))
	#COL <- rgb(t(col2rgb(col.shade)),alpha=100,maxColorValue=255)
	#polygon(forPolygon[,1],forPolygon[,2],col=COL,border=NA)
	#hDens <- cbind(disSite=ftpts,hDens=hDens.dis)
	#hDens
}
extractChrRoh <- function(roh) {
	roh <- roh[,2:3]
	row2vec4 <- function(row) {
		a <- row[1]
		b <- row[2]
		vec4 <- c(a,a,b,b)
	}
	vec4 <- as.vector(apply(roh,1,row2vec4))
	forPolygon <- cbind(vec4,rep(c(0,1,1,0),nrow(roh)))
	#res <- list(forPolygon=forPolygon,roh=roh)
}
roh <- function(udnID='UDN000000',affStatus=1,wsz=1,sep=0.5,minROH=2,qual=0,col.shade='orange') {
	rohFile=paste0(udnID,'.roh')
	dosFile=paste0(udnID,'.dos')
	dos <- read.delim(dosFile,header=F,sep=' ')
	chr.dens <- hDens(dos,wsz,sep)
	roh <- read.delim(rohFile,header=T)
	roh <- roh[roh[,5]>=qual&roh[,6]>=minROH&roh[,4]==1,c(1:3,5:6)]
	roh.fct <- factor(roh[,1],levels=unique(roh[,1]),ordered=T)
	chr.roh <- by(roh,roh.fct,extractChrRoh,simplify=F)	
	chrs <- names(chr.roh)[!grepl('M',names(chr.roh))] #allowing X & Y ROHs in result
	#layout(matrix(1:length(chrs),nc=2))
	COL <- rgb(t(col2rgb(col.shade)),alpha=100,maxColorValue=255)
	pdf(paste0(udnID,'_roh.pdf'))
	for  (chr in chrs) {
		chrI.dens <- chr.dens[[chr]]
		chrI.roh <- chr.roh[[chr]]
		layout(matrix(1:2,nr=2))
		plot(chrI.dens[,1],chrI.dens[,2],col='blue',type='l',ylim=c(0,1),xlim=c(0,max(chrI.dens[,1])),lwd=2,las=2,xlab='',xaxt='n',ylab='Homozygosity',main=paste0('chr',chr)) #11/10/2016: xlim max border changed to be max density coordinate. #max(chrI.roh[,1])
		xtag = seq(from=0,to=240,by=20)
		axis(side=1,at=xtag*1e6,label=xtag)
		#par(adj=1)
		title(xlab='Genomic coordinate  (x 1e6)')
		#resetPar()
		for (i in 1:length(xtag)) {
			abline(v=xtag[i]*1e6,col='gray',lty='dotted')
		}
		polygon(chrI.roh[,1],chrI.roh[,2],col=COL,border=NA)
		plot(0, type='n', xlim=c(0,max(chrI.roh[,1])), ylim=c(0,1), xaxt='n', yaxt='n', xlab='', ylab='', main='ROH of longer length (>2M)' )
		roh.i <- roh[roh[,1]==chr,]
		roh.i.text <- apply(roh.i,1,paste,collapse='  ')
		legend('center',legend=roh.i.text,col=COL,pch=15,bty='n')
	}
	dev.off()
	if (nrow(roh)==0) {
		for (i in 1:ncol(roh)) 
			roh[1,i] <- NA
	}
	roh <- data.frame(udnID=udnID,affStatus=affStatus,roh)
	colnames(roh) <- c('udnID','affStatus','chr','start','end','minQual','runLen_M')
	write.table(roh,paste0(udnID,'_roh.txt'),row.names=F,sep='\t',quote=F)
	#write.table(roh,paste(filePre,'minROH',minROH,'qual',qual,'txt',sep='.'),row.names=F,sep='\t',quote=F)
	res <- list(roh=roh,dens=chr.dens)	
}

resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}

#if (FALSE) {" ### begin of annotation

args=unlist((commandArgs(TRUE)))
udnID=args[1]
affStatus=args[2]
#### customize filePre
if (FALSE) {"
if (sum(grepl('filePre=',args))>0) {
	eval( parse( text=args[grepl('filePre=',args)] ) )
} else {
	filePre='hDensity'
}
"}
#### customize wsz
if (sum(grepl('wsz=',args))>0) {
	eval( parse( text=args[grepl('wsz=',args)] ) )
} else {
	wsz=1
}
#### customize sep
if (sum(grepl('sep=',args))>0) {
	eval( parse( text=args[grepl('sep=',args)] ) )
} else {
	sep=0.5
}
#### customize minROH
if (sum(grepl('minROH=',args))>0) { 
	eval( parse( text=args[grepl('minROH=',args)] ) )
} else {
	minROH=2
}
#### customize qual
if (sum(grepl('qual=',args))>0) {
	eval( parse( text=args[grepl('qual=',args)] ) )
} else {
	qual=0
}
#### customize col.shade
if (sum(grepl('col.shade=',args))>0) { 
	eval( parse( text=args[grepl('col.shade=',args)] ) )
} else {
	col.shade='orange'
}
roh(udnID,affStatus,wsz,sep,minROH,qual,col.shade)
#} ### end of annotation
