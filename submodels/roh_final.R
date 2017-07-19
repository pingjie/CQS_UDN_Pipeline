options(echo=F)
### roh_summ(): check if the roh segments of proband are shared across all affected relatives and absent in any non-affected relative.
### OUTPUT proTbl: a matrix of three columns (chr, start, end). Can be a matrix of 0 row.
####### Note: "sharing" is defined as overlapping by a minimum of 1 bp.
roh_summ <- function(proID='UDN688025') {#,udnIDs=c('UDN024638','UDN581290','UDN483335','UDN132819'),affS=c(0,0,1,0)) {
	proTbl <- read.delim(paste0(proID,'_roh.txt'),as.is=T)[,3:5]
	rohFiles <- dir(pattern='_roh.txt')
	udnIDs <- gsub('_roh.txt','',rohFiles)
	affS <- rep(0,length(udnIDs))
	for (i in 1:length(udnIDs)) {
		tbl <- read.delim(paste0(udnIDs[i],'_roh.txt'),as.is=T)[,1:2]
		if (tbl[1,'affStatus']==1)
			affS[i] <- 1 
	}	
	rels.pos <- udnIDs[affS==1]
	rels.neg <- udnIDs[affS==0]
	nullRes <- matrix(nr=0,nc=ncol(proTbl))
	colnames(nullRes) <- c('chr','start','end')
	#if (length(rels.pos)>0) {
		for (i in 1:length(rels.pos)) { # there will be at least one rel.pos, which is proID itself [self-to-self ovalapping check]
			rel.i <- rels.pos[i]
			relTbl <- read.delim(paste0(rel.i,'_roh.txt'),as.is=T)[,3:5]
			if (length(unlist(relTbl))==3) ### NEWLY ADDED
				relTbl <- data.frame(t(relTbl))
			cmmChrs <- setdiff(intersect(unique(proTbl[,1]),unique(relTbl[,1])),NA)
                	PROTBL <- matrix(nr=0,nc=ncol(proTbl))
			if (length(cmmChrs)>0) {
				proTbl <- proTbl[proTbl[,1]%in%cmmChrs,] # proTbl is shrinking, but a new variable to substitute it is in growth which absorbs useful info from disappearing proTbl lines.
				relTbl <- relTbl[relTbl[,1]%in%cmmChrs,]
				proList <- break_by_chr(proTbl)
				relList <- break_by_chr(relTbl)
				for (j in 1:length(cmmChrs)) {
					chr <- cmmChrs[j]
					shareMat <- mapply(overlap_or_not,rep(proList[[chr]],times=length(relList[[chr]])),rep(relList[[chr]],each=length(proList[[chr]]) ))
					shareMat <- matrix(shareMat,nr=length(proList[[chr]]))
					idx <- rowSums(shareMat,na.rm=T)>=1 # taking the overlapping rows.
					if (sum(idx)>0) {
						interim <- matrix(cbind(chr=chr,do.call(rbind,proList[[chr]][idx])),nc=3)				
						PROTBL <- rbind(PROTBL,interim)
					}
				}			

			}		
			colnames(PROTBL) <- c('chr','start','end')	
			proTbl <- data.frame(PROTBL)
			if (nrow(PROTBL)==0) 
				break
		}
	#}
	if ( nrow(proTbl)>0 & length(rels.neg)>0) {
		for (i in 1:length(rels.neg)) {
	        	rel.i <- rels.neg[i]
        		relTbl <- read.delim(paste0(rel.i,'_roh.txt'),as.is=T)[,3:5]
                	cmmChrs <- intersect(unique(proTbl[,1]),unique(relTbl[,1]))
                	PROTBL <- matrix(nr=0,nc=ncol(proTbl))
	        	if (length(cmmChrs)>0) {
				proTbl.0 <- proTbl[!proTbl[,1]%in%cmmChrs,] # the complement set proTbl.0 is kept aside, and will be combined with non-overlapping rohs of cmmChrs, and then collectively will replace proTbl. 
        	        	proTbl <- proTbl[proTbl[,1]%in%cmmChrs,]
                        	relTbl <- relTbl[relTbl[,1]%in%cmmChrs,]
	              		proList <- break_by_chr(proTbl)
                       		relList <- break_by_chr(relTbl)
                       		for (j in 1:length(cmmChrs)) {
                               		chr <- cmmChrs[j]
                               		shareMat <- mapply(overlap_or_not,rep(proList[[chr]],times=length(relList[[chr]])),rep(relList[[chr]],each=length(proList[[chr]]) ))
                               		shareMat <- matrix(shareMat,nr=length(proList[[chr]]))
                               		idx <- rowSums(shareMat,na.rm=T)==0 # taking the non-overlapping rows.
					if (sum(idx)>0) {
	                               		interim <- cbind(chr=chr,do.call(rbind,proList[[chr]][idx]))
        	                       		PROTBL <- rbind(PROTBL,interim)
					}
                       		}
				colnames(PROTBL) <- colnames(proTbl.0)
				proTbl <- rbind(proTbl.0,PROTBL)
				colnames(proTbl) <- c('chr','start','end')
                	}
			if (nrow(proTbl)==0) 
				break
		}
	}
	if (nrow(proTbl)>0) #Add column runLen_M if there is at least one row in proTbl. 
		proTbl <- data.frame(proTbl,runLen_M=round((as.numeric(as.character(proTbl[,3]))-as.numeric(as.character(proTbl[,2])))/1000000,1))
	write.table(proTbl,paste0(proID,'_roh_final.txt'),row.names=F,sep='\t',quote=F)
	#### processing step4 tbl ####
        step4 <- read.delim(paste0(proID,'_step4.txt'),as.is=T)
        varID <- step4$varID
        relTbl <- matrix(unlist(strsplit(varID,'_')),byrow=T,nc=6)[,2:5]
	if (length(unlist(relTbl))<=4) # NEWLY ADDED. this if sentence is added 3/20/2017.
		relTbl <- data.frame(t(relTbl))
	rohCol <- rep(0,nrow(step4))
        cmmChrs <- intersect(unique(relTbl[,1]),unique(proTbl[,1]))
        if (nrow(proTbl)>0 & length(cmmChrs)>0) {
                        rohCol <- rep(0,nrow(step4))
                        #proTbl <- proTbl[proTbl[,1]%in%cmmChrs,] # proTbl is shrinking, but a new variable to substitute it is in growth which absorbs useful info from disappearing proTbl lines.
                        #relTbl <- relTbl[relTbl[,1]%in%cmmChrs,]
                        proList <- break_by_chr(proTbl)
                        relList <- break_by_chr(relTbl)
                        #browser()
                        for (j in 1:length(cmmChrs)) {
                                chr <- cmmChrs[j]
                                shareMat <- mapply(overlap_or_not,rep(proList[[chr]],times=length(relList[[chr]])),rep(relList[[chr]],each=length(proList[[chr]]) ))
                                shareMat <- matrix(shareMat,nr=length(proList[[chr]]))
                                idx <- colSums(shareMat,na.rm=T)>=1 # taking the overlapping rows.
                                if (sum(idx)>0)
                                        rohCol[which(relTbl[,1]==chr)[idx]]=1
                        }
        }
	relTbl <- data.frame(inROH=rohCol,step4)
	write.table(relTbl,paste0(proID,'_step4.txt'),quote=F,row.names=F,sep='\t')
	proTbl

}
roh_final <- function(proID='UDN000000',col.shade='red') {
	rohFile=paste0(proID,'_roh_final.txt')
	dosFile=paste0(proID,'.dos')
	dos <- read.delim(dosFile,header=F,sep=' ')
	roh <- read.delim(rohFile,header=T)
	if (nrow(roh)==0) {
                return(NULL)
        }
	chr.dens <- hDens(dos) # using default wsz & sep specified in hDens() function.
        roh.fct <- factor(roh[,1],levels=unique(roh[,1]),ordered=T)
        chr.roh <- by(roh,roh.fct,extractChrRoh,simplify=F)
	chrs <- names(chr.roh)
	COL <- rgb(t(col2rgb(col.shade)),alpha=100,maxColorValue=255)

        pdf(paste0(proID,'_roh_final.pdf')) 
	### end of hurriedly
	for (chr in chrs) {
		chrI.dens <- chr.dens[[chr]]
		chrI.roh <- chr.roh[[chr]]
                layout(matrix(1:2,nr=2))
                plot(chrI.dens[,1],chrI.dens[,2],col='blue',type='l',ylim=c(0,1),xlim=c(0,max(chrI.dens[,1])),lwd=2,las=2,xlab='',xaxt='n',ylab='Homozygosity',main=paste0('chr',chr))
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
}
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
        #roh <- roh[order(roh[,2]),]
        #roh <- roh[roh[,5]>=qual&roh[,6]>=minROH,]
        #roh <- roh[roh[,4]==1,c(2:3,5:6)]
        #colnames(roh) <- c('start','end','qualScore','runLen_M')
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

### rohtbl: roh rows that are concerned in memember-to-memember comparisons.
break_by_chr <- function(rohtbl) {
        chrs <- factor(rohtbl[,1],levels=unique(rohtbl[,1]),ordered=T)
        res <- by(rohtbl[,2:3],chrs,function(x) split(t(x),rep(1:nrow(x),each=ncol(x)) ),simplify=FALSE ) #res will be a two-layered list: 1st, chrs; 2nd, each interval as a list component.
        #names(res) <- unique(chrs)
        res
}
### overlap_or_not(): judge if two intervals (of the same chromosome) are overlapping with each other or not.
### EX: itv1 = c(186108812,188663383), itv2 = c(188563383,188863383)
overlap_or_not <- function(itv1,itv2) {
        itv1 <- as.numeric(unlist(itv1)) # may be unnecessary.
        itv2 <- as.numeric(unlist(itv2)) # may be unnecessary.
        if (itv1[1]>=itv2[1]) {
                tmp <- itv1
                itv1 <- itv2
                itv2 <- tmp
        }
        if (itv2[1]<=itv1[2]) {
                return(TRUE)
        } else {
                return(FALSE)
        }
}

args=unlist((commandArgs(TRUE)))
proID=args[1]
temp = roh_summ(proID)
roh_final(proID)
