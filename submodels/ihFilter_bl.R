options(echo=F)
####### ARch_common(): examine which ARch variants are shared between one proband and two additional members.
### OUTPUT: ARch variants in proband that are shared with two additional members.
### GT3 contains rows for only ARch variants; gene is a vector of the corresponding gene names of these ARch vars.
### GT3 must be row-named with variant names.
ARch_common <- function(GT3,gene) { #combined: var table for one individual
#to identify genes with multiple heterozygous variants, and collapse each set of multiple varIDs for each such gene.
	if (unique(GT3[,1])!='het') warning('proband column should contains purely het variants!')
	var1 <- rownames(GT3)[GT3[,2]=='het']
	var2 <- rownames(GT3)[GT3[,3]=='het']
	gene1 <- gene[GT3[,2]=='het']
	gene2 <- gene[GT3[,3]=='het']
	fct1 <- factor(gene1,levels=unique(gene1),ordered=T)
	fct2 <- factor(gene2,levels=unique(gene2),ordered=T)
	varsInGene.1 <- tapply(var1,fct1,function(x) x)
	varsInGene.2 <- tapply(var2,fct2,function(x) x)
	gene.1 <- names(varsInGene.1)[sapply(varsInGene.1,length)>1]
	gene.2 <- names(varsInGene.2)[sapply(varsInGene.2,length)>1]
	sharedGenes <- intersect(intersect(gene,gene.1),gene.2)
	ix <- gene%in%sharedGenes
	reservVar <- rownames(GT3)[ix]
	reservVar
}

##### ARch_from_GT3(): identify ARch variants from trio-collated Genotype matrix and the proband data matrix.
### INPUT GT3: three columns indicating GT in offspring, dad, and mom.
### INPUT pro: step-4 format data of proband,with a column "geneName".
### OUTPUT: a logic vector of length nrow(pro), incicating which rows of pro should be kept as ARch variants.
### NOTE: GT3 and pro have the same set of rownames, in the form of chr_start_end_ref_alt.
### NOTE: GT3 take such values: {hom, het, wdt} 
ARch_from_GT3 <- function(GT3,pro) {
	f1occ <- (GT3[,1]=='het') #f1occ means Occurence in F1 offspring.
        DadOcc <- (GT3[,2]=='het')
        MomOcc <- (GT3[,3]=='het')
        dadTrans <- f1occ & DadOcc & (GT3[,3]=='wdt')
        momTrans <- f1occ & MomOcc & (GT3[,2]=='wdt')
        middle <- pro[dadTrans|momTrans,] #middle stage variants, still contains false positives as some genes have variants all from single parent.
        geneCnt <- table(middle$geneName)
        middle <- middle[middle$geneName%in%names(geneCnt)[geneCnt>1],]
        middle <- middle[order(middle$geneName),]
        middle.var <- middle$varID
        GT3ch <- GT3[middle.var,]
        fct.gene <- factor(middle$geneName,levels=unique(middle$geneName),ordered=T)
        parent_each_ctrb <- function(x) { #within a gene scope, each parent must transmit at least 1 het variant.
                logic2 <- cbind(x[,1]=='het', x[,2]=='het')
                ctrbSum <- colSums(logic2)
                qualify= (ctrbSum[1]>0 & ctrbSum[2]>0)
                qualify # one logic value indicating if this gene suffice parent-each-ctrb-1 criterion.
        }
        qualify_middleGenes = by(GT3ch[,2:3],fct.gene,parent_each_ctrb)
        qualify.genes <- levels(fct.gene)[qualify_middleGenes]
        qualify <- middle$geneName%in%qualify.genes
        ix.arch <- rownames(pro)%in%rownames(middle[qualify,])
	ix.arch
}

### ihFilter_parImpute(): in case of missing complete parental data, use this function to format the output in compatible with those out of the three formal functions.
######################### first written to address the incomplete family of case GB.
######################### Mostly identical to the counterpart for HA, but we apply frequency filter only to gnomad_exome but no gnomad_genome.
#### Task 1: remove leftmost column "inROH".
#### Task 2: append a rightmost column "ihMode". Basically keep all variants, labeling hom vars as ARs, het vars as deNovo. Not considering ARch & X-Linked.
#### Task 3: remove deNovo (het) variants with freq exceeding thStrict.
ihFilter_parImpute <- function(previous,thStrict=1E-4) {
        thStrict <- as.numeric(thStrict)
        proID.start <- regexpr('UDN\\d+',previous,perl=T,ignore.case=T)
        proID <- substr(previous,proID.start,proID.start+8)
        pro <- read.delim(previous,as.is=T,header=T)
	pro$gnomADexomeFreq <- as.numeric(pro$gnomADexomeFreq)
	pro$gnomADgenomeFreq <- as.numeric(pro$gnomADgenomeFreq)
        inROH <- pro[,'inROH']
        pro <- pro[,setdiff(colnames(pro),'inROH')]
        rownames(pro) <- pro[,'varID']
        colnames(pro)[3] <- 'geneName' ### GeneDetail.refGene
        ihMode <- character(nrow(pro))
        ix.denovo.a <- (is.na(pro[,'exacFreq']) | pro[,'exacFreq']<=thStrict)
        ix.denovo.b <- (is.na(pro[,'espFreq']) | pro[,'espFreq']<=thStrict)
        ix.denovo.c <- (is.na(pro[,'g1kFreq']) | pro[,'g1kFreq']<=thStrict)
	ix.denovo.d <- (is.na(pro[,'gnomADexomeFreq']) | pro[,'gnomADexomeFreq']<=thStrict) 
	#ix.denovo.d <- (is.na(pro[,'gnomADexomeFreq']) | pro[,'gnomADexomeFreq']<=thStrict) & (is.na(pro[,'gnomADgenomeFreq']) | pro[,'gnomADgenomeFreq']<=thStrict) 
        ihMode[pro[,'GT']=='het'&ix.denovo.a&ix.denovo.b&ix.denovo.c&ix.denovo.d] <- 'deNovo'
        ihMode[pro[,'GT']=='hom'] <- 'ARs'
        ihMode[inROH==1] <- paste('inROH', ihMode[inROH==1])
        cat(sum(pro[,'GT']=='het'&ix.denovo.a&ix.denovo.b&ix.denovo.c),' deNovo variants\n')
        cat(sum(pro[,'GT']=='hom'),' AR-simple variants\n')
        cat(sum(inROH==1),' inROH variants (sex-chr variants are imprecise)\n')
        pro <- data.frame(pro,ihMode=ihMode)
        pro <- pro[nchar(ihMode)>0,]
        write.table(pro,paste0(proID,'_parFiltered.txt'),row.names=F,sep='\t',quote=F)
}


#### ihFilter_par2pro(): given proband data matrix (usually outputted from step4) and two step2 files of both parents, 
#######################  identify three types of variants: ARs, ARch, and deNovo.
### INPUT previous: proID_step4.txt. Indicating the centering subject (proband) we are primarily consider. Contains 10 columns.
### INPUT dad & mom: udnIDs for the two parents. Corresponding udnID_step0.txt must be present in the wd.
### OUTPUT: a data matrix with reduced rows, but 4 more columns than "previous" data matrix. 
### NOTE: as we apply variant filtration criteria to affected and unaffected members, we regard parents as unaffected (as in general they are unaffected in our trios).
########### last (14th) column specifies the inheritance mode the variant abides by. Values are {'deNovo','ARs','ARch'}
ihFilter_par2pro <- function(previous,dadID,momID,thStrict=1E-5,Gender=1) {
        thStrict <- as.numeric(thStrict)
	proID.start <- regexpr('UDN\\d+',previous,perl=T,ignore.case=T)
	proID <- substr(previous,proID.start,proID.start+8)
	pro <- read.delim(previous,as.is=T,header=T)
        pro$gnomADexomeFreq <- as.numeric(pro$gnomADexomeFreq)
        pro$gnomADgenomeFreq <- as.numeric(pro$gnomADgenomeFreq)
	inROH <- pro[,'inROH'] #updated
	pro <- pro[,setdiff(colnames(pro),'inROH')] #updated
	rownames(pro) <- pro[,'varID']
	colnames(pro)[3] <- 'geneName' ### GeneDetail.refGene
	dad <- paste0(dadID,'_step0.txt')
	dad <- read.delim(dad,as.is=T,header=F,quote='')
	mom <- paste0(momID,'_step0.txt')
	mom <- read.delim(mom,as.is=T,header=F,quote='')
	dad.rnms <- paste('var',dad[,1],dad[,2],dad[,3],dad[,4],dad[,5],sep='_') 
        mom.rnms <- paste('var',mom[,1],mom[,2],mom[,3],mom[,4],mom[,5],sep='_') 
	dad <- cbind(dad.rnms,dad[,6])
	mom <- cbind(mom.rnms,mom[,6])
	proGT <- pro[,c('varID','GT')]
	GT3 <- merge(proGT,dad,by.x=1,by.y=1,all.x=T)
	GT3 <- merge(GT3,mom,by.x=1,by.y=1,all.x=T)
	GT3 <- cbind(as.character(GT3[,1]),as.character(GT3[,2]),as.character(GT3[,3]),as.character(GT3[,4]))
	GT3[is.na(GT3)] <- 'wdt' # for wild-type
	colnames(GT3) <- c('varID','proGT',dadID,momID)
	rownames(GT3) <- GT3[,1]	
	GT3 <- GT3[,-1]
	GT3 <- GT3[rownames(pro),]
       ###### deNovo #########
	ix.denovo <- (GT3[,1]=='het') & GT3[,2]=='wdt' & GT3[,3]=='wdt' #GT3[,1]=='hom'
        ix.denovo.a <- (is.na(pro[,'exacFreq']) | pro[,'exacFreq']<=thStrict)
        ix.denovo.b <- (is.na(pro[,'espFreq']) | pro[,'espFreq']<=thStrict)
        ix.denovo.c <- (is.na(pro[,'g1kFreq']) | pro[,'g1kFreq']<=thStrict)
	#ix.denovo.d <- (is.na(pro[,'gnomADexomeFreq']) | pro[,'gnomADexomeFreq']<=thStrict) & (is.na(pro[,'gnomADgenomeFreq']) | pro[,'gnomADgenomeFreq']<=thStrict)
	ix.denovo.d <- (is.na(pro[,'gnomADexomeFreq']) | pro[,'gnomADexomeFreq']<=thStrict)
        ix.denovo <- ix.denovo & (ix.denovo.a & ix.denovo.b & ix.denovo.c & ix.denovo.d)
        cat(sum(ix.denovo),' deNovo variants\n')
	#### ARs #########
     	ix.ars <- (GT3[,1]=='hom') & (GT3[,2]=='het') & (GT3[,3]=='het')
        cat(sum(ix.ars),' AR-simple variants\n')
        ####### ARch ##########
	ix.arch <- ARch_from_GT3(GT3,pro)
        cat(sum(ix.arch),' AR-ch variants for ',length(unique(pro[ix.arch,'geneName'])),' genes\n')
	#ix.keep <- (ix.denovo|ix.ars|ix.arch)
	#ix.drop <- !ix.keep
	ihMode <- character(nrow(pro))
	ihMode[ix.denovo] <- 'deNovo'
	ihMode[ix.ars] <- 'ARs'
	ihMode[ix.arch] <- 'ARch'
	ihMode[inROH==1] <- paste('inROH', ihMode[inROH==1])
	cat(sum(inROH==1),' inROH variants (sex-chr variants are imprecise)\n')
        ##### Impose X-linked filtration if it is a male proband ############
        if (Gender==1) {
                ix.xlinked <- grepl('var_X_',rownames(GT3)) & GT3[,1]!='wdt' & GT3[,3]!='wdt' & GT3[,2]!=GT3[,1] # Mom and son must be non-wdt; Dad must be different from son (usually wdt, or het when son is hom).
                cat(sum(ix.xlinked),' X-linked variants\n')
                ihMode[ix.xlinked] <- 'Xlinked'
        }
	pro.expanded <- data.frame(pro,GT3[,-1],ihMode=ihMode)
	kept <- pro.expanded[!is.na(ihMode),]
	kept <- kept[order(kept$ihMode),]
	dropped <- pro.expanded[is.na(ihMode),]
        write.table(kept,paste0(proID,'_parFiltered.txt'),row.names=F,sep='\t',quote=F)
	write.table(dropped,'parFiltered.txt',row.names=F,sep='\t',quote=F)
}

##### ihFilter_common(): Identify which variants (rows) are shared among all members (proband and 1 or 2 additional members).
### INPUT previous: proband data matrix (output of any of the three ihFilter functions), must have a last column for ihMode. 
################### 3rd col must be gene names; rightmost column must assume values {'ARs','deNovo','AD','ARch'}.
### INPUT sib1ID & sib2ID: udnIDs of two relatives. Initial case involved proband and his two siblings, so sib# is used to denote relatives.
### OUTPUT: a data matrix derived from "previous", with several GT columns inserted before the last column of ihMode.
ihFilter_common <- function(previous,sib1ID,sib2ID.0=NULL) { # can considering sharing situation relative to one or two additional members.
        proID.start <- regexpr('UDN\\d+',previous,perl=T,ignore.case=T)
        proID <- substr(previous,proID.start,proID.start+8)
        pro <- read.delim(previous,as.is=T,header=T)
        rownames(pro) <- pro[,'varID']
        colnames(pro)[3] <- 'geneName' #Gene.refGene
	sib1 <- paste0(sib1ID,'_step2.txt')
        sib1 <- read.delim(sib1,as.is=T,header=F,quote='')
	if ( is.null(sib2ID.0) ) {
		sib2ID <- sib1ID
	} else {
		sib2ID <- sib2ID.0
	}
	sib2 <- paste0(sib2ID,'_step2.txt')
	sib2 <- read.delim(sib2,as.is=T,header=F,quote='')
        sib1.rnms <- paste('var',sib1[,1],sib1[,2],sib1[,3],sib1[,4],sib1[,5],sep='_') 
        sib2.rnms <- paste('var',sib2[,1],sib2[,2],sib2[,3],sib2[,4],sib2[,5],sep='_') 
        sib1 <- cbind(sib1.rnms,sib1[,6])
        sib2 <- cbind(sib2.rnms,sib2[,6])
	proGT <- pro[,c('varID','GT')]
        GT3 <- merge(proGT,sib1,by.x=1,by.y=1,all.x=T)
        GT3 <- merge(GT3,sib2,by.x=1,by.y=1,all.x=T)
        GT3 <- cbind(as.character(GT3[,1]),as.character(GT3[,2]),as.character(GT3[,3]),as.character(GT3[,4]))
        GT3[is.na(GT3)] <- 'wdt' # for wild-type
        colnames(GT3) <- c('varID','proGT',sib1ID,sib2ID)
        rownames(GT3) <- GT3[,1]
        GT3 <- GT3[,-1]
        GT3 <- GT3[rownames(pro),]
	
	ix.ARs <- (pro[,ncol(pro)]=='ARs' & GT3[,2]=='hom' & GT3[,3]=='hom')
	ix.AD <- ( pro[,ncol(pro)]%in% c('deNovo','AD','Xlinked' ) & (GT3[,2]!='wdt') & (GT3[,3]!='wdt'))
	oldARch.ix <- (pro[,ncol(pro)]=='ARch')
	if (sum(oldARch.ix>0)) {
		var.ARch <- ARch_common(GT3[oldARch.ix,],pro[oldARch.ix,'geneName'])
		ix.ARch <- rownames(pro)%in%var.ARch
	} else {
		ix.ARch <- rep(FALSE,nrow(pro))
	}
	ix.ROH <- grepl('inROH',pro[,ncol(pro)])
	ix.keep <- (ix.ARs|ix.AD|ix.ARch|ix.ROH)
	ix.drop <- !ix.keep
	if (is.null(sib2ID.0)) {
		GT3 <- GT3[,1:2]
		colnames(GT3) <- c('proGT',sib1ID)
	}
	if (!is.null(sib2ID.0)) {
		pro.expanded <- data.frame(pro[,1:(ncol(pro)-1)],GT3[,-1],ihMode=pro[,ncol(pro)])
	} else {
		pro.expanded <- data.frame(pro[,1:(ncol(pro)-1)],GT3[,2],ihMode=pro[,ncol(pro)])
		colnames(pro.expanded)[ncol(pro.expanded)-1] = sib1ID
	}
	kept <- pro.expanded[ix.keep,]
	dropped <- pro.expanded[ix.drop,]
	write.table(kept,paste0(proID,'_commonFiltered.txt'),row.names=F,sep='\t',quote=F)
	write.table(dropped,'commonFiltered.txt',row.names=F,sep='\t',quote=F)
}

##### ihFilter_uniq(): Identify which variants (rows) are uniquely present in proband (not observed in additional members).
### INPUT previous: proband data matrix (output of any of the three ihFilter functions).
### INPUT negID, negID2: udnIDs of two relatives to exclude from.
### OUTPUT: a data matrix of the same colnames as "previous".
ihFilter_uniq <- function(previous,negID,negID2=NULL) { #negID corresponds to an unaffected subject, usually a sibling or an offspring of the proband. 
        proID.start <- regexpr('UDN\\d+',previous,perl=T,ignore.case=T)
        proID <- substr(previous,proID.start,proID.start+8)
        pro <- read.delim(previous,as.is=T,header=T)
        rownames(pro) <- pro[,'varID']
        colnames(pro)[3] <- 'geneName' ### ENSGN
	if (is.null(negID2)) {
		nRel=1 
		Rels = c(negID)
	} else {
		nRel=2
		Rels = c(negID,negID2)
	}
	for (i in 1:nRel) {
		negID = Rels[i]
		neg <- paste0(negID,'_step0.txt')
        	neg <- read.delim(neg,as.is=T,header=F,quote='')
        	neg.rnms <- paste('var',neg[,1],neg[,2],neg[,3],neg[,4],neg[,5],sep='_') ### varID: using 2rd VCF column twice.i
		neg <- cbind(neg.rnms,neg[,6])
		proGT <- pro[,c('varID','GT')]
        	GT2 <- merge(proGT,neg,by.x=1,by.y=1,all.x=T)
        	GT2 <- cbind(as.character(GT2[,1]),as.character(GT2[,2]),as.character(GT2[,3]))
        	GT2[is.na(GT2)] <- 'wdt' # for wild-type
        	colnames(GT2) <- c('varID','proGT',negID)
        	rownames(GT2) <- GT2[,1]
        	GT2 <- GT2[,-1]
        	GT2 <- GT2[rownames(pro),]
	
		ix.ARs <- (pro[,ncol(pro)]=='ARs' & GT2[,2]!='hom')
		ix.AD <- (pro[,ncol(pro)]%in% c('deNovo','AD','Xlinked' ) & (GT2[,2]=='wdt'))
			
		GT3 <- cbind(GT2,GT2[,2]) # generate a pseudo 3rd column of GT3 for the use by ARch_common().
		rownames(GT3) <- rownames(GT2)
		oldARch.ix <- (pro[,ncol(pro)]=='ARch')
	        if (sum(oldARch.ix>0)) {
	                sharedVar.ARch <- ARch_common(GT3[oldARch.ix,],pro[oldARch.ix,'geneName'])
        	        ix.ARch <- oldARch.ix & !(rownames(pro)%in%sharedVar.ARch)
        	} else {
                	ix.ARch <- rep(FALSE,nrow(pro))
        	}

		ix.ROH <- grepl('inROH',pro[,ncol(pro)])
		ix.keep <- (ix.ARs|ix.AD|ix.ARch|ix.ROH)
		ix.drop <- !ix.keep
		pro.expanded <- data.frame(pro[,1:(ncol(pro)-1)],GT2,ihMode=pro[,ncol(pro)])
		pro.expanded <- pro.expanded[,-(ncol(pro.expanded)-2)] #remove the second to last column: the redundant proband GT column.
	        pro <- pro.expanded[ix.keep,]
		if (!exists('dropped')) {
	        	dropped <- pro.expanded[ix.drop,]
		} else {
			dropped <- data.frame(dropped[,1:(ncol(dropped)-1)],temp=rep(NA,nrow(dropped)),ihMode=dropped[,ncol(dropped)])
			new <- pro.expanded[ix.drop,]
			colnames(dropped) <- colnames(new)
			dropped <- rbind(dropped,new)
		}
	}
	write.table(pro,paste0(proID,'_uniqFiltered.txt'),row.names=F,sep='\t',quote=F)
	write.table(dropped,'uniqFiltered.txt',row.names=F,sep='\t',quote=F)
}

args=(commandArgs(TRUE))
########### Display main input: previous file & additional udnIDs ###########
for (i in 1:length(args)) {
        if (i==1) {
                cat('method choice: ',args[[i]],'\n')
        } else {
                if (i==2)
                        cat('previous data file: ',args[[i]],'\n','additional relative IDs: ')
                else {
                        if (i==3 | i==4)
                                cat(args[[i]],' ')
                }
        }
}
cat('\n')

################ compile command and execute the command ############
method=switch(args[[1]],
        parTrio='ihFilter_par2pro',
        common='ihFilter_common',
        unique='ihFilter_uniq',
	impute='ihFilter_parImpute'
)
previous=args[[2]]
addID=args[[3]]
if (method=='ihFilter_par2pro') {
        Gender.arg <- grep('Gender=',unlist(args),value=T)
        if (length(Gender.arg)==1)
                eval(parse(text=Gender.arg)) #suppose the argument is like 'Xlinked=TRUE'
        else
                Gender=1
}

if (length(args)==3) {
        cmdTxt <- paste0(method,'("',previous,'","',addID,'")')
} else { ### when we have more than 3 arguments. Two additional udnIDs are involved.
        addID2 <- args[[4]]
        if (length(args)==4)
                cmdTxt <- paste0(method,'("',previous,'","',addID,'","',addID2,'")')
        else { ### now the method must be parTrio, and user has supplied custom thStrict & Gender.
                Gender.arg <- grep('Gender=',unlist(args),value=T)
                thStrict.arg <- grep('thStrict=',unlist(args),value=T)
                if (length(Gender.arg)==1)
                        eval(parse(text=Gender.arg))
                else
                        Gender=1
                if (length(thStrict.arg)==1)
                        eval(parse(text=thStrict.arg))
                else
                        thStrict=1E-5
                cmdTxt <- paste0(method,'("',previous,'","',addID,'","',addID2,'",',thStrict,',',Gender,')')
        }
}
cat('COMMAND: ',cmdTxt,'\n')
eval(parse(text=cmdTxt))




