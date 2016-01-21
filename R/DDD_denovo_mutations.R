DDD_mis <- function(){
    
    options(stringsAsFactors=FALSE)
    ## DDD
    filename <- "DDD_mutations/nature14135-s2/Table s2.csv" ## SNV and indels
    deddd <- read.csv(filename,stringsAsFactors=FALSE)
    deddda <- deddd[,c("Chr","Pos","Pos")]
    deddda[,"Ref"] <- sapply(1:dim(deddd)[1],function(i) unlist(strsplit(deddd[i,"Ref.Alt"],"/"))[1])
    deddda[,"Alt"] <- sapply(1:dim(deddd)[1],function(i) unlist(strsplit(deddd[i,"Ref.Alt"],"/"))[2])
    deddda[,"Gene"] <- deddd[,"Gene"]
    deddda[,"Category"] <- deddd[,"Consequence"]
    deddda[,"AAchanges"] <- deddd[,"AAchange"]
    colnames(deddda) <- c("Chr","Start","End","Ref","Alt","Gene","Category","AAchanges")
    deddda[,"sampleID"] <- ""
    deddda[,"Disease"] <- "DDD"
    deddda[,"From"] <- "nature14135_DDD"
    
    ### ID
    filename <- "DDD_mutations/ID/plosone/journal.pgen.1004772.s002.csv"
    deidplos <- read.csv(filename,skip=1)
    deidplosa <- deidplos[,c("Chr","Position","Position","Reference.Allele","Mutant.Allele","Gene.symbol","variant.type","Family.ID","Detailed.annotation.of.the.variant")]
    colnames(deidplosa) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges")
    deidplosa[,"Disease"] <- "ID"
    deidplosa[,"From"] <- "Hamdan et al., PLOS genetics"
    
    
    idsampleMap <- read.csv("DDD_mutations/ID/id100-idna13394.csv",skip=2,header=FALSE)
    filename <- "DDD_mutations/ID/100/all.csv"
    deid100 <- read.csv(filename)
    deid100 <- deid100[!(deid100[,"V1"] %in% idsampleMap[,2]),]
    
    deid100[is.na(deid100[,8]),8] <- ""
    tmp <- gsub("Chr","",deid100[,"V3"])
    deid100a <- sapply(1:dim(deid100)[1], function(i) {
        
        sub1 <- grepl("del",tmp[i])
        sub2 <- grepl("-",tmp[i])
        sub3 <- grepl("_",tmp[i])
        a1 <- unlist(strsplit(tmp[i],"\\("))[1]
        subtmp <- unlist(strsplit(tmp[i],"g\\."))[2]
        
        if(sub1 | sub2 | sub3){
            
                subtmp <-  gsub("\\D*$","",subtmp)
                if(sub2){
                    a2 <- as.numeric(unlist(strsplit(subtmp,"-"))[1])
                    a3 <- as.numeric(unlist(strsplit(subtmp,"-"))[2])
                }else if(sub3){
                    a2 <- as.numeric(unlist(strsplit(subtmp,"_"))[1])
                    a3 <- as.numeric(unlist(strsplit(subtmp,"_"))[2])
                }else{
                    a2 <- as.numeric(subtmp)
                    a3 <- a2
                }
            
        }else{
            a2 <- as.numeric(substr(subtmp,1,nchar(subtmp)-3))
            a3 <- a2        
        }
        
        if(grepl(">",subtmp)){
            a4 <- substr(subtmp,nchar(subtmp)-2,nchar(subtmp)-2)
            a5 <- substr(subtmp,nchar(subtmp),nchar(subtmp))
        }else{a4 <- "-";a5 <- "-";}
        
        c(a1,a2,a3,a4,a5)
        
        })
    deid100a <- t(deid100a)
    deid100a <- cbind(deid100a,deid100[,"V2"])
    deid100a <- cbind(deid100a,"synonmous")
    
    subs <- (deid100[,7]=="" & deid100[,8]=="" & deid100[,9]=="" & deid100[,10]=="" & deid100[,11]=="")
    deid100a[deid100[,9]=="D",7] <- "LOF"
    deid100a[deid100[,9]!="D" & !subs ,7] <- "MIS"
    deid100a <- cbind(deid100a,deid100[,"V5"])
    deid100a <- cbind(deid100a,deid100[,"V1"])
    deid100a <- cbind(deid100a,"ID")
    deid100a <- cbind(deid100a,"deLigt_et_al.,_NEJM")
    colnames(deid100a) <- c("Chr","Start","End","Ref","Alt","Gene","Category","AAchanges","sampleID","Disease","From")

    
    filename <- "DDD_mutations/ID/nature13394/all.csv"
    deidna <- read.csv(filename,header=FALSE)
    for(i in 1:16){
        deidna[,i] <- gsub(" $", "", deidna[,i])
    }
    tmp <- gsub("Chr","",deidna[,"V3"])
    deidnaa <- sapply(1:dim(deidna)[1], function(i) {
        
        sub1 <- grepl("del",tmp[i])
        sub2 <- grepl("-",tmp[i])
        sub3 <- grepl("_",tmp[i])
        sub4 <- grepl("dup",tmp[i])
        sub5 <- grepl("ins",tmp[i])
        a1 <- unlist(strsplit(tmp[i],"\\("))[1]
        subtmp <- unlist(strsplit(tmp[i],"g\\."))[2]
        
        if(sub1 | sub2 | sub3 | sub4 | sub5){
            subtmp <-  gsub("\\D*$","",subtmp)
            if(sub2){
                a2 <- as.numeric(unlist(strsplit(subtmp,"-"))[1])
                a3 <- as.numeric(unlist(strsplit(subtmp,"-"))[2])
            }else if(sub3){
                a2 <- as.numeric(unlist(strsplit(subtmp,"_"))[1])
                a3 <- as.numeric(unlist(strsplit(subtmp,"_"))[2])
            }else{
                a2 <- as.numeric(subtmp)
                a3 <- a2
            }
        }else{
            a2 <- as.numeric(substr(subtmp,1,nchar(subtmp)-3))
            a3 <- a2        
        }
        
        if(grepl(">",subtmp)){
            a4 <- substr(subtmp,nchar(subtmp)-2,nchar(subtmp)-2)
            a5 <- substr(subtmp,nchar(subtmp),nchar(subtmp))
        }else{a4 <- "-";a5 <- "-";}
        
        c(a1,a2,a3,a4,a5)
        
    })
    deidnaa <- t(deidnaa)
    deidnaa <- cbind(deidnaa,deidna[,"V2"])
    deidnaa <- cbind(deidnaa,deidna[,"V6"])
    deidnaa <- cbind(deidnaa,deidna[,"V1"])
    deidnaa <- cbind(deidnaa,deidna[,"V4"])
    deidnaa <- cbind(deidnaa,"ID")
    deidnaa <- cbind(deidnaa,"nature13394")
    colnames(deidnaa) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From")
    
    ## Lancet
    deidlan <- read.csv("DDD_mutations/ID/Lancet/all.csv",header=FALSE)
    tmp <- gsub("chr","",deidlan[,"V4"])
    tmp <- gsub(" $", "", tmp)
    deidlana <- sapply(1:dim(deidlan)[1], function(i) {
        sub1 <- grepl("del",tmp[i])
        sub2 <- grepl("-",tmp[i])
        sub3 <- grepl("_",tmp[i])
        a1 <- unlist(strsplit(tmp[i],":"))[1]
        subtmp <- unlist(strsplit(tmp[i],"g\\."))[2]
        
        if(sub1 | sub2 | sub3){
                subtmp <-  gsub("\\D*$","",subtmp)
                if(sub2){
                    a2 <- as.numeric(unlist(strsplit(subtmp,"-"))[1])
                    a3 <- as.numeric(unlist(strsplit(subtmp,"-"))[2])
                }else if(sub3){
                    a2 <- as.numeric(unlist(strsplit(subtmp,"_"))[1])
                    a3 <- as.numeric(unlist(strsplit(subtmp,"_"))[2])
                }else{
                    a2 <- as.numeric(subtmp)
                    a3 <- a2
                }
        }else{
            a2 <- as.numeric(substr(subtmp,1,nchar(subtmp)-3))
            a3 <- a2        
        }
        
        if(grepl(">",subtmp)){
            a4 <- substr(subtmp,nchar(subtmp)-2,nchar(subtmp)-2)
            a5 <- substr(subtmp,nchar(subtmp),nchar(subtmp))
        }else{a4 <- "-";a5 <- "-";}
        
        c(a1,a2,a3,a4,a5)
        
    })
    deidlana <- t(deidlana)
    deidlana <- cbind(deidlana,deidlan[,"V2"])
    deidlana <- cbind(deidlana,deidlan[,"V3"])
    deidlana <- cbind(deidlana,deidlan[,"V6"])
    deidlana <- cbind(deidlana,deidlan[,"V1"])
    deidlana <- cbind(deidlana,"ID")
    deidlana <- cbind(deidlana,"Rauch_et_al.,_Lancet")
    colnames(deidlana) <- c("Chr","Start","End","Ref","Alt","Gene","Category","AAchanges","sampleID","Disease","From")
    
    deIDa <- rbind(deidplosa,deid100a,deidnaa,deidlana)


    ### EE
    filename <- "DDD_mutations/epileptic/nature12439/all.csv" 
    deeena <- read.csv(filename,header=FALSE)
    for(i in 1:9){
        deeena[,i] <- gsub(" $", "", deeena[,i])
    }
    deeenaa <- sapply(1:dim(deeena)[1],function(i) unlist(strsplit(deeena[i,"V2"],":"))[1])
    deeenaa <- as.data.frame(deeenaa)
    deeenaa[,"Start"] <- sapply(1:dim(deeena)[1],function(i) unlist(strsplit(deeena[i,"V2"],":"))[2])
    deeenaa[,"End"] <- deeenaa[,"Start"]
    deeenaa[,"Ref"] <- sapply(1:dim(deeena)[1],function(i) unlist(strsplit(deeena[i,"V3"],"/"))[1])
    deeenaa[,"Alt"] <- sapply(1:dim(deeena)[1],function(i) unlist(strsplit(deeena[i,"V3"],"/"))[2])
   # deeenaa[,"Alt"] <- gsub(" $", "", deeenaa[,"Alt"])
    deeenaa[,"Gene"] <- deeena[,"V5"]
    deeenaa[,"AAchanges"] <- ""
    deeenaa[,"Category"] <- deeena[,"V7"]
    deeenaa[,"sampleID"] <- deeena[,"V1"]
    deeenaa[,"Disease"] <- "EE"
    colnames(deeenaa) <- c("Chr","Start","End","Ref","Alt","Gene","AAchanges","Category","sampleID","Disease")
    deeenaa[,"From"] <- "nature12439"
   
    ## ASD
    asddataa <- read.csv("ASD/ASDmutationlists_annovar.csv")
    
    ALLde <- rbind(deIDa,deeenaa,asddataa)
    source("Network_analysis.R")
    ALLde[,"Gene"] <- mapping_to(ALLde[,"Gene"])
    write.csv(ALLde,file="DDD_mutations/datasheet/ALLdemutationlist_annovar.csv",row.names=FALSE)
    write.table(ALLde[,1:5],file="DDD_mutations/datasheet/ALLdemutationlist_annovar.txt",row.names=FALSE,sep="\t",quote=FALSE,col.names=FALSE)


    ALLde <- deddda
    source("Network_analysis.R")
    ALLde[,"Gene"] <- mapping_to(ALLde[,"Gene"])
    write.csv(ALLde,file="DDD_mutations/datasheet/DDDdemutationlist_annovar.csv",row.names=FALSE)
    write.table(ALLde[,1:5],file="DDD_mutations/datasheet/DDDdemutationlist_annovar.txt",row.names=FALSE,sep="\t",quote=FALSE,col.names=FALSE)

}

Clear_mis_metaSVM <- function(){
    
    reann <- read.csv("DDD_mutations/datasheet/Annotated_mutationlist_4_15_metaSVM.csv")
    reall <- read.csv("DDD_mutations/datasheet/mutationlist_4_15_all.csv") 
    reall <- cbind(reall,reann[,c("Func.refgene","Gene.refgene","ExonicFunc.refgene","RadialSVM_pred","AAChange.refgene")])
    reall <- cbind(reall,"")
    colnames(reall) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.annovar","Gene.annovar","Category_annovar","Polyphen2_HDIV_pred_annovar","AAChange.annovar","mutation_type_annovar","Func.metaSVM","Gene.metaSVM","Category_metaSVM","RadialSVM_pred_metaSVM","AAChange.metaSVM","mutation_type_metaSVM")
    
    ## table(reann[,"ExonicFunc.refgene"])
    # c("","frameshift deletion","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV;nonsynonymous SNV","stopgain","stopgain;stopgain","stoploss","synonymous SNV","synonymous SNV;synonymous SNV","synonymous SNV;synonymous SNV;synonymous SNV","unknown")

    #subs <- reann[,"ExonicFunc.refgene"]==""
    #table(reann[subs,"Func.refgene"])
    ##c("","downstream","intergenic","intronic","ncRNA_exonic","ncRNA_exonic;ncRNA_exonic;ncRNA_exonic","ncRNA_exonic;ncRNA_exonic;ncRNA_exonic;ncRNA_exonic;ncRNA_exonic","ncRNA_exonic;splicing","ncRNA_intronic","splicing","splicing;splicing","upstream","UTR3","UTR3;UTR3;UTR3","UTR5") 
    annlof <- c("frameshift deletion","frameshift insertion","stopgain","stopgain;stopgain","stoploss")
    subslof <- reann[,"ExonicFunc.refgene"] %in% annlof | (reann[,"ExonicFunc.refgene"]=="" & (reann[,"Func.refgene"] %in% c("splicing","splicing;splicing","")))
    #subslof1  <- reann[,"ExonicFunc.refgene"]=="" & reann[,"Func.refgene"]=="" ## table(reall[subslof1,"Category"])
    subsmisd <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift deletion",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift insertion",reann[,"ExonicFunc.refgene"])) & reann[,"RadialSVM_pred"]=="D"
    subsmis <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift deletion",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift insertion",reann[,"ExonicFunc.refgene"])) & reann[,"RadialSVM_pred"]!="D"
    subssyn <- grepl("synonymous SNV",reann[,"ExonicFunc.refgene"]) & !grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"])
    subsintron <- reann[,"ExonicFunc.refgene"]=="" & !(reann[,"Func.refgene"] %in% c("splicing","splicing;splicing",""))
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    
    reall[subslof,"mutation_type_metaSVM"] <- "LOF"
    reall[subsmisd,"mutation_type_metaSVM"] <- "dMIS"
    reall[subsmis,"mutation_type_metaSVM"] <- "MIS"
    reall[subssyn,"mutation_type_metaSVM"] <- "SYNONMOUS"
    reall[subsintron,"mutation_type_metaSVM"] <- "Other"
    
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    subsunlof <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"] %in% c("frame-shift")
    subsunsyn <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"] %in% c("synonmous","synonymous")
    subsundmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"]== "missense" & reann[,"RadialSVM_pred"]=="D"
    subsunmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"]== "missense" & reann[,"RadialSVM_pred"]!="D"
    
    reall[subsunlof,"mutation_type_metaSVM"] <- "LOF"
    reall[subsunsyn,"mutation_type_metaSVM"] <- "SYNONMOUS"
    reall[subsunmis,"mutation_type_metaSVM"] <- "MIS"
    reall[subsundmis,"mutation_type_metaSVM"] <- "dMIS"
    reall <- reall[,c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.annovar","AAChange.annovar","Gene.annovar","Category_annovar","Polyphen2_HDIV_pred_annovar","mutation_type_annovar","Func.metaSVM","AAChange.metaSVM","Gene.metaSVM","Category_metaSVM","RadialSVM_pred_metaSVM","mutation_type_metaSVM")]
    
    write.csv(reall,file="DDD_mutations/datasheet/anno_mutationlist_4_15_all.csv",row.names=FALSE)
    
    reallnew1 <- reall[,c("Chr","Start","Ref","Alt","sampleID","Disease","Gene","AAChange.metaSVM","Category","Category_metaSVM","RadialSVM_pred_metaSVM","mutation_type_metaSVM")]
    write.csv(reallnew1,file="DDD_mutations/datasheet/anno_mutationlist_4_15.csv",row.names=FALSE)
        
}

Clear_misDDD_metaSVM <- function(){
    
    reann <- read.csv("DDD_mutations/datasheet/Annotated_DDDmutationlist_4_15_metaSVM.csv")
    reall <- read.csv("DDD_mutations/datasheet/DDDmutationlist_4_15_all.csv") 
    
    reall <- cbind(reall,reann[,c("Func.refgene","Gene.refgene","ExonicFunc.refgene","RadialSVM_pred","AAChange.refgene")])
    reall <- cbind(reall,"")
    colnames(reall) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.annovar","Gene.annovar","Category_annovar","Polyphen2_HDIV_pred_annovar","AAChange.annovar","mutation_type_annovar","Func.metaSVM","Gene.metaSVM","Category_metaSVM","RadialSVM_pred_metaSVM","AAChange.metaSVM","mutation_type_metaSVM")
    
    ## table(reann[,"ExonicFunc.refgene"])
    # c("","frameshift deletion","frameshift deletion;frameshift deletion","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV;nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV;nonsynonymous SNV;nonsynonymous SNV","stopgain","stoploss","synonymous SNV","synonymous SNV;synonymous SNV","unknown")
    
    #subs <- reann[,"ExonicFunc.refgene"]==""
    #table(reann[subs,"Func.refgene"])
    ##c("","downstream","intergenic","intronic","ncRNA_exonic","ncRNA_intronic","splicing","UTR3") 
    
    annlof <- c("frameshift deletion","frameshift deletion;frameshift deletion","frameshift insertion","stopgain","stoploss")
    subslof <- reann[,"ExonicFunc.refgene"] %in% annlof | (reann[,"ExonicFunc.refgene"]=="" & (reann[,"Func.refgene"] %in% c("splicing")))
    #subslof1  <- reann[,"ExonicFunc.refgene"]=="" & reann[,"Func.refgene"]=="" ## table(reall0[subslof1,"Category"])
    subsmisd <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift deletion",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift insertion",reann[,"ExonicFunc.refgene"])) & reann[,"RadialSVM_pred"]=="D"
    subsmis <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift deletion",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift insertion",reann[,"ExonicFunc.refgene"])) & reann[,"RadialSVM_pred"]!="D"
    subssyn <- grepl("synonymous SNV",reann[,"ExonicFunc.refgene"]) & !grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"])
    subsintron <- reann[,"ExonicFunc.refgene"]=="" & !(reann[,"Func.refgene"] %in% c("splicing"))
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    
    reall[subslof,"mutation_type_metaSVM"] <- "LOF"
    reall[subsmisd,"mutation_type_metaSVM"] <- "dMIS"
    reall[subsmis,"mutation_type_metaSVM"] <- "MIS"
    reall[subssyn,"mutation_type_metaSVM"] <- "SYNONMOUS"
    reall[subsintron,"mutation_type_metaSVM"] <- "Other"
    
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    subsunsyn <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"] %in% c("synonymous_variant")
    subsundmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"]== "missense_variant" & reann[,"RadialSVM_pred"]=="D"
    subsunmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"]== "missense_variant" & reann[,"RadialSVM_pred"]!="D"
    
    reall[subsunsyn,"mutation_type_metaSVM"] <- "SYNONMOUS"
    reall[subsunmis,"mutation_type_metaSVM"] <- "MIS"
    reall[subsundmis,"mutation_type_metaSVM"] <- "dMIS"
    reall <- reall[,c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.annovar","AAChange.annovar","Gene.annovar","Category_annovar","Polyphen2_HDIV_pred_annovar","mutation_type_annovar","Func.metaSVM","AAChange.metaSVM","Gene.metaSVM","Category_metaSVM","RadialSVM_pred_metaSVM","mutation_type_metaSVM")]
    
    write.csv(reall,file="DDD_mutations/datasheet/anno_DDDmutationlist_4_15_all.csv",row.names=FALSE)
    
    reallnew1 <- reall[,c("Chr","Start","Ref","Alt","sampleID","Disease","Gene","AAChange.metaSVM","Category","Category_metaSVM","RadialSVM_pred_metaSVM","mutation_type_metaSVM")]
    write.csv(reallnew1,file="DDD_mutations/datasheet/anno_DDDmutationlist_4_15.csv",row.names=FALSE)
    
}

clear_ncounttable <- function(){
    source("DDD.R")
    tmp <- clear_number_metaSVM("DDD_mutations/datasheet/anno_mutationlist_4_15.csv")
    write.table(tmp,file="DDD_mutations/datasheet/mutation_4_15_table1.txt",quote=FALSE,sep="\t")
    
    tmp <- clear_number_metaSVM("DDD_mutations/datasheet/anno_DDDmutationlist_4_15.csv")
    write.table(tmp,file="DDD_mutations/datasheet/DDDmutation_4_15_table1.txt",quote=FALSE,sep="\t")
}

clear_number_metaSVM <- function(filename){
    tmp <- read.csv(filename)
    
    genes <- unique(tmp[,"Gene"])
    dis <- unique(tmp[,"Disease"])
    II <- length(dis)
    Mut <- c("LOF","dMIS","MIS","SYNONMOUS")
    JJ <- length(Mut)
    
    numberT <- matrix(0,length(genes),II*JJ)
    rownames(numberT) <- genes
    colnames(numberT) <- rep(c("number_of_LGD","number_of_D-mis","number_of_B-mis","number_of_silent"),II)
    for(i in 1:II){
        for(j in 1:JJ){
            tmpdata <- tmp[tmp[,"Disease"]==dis[i] & tmp[,"mutation_type_metaSVM"]==Mut[j],]
            tmpg <- names(table(tmpdata[,"Gene"]))
            numberT[match(tmpg,genes),((i-1)*JJ + j)] <- table(tmpdata[,"Gene"])
        }
    }
    
    numberT
    
}

clear_number <- function(filename){
    tmp <- read.csv(filename)
    
    genes <- unique(tmp[,"Gene"])
    dis <- unique(tmp[,"Disease"])
    II <- length(dis)
    Mut <- c("LOF","dMIS","MIS","SYNONMOUS")
    JJ <- length(Mut)
    
    numberT <- matrix(0,length(genes),II*JJ)
    rownames(numberT) <- genes
    colnames(numberT) <- rep(c("number_of_LGD","number_of_D-mis","number_of_B-mis","number_of_silent"),II)
    for(i in 1:II){
        for(j in 1:JJ){
            tmpdata <- tmp[tmp[,"Disease"]==dis[i] & tmp[,"mutation_type"]==Mut[j],]
            tmpg <- names(table(tmpdata[,"Gene"]))
            numberT[match(tmpg,genes),((i-1)*JJ + j)] <- table(tmpdata[,"Gene"])
        }
    }
    
    numberT

}

