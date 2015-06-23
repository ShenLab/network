Multi_net <- function(netflag,strn="LBP",filenames,fileexp,modulefile,cutn=15000,idm=FALSE,beta,net_e=4){
    #netflag <- 1 #0: mix net ??; 1: expression; 2: regulatory ARACNE; 3: PPI; 4: Phenotype based; 5: GRN from CellNet;
    #filenames <- c("TADAinfo.txt","CNVinfo.txt","otherinfo.txt")
    options(expressions=500000)
    
    allnet <- build_net(netflag,fileexp,net_e)
    load(modulefile)
    module <- module[module[,1] %in% allnet$node,]
    #module <- build_module(allnet$matrix)
    nodeinfo <- build_node(module,1,filenames)
    POSpro <- build_CRF(nodeinfo,allnet,module,netflag,strn,cutn,beta,net_e)
    if(idm){
        mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
        pros1 <- POSpro
        pros1[,1] <- mapT[match(pros1[,1],mapT[,1]),2]
        pros1 <- pros1[!is.na(pros1[,1]),]
        write.table(pros1,file=paste(strn,"LBP_",netflag,"1.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    }
    
}

batch_Multi_net <- function(){
    
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 

    modulefiles <- c("module1","module3","iRefmodule","DAWNmodule","data/expressiondata/GTExbatch1modulem","data/expressiondata/GTExbatch2modulem","data/expressiondata/PCGCbatch1module","data/expressiondata/PCGCbatch2module","data/expressiondata/PCGCnewbornmodule","data/expressiondata/PCGCbatch1module_t","data/expressiondata/PCGCbatch2module_t")
    netstr <- c("coexp","STRING","iRef","Infer","coexp","coexp","coexp","coexp","coexp","coexp","coexp")
    
    for(netflag in c(3,6,11,12)){
        if(netflag==3) k=2;
        if(netflag==6) k=3;
        if(netflag==11) k=10;
        if(netflag==12) k=11;
        if(netflag==3 | netflag==1 | netflag>10){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        modulefile <- modulefiles[k]
        dirstr <- "result/PCGC/" 
        strname <- "PCGC"
        filename <- paste(dirstr,instr,strname,".txt",sep="")
        strn <- paste("result/PCGC/",netstr[k],"CRFresult_",strname,sep="")
        Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=15000,idm);          
        
    } # end for netflag
    filenames <- c("result/PCGC/STRINGCRFresult_PCGCLBP_31.txt","result/PCGC/iRefCRFresult_PCGCLBP_6.txt","result/PCGC/coexpCRFresult_PCGCLBP_111.txt","result/PCGC/coexpCRFresult_PCGCLBP_121.txt")
    
    rank_combine2(filenames[c(1,2,4)],"Com",flag=8,0.1)
    
    #flag=9
    #rank_combine(filenames,strn0,flag,0.1)
    #rank_combine1(paste(strn,"LBP_",netflag,"1.txt",sep=""),strn,netflag,0.1)
    #filenames <- union(filenames,paste(strn,"LBP_",netflag,".txt",sep=""))
    for(netflag in c(9,10,11,12,13)){
        a <- read.table(paste(strn,"LBP_",netflag,".txt",sep=""))
    if(netflag >10){
        a <- read.table(paste(strn,"LBP_",netflag,"1.txt",sep=""))
    }
    testgenes <- c("RBFOX2","KMT2D","EP300","NAA16","KDM5B","NCOA3","SMAD2","PTPN11","KDM5B","NAA15","POGZ","AHNAK","MYH6","NOTCH1","SUV420H1","MASTL","CHD7","BRAF","SMARCC1","CDKL2","JAG1","FLT4","PCNXL2","MKRN2","MROH7","TLK2","ARID1B","CTNNB1","ETS1","CUL3")
    print(match(testgenes,a[,1]))
    a1 <- match(testgenes,a[,1])
    print(sum(a1[!is.na(a1)]<200))
    }
}

batch_control <- function(){

    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 
    
    netstr <- c("STRING","iRef","Infer","coexp")
    k = 1 
    modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule")
    
    for(netflag in c(3,6,8,7)){
        
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        modulefile <- modulefiles[k]
        
            dirstr <- "result/"
            strname <- "control"
            filename <- paste(dirstr,instr,strname,".txt",sep="")
            strn <- paste("result/",netstr[k],"CRFresult_",strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm);
        k=k+1
    }
    
}

batch_DDD <- function(){
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "" 
    
    modulefiles <- "data/brain_expression/expressedmodule"
    netstr <- "coexp"
    netflag=14
    instr <- "hotnet_input"; idm=FALSE;
    k=1
    
    filenames <- c("DDD_mutations/datasheet/TADAdenovo_DDD.csv","DDD_mutations/datasheet/TADAdenovo_ID.csv","DDD_mutations/datasheet/TADAdenovo_EE.csv","DDD_mutations/datasheet/TADAdenovo_ASD.csv","DDD_mutations/datasheet/TADAdenovo_ALLde.csv")
    strname <- c("DDD","ID","EE","ASD","ALL")
    modulefile <- modulefiles[k]
    dirstr <- "DDD_mutations/netinput/" 
    for(i in 1:5){
        filename <- paste(dirstr,instr,strname[i],".txt",sep="")
        strn <- paste("result/DDD/","CRFresult",strname[i],sep="")
        Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=15000,idm);          
    }   
      
}

batch_ASD <- function(){
    
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 
    
    netstr <- c("STRING/","iRef/","Infer/","coexp/","GNAT/")
    k = 1 
    modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule","data/GNATnet/module_brain")
    
for(netflag in c(3,6,8,7)){
    
    if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    modulefile <- modulefiles[k]
    
    dirstr <- "result/randset4/"
    strname <- "ASD4_16"
    filename <- paste(dirstr,instr,strname,".txt",sep="")
    strn <- paste("result/randresult4_1/",netstr[k],"CRFresult_",strname,sep="")
    Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=0.5,net_e=4);
    
    ## random set 1: different sample size and exclude sample sets
#     for(j in 2:9){
#         for(i in 1:10){
#             strname <- paste("part",j,"_",i,sep="")
#             filename <- paste(dirstr,instr,strname,".txt",sep="")
#             strn <- paste("result/randresult3/",netstr[k],"CRFresult_",strname,sep="")
#             Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm);
            
#             strname <- paste("rest",j,"_",i,sep="")
#             filename <- paste(dirstr,instr,strname,".txt",sep="")
#             strn <- paste("result/randresult/",netstr[k],"CRFresult_",strname,sep="")
#             Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm);         
#        }
#    }

#     dirstr <- "result/randset2/" 
#     ## random set 3: different sample size and more samples
#     for(j in 1:3){
#         for(i in 1:5){
#             strname <- paste("part",j,"_",i,sep="")
#             filename <- paste(dirstr,instr,strname,".txt",sep="")
#             strn <- paste("result/randresult2/",netstr[k],"CRFresult_",strname,sep="")
#             Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm);       
#         }
#     }

    dirstr <- "result/randset4/" 
    ## random set 4: different sample size and more samples
    for(j in 3){
        for(i in 1:20){
            strname <- paste("part",j,"_",i,sep="")
            filename <- paste(dirstr,instr,strname,".txt",sep="")
            strn <- paste("result/randresult4/",netstr[k],"CRFresult_",strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=0);       
        }
    }
    
#     dirstr <- "result/leaveone1/"
#     ## random set 2: one missense mutations
#     for(i in 1:162){
#         strname <- paste("rand2_1_",i,sep="")
#         filename <- paste(dirstr,instr,"rand2_1_",i,".txt",sep="")
#         strn <- paste("result/leaveone1result/",netstr[k],"CRFresult_",strname,sep="")
#         Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm); 
#     }
    
    
#     dirstr <- "result/leaveone3/"
#     ## random set 3: one denovo mutation
#     for(i in 1:251){
#         strname <- paste("rand2_de1_",i,sep="")
#         filename <- paste(dirstr,instr,"rand2_de1_",i,".txt",sep="")
#         strn <- paste("result/leaveone3result/",netstr[k],"CRFresult_",strname,sep="")
#         Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm); 
#     }

    k <- k + 1
} # end for netflag


#     fileexp <- "PCGCall.txt" 
#     strn0 <- c("result/CRFresult_ASDall","result/CRFresult_ASD13772","result/CRFresult_ASD13908","result/CRFresult_ASD932")
#     fnode <- c("ASD/TADAresult/CRF_inputall.txt","ASD/TADAresult/CRF_input13772.txt","ASD/TADAresult/CRF_input13908.txt","ASD/TADAresult/CRF_input932.txt")
#     for(k in 1:length(fnode)){
#         for(netflag in c(6)){
#             strn <- paste(strn0[k],netflag,sep="_")
#             Multi_net(netflag,strn,fnode[k],fileexp,cutn=20000);
#         }
#     }

}

batch_test <- function(){

    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 
    
    netstr <- c("STRING/","iRef/","Infer/","coexp/")
    k = 4 
    modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule")
    
    for(netflag in c(7)){
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        modulefile <- modulefiles[k]
        
            dirstr <- "result/randset3/"
            strname <- "ASDall"
            filename <- paste(dirstr,instr,strname,".txt",sep="")
            strn <- paste("result/randresult3/",netstr[k],"CRFresult1_",strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm);
            filename <- paste(strn,"LBP_7.txt",sep="")
            totauc <- test_auc(filename)
    }
    
    auc1 <- matrix(0,10,3)
    netflag=7
    # random set 1: different sample size and exclude sample sets
            for(j in 2){
                for(i in 1:10){
                    strname <- paste("part",j,"_",i,sep="")
                    filename <- paste(dirstr,instr,strname,".txt",sep="")
                    strn <- paste("result/randresult3/",netstr[k],"CRFresult_",strname,sep="")
                    Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm);
                    
                    filename <- paste(strn,"LBP_7.txt",sep="")
                    auc1[i,] <- test_auc(filename)
               }
           }
    colMeans(auc1)
    totauc

    auc2 <- matrix(0,10,3)
    for(i in 1:10){
        filename <- paste("ASD/TADAresult/randset3/TADAdenovo_part2_",i,".csv",sep="")
        auc2[i,] <- test_auc(filename)
    }
    colMeans(auc2)
    
}

test_auc <- function(filename){
    library(ROCR)
    source("enrichana.R")
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= 0.1,1]
    
    if(grepl(".csv",filename)){ 
        result <- read.csv(filename); 
        predictions <- 1- result[,"qvalue.dn"];
    }else if(grepl(".txt",filename)){
        result <- read.table(filename);
        predictions <- result[,2]/max(result[,2]);
    }
    
    labels <- rep(-1,length(predictions))
    labels[match(Tset,result[,1])] <- 1
    pred <- prediction(predictions, labels)
    auc <- unlist(performance(pred, "auc")@y.values)
    subs1 <- match(Tset,result[,1])
    subn <- max(subs1[!is.na(subs1)])
    
    pred1 <- prediction(predictions[1:min(300,(2*subn))], labels[1:min(300,(2*subn))])
    auc1 <- unlist(performance(pred1, "auc")@y.values)
    
    
    c(auc,auc1,subn)
}

batch_subset_PCGC <- function(){

    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    
    fileexp <- "PCGCall.txt" 
    strn0 <- c("result/CRFresult_rand1","result/CRFresult_rand2","result/CRFresult_rand3")
    fnode <- c("random_samples/CRF_rand1.txt","random_samples/CRF_rand2.txt","random_samples/CRF_rand3.txt")
    for(k in 1:length(fnode)){
        for(netflag in c(6)){
            strn <- paste(strn0[k],netflag,sep="_")
            Multi_net(netflag,strn,fnode[k],fileexp,cutn=20000);
        }
    }
    
}

build_node <- function(module,flag=1,filenames=c("TADAinfo.txt","CNVinfo.txt","otherinfo.txt")){

    #mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t")) ###!!! 可以进一步精细
    genes <- module[,1]
    #mapT <- mapT[mapT[,1] %in% genes,]
    
    if(flag==1){
        nodeinfo <- matrix(0,length(genes),1,dimnames=list(genes,1))
        #nodeinfo[,1] <- 0.1 ## related genes percent estimate
        tmp <- read.table(filenames[1],stringsAsFactors=FALSE)
        tmp <- tmp[tmp[,1] %in% genes,]
        if(FALSE){
            tmp[,2] <- fit_mix_nor(tmp[,2]) ##!!!!!!!!!!!!!!!!!!!
        }
        nodeinfo[match(tmp[,1],rownames(nodeinfo)),1] <- tmp[,2]
    }else{
        n.info <- length(filenames)
        nodeinfo <- matrix(0,length(genes),n.info,dimnames=list(genes,1:n.info))
        #nodeinfo[,1] <- 0.1 ## related genes percent estimate
        #nodeinfo[,2:n.info] <- 0
        #nodeinfo[match(names(nodesim),rownames(nodeinfo)),1] <- nodesim
        for(i in 1:n.info){
            tmp <- read.table(filenames[i],stringsAsFactors=FALSE)
            tmp <- tmp[tmp[,1] %in% genes,]
            nodeinfo[match(tmp[,1],rownames(nodeinfo)),i] <- tmp[,2]
        }
    }
    nodeinfo
}

build_node1 <- function(module,net,filenames,netflag){
    
    genes <- module[,1]
    
    nodeinfo <- matrix(0,length(genes),1,dimnames=list(genes,1))
    #nodeinfo[,1] <- 0.1 ## related genes percent estimate
    tmp <- read.table(filenames[1],stringsAsFactors=FALSE)
    tmp <- tmp[tmp[,1] %in% genes,]
    nodeinfo[match(tmp[,1],rownames(nodeinfo)),1] <- tmp[,2]
    
    if(netflag==3){
        filename <- "data/Network_betweenness/Betweenness_node_STRING.txt"
    }else if(netflag==6){
        filename <- "data/Network_betweenness/Betweenness_node_iRef.txt"
    }else if(netflag==7){
        filename <- "data/Network_betweenness/Betweenness_node_coexp.txt"
    }else if(netflag==8){
        filename <- "data/Network_betweenness/Betweenness_node_Infer.txt"
    }
    
    TADAFile <- filenames[2]
    mutaT <- read.csv(TADAFile)
    if(netflag==3){
        mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
        mutaT[,1] <- mapT[match(mutaT[,1],mapT[,2]),1]
    }
    mutg <- mutaT[rowSums(mutaT[,c("dn.LoF","dn.mis3")]) > 0,1]
    
    tmp <- read.table(filename)
    mutg <- intersect(mutg,tmp[,1])
    theta <- sum(tmp[match(mutg,tmp[,1]),2]) / sum(tmp[,2])
    
    theta_1 <- pbinom(0, ceiling(tmp[,2]), theta, lower.tail = TRUE, log.p = FALSE)
    net$matrix[net$matrix>0] <- 1
    diag(net$matrix) <- 1
    P1 <- length(mutg)/sum(colSums(net$matrix[mutg,])>0)
    
    theta_2 <- rep(1,length(nodeinfo))
    genes <- intersect(genes,tmp[,1])
    theta_2[match(genes,rownames(nodeinfo))] <- theta_1[match(genes,tmp[,1])]
    
    nodeinfo <- (1-P1)*nodeinfo + (1-theta_2) * P1
    
    nodeinfo
}

build_net <- function(netflag,fileexp="PCGCall.txt",net_e=FALSE){
    # net: a list: size node matrix (systemtic and zero diag)
    net <- list()
    if(netflag==1 | netflag==2){
        datExpr <- as.matrix(read.delim(fileexp,sep="\t",row.names=1))
        datExpr <- t(datExpr)
        if(netflag==1){library(WGCNA); distM = adjacency(datExpr, type = "unsigned", power=1);}
        if(netflag==2){library(minet); distM <- build.mim(datExpr,estimator="spearman"); distM <- distM/max(distM);}
        genes <- colnames(datExpr)
        net$matrix <- distM
        diag(net$matrix) <- 0
        net$size <- length(genes)
        net$node <- genes
    }else if(netflag==3 | netflag==4 | netflag>=6){
        if(netflag==3){ filename <- "data/STRINGnetmap.txt";}
        if(netflag==4){ filename <- "data/phenotype_net.txt";}
        if(netflag==6){ filename <- "data/hotnet/iRefIndexm.txt";}
        if(netflag==7){ filename <- "data/network_inference/ASDexp_net_top5.txt";} ## ASD co-expression network
        if(netflag==8){ filename <- "data/network_inference/Infer_net_O4m.txt";} ## ASD inferred network
        if(netflag==9){ filename <- "data/expressiondata/GTExbatch1net_truned.txt";} ## GTEx batch 1
        if(netflag==10){ filename <- "data/expressiondata/GTExbatch2net_truned.txt";} ## GTEx batch 2
        if(netflag==11){ filename <- "data/expressiondata/PCGCbatch1net_t.txt";} ## PCGC batch 1
        if(netflag==12){ filename <- "data/expressiondata/PCGCbatch2net_t.txt";} ## PCGC batch 2
        if(netflag==13){ filename <- "data/expressiondata/PCGCnewbornnet_truned.txt";} ## PCGC new born
        if(netflag==14){ filename <- "data/brain_expression/expressednet.txt";} ## brain expressed net
        
        if(netflag==15){ filename <- "data/network_inference/Overlap_O4_top5.txt";} ## overlap coexp and infer
        if(netflag==16){ filename <- "data/network_inference/Joint_O4_top5.txt";} ## joint coexp and infer
        if(netflag==17){ filename <- "data/network_inference/Overlap_O4_top5_PPI.txt";} ## overlap coexp and infer + overlap PPI
        if(netflag==18){ filename <- "data/network_inference/Joint_O4_top5_PPI.txt";} ## joint coexp and infer + overlap PPI
        if(netflag==19){ filename <- "data/network_inference/OverlapPPI.txt";} ## overlap PPI
        
        if(netflag==20){ filename <- "data/StringNew_HPRDnet.txt";} ## MAGI
        if(netflag==21){ filename <- "data/GNATnet/network_8_brain.txt";} ## GNAT network brain
        if(netflag==22){ filename <- "data/GNATnet/network_26_Atrial.txt";} ## GNAT network Atrial
        if(netflag==23){ filename <- "data/GNATnet/network_30_left.txt";} ## GNAT network Left
        if(netflag==24){ filename <- "data/GNATnet/GMRFNet.txt";} ## Graphical lasso method 
        
        net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
        net.text <- rbind(net.text,net.text[,c(2,1,3)])
        net.text[,3] <- as.numeric(net.text[,3])/max(as.numeric(net.text[,3]))
        net <- read_net(net.text)
        diag(net$matrix) <- 0
    }else if(netflag==0){
        load("net0")
        net <- net0
    }else if(netflag==5){
        filename <- "GRN_net.txt";
        net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
        net.text[,3] <- abs(as.numeric(net.text[,3]))
        net.text <- rbind(net.text,net.text[,c(2,1,3)])
        net.text[,3] <- as.numeric(net.text[,3])/max(as.numeric(net.text[,3]))
        net <- read_net(net.text)
        diag(net$matrix) <- 0
        net$matrix <- (net$matrix + t(net$matrix))/2
    }
    # delete isolated node
    subs <- colSums(abs(net$matrix)) > 0
    net$node <- net$node[subs]
    net$matrix <- net$matrix[subs,subs]
    net$size <- sum(subs)

    ### add network edge betweenness 
    if(net_e){
        if(netflag==3){
            filename <- "data/Network_betweenness/Betweenness_edge_STRING.txt"
        }else if(netflag==6){
            filename <- "data/Network_betweenness/Betweenness_edge_iRef.txt"
        }else if(netflag==7){
            filename <- "data/Network_betweenness/Betweenness_edge_coexp.txt"
        }else if(netflag==8){
            filename <- "data/Network_betweenness/Betweenness_edge_Infer.txt"
        }else if(netflag==20){
            filename <- "data/Network_betweenness/Betweenness_edge_MAGI.txt"
        }else if(netflag==21){
            filename <- "data/Network_betweenness/Betweenness_edge_GNAT_brain.txt"
        }else if(netflag==22){
            filename <- "data/Network_betweenness/Betweenness_edge_GNAT_Atrial.txt"
        }else if(netflag==23){
            filename <- "data/Network_betweenness/Betweenness_edge_GNAT_Left.txt"
        }
        
        edgebe <- read.table(filename)
        pedge <- ecdf(edgebe[,4])
        
        edgebe <- edgebe[edgebe[,1] %in% net$node & edgebe[,3] %in% net$node,]
        we <- matrix(0,net$size,net$size, dimnames= list(net$node,net$node))
        subs1 <- match(edgebe[,1],net$node)
        subs2 <- match(edgebe[,3],net$node)
        we[cbind(subs1,subs2)] <- pedge(edgebe[,4])
        we[cbind(subs2,subs1)] <- pedge(edgebe[,4])
        net$we <- we
    }
    
    net
}

build_module <- function(distgene){
    library(WGCNA)
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft = pickSoftThreshold.fromSimilarity(distgene, powerVector = powers, verbose = 5)
    softPower <- sft$fitIndices[sft$fitIndices[,2]==max(sft$fitIndices[,2]),1] #sft$fitIndices[,1],-sign(sft$fitIndices[,3]*sft$fitIndices[,2]), #sft$fitIndices[,5]

    # WGCNA to module relationships
    #MEDissThres = quantile(distgene[upper.tri(distgene)],0.2)
    minModuleSize = 200; # set the minimum module size relatively high:
    diag(distgene) <- 1
    adjacency = adjacency.fromSimilarity(distgene, type = "unsigned", power=softPower)
    TOM = TOMsimilarity(adjacency);
    dissTOM = 1-TOM
    geneTree = flashClust(as.dist(dissTOM), method = "average"); # Call the hierarchical clustering function
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize); # Module identification using dynamic tree cut
    dynamicColors = labels2colors(dynamicMods) # Convert numeric lables into colors
    #merge = mergeCloseModules(distgene, dynamicColors, cutHeight = MEDissThres, verbose = 3) # Call an automatic merging function
    #moduleColors = merge$colors; # The merged module colors
    moduleColors = dynamicColors
    colorOrder = c("grey", standardColors());
    Labels = match(moduleColors, colorOrder)-1;
    module <- cbind(rownames(distgene),Labels)
    module
}

build_CRF <- function(nodeinfo,allnet,module,netflag,strn,cutn=100000,beta=0,net_e=4){
    library(CRF)
    library(Corbi)
    
    modLab <- setdiff(unique(module[,2]),0)
    n.module <- length(modLab)
    n.genes <- length(module[module[,2]!=0,1])
    
    POSpro <- matrix(0,n.genes,3)
    POSpro[,1] <- as.vector(module[module[,2]!=0,1]) # gene, probability; FDR, node score, edge score
    for (i in 1:n.module){
        modgenes <- as.vector(module[module[,2]==modLab[i],1])
        if(length(netflag)==1){
            net <- discrete_net(allnet,modgenes,netflag,cutn,nodeinfo,net_e)
        }else{
            net <- combine_net(allnet,modgenes,netflag)
        }
        if(sum(net$matrix)/2==0){break;}
        print(sum(net$matrix)/2)
        
        modnodeinfo <- nodeinfo[match(modgenes,rownames(nodeinfo)),]
        ### change 6_3
        model <- build_model_0(modnodeinfo,net,net_e,beta)
           
        crfresult <- infer_crf(model, query.type=4)
        POSpro[match(modgenes,POSpro[,1]),2:3] <- crfresult
        ###print(min(crfresult))
    }
    
    POSpro <- est.fdr(POSpro,flag=5)
    write.table(POSpro,file=paste(strn,"LBP_",netflag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
    POSpro
}

build_model_0 <- function(modnodesim,net,net_e,beta=0){   
    
    n.nf <- dim(as.matrix(modnodesim))[2]
    S <- modnodesim
    S <- cbind(S,1-S)
    colnames(S) <- 1:2
    
    n.states <- 2
    n.node <- net$size
    query.net <- net$matrix
    crf <- make.crf(net$matrix, rep(n.states,net$size))
    
    crf$state.map <- matrix(n.node, nrow=crf$n.nodes, ncol=crf$max.state)
    for (i in 1:crf$n.nodes)
    {
        crf$state.map[i, 1:crf$n.states[i]] <- 1:n.states
        crf$node.pot[i,] <- exp(S[i, crf$state.map[i,]])
    }   
    
    if(net_e>=3){
        we <- net$we
    }
    
    deg <- colSums(net$matrix)
    W1 <- matrix(1,n.states,n.states)
    for (e in 1:crf$n.edges){
        n1 <- crf$edges[e, 1]
        n2 <- crf$edges[e, 2]
        m1 <- 1:crf$n.states[n1]
        m2 <- 1:crf$n.states[n2]
        S1 <- matrix(log(crf$node.pot[n1, m1]), crf$n.states[n1], crf$n.states[n2])
        S2 <- matrix(log(crf$node.pot[n2, m2]), crf$n.states[n1], crf$n.states[n2], byrow=T)
        
        if(net_e==1){
            W1[1,1] <- (max(S1[1,1],S2[1,1])+beta)/sqrt(deg[n1]*deg[n2]) ## version 1
        }else if(net_e==2){
            W1[1,1] <- (max(S1[1,1],S2[1,1])*beta)/sqrt(deg[n1]*deg[n2]) ## version 2
        }else if(net_e==3){
            W1[1,1] <- (beta+max(S1[1,1],S2[1,1]))*we[n1,n2] ## version 3
        }else if(net_e==4){
            W1[1,1] <- beta*max(S1[1,1],S2[1,1])*we[n1,n2] ## version 4
        }
        
        W1[2,2] <- 0
        W1[1,2] <- 0
        W1[2,1] <- 0
        crf$edge.pot[[e]] <- exp(W1)
    }
    crf
}

rank_combine <- function(filenames,strn,flag=1,n){
    source("CRF_build.R")
    source("Network_analysis.R")
    
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    n.net <- length(filenames)
    genes <- c()
    for(i in 1:n.net){
        tmp <- read.table(filenames[i])
        genes <- union(genes,tmp[,1])
    }
    n.genes <- length(genes)
    
    scoreM <- matrix(0,n.genes,n.net+1,dimnames=list(genes,1:(n.net+1)))
    pM <- matrix(1,n.genes,n.net+2,dimnames=list(genes,1:(n.net+2)))
    for(i in 1:n.net){
        tmp <- read.table(filenames[i])
        scoreM[match(tmp[,1],genes),i] <- as.numeric(tmp[,2])
        scoreM[is.na(scoreM[,i]),i] <- 0
        
        pM[match(tmp[,1],genes),i] <- as.numeric(tmp[,3])
        pM[pM[,i]==0,i] <- 1e-50
        pM[is.na(pM[,i]),i] <- 1
    }
    
    if(n.net==4){
    csubs <- matrix(c(1,1,1,1,
                      1,1,1,2,
                      1,1,2,1,
                      1,2,1,1,
                      2,1,1,1,
                      1,1,2,2,
                      1,2,1,2,
                      2,1,2,1,
                      2,2,1,1,
                      2,1,1,2,
                      1,2,2,1,
                      2,2,1,2,
                      2,1,2,2,
                      2,2,2,1,
                      1,2,2,2,
                      2,2,2,2),64,1,byrow=TRUE)
    subs <- cbind(csubs,1:4)
    
    esubsn1 <- matrix(c(1,1,1,2,2,3),6,16)
    esubsn2 <- sapply(1:16,function(i) {tmp <- subs[((i-1)*4+1):(i*4),1]; tmp[esubsn1[,i]]})
    esubs1 <- cbind(matrix(esubsn2,96,1),matrix(esubsn1,96,1))
    esubsn1 <- matrix(c(2,3,4,3,4,4),6,16)
    esubsn2 <- sapply(1:16,function(i) {tmp <- subs[((i-1)*4+1):(i*4),1]; tmp[esubsn1[,i]]})
    esubs2 <- cbind(matrix(esubsn2,96,1),matrix(esubsn1,96,1))
    }
    if(n.net==3){
        csubs <- matrix(c(1,1,1,
                          1,1,2,
                          1,2,1,
                          2,1,1,
                          1,2,2,
                          2,1,2,
                          2,2,1,
                          2,2,2),24,1,byrow=TRUE)
        subs <- cbind(csubs,1:3)    
        
        esubsn1 <- matrix(c(1,1,2),3,8)
        esubsn2 <- sapply(1:8,function(i) {tmp <- subs[((i-1)*n.net+1):(i*n.net),1]; tmp[esubsn1[,i]]})
        esubs1 <- cbind(matrix(esubsn2,24,1),matrix(esubsn1,24,1))
        esubsn1 <- matrix(c(2,3,3),3,8)
        esubsn2 <- sapply(1:8,function(i) {tmp <- subs[((i-1)*n.net+1):(i*n.net),1]; tmp[esubsn1[,i]]})
        esubs2 <- cbind(matrix(esubsn2,24,1),matrix(esubsn1,24,1))
    }
    
    if(n.net==2){
        csubs <- matrix(c(1,1,
                          1,2,
                          2,1,
                          2,2),8,1,byrow=TRUE)
        subs <- cbind(csubs,1:2)    
        
        esubsn1 <- matrix(1,1,2^n.net)
        esubsn2 <- sapply(1:(2^n.net),function(i) {tmp <- subs[((i-1)*n.net+1):(i*n.net),1]; tmp[esubsn1[,i]]})
        esubs1 <- cbind(matrix(esubsn2,4,1),matrix(esubsn1,4,1))
        esubsn1 <- matrix(2,1,2^n.net)
        esubsn2 <- sapply(1:(2^n.net),function(i) {tmp <- subs[((i-1)*n.net+1):(i*n.net),1]; tmp[esubsn1[,i]]})
        esubs2 <- cbind(matrix(esubsn2,4,1),matrix(esubsn1,4,1))
    }

    library(matrixStats)
    
    if(flag==6){
        y1 <- sapply(1:n.genes,function(i){
                am <- rbind(scoreM[i,1:n.net],pM[i,1:n.net])
                subst <- am[1,]==0
                am[1,subst] <- ifelse(any(am[1,]>0),0.5,0)
                am[2,subst] <- 1 - am[1,subst]
                bm <- matrix(am[subs],2^n.net,n.net,byrow=TRUE)
                #exp(sum(log(bm[1,bm[1,]>0]))) # /sum(rowProds(bm, na.rm=FALSE, method="expSumLog"))
                #1 - sum(rowProds(bm[1,], na.rm=FALSE, method="expSumLog"))/sum(rowProds(bm, na.rm=FALSE, method="expSumLog"))
                #prod(bm[2^n.net,])/sum(rowProds(bm, na.rm=FALSE, method="expSumLog"))
                prod(bm[1,])/sum(rowProds(bm, na.rm=FALSE, method="expSumLog"))
            })
        y <- 1 - y1^(0.25)
    }else if(flag==7){
        pi0=0.94
        pi= 1 - pi0
        atmp <- pi0 * pi^n.net
        y1 <- sapply(1:n.genes,function(i){
            am <- scoreM[i,1:n.net]
            subst <- am==0
            am[subst] <- ifelse(any(am>0),0.5,0)
            (pi*prod(am))/(pi*prod(am) + atmp)
        })
        y <- 1 - y1
    }else if(flag==8){
        y1 <- sapply(1:n.genes,function(i){
        am <- rbind(scoreM[i,1:n.net],pM[i,1:n.net])
        subst <- am[1,]==0
        am[1,subst] <- ifelse(any(am[1,]>0),0.5,0)
        am[2,subst] <- 1 - am[1,subst]
        bm <- matrix(am[subs],2^n.net,n.net,byrow=TRUE)
        #cm <- matrix(am[esubs1],2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
        #dm <- matrix(am[esubs2],2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
        #fm <- matrix(esubs1[,1]==esubs2[,1],2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
        #em <- abs(abs(cm-dm) - fm)
        cm <- matrix(am[1,esubs1[,2]],2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
        dm <- matrix(am[1,esubs2[,2]],2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
        fm <- matrix((esubs1[,1]!=esubs2[,1])+1,2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
        em <- c(1,-1)[fm]
        gm <- 0.5 + em*sqrt(cm*dm)*(1-abs(cm-dm))/(n.net-1)
        sm <- cbind(bm,gm)
        prod(sm[1,])/sum(rowProds(sm, na.rm=FALSE, method="expSumLog"))    
        })
        y <- 1 - y1
    }
    
    pM[,n.net+1] <- y
    #y[y==0] <- 10^-20
    BF <- y1/y
    BF[BF==Inf] <- 2*max(BF[BF<Inf])
    #BF[is.na(BF)] <- 10^50
    pM[,n.net+2] <- BF
    scoreM[,n.net+1] <- y1
    
    result <- as.matrix(cbind(scoreM,pM))
    result <- result[order(as.numeric(pM[,n.net+2]),decreasing=TRUE),]
    
    FDR <- Bayesian.FDR(result[,2*n.net+3])$FDR
    result <- cbind(result,FDR)
    result <- result[order(FDR,decreasing=FALSE),]
    result0 <- result
    
    rownames(result) <- mapT[match(rownames(result),mapT[,1]),2]
    result <- result[!is.na(rownames(scoreM)),]

    filename <- paste(strn,flag,"resultall.txt",sep="_")
    filename0 <- paste(strn,flag,"resultall0.txt",sep="_")
    testname <- paste(strn,flag,"enrichall","1.txt",sep="_")    
    filetopgene <- paste(strn,flag,"topgeneall.txt",sep="_")
    write.table(result0,file=filename0,quote=FALSE,col.names=FALSE,sep="\t")
    write.table(result,file=filename,quote=FALSE,col.names=FALSE,sep="\t")
    
   # Output_file1(filename,mapT,testname,netflag,n)
    #testgenes <- rownames(result)[result[,2*n.net+2] < 0.1]
    #write.table(testgenes,file=filetopgene,quote=FALSE,col.names=FALSE,row.names=FALSE)
}

rank_combine1 <- function(filename,strn,netflag,n){
  
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    testname <- paste(strn,"enrich",netflag,"1.txt",sep="_") 
    Output_file1(filename,mapT,testname,netflag,n)

}

rank_combine2 <- function(filenames,strn,flag=1,n){
    source("CRF_build.R")
    source("Network_analysis.R")
    
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    n.net <- length(filenames)
    genes <- c()
    for(i in 1:n.net){
        tmp <- read.table(filenames[i])
        genes <- union(genes,tmp[,1])
    }
    n.genes <- length(genes)
    
    scoreM <- matrix(0,n.genes,n.net+1,dimnames=list(genes,1:(n.net+1)))
    pM <- matrix(1,n.genes,n.net+1,dimnames=list(genes,1:(n.net+1)))
    for(i in 1:n.net){
        tmp <- read.table(filenames[i])
        scoreM[match(tmp[,1],genes),i] <- as.numeric(tmp[,2])
        scoreM[is.na(scoreM[,i]),i] <- 0
        
        pM[match(tmp[,1],genes),i] <- as.numeric(tmp[,3])
        pM[pM[,i]==0,i] <- 1e-50
        pM[is.na(pM[,i]),i] <- 1
    }
    
    if(n.net==4){
        csubs <- matrix(c(1,1,1,1,
                          1,1,1,2,
                          1,1,2,1,
                          1,2,1,1,
                          2,1,1,1,
                          1,1,2,2,
                          1,2,1,2,
                          2,1,2,1,
                          2,2,1,1,
                          2,1,1,2,
                          1,2,2,1,
                          2,2,1,2,
                          2,1,2,2,
                          2,2,2,1,
                          1,2,2,2,
                          2,2,2,2),64,1,byrow=TRUE)
        subs <- cbind(csubs,1:4)
        
        esubsn1 <- matrix(c(1,1,1,2,2,3),6,16)
        esubsn2 <- sapply(1:16,function(i) {tmp <- subs[((i-1)*4+1):(i*4),1]; tmp[esubsn1[,i]]})
        esubs1 <- cbind(matrix(esubsn2,96,1),matrix(esubsn1,96,1))
        esubsn1 <- matrix(c(2,3,4,3,4,4),6,16)
        esubsn2 <- sapply(1:16,function(i) {tmp <- subs[((i-1)*4+1):(i*4),1]; tmp[esubsn1[,i]]})
        esubs2 <- cbind(matrix(esubsn2,96,1),matrix(esubsn1,96,1))
    }
    if(n.net==3){
        csubs <- matrix(c(1,1,1,
                          1,1,2,
                          1,2,1,
                          2,1,1,
                          1,2,2,
                          2,1,2,
                          2,2,1,
                          2,2,2),24,1,byrow=TRUE)
        subs <- cbind(csubs,1:3)    
        
        esubsn1 <- matrix(c(1,1,2),3,8)
        esubsn2 <- sapply(1:8,function(i) {tmp <- subs[((i-1)*n.net+1):(i*n.net),1]; tmp[esubsn1[,i]]})
        esubs1 <- cbind(matrix(esubsn2,24,1),matrix(esubsn1,24,1))
        esubsn1 <- matrix(c(2,3,3),3,8)
        esubsn2 <- sapply(1:8,function(i) {tmp <- subs[((i-1)*n.net+1):(i*n.net),1]; tmp[esubsn1[,i]]})
        esubs2 <- cbind(matrix(esubsn2,24,1),matrix(esubsn1,24,1))
    }
    
    if(n.net==2){
        csubs <- matrix(c(1,1,
                          1,2,
                          2,1,
                          2,2),8,1,byrow=TRUE)
        subs <- cbind(csubs,1:2)    
        
        esubsn1 <- matrix(1,1,2^n.net)
        esubsn2 <- sapply(1:(2^n.net),function(i) {tmp <- subs[((i-1)*n.net+1):(i*n.net),1]; tmp[esubsn1[,i]]})
        esubs1 <- cbind(matrix(esubsn2,4,1),matrix(esubsn1,4,1))
        esubsn1 <- matrix(2,1,2^n.net)
        esubsn2 <- sapply(1:(2^n.net),function(i) {tmp <- subs[((i-1)*n.net+1):(i*n.net),1]; tmp[esubsn1[,i]]})
        esubs2 <- cbind(matrix(esubsn2,4,1),matrix(esubsn1,4,1))
    }
    
    library(matrixStats)
    
    if(flag==6){
        y1 <- sapply(1:n.genes,function(i){
            am <- rbind(scoreM[i,1:n.net],pM[i,1:n.net])
            subst <- am[1,]==0
            am[1,subst] <- ifelse(any(am[1,]>0),0.5,0)
            am[2,subst] <- 1 - am[1,subst]
            bm <- matrix(am[subs],2^n.net,n.net,byrow=TRUE)
            #exp(sum(log(bm[1,bm[1,]>0]))) # /sum(rowProds(bm, na.rm=FALSE, method="expSumLog"))
            #1 - sum(rowProds(bm[1,], na.rm=FALSE, method="expSumLog"))/sum(rowProds(bm, na.rm=FALSE, method="expSumLog"))
            #prod(bm[2^n.net,])/sum(rowProds(bm, na.rm=FALSE, method="expSumLog"))
            prod(bm[1,])/sum(rowProds(bm, na.rm=FALSE, method="expSumLog"))
        })
        y <- 1 - y1^(0.25)
    }else if(flag==7){
        pi0=0.94
        pi= 1 - pi0
        atmp <- pi0 * pi^n.net
        y1 <- sapply(1:n.genes,function(i){
            am <- scoreM[i,1:n.net]
            subst <- am==0
            am[subst] <- ifelse(any(am>0),0.5,0)
            (pi*prod(am))/(pi*prod(am) + atmp)
        })
        y <- 1 - y1
    }else if(flag==8){
        y1 <- sapply(1:n.genes,function(i){
            am <- rbind(scoreM[i,1:n.net],pM[i,1:n.net])
            subst <- am[1,]==0
            am[1,subst] <- ifelse(any(am[1,]>0),0.5,0)
            am[2,subst] <- 1 - am[1,subst]
            bm <- matrix(am[subs],2^n.net,n.net,byrow=TRUE)
            #cm <- matrix(am[esubs1],2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
            #dm <- matrix(am[esubs2],2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
            #fm <- matrix(esubs1[,1]==esubs2[,1],2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
            #em <- abs(abs(cm-dm) - fm)
            cm <- matrix(am[1,esubs1[,2]],2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
            dm <- matrix(am[1,esubs2[,2]],2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
            fm <- matrix((esubs1[,1]!=esubs2[,1])+1,2^n.net,(n.net*(n.net-1))/2,byrow=TRUE)
            em <- c(1,-1)[fm]
            gm <- 0.5 + em*sqrt(cm*dm)*(1-abs(cm-dm))/(n.net-1)
            sm <- cbind(bm,gm)
            prod(sm[1,])/sum(rowProds(sm, na.rm=FALSE, method="expSumLog"))    
        })
        y <- 1 - y1
    }
    
    pM[,n.net+1] <- y
    scoreM[,n.net+1] <- y1
    
    result <- as.matrix(cbind(scoreM,pM))
    result <- result[order(as.numeric(scoreM[,n.net+1]),decreasing=TRUE),]
    
    pnew <- fit_mix_nor(scoreM[,n.net+1],k=1)
    pnew <- p.adjust(pnew,method="fdr")
    result <- cbind(result,pnew)

    filename0 <- paste(strn,flag,"resultall0.txt",sep="_")
    write.table(result,file=filename0,quote=FALSE,col.names=FALSE,sep="\t")

}
#====================
# auxiliary function
#====================
# tmp used from CRF_build.R

combine_net <- function(allnet,modgenes,netflag){


}

discrete_net <- function(allnet,modgenes,netflag,cutn=20000,nodeinfo,net_e){
    net <- list()
    net$node <- modgenes
    net$size <- length(modgenes)
    net$matrix <- allnet$matrix[modgenes,modgenes]
    if(net_e){net$we <- allnet$we[modgenes,modgenes];}
    tmp <- net$matrix
    
    if(netflag==1){
        cutf <- quantile(net$matrix[upper.tri(net$matrix)],0.95)
        net$matrix[net$matrix >= cutf] <- 1
        net$matrix[net$matrix < cutf] <- 0
    }
    
    if(netflag==2){
        neta <- minet::aracne(net$matrix)
        net$matrix <- neta[modgenes,modgenes]
        net$matrix[net$matrix>0] <- 1
    } 
    
    if(netflag==3 | netflag==4 | netflag==0 | netflag>=5){
        net$matrix[net$matrix > 0] <- 1
        #net$matrix[net$matrix != 0] <- 1
    }
    
    nodset <- nodeinfo[match(modgenes,rownames(nodeinfo))]
    if(sum(net$matrix==1) > 2*cutn){ # memory limit
        net$matrix[lower.tri(net$matrix,diag=TRUE)] <- 0
        edges <- which(net$matrix==1,arr.ind=TRUE)
        tmpS <- min(nodset[edges[,1]],nodset[edges[,2]])
        edges <- cbind(edges,tmpS)
        #edges <- as.matrix(cbind(edges,tmp[edges]))
        edges <- edges[order(as.numeric(edges[,3]),decreasing=TRUE),]
        edges <- edges[1:cutn,1:2]
        net$matrix[upper.tri(net$matrix,diag=TRUE)] <- 0
        net$matrix[edges] <- 1
        net$matrix[edges[,c(2,1)]] <- 1        
    }
    
#     #if(netflag==1){
#     mindeg <- 2 # eliminate isolated nodes
#     subs <- which(rowSums(net$matrix) == 0)
#     if(length(subs)>0){
#     #tmp <- allnet$matrix[modgenes,modgenes]
#         for (i in 1:length(subs)){
#             net$matrix[subs[i],sort(tmp[subs[i],],decreasing=TRUE,index.return=TRUE)$ix[1:mindeg]] <- 2
#             net$matrix[sort(tmp[subs[i],],decreasing=TRUE,index.return=TRUE)$ix[1:mindeg],subs[i]] <- 2
#         }
#     }
#     
#     if(sum(net$matrix >= 1) > 2*cutn){ # memory limit
#         edges <- which(net$matrix==1,arr.ind=TRUE)
#         edges <- edges[edges[,1] < edges[,2],]
#         edges <- as.matrix(cbind(edges,tmp[edges]))
#         edges <- edges[order(as.numeric(edges[,3]),decreasing=FALSE),]
#         ntmp <- min( (sum(net$matrix >= 1))/2 - cutn, dim(edges)[1])
#         edges <- edges[1:ntmp,1:2]
#         net$matrix[edges] <- 0
#         net$matrix[edges[,c(2,1)]] <- 0       
#     }
#     
#     net$matrix[net$matrix > 1] <- 1
    #}
    
    diag(net$matrix) <- 0
    
    net
}

combine_allnet <- function(){
    source("Multi_net.R")
    library(WGCNA)
    library(minet)
    fileexp <- c("PCGCall.txt","Brain_expression_Ref.txt")
    k=1
    # net: a list: size node matrix (systemtic and zero diag)
    netlist <- list()
    nodes <- c()
   
    netlist[[1]] <- build_net(1,fileexp[k])
    cutf <- quantile(netlist[[1]]$matrix[upper.tri(netlist[[1]]$matrix)],0.999)
    netlist[[1]]$matrix[netlist[[1]]$matrix >= cutf] <- 1
    netlist[[1]]$matrix[netlist[[1]]$matrix < cutf] <- 0
    print("Network1")
    
    netlist[[2]] <- build_net(2,fileexp[k])
    neta <- minet::aracne(netlist[[2]]$matrix)
    netlist[[2]]$matrix <- neta
    netlist[[2]]$matrix[netlist[[2]]$matrix > 0] <- 1
    print("Network2")
    
    netlist[[3]] <- build_net(3,fileexp[k])
    netlist[[3]]$matrix[netlist[[3]]$matrix > 0] <- 1
    netlist[[4]] <- build_net(4,fileexp[k])
    netlist[[4]]$matrix[netlist[[4]]$matrix > 0] <- 1
    print("Network4")
    
    nodes <- union(nodes,netlist[[1]]$node)
    nodes <- union(nodes,netlist[[2]]$node)
    nodes <- union(nodes,netlist[[3]]$node)
    nodes <- union(nodes,netlist[[4]]$node)
      

    net0 <- list()
    net0$node <- nodes
    net0$size <- length(nodes)
    net0$matrix <- matrix(0,net0$size,net0$size,dimnames=list(nodes,nodes))
    for(i in 1:4){
        net0$matrix[netlist[[i]]$node,netlist[[i]]$node] <- net0$matrix[netlist[[i]]$node,netlist[[i]]$node] + netlist[[i]]$matrix
    }
    rm(netlist)
    print("RM") 
    net0$matrix[net0$matrix < 2] <- 0
    subs <- rowSums(net0$matrix) > 0
    net0$node <- net0$node[subs]
    net0$size <- sum(subs)
    net0$matrix <- net0$matrix[subs,subs]
    net0$matrix <- net0$matrix/max(net0$matrix)
    print("Saving")
    save(net0,file="net0");
    
    
    Labels <- build_module(net0$matrix)
    module <- cbind(net0$node,Labels)
    print("Module done")
    save(module,file="module0");
    
}

Bayesian.FDR <- function(BF, pi0=0.94, alpha=0.05) {
    # Bayesian FDR control (PMID:19822692, Section2.3)
    # BF: a sorted vector of BFs (in decreasing order)
    # pi0: the prior probability that the null model is true
    # alpha: the FDR target
    # Return: the q-value of each BF, and the number of findings with q below alpha. 
    # convert BFs to PPA (posterior probability of alternative model)
    
    BF <- as.numeric(BF)
    BF0 <- BF
    BF <- sort(BF,decreasing = TRUE)
    
    pi <- 1-pi0
    q <- pi*BF/(1-pi+pi*BF) # PPA
    q0 <- 1 - q # posterior probability of null model
    
    # the FDR at each PPA cutoff
    n <- length(BF)
    FDR <- numeric(n)
    for (i in 1:n) FDR[i] <- sum(q0[1:i]) / i 

    # the cutoff
    #t <- 1
    #while (t <= length(q0) & mean(q0[1:t]) <= alpha) { t <- t+1 }
    #return (list(FDR=FDR, ND=t))

    rankp=rank(BF0,ties.method="random")
    FDR=FDR[rankp]
    
    list(FDR=FDR)
}

Posterior.FDR <- function(P0){
    #P0 is the posterior probability for null model
    P0 <- as.numeric(P0)
    rankpost=sort(P0)
    n=length(P0)
    localfdr=rep(0,n)
    for (i in 1:n){
        localfdr[i]=mean(rankpost[1:i])
    }
    
    fdr=rep(0,n)
    rankp=rank(P0,ties.method="random")
    fdr=localfdr[rankp]
    
    fdr
}

estimate.FDR <- function(P0){
    #P0 is the p value for each test
    P0 <- as.numeric(P0)
    Ppost=sort(P0)
    n=length(P0)
    Cp <- n*Ppost*1/(1:n)
    localfdr=rep(0,n)
    for (i in 1:n){
        localfdr[i]=min(Cp[i:n])
    }
    
    fdr=rep(0,n)
    rankp=rank(P0,ties.method="random")
    fdr=localfdr[rankp]
    
    fdr
}

TADAinput <- function(){
    library(mixtools)
    source("Network_analysis.R")
    datdir <- "random_samples/TADA_results"
    strname <- c("NDD","nonNDD","other","CHD","chdother","P2","T2","P3","T3")
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    for(i in 1:length(strname)){
        TADAFile <- paste(datdir,strname[i],".csv",sep="")
        geneinfo <- as.matrix(read.csv(TADAFile))
        genelist <- geneinfo[,1]
        genelist <- mapping_to(genelist)
        genelist <- mapT[match(genelist,mapT[,2]),1]
        pi0 <- 0.85 #!!!!!!
        pi <- 1-pi0
        nodesim <- as.matrix((as.numeric(geneinfo[,"BF"])*pi)/(as.numeric(geneinfo[,"BF"])*pi + 1-pi))
        names(nodesim) <- genelist
        subs <- !is.na(genelist)
        nodesim <- nodesim[subs]
        infofile <- paste("TADAinfo0210",strname[i],".txt",sep="")
        write.table(cbind(names(nodesim),nodesim),file=infofile,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
        nodes  <- as.matrix(read.table(infofile))
        nodesV <- as.numeric(nodes[,2])
        nodesV <- fit_mix_nor(nodesV)
        nodes2 <- nodes
        nodes2[,2] <- nodesV
        infofile <- paste("TADAinfo0210",strname[i],"_c.txt",sep="")
        write.table(nodes2,file=infofile,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    }

    
    ## mutation rand set 1 2 3    
    source("Network_analysis.R")
    source("Multi_net.R")
    datdir <- "random_samples/"
    filenames <- c("TADA_resultsrand1.csv","TADA_resultsrand2.csv","TADA_resultsrand3.csv")
    filenames <- paste(datdir,filenames,sep="")
    strname <- c("rand1","rand2","rand3")
    TADAinput1(filenames,strname)
    
    filenames <- c("TADA_resultsPCGC_2_4.csv")
    strname <- "PCGC"
    dirstr <- ""
    TADAinput1(filenames,strname,dirstr)
}

TADAinput1 <- function(filenames,strname,dirstr="random_samples/",pi0=0.94,genelist0=""){
    
    if(genelist0==""){
        bfs <- "BF"
    }else{
        bfs <- "BF.dn"
    }    
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    for(i in 1:length(filenames)){
        TADAFile <- filenames[i]
        geneinfo <- read.csv(TADAFile,stringsAsFactors=FALSE)
        pi <- 1-pi0
 
        posP <- (geneinfo[,bfs]*pi)/(geneinfo[,bfs]*pi + 1-pi)
        infofile <- paste(dirstr,"hotnet_input",strname[i],".txt",sep="")
        write.table(cbind(geneinfo[,1],posP),file=infofile,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
        
        genelist <- geneinfo[,1]
        if(genelist0==""){
            genelist <- mapping_to(genelist)
        }
        genelist <- mapT[match(genelist,mapT[,2]),1]
        nodesim <- (geneinfo[,bfs]*pi)/(geneinfo[,bfs]*pi + 1-pi)
        names(nodesim) <- genelist
        subs <- !is.na(genelist)
        nodesim <- nodesim[subs]
        infofile <- paste(dirstr,"CRF_input",strname[i],".txt",sep="")
        write.table(cbind(names(nodesim),nodesim),file=infofile,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    }
    
}

fit_mix_nor <- function(p1,k=2,plot=FALSE){
    #library(nor1mix) #re <- norMixEM(as.vector(p1),k)
    p1 <- as.numeric(p1)
    if(k>1){
        library(mixtools)
        re <- normalmixEM(p1,k) 
        paM <- matrix(c(re$mu,re$sigma,re$lambda),k,3)
        paM <- paM[order(paM[,1]),]
        vaM <- c()
        for(i in 1:k){
            vaM <- cbind(vaM,dnorm(p1,mean=paM[k,1],sd=paM[k,2]) * paM[k,3])
        }
        pnew <- pnorm(p1,mean=paM[k,1],sd=paM[k,2])*exp(log(vaM[,k])-log(rowSums(vaM)))
    }else if(k==1){
        pnew <- pnorm(p1,mean=mean(p1),sd=sd(p1),lower.tail=FALSE) # + pnorm(2*mean(p1) - p1,mean=mean(p1),sd=sd(p1),lower.tail=TRUE)
    }
    
    if(plot){
        library(vioplot); vioplot(p1,pnew,names=c("p0","p1"))
    }
    pnew

}

est.fdr <- function(POSpro,flag=4){
    
    if(flag==1){
        library(fdrtool)
        p1 <- as.numeric(POSpro[,2])
        fdrest <- fdrtool((p1-mean(p1))/sd(p1),statistic="normal",plot=FALSE)
        POSpro <- cbind(POSpro,cbind(fdrest$pval,fdrest$qval))
    }else if(flag==2){
        pnew <- fit_mix_nor(POSpro[,2],k=1)
        POSpro <- as.matrix(cbind(POSpro,pmin(pnew * dim(POSpro)[1],1)))
        POSpro[,3] <- 1 - as.numeric(POSpro[,4])
        POSpro <- POSpro[order(as.numeric(POSpro[,2]),decreasing = TRUE),]
        b <- as.numeric(POSpro[,4])
        BF <- as.numeric(POSpro[,3])/b
        BF[BF==Inf] <- 2*max(BF[BF<Inf])
        FDR <- Bayesian.FDR(BF)$FDR 
        POSpro <- cbind(POSpro,FDR)
    }else if(flag==3){
        pnew <- fit_mix_nor(POSpro[,2],k=2)
        POSpro[,2] <- pnew
        POSpro[,3] <- 1-pnew
        POSpro <- POSpro[order(as.numeric(POSpro[,2]),decreasing = TRUE),]
        POSpro <- as.matrix(cbind(POSpro,as.numeric(POSpro[,3]) * dim(POSpro)[1]))
        POSpro[as.numeric(POSpro[,4])>1,4] <- 1
        b <- as.numeric(POSpro[,3])
        BF <- as.numeric(POSpro[,2])/b
        BF[BF==Inf] <- 2*max(BF[BF<Inf])
        FDR <- Bayesian.FDR(BF)$FDR 
        POSpro <- cbind(POSpro,FDR)    
    }else if(flag==4){
        FDR <- Posterior.FDR(POSpro[,3])
        FDR[as.numeric(POSpro[,3])==1] <- 1 ##!!!!!
        pnew <- p.adjust(as.numeric(POSpro[,3]),method="fdr")
        POSpro <- cbind(POSpro,pnew)
        POSpro <- cbind(POSpro,FDR)
    }else if(flag==5){
        p1 <- as.numeric(POSpro[,2])
        pnew <- fit_mix_nor(p1,k=1)
        POSpro <- cbind(POSpro,pnew)
        pnew <- p.adjust(pnew,method="fdr")
        POSpro <- cbind(POSpro,pnew)
    }else if(flag==6){
        p1 <- as.numeric(POSpro[,2])
        pnew <- sampling_p(p1)
        POSpro <- cbind(POSpro,pnew)
        pnew <- p.adjust(pnew,method="fdr")
        POSpro <- cbind(POSpro,pnew)    
    }
    POSpro <- POSpro[order(as.numeric(POSpro[,2]),decreasing=TRUE),]
    POSpro
    
}

sampling_p <- function(p,si=10^8){
    
#     y <- spline(p,n=10*length(p),xmin=0,xmax=1)$y
#     edd <- ecdf(p)
#     pnew <- 1 - edd(p)
    
#     den <- density(p,n=1024,from=0,to=1)
#     sam <- sample(den$x, size=si, replace = TRUE, prob = den$y)
#     n <- length(den$x)-1
#     inv <- den$x[2:(n+1)] - den$x[1:n]
#     linv <- max(inv)/2
#     sam1 <- sam + runif(si,min=-linv,max=linv)
    
#     sam1 <- c()
#     for(i in 1:n){
#         subs <- sam==den$x[i]
#         subn <- sum(subs)
#         tmp <- runif(subn,min=den$x[i],max=den$x[i+1])
#         sam1 <- c(sam1,tmp)
#         print(i)
#     }
#     subs <- sam==den$x[n+1]
#     subn <- sum(subs)
#     tmp <- runif(subs,min=den$x[n],max=den$x[n+1])
#     sam1 <- c(sam1,tmp)

#     pnew <- sapply(1:length(p), function(i) sum(sam1>p[i]))
#     pnew <- pnew/si
    pnew
}

net_combine <- function(){
        
    source("Multi_net.R")
    source("CRF_build.R")
    fileexp=""
    ppistring <- build_net(netflag=3,fileexp) #237981
    ppiiref <- build_net(netflag=6,fileexp) #91807
    ## two ppi network overlap
    netpair_overlap(ppistring,ppiiref,idmap=TRUE,plot=TRUE)
    # overlap nodes: 10284; overlap edges: 38152; smaller net edges: 75749
    
    asdcoexp <- build_net(netflag=7,fileexp) #39966
    asdinfer <- build_net(netflag=8,fileexp) #65011
    ## two co-expression network overlap
    netpair_overlap(asdcoexp,asdinfer,idmap=FALSE,plot=TRUE)
    # overlap nodes: 9157; overlap edges: 20798; smaller net edges: 38476
    
    #### different co-expression networks
    pcgccoexp <- build_net(netflag=12,fileexp)
    netpair_overlap(pcgccoexp,asdcoexp,idmap=TRUE,plot=TRUE)
    # overlap nodes: 6985; overlap edges: 159; smaller net edges: 20627
    
    braincoexp <- build_net(netflag=14,fileexp)
    netpair_overlap(asdcoexp,braincoexp,idmap=FALSE,plot=TRUE)
    # overlap nodes: 6813; overlap edges: 451; smaller net edges: 17904
       
    netpair_overlap(pcgccoexp,braincoexp,idmap=TRUE,plot=TRUE)
    # overlap nodes: 8962; overlap edges: 246; smaller net edges: 29240
    
    ### PPI and co-expression networks
    netpair_overlap(ppiiref,asdcoexp,idmap=FALSE,plot=TRUE)
    # overlap nodes: 6918; overlap edges: 111; smaller net edges: 22682
    netpair_overlap(ppiiref,braincoexp,idmap=FALSE,plot=TRUE)
    # overlap nodes: 8289; overlap edges: 145; smaller net edges: 25378
    netpair_overlap(pcgccoexp,ppiiref,idmap=TRUE,plot=FALSE)
    # overlap nodes: 8946; overlap edges: 204; smaller net edges: 33274
    netpair_overlap(ppistring,asdcoexp,idmap=TRUE,plot=FALSE)
    # overlap nodes: 7331; overlap edges: 586; smaller net edges: 24666
    netpair_overlap(ppistring,braincoexp,idmap=TRUE,plot=FALSE)
    # overlap nodes: 8331; overlap edges: 830; smaller net edges: 24983

    ###===============================================================###
    ### four combined methods
    ## 1) overlap co-exp and infer
    asdcoexp <- build_net(netflag=7,fileexp) #39966
    asdinfer <- build_net(netflag=8,fileexp) #65011
    netpair_merge(asdcoexp,asdinfer,flag=1,filename="data/network_inference/Overlap_O4_top5.txt",idmap=FALSE)
    ## 2) joint co-exp and infer
    netpair_merge(asdcoexp,asdinfer,flag=2,filename="data/network_inference/Joint_O4_top5.txt",idmap=FALSE)
    ## 3) overlap PPI + overlap co-exp and infer::  expected
    source("Multi_net.R")
    source("CRF_build.R")
    fileexp=""
    ppistring <- build_net(netflag=3,fileexp) 
    ppiiref <- build_net(netflag=6,fileexp) 
    netpair_merge(ppistring,ppiiref,flag=1,filename="data/network_inference/OverlapPPI.txt",idmap=TRUE)
   
    overcoexp <- build_net(netflag=15,fileexp) 
    overppi <- build_net(netflag=19,fileexp) 
    netpair_merge(overcoexp,overppi,flag=2,filename="data/network_inference/Overlap_O4_top5_PPI.txt",idmap=FALSE)
    ## 4) overlap PPI + joint co-exp and infer
    jointcoexp <- build_net(netflag=16,fileexp) 
    netpair_merge(jointcoexp,overppi,flag=2,filename="data/network_inference/Joint_O4_top5_PPI.txt",idmap=FALSE)
    
    a <- read.table("data/StringNew_HPRD.txt")
    a <- cbind(a,1)
    write.table(a,file="data/StringNew_HPRDnet.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
    ppimagi <- build_net(netflag=20,fileexp) #52801
    
    genes <- intersect(rownames(distM),ppimagi$node)
    tmp <- distM[genes,genes]
    cormagi <- tmp[ppimagi$matrix[genes,genes]==1]
    
    a <- unlist(read.table("plot/test_output.txt",sep=""))
    a <- a[2:51]
    genem <- intersect(genes,a)
    tmp <- distM[genem,genem]
    cormod <- tmp[ppimagi$matrix[genem,genem]==1]
    plot(density(cormagi),xlim=c(0,1),ylim=c(0,3),xlab="Pearson Correlation",main="MAGI PPI and its best module",col=1,lwd=2)
    lines(density(cormod),col=2,lwd=2)
    legend("topright",col=1:2,lty=c(1,1),legend=c("PPI","Best_module"))
    
    
    ### correlation distribution on overlapped edges    
    load("geneexp.Rdata")
    library(WGCNA)
    distM = adjacency(t(data), type = "unsigned", power=1);
    cors1 <- netpair_correlation(asdcoexp,asdinfer,distM,idmap=FALSE,plot=TRUE)
    cors <- netpair_correlation(ppistring,ppiiref,distM,idmap=TRUE,plot=FALSE)
    plot(density(cors1),xlim=c(0,1),ylim=c(0,5),xlab="Pearson Correlation",main="",col=2,lwd=2)
    lines(density(cors),col=1,lwd=2)
    legend("topleft",col=1:2,lty=c(1,1),legend=c("STRING & iRef","Co-exp & Infer"))
    
}

netpair_correlation <- function(net1,net2,distM,idmap=FALSE,plot=TRUE){

    if(idmap){
        mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
        ppinode <- mapT[match(net1$node,mapT[,1]),2]
        net1$node <- ppinode
        colnames(net1$matrix) <- ppinode
        rownames(net1$matrix) <- ppinode
    }
    
    nodes  <- intersect(net1$node,net2$node)
    
    M1 <- net1$matrix[nodes,nodes]
    M2 <- net2$matrix[nodes,nodes]
    
    M1[M1>0] <- 1
    M2[M2>0] <- 1
    M <- M1 + M2
    
    genes <- intersect(rownames(distM),nodes)
    tmp <- distM[genes,genes]
    cors <- tmp[M[genes,genes]==2]

    if(plot){
        plot(density(cors),col=1)
    }else{
        lines(density(cors),col=2)
    }
    
    cors

}

netpair_overlap <- function(net1,net2,idmap=FALSE,plot=FALSE){
    
    if(idmap){
        mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
        ppinode <- mapT[match(net1$node,mapT[,1]),2]
        net1$node <- ppinode
        colnames(net1$matrix) <- ppinode
        rownames(net1$matrix) <- ppinode
    }
   
    nodes  <- intersect(net1$node,net2$node)
  
    M1 <- net1$matrix[nodes,nodes]
    M2 <- net2$matrix[nodes,nodes]

    M1[M1>0] <- 1
    M2[M2>0] <- 1
    M <- M1 + M2
    
    mm <- min(sum(M1)/2,sum(M2)/2)
    
    if(plot){
        a <- rowSums(M1)
        b <- rowSums(M2)
        plot(density(a),col=1)
        lines(density(b),col=2)    
    }
    
    c(sum(M==2)/2,length(nodes),mm)
}

netpair_merge <- function(net1,net2,flag=1,filename,idmap=FALSE){
    ## flag =1; overlap edges for two networks
    ## flag =2; joint edges for two networks
    
    if(idmap){
        mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
        ppinode <- mapT[match(net1$node,mapT[,1]),2]
        net1$node <- ppinode
        colnames(net1$matrix) <- ppinode
        rownames(net1$matrix) <- ppinode
    }
    
    if(flag==1){
        nodes  <- intersect(net1$node,net2$node)
        M1 <- net1$matrix[nodes,nodes]
        M2 <- net2$matrix[nodes,nodes]
        M1[M1>0] <- 1
        M2[M2>0] <- 1
        M <- M1 + M2
        subedges <- which(M==2,arr.ind=TRUE)
    }else if(flag==2){
        nodes  <- union(net1$node,net2$node)
        M <- matrix(0,length(nodes),length(nodes),dimnames=list(nodes,nodes))
        M[net1$node,net1$node] <- net1$matrix
        M[net2$node,net2$node] <- M[net2$node,net2$node] + net2$matrix
        subedges <- which(M>0,arr.ind=TRUE)
    }
    
    net <- cbind(nodes[subedges[,1]],nodes[subedges[,2]],1)
    print(dim(net))
    write.table(net,filename,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

}

FDR_select <- function(){
    
filenames <- c("result/DDD/CRFresultALLde.txt","result/DDD/CRFresultASD.txt","result/DDD/CRFresultDDD.txt","result/DDD/CRFresultEE.txt","result/DDD/CRFresultID.txt","result/PCGC/coexpCRFresult_PCGCLBP_121.txt","result/randresult/coexp/CRFresult_ASDallLBP_7.txt")

for(i in 1:length(filenames)){
    tmp <- read.table(filenames[i])
    if(i==1){
        plot(density(tmp[,2]),col=1)
    }else{
        lines(density(tmp[,2]),col=i)
    }
}

for(i in c(2,7)){
    tmp <- read.table(filenames[i])
    if(i==2){
        plot(density(tmp[,2]),col=1)
    }else{
        lines(density(tmp[,2]),col=i)
    }
}

filenames <- c("DDD_mutations/netinput/CRF_inputASD.txt","result/randset/CRF_inputASDall.txt")
for(i in 1:length(filenames)){
    tmp <- read.table(filenames[i])
    if(i==1){
        plot(density(tmp[,2]),col=1)
    }else{
        lines(density(tmp[,2]),col=i)
    }
}

}