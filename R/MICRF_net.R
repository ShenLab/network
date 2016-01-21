source("misc.R")

MICRF_net <- function(netflag,strn="LBP",filenames,fileexp,modulefile,cutn=15000,idm=FALSE,beta,net_e=4){
    
    allnet <- build_net(netflag,fileexp,net_e)
    load(modulefile)
    module <- module[module[,1] %in% allnet$node,]
    #module <- build_module(allnet$matrix)
    
    if(net_e <= 4){
        nodeinfo <- build_node(module,1,filenames)
        ## deal with NA
        nodeinfo[is.na(nodeinfo[,1]),1] <- 0
    }else{
        nodeinfo <- build_node_1(module,filenames)
        mode(nodeinfo) <- "numeric"
        if(any(is.na(nodeinfo))) print("NA error!!")
    }
    
    POSpro <- build_CRF(nodeinfo,allnet,module,netflag,strn,cutn,beta,net_e)
    if(idm){
        mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
        pros1 <- POSpro
        pros1[,1] <- mapT[match(pros1[,1],mapT[,1]),2]
        pros1 <- pros1[!is.na(pros1[,1]),]
        write.table(pros1,file=paste(strn,"LBP_",netflag,"1.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
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

build_node_1 <- function(module,filenames){
    
    genes <- module[,1]
    nodeinfo <- matrix(0,length(genes),2,dimnames=list(genes,1:2))
    nodeinfo[,1] <- 0.06
    nodeinfo[,2] <- 0.94
    tmp <- as.matrix(read.table(filenames[1]))
    tmp <- tmp[tmp[,1] %in% genes,]
    nodeinfo[match(tmp[,1],rownames(nodeinfo)),] <- tmp[,2:3]
    
    nodeinfo
}

build_net <- function(netflag,fileexp="PCGCall.txt",net_e=FALSE){
    # net: a list: size node matrix (systemtic and zero diag)
    net <- list()
    if(netflag==3 | netflag==4 | netflag>=6){
        if(netflag==3){ filename <- "data/STRINGnetmap.txt";}
        if(netflag==6){ filename <- "data/hotnet/iRefIndexm.txt";}
        if(netflag==7){ filename <- "data/network_inference/ASDexp_net_top5.txt";} ## ASD co-expression network
        if(netflag==8){ filename <- "data/network_inference/Infer_net_O4m.txt";} ## ASD inferred network
        if(netflag==20){ filename <- "data/StringNew_HPRD_mnet.txt";} ## MAGI
        if(netflag==21){ filename <- "data/network_inference/ASDexp_net.txt";} ## ASD co-expression network 0.7
        
        if(netflag==26){ filename <- "data/network_inference/brainspan_net_cor.txt";} ## brain corr 0.7
        if(netflag==27){ filename <- "data/network_inference/brainspan_net_top5.txt";} ## brain coexp 5
        if(netflag==28){ filename <- "data/network_inference/brainspan_net_top7.txt";} ## brain coexp 7
        if(netflag==29){ filename <- "data/network_inference/ComNet.txt";} ## combined network
        if(netflag==34){ filename <- "data/PrePPI.txt";}
        if(netflag==36){ filename <- "data/network_inference/ComCo_PrePPI.txt";}
        
        net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
        net.text <- rbind(net.text,net.text[,c(2,1,3)])
        net.text[,3] <- as.numeric(net.text[,3])/max(as.numeric(net.text[,3]))
        net <- read_net(net.text)
        diag(net$matrix) <- 0
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
            filename <- "data/Network_betweenness/Betweenness_edge_HPRD.txt"
        }else if(netflag==21){
            filename <- "data/Network_betweenness/Betweenness_edge_coexp1.txt"
        }else if(netflag==26){
            filename <- "data/Network_betweenness/Betweenness_edge_corr1.txt"
        }else if(netflag==27){
            filename <- "data/Network_betweenness/Betweenness_edge_coexp5.txt"
        }else if(netflag==28){
            filename <- "data/Network_betweenness/Betweenness_edge_coexp7.txt"
        }else if(netflag==29){
            filename <- "data/Network_betweenness/Betweenness_edge_ComNet.txt"
        }else if(netflag==34){
            filename <- "data/Network_betweenness/Betweenness_edge_PrePPI.txt"
        }else if(netflag==36){
            filename <- "data/Network_betweenness/Betweenness_edge_Co_PrePPI.txt"
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
            net <- discrete_net(allnet,modgenes,netflag,cutn,nodeinfo[,1],net_e)
        }else{
            net <- combine_net(allnet,modgenes,netflag)
        }
        if(sum(net$matrix)/2==0){break;}
        print(sum(net$matrix)/2)
        
        modnodeinfo <- nodeinfo[match(modgenes,rownames(nodeinfo)),]
        ### change 6_3
        if(net_e < 5){
            model <- build_model_0(modnodeinfo,net,net_e,beta)
        }else{
            model <- build_model_1(modnodeinfo,net,net_e,beta)
        }
        crfresult <- infer_crf(model, query.type=4)
        labels <- decode.lbp(model)
        POSpro[match(modgenes,POSpro[,1]),2:3] <- cbind(crfresult[,1],labels)
        
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

build_model_1 <- function(modnodesim,net,net_e=5,beta){   
    ## net_e discard parameter for old version
    ## beta: weights for two features
    
    n.nf <- dim(as.matrix(modnodesim))[2]
    S <- modnodesim
    colnames(S) <- 1:2
    subs <- !(S[,1]==0.06 & S[,2]==0.94)
    S[subs,1] <- 10^S[subs,1]
    S[subs,2] <- 1
    S[subs,] <- S[subs,]/matrix(apply(S[subs,],1,max),sum(subs),2,byrow=FALSE)
    
    n.states <- 2
    n.node <- net$size
    query.net <- net$matrix
    crf <- make.crf(net$matrix, rep(n.states,net$size))
    
    crf$state.map <- matrix(n.node, nrow=crf$n.nodes, ncol=crf$max.state)
    for (i in 1:crf$n.nodes)
    {
        crf$state.map[i, 1:crf$n.states[i]] <- 1:n.states
        crf$node.pot[i,] <- exp(beta[1] * S[i, crf$state.map[i,]])
    }   
    
    we <- net$we
    
    W1 <- matrix(1,n.states,n.states)
    for (e in 1:crf$n.edges){
        n1 <- crf$edges[e, 1]
        n2 <- crf$edges[e, 2]
        
        W1[1,1] = 0.5 * (S[n1,1] + S[n2,1]) * we[n1,n2]
        W1[1,2] = ( 1-0.5 * (S[n1,1] + S[n2,1]) ) * we[n1,n2]
        W1[2,1] = ( 1-0.5 * (S[n1,2] + S[n2,2]) ) * we[n1,n2]
        W1[2,2] = 0.5 * (S[n1,2] + S[n2,2]) * we[n1,n2]
        
        W1 <- beta[2] * W1
        crf$edge.pot[[e]] <- exp(W1)
    }
    
    crf
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
    
    diag(net$matrix) <- 0
    
    net
}

infer_crf <- function(model, query.type){
    query.type=4
    result <- infer.lbp(model)
    if(is.na(max(result$node.bel))){
        result <- infer.tree(model)
    }
    result <- result$node.bel
    result
}
