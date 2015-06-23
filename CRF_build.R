CRF_build <- function(nodesim,module,datExpr,strn,nodeflag){
    options(stringsAsFactors=FALSE)
    library(CRF)
    library(Corbi)
    library(bnlearn)
    library(minet)
    # deal with 0 express
    modLab <- setdiff(unique(module[,2]),0)
    for(i in 1:length(modLab)){
        modgenes <- module[module[,2]==modLab[i],1]
        if(2*sum(datExpr[,modgenes]==0) > dim(datExpr[,modgenes])[1] * dim(datExpr[,modgenes])[2]){
            module[module[,1]%in%modgenes,2] <- 0
        }
    }
    
    
    modLab <- setdiff(unique(module[,2]),0)
    n.module <- length(modLab)
    n.genes <- length(module[module[,2]!=0,1])

    rankM <- matrix(0,7,10,dimnames=list(c("RBFOX2","KMT2D","EP300","NAA16","KDM5B","NCOA3","SMAD2"),1:10))
    for(flag in 8){
    POSpro <- matrix(0,n.genes,3)
    POSpro[,1] <- as.vector(module[module[,2]!=0,1]) # gene, probability; FDR, node score, edge score
    for (i in 1:n.module){
        #print(i)
        modgenes <- as.vector(module[module[,2]==modLab[i],1])
        #print(length(modgenes))
        if(flag==1){  net <- build_net_exp(datExpr,modgenes,times=4);}
        if(flag==2){  net <- build_net_aracne(datExpr,modgenes);}
        if(flag==3){
            neta <- build_net_exp(datExpr,modgenes,times=4);
            net <- build_net_aracne(datExpr,modgenes);
            net$matrix[modgenes,modgenes] <- neta$matrix[modgenes,modgenes] + net$matrix[modgenes,modgenes]
            net$matrix[net$matrix > 0] <- 1
        }
        if(flag==4){
            neta <- build_net_exp(datExpr,modgenes,times=4);
            net <- build_net_PPI(datExpr,modgenes,1);
            net$matrix[modgenes,modgenes] <- neta$matrix[modgenes,modgenes] + net$matrix[modgenes,modgenes]
            net$matrix[net$matrix > 0] <- 1
        }
        
        if(flag==5){
            neta <- build_net_aracne(datExpr,modgenes);
            net <- build_net_PPI(datExpr,modgenes,1);
            net$matrix[modgenes,modgenes] <- neta$matrix[modgenes,modgenes] + net$matrix[modgenes,modgenes]
            net$matrix[net$matrix > 0] <- 1
            # deal with memory limit
            if((sum(net$matrix)/2) > 20000){
                neta <- build_net_aracne(datExpr,modgenes);
                net <- build_net_PPI(datExpr,modgenes,1);
                tmpnodes <- net$node[which(rowSums(net$matrix)==0)]
                net$matrix[tmpnodes,tmpnodes] <- neta$matrix[tmpnodes,tmpnodes] + net$matrix[tmpnodes,tmpnodes]
                net$matrix[net$matrix > 0] <- 1
            }
        }
        
        if(flag==6){
            neta <- build_net_aracne(datExpr,modgenes)
            netb <- build_net_exp(datExpr,modgenes,times=4);
            net <- build_net_PPI(datExpr,modgenes,1)
            net$matrix[modgenes,modgenes] <- neta$matrix[modgenes,modgenes] + netb$matrix[modgenes,modgenes] + net$matrix[modgenes,modgenes]
            net$matrix[net$matrix > 0] <- 1  
        }
        
        if(flag==7){
            neta <- build_net_aracne(datExpr,modgenes);
            net <- build_net_phe(datExpr,modgenes,0.3);
            net$matrix[modgenes,modgenes] <- neta$matrix[modgenes,modgenes] + net$matrix[modgenes,modgenes]
            net$matrix[net$matrix > 0] <- 1
            if((sum(net$matrix)/2) > 20000){     # deal with memory limit
                neta <- build_net_aracne(datExpr,modgenes);
                net <- build_net_phe(datExpr,modgenes,0.2);
                tmpnodes <- net$node[which(rowSums(net$matrix)==0)]
                net$matrix[tmpnodes,tmpnodes] <- neta$matrix[tmpnodes,tmpnodes] + net$matrix[tmpnodes,tmpnodes]
                net$matrix[net$matrix > 0] <- 1
            }        
        }
        
        if(flag==8){
            net <- build_net_phe(datExpr,modgenes,0.3);
            if((sum(net$matrix)/2) > 20000){     # deal with memory limit
                net <- build_net_phe(datExpr,modgenes,0.5);
                if((sum(net$matrix)/2) > 20000){
                    net <- build_net_PPI(datExpr,modgenes,1);
                }
            } 
        }
        
        print(sum(net$matrix)/2)
        
        n.modgenes <- length(modgenes)
        modnodesim <- nodesim[match(modgenes,rownames(nodesim)),]
        #modnodesim[is.na(modnodesim[,1]),] <- 0.5 #???
        #rownames(modnodesim) <- modgenes
        
        model <- build_model(modnodesim,net,bflag=1)
        #crfresult <- solve_crf(model, query.type=4)
        crfresult <- infer_crf(model, query.type=4)
        POSpro[match(modgenes,POSpro[,1]),2:3] <- crfresult
        
        #POSpro[match(modgenes,POSpro[,1]),3] <- ifelse(crfresult[,2] * net$size >=1, 1, crfresult[,2] * net$size)
        #POSpro[match(modgenes,POSpro[,1]),2] <- 1 - as.numeric(POSpro[match(modgenes,POSpro[,1]),3])
    }

    POSpro <- POSpro[order(as.numeric(POSpro[,2]),decreasing = TRUE),]
    write.table(POSpro,file=paste("LBP_",strn,"_",flag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    pros1 <- mapT[match(POSpro[,1],mapT[,1]),2]
    pros1 <- cbind(pros1,POSpro[,2])
    write.table(pros1,file=paste("LBP",strn,flag,"1.txt",sep="_"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

    rankM[,flag] <- match(rownames(rankM),pros1[,1])
    testname <- paste("enrich",strn,flag,"1.txt",sep="_")
    Output_file1(paste("LBP",strn,flag,"1.txt",sep="_"),mapT,testname,flag)
    }
    write.table(rankM,file=paste(strn,"result.txt",sep=""),quote=F,col.names=F,sep="\t")
    #Result
}

build_net_exp <- function(datExpr,modgenes,times=4){
    library(WGCNA)
    modExpr <- datExpr[,modgenes]
    distgene = adjacency(modExpr, type = "unsigned", power=1) # co-expression network
    adj <- distgene
    diag(adj) <- 0
    cutf <- quantile(adj[upper.tri(adj)],0.95)
    adj[adj <= cutf] <- 0
    adj[adj > cutf] <- 1
    
    flag=FALSE
    if(dim(adj)[1]>3000){
        cute <- ifelse(dim(adj)[1] * times > 6000,dim(adj)[1] * times/2,dim(adj)[1] * times)
        if(sum(adj) > 2*cute) flag=TRUE;
    }else if(dim(adj)[1] < 500){
        cute <- times*dim(adj)[1]
        if(sum(adj) < 2*cute) flag=TRUE;
    }else if(dim(adj)[1] >= 500 & dim(adj)[1] <= 3000){
        cute <- times*dim(adj)[1]
        if(sum(adj) > 2*cute) flag=TRUE;    
    }
    
    if(flag==TRUE){
        adj <- distgene
        diag(adj) <- 0
        cutf <- sort(adj[upper.tri(adj)],decreasing=TRUE)[cute]
        adj[adj <= cutf] <- 0
        adj[adj > cutf] <- 1    
    } 
    
    net <- list()
    net$matrix <- adj
    net$size <- dim(adj)[1]
    net$node <- colnames(adj)
    net
}

build_net_aracne <- function(datExpr,modgenes){

    modExpr <- datExpr[,modgenes]
    modExpr <- as.data.frame(modExpr)
    neta <- bnlearn::aracne(modExpr)
    net <- list()
    net$size <- length(neta$nodes)
    net$node <- names(neta$nodes)
    net$matrix <- matrix(0,net$size,net$size,dimnames=list(net$node,net$node))
    
    for(i in 1:net$size){
        net$matrix[net$node[i],neta$nodes[[i]]$nbr] <- 1
        net$matrix[neta$nodes[[i]]$nbr, net$node[i]] <- 1
    }
    diag(net$matrix) <- 0
    
    net
    
}

build_net_aracne1 <- function(datExpr,modgenes){

    modExpr <- datExpr[,modgenes]
    mim <-  build.mim(modExpr,estimator="spearman")
    neta <- minet::aracne(mim)

    net <- list()
    net$size <- dim(neta)[1]
    net$node <- colnames(neta)
    net$matrix <- neta
    net$matrix[net$matrix>0] <- 1
    
    diag(net$matrix) <- 0
    
    net
}

build_net_PPI <- function(datExpr,modgenes,flag=1){
    
    if(flag==1){filename <- "STRINGnetmap.txt";}else{filename <- "FUNCTIONnetmap.txt";}
    net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
    net.text <- cbind(net.text,1)
    net.text <- rbind(net.text,net.text[,c(2,1,3)])
    allnet <- read_net(net.text)
    
    net <- list()
    net$node <- modgenes
    net$size <- length(modgenes)
    net$matrix <- matrix(0,net$size,net$size,dimnames=list(modgenes,modgenes))
    tmp <- intersect(modgenes,allnet$node)
    net$matrix[tmp,tmp] <- allnet$matrix[tmp,tmp]
    diag(net$matrix) <- 0
    net
    
    
    # GRN known
    # modgenes <- intersect(modgenes,netnodes)
    # net.edges <- cbind(GRNTable[GRNTable[,1] %in% modgenes & GRNTable[,2] %in% modgenes,1:2],1)
    # net <- read_net(net.edges)
    # modgenes <- net$node

}

build_net_phe <- function(datExpr,modgenes,cut=0.2){
    
    filename <- "phenotype_net.txt"
    net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
    net.text <- net.text[as.numeric(net.text[,3]) > cut,]
    net.text[,3] <- 1
    net.text <- rbind(net.text,net.text[,c(2,1,3)])
    allnet <- read_net(net.text)
    
    net <- list()
    net$node <- modgenes
    net$size <- length(modgenes)
    net$matrix <- matrix(0,net$size,net$size,dimnames=list(modgenes,modgenes))
    tmp <- intersect(modgenes,allnet$node)
    net$matrix[tmp,tmp] <- allnet$matrix[tmp,tmp]
    diag(net$matrix) <- 0
    net
}

build_model <- function(modnodesim,net,bflag=1,beta=0){   

    n.nf <- dim(as.matrix(modnodesim))[2]
    S <- modnodesim
    # S <- node_information(modnodesim,n.nf,bflag)
    # bflag==1 # different node features as different features
    # bflag==2 # integrate different node features first
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
    
    deg <- colSums(net$matrix)
    W <- matrix(1,n.states,n.states)
    W1 <- W
    W2 <- W
    for (e in 1:crf$n.edges){
        n1 <- crf$edges[e, 1]
        n2 <- crf$edges[e, 2]
        m1 <- 1:crf$n.states[n1]
        m2 <- 1:crf$n.states[n2]
        S1 <- matrix(log(crf$node.pot[n1, m1]), crf$n.states[n1], crf$n.states[n2])
        S2 <- matrix(log(crf$node.pot[n2, m2]), crf$n.states[n1], crf$n.states[n2], byrow=T)
        
#         W[1,1] <- exp(-abs(S1[1,1]-S2[1,1])^2)
#         W[2,2] <- exp(-abs(S1[2,2]-S2[2,2])^2)
#         W[1,2] <- min(exp(-abs(S1[1,1]-S2[1,2])^2),W[1,1])
#         W[2,1] <- min(exp(-abs(S1[2,1]-S2[1,1])^2),W[2,2])
#         crf$edge.pot[[e]] <- W*exp(S1+S2)/2
        
        ####W1 <- beta*S1*S2/sqrt(deg[n1]*deg[n2]) ## version 4 del
        #W1[1,1] <- (max(S1[1,1],S2[1,1])+beta)/sqrt(deg[n1]*deg[n2]) ## version 1
        W1[1,1] <- (max(S1[1,1],S2[1,1])*beta)/sqrt(deg[n1]*deg[n2]) ## version 2
        W1[2,2] <- 0
        W1[1,2] <- 0
        W1[2,1] <- 0
        crf$edge.pot[[e]] <- exp(W1)
    }
    crf
}

build_model1 <- function(modnodesim,net,bflag=1,beta=0){   
    
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
    
    we <- net$we
    
    deg <- colSums(net$matrix)
    W <- matrix(1,n.states,n.states)
    W1 <- W
    W2 <- W
    for (e in 1:crf$n.edges){
        n1 <- crf$edges[e, 1]
        n2 <- crf$edges[e, 2]
        m1 <- 1:crf$n.states[n1]
        m2 <- 1:crf$n.states[n2]
        S1 <- matrix(log(crf$node.pot[n1, m1]), crf$n.states[n1], crf$n.states[n2])
        S2 <- matrix(log(crf$node.pot[n2, m2]), crf$n.states[n1], crf$n.states[n2], byrow=T)
                
        #####W1 <- beta*S1*S2*we[n1,n2] # a little better in rand set 3 ## version 1 del
        ##W1[1,1] <- (beta+max(S1[1,1],S2[1,1]))*we[n1,n2] ## version 3
        W1[1,1] <- beta*max(S1[1,1],S2[1,1])*we[n1,n2] ## version 4
        W1[2,2] <- 0
        W1[1,2] <- 0
        W1[2,1] <- 0
        crf$edge.pot[[e]] <- exp(W1)
    }
    crf
}

decode_heuristic <- function(crf){
    result <- try(decode.junction(crf), T)
    if (class(result) == "try-error")
    {
        result <- decode.lbp(crf)
    }
    result
}

solve_crf <- function(model, query.type){
    query.type = 4
    result <- decode.lbp(model)
    #result <- decode.tree(model)
    result <- model$state.map[cbind(1:model$n.nodes, result)]
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
    
read_net <- function(net.text){
    
    net.node <- unique(union(net.text[,1],net.text[,2]))
    net.node <- net.node[net.node != ""]
    net.size <- length(net.node)
    net.edge <- cbind(as.character(net.text[,1]), as.character(net.text[,2]))
    net.edge <- net.edge[net.edge[,2] != "", ]
    net.edge <- net.edge[net.edge[,1] != "", ]
    net.matrix <- matrix(0, net.size, net.size, dimnames=list(net.node, net.node))
    net.matrix[net.edge] <- as.numeric(net.text[,3])
    list(size=net.size, node=net.node, matrix=net.matrix)
    
}

Output_file1 <- function(filename,mapT,testname,fileflag,flag=0){
    
    HMGFile <- "HMGs.csv"
    TGFFile <- "TGF-beta.txt"
    HMGgene <- read.csv(HMGFile,header=FALSE)[,2] # there is an unkonwn in 152 HMGgene
    TGFgene <- read.delim(TGFFile,header=FALSE,sep="\t")[,2]
    HMGgene <- mapping_to(HMGgene)
    TGFgene <- mapping_to(TGFgene)
    
    if(grepl(".csv",filename)){
        result <- read.csv(filename)
    }else{
        result <- as.matrix(read.table(filename,sep="\t"))
    }
    ntmp <- dim(result)[2]
    
    result[is.na(result[,ntmp]),ntmp] <- 1
    if(flag<=1){
        rASDgene <- result[as.numeric(result[,ntmp]) < flag,1] # cutoff
    }else if(flag >1 ){
        rASDgene <- result[1:flag,1] # number 
    }
    print(length(rASDgene))
    
    con <- file(testname,"w")
    tmp <- "RBFOX2" %in% toupper(rASDgene)
    writeLines("RBFOX2 test:",con)
    writeLines(as.character(tmp+0),con)
    
    writeLines("HMG genes test p value is:",con)
    tmp <- hyper_test(HMGgene,rASDgene,result[,1])
    writeLines(as.character(tmp),con)
    writeLines(intersect(HMGgene,rASDgene),con)
    writeLines("TGF-beta gene test p value is:",con)
    tmp <- hyper_test(TGFgene,rASDgene,result[,1])
    writeLines(intersect(TGFgene,rASDgene),con)
    writeLines(as.character(tmp),con)
    
    knownFile <- "known_genes.csv"
    knowngenes <- as.matrix(read.csv(knownFile,header=FALSE))[,1]
    knowngenes <- mapping_to(knowngenes)
    writeLines("Total known genes and in rASD genes: ",con)  
    writeLines(knowngenes[which(toupper(knowngenes) %in% toupper(rASDgene))],con)
    close(con)     
    
}

hyper_test <- function(geneset,rgeneset,genes){
    m <- length(intersect(genes,geneset))
    k <- length(rgeneset)
    n <- length(genes) - m
    x <- length(intersect(rgeneset,geneset))
    P <- phyper(x-1, m, n, k, lower.tail = FALSE, log.p = FALSE) 
    P
}

to_result <- function(){

filename <- "Bayes_Posterior_5_1.txt"
TADAFile <- "TADA_lofmis1202.csv"
TADAresult <- read.csv(TADAFile)
TADAresult[,1] <- mapping_to(TADAresult[,1])
result <- as.matrix(read.table(filename))
result[,2] <- 1:dim(result)[1]

subs <- match(result[,1],TADAresult[,1])
result <- cbind(result,TADAresult[subs,2:dim(TADAresult)[2]])
colnames(result) <- c("Gene","Rank",colnames(TADAresult)[2:dim(TADAresult)[2]])
write.csv(result,file=gsub(".txt",".csv",filename),row.names=F,quote=F)

}

node_information <- function(nodesim,n.nf,bflag=1){
    
    if(n.nf ==1 ){
        n.pot <- nodesim
    }else if(n.nf > 1){
        if(bflag==1){
            x <- sapply(1:dim(nodesim)[1],function(i) mean(nodesim[i, nodesim[i,] > 0]))        
            x[is.na(x)] <- nodesim[is.na(x),1]
            n.pot <- plogis(x, location = 0.5, scale = 1, lower.tail = TRUE, log.p = FALSE)
        }else if(bflag==2){
            W <- matrix(1,n.nf,1); W[1,1] <- 2; #####weights for different features
            x <- nodesim %*% W
            n.pot <- x/max(x)
        }
        names(n.pot) <- rownames(nodesim)
    }
    
    n.pot
}