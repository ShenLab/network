comNet <- function(){
    source("misc.R")
    source("comNet.R")
    
    ## combined the two best performance networks: iRefIndex and co-expression networks
    filename <- "data/hotnet/iRefIndexm.txt";
    netiRef <- net_matrix(filename)
    netiRef$matrix[netiRef$matrix > 0] <- 1
    
    filename <- "data/network_inference/brainspan_net_top5.txt"; ## brain coexp 5
    netCo <- net_matrix(filename)
    netCo$matrix[netCo$matrix > 0] <- 1
    
    ## read the edge betweenness 
    filename <- "data/Network_betweenness/Betweenness_edge_iRef.txt"
    netiRef <- net_between(filename,netiRef)
    filename <- "data/Network_betweenness/Betweenness_edge_coexp5.txt"
    netCo <- net_between(filename,netCo)
    
    genes <- union(netiRef$node,netCo$node)
    ig <- intersect(netiRef$node,netCo$node)
    tmp <- netiRef$matrix[ig,ig] + netCo$matrix[ig,ig]
    tmp[tmp>0] <- 1
    uniRef <- setdiff(netiRef$node,ig)
    uniCo  <- setdiff(netCo$node,ig)
    
    net <- list()
    net$size <- length(genes)
    net$node <- genes
    net$matrix <- matrix(0,net$size,net$size,dimnames=list(genes,genes))
    net$matrix[uniRef,uniRef] <- netiRef$matrix[uniRef,uniRef]
    net$matrix[uniRef,ig] <- netiRef$matrix[uniRef,ig]
    net$matrix[ig,uniRef] <- netiRef$matrix[ig,uniRef]
    net$matrix[uniCo,uniCo] <- netCo$matrix[uniCo,uniCo]
    net$matrix[uniCo,ig] <- netCo$matrix[uniCo,ig]
    net$matrix[ig,uniCo] <- netCo$matrix[ig,uniCo]
    net$matrix[ig,ig] <- tmp
    
    rm(tmp)
    
    net$we <- matrix(0,net$size,net$size,dimnames=list(genes,genes))
    net$we[uniRef,uniRef] <- netiRef$we[uniRef,uniRef]
    net$we[uniCo,uniCo] <- netCo$we[uniCo,uniCo]
    net$we[uniRef,uniRef] <- netiRef$we[uniRef,uniRef]
    net$we[uniRef,ig] <- netiRef$we[uniRef,ig]
    net$we[ig,uniRef] <- netiRef$we[ig,uniRef]
    net$we[uniCo,uniCo] <- netCo$we[uniCo,uniCo]
    net$we[uniCo,ig] <- netCo$we[uniCo,ig]
    net$we[ig,uniCo] <- netCo$we[ig,uniCo]
    net$we[ig,ig] <- pmax(netiRef$we[ig,ig],netCo$we[ig,ig])## this the method to combine network betweenness
    ## maybe recalculte the network betweenness for combined networks
    
    #net
    
    ## write the net file
    net$matrix[upper.tri(net$matrix,diag=TRUE)] <- 0
    edges <- which(net$matrix>0,arr.ind=TRUE)
    net.text <- cbind(genes[edges[,1]],genes[edges[,2]],1)
    qwt(net.text,file="data/network_inference/ComNet.txt")
    
    ## write the net betweenness
    net.we <- cbind(genes[edges[,1]],"--",genes[edges[,2]],net$we[edges])
    qwt(net.we,file="data/Network_betweenness/Betweenness_edge_ComNet.txt")
    
}

comNet_glasso <- function(){
    source("misc.R")
    source("comNet.R")
    
    ## combined the two best performance networks: iRefIndex and co-expression networks
    filename <- "data/hotnet/iRefIndexm.txt";
    netiRef <- net_matrix(filename)
    netiRef$matrix[netiRef$matrix > 0] <- 1
    
    filename <- "data/network_inference/brainspan_net_top5.txt"; ## brain coexp 5
    netCo <- net_matrix(filename)
    netCo$matrix[netCo$matrix > 0] <- 1
  
    filename <- "data/network_inference/glassoBrainNet_ab.txt"; ##
    netLa <- net_matrix(filename)
    netLa$matrix <- netLa$matrix + t(netLa$matrix)
    netLa$matrix[netLa$matrix > 0] <- 1
   
    
    genes <- union(netiRef$node,union(netLa$node,netCo$node))
    net <- list()
    net$size <- length(genes)
    net$node <- genes
    net$matrix <- matrix(0,net$size,net$size,dimnames=list(genes,genes))
    net$matrix[netiRef$node,netiRef$node] <- net$matrix[netiRef$node,netiRef$node] + netiRef$matrix[netiRef$node,netiRef$node]
    net$matrix[netCo$node,netCo$node] <- net$matrix[netCo$node,netCo$node] + netCo$matrix[netCo$node,netCo$node]
    net$matrix[netLa$node,netLa$node] <- net$matrix[netLa$node,netLa$node] + netLa$matrix[netLa$node,netLa$node]
    net$matrix[net$matrix > 0] <- 1
    
    ## write the net file
    net$matrix[upper.tri(net$matrix,diag=TRUE)] <- 0
    edges <- which(net$matrix>0,arr.ind=TRUE)
    net.text <- cbind(genes[edges[,1]],genes[edges[,2]],1)
    qwt(net.text,file="data/network_inference/ComNet3.txt")
    
}

comNet2iRef <- function(){
    source("misc.R")
    source("comNet.R")
    
    ## combined the two best performance networks: iRefIndex and co-expression networks
    filename <- "data/hotnet/iRefIndexm.txt";
    netiRef <- net_matrix(filename)
    diag(netiRef$matrix) <- 0
    netiRef$matrix[netiRef$matrix > 0] <- 1
       
    filename <- "data/network_inference/glassoBrainNet_ab.txt"; ##
    netLa <- net_matrix(filename)
    diag(netLa$matrix) <- 0
    netLa$matrix <- netLa$matrix + t(netLa$matrix)  
    netLa$matrix[netLa$matrix > 0] <- 1
  
    genes <- union(netiRef$node,netLa$node)
    net <- list()
    net$size <- length(genes)
    net$node <- genes
    net$matrix <- matrix(0,net$size,net$size,dimnames=list(genes,genes))
    net$matrix[netiRef$node,netiRef$node] <- net$matrix[netiRef$node,netiRef$node] + netiRef$matrix[netiRef$node,netiRef$node]
    net$matrix[netLa$node,netLa$node] <- net$matrix[netLa$node,netLa$node] + netLa$matrix[netLa$node,netLa$node]
    net$matrix[net$matrix > 0] <- 1
    
    ## write the net file
    net$matrix[upper.tri(net$matrix,diag=TRUE)] <- 0
    edges <- which(net$matrix>0,arr.ind=TRUE)
    net.text <- cbind(genes[edges[,1]],genes[edges[,2]],1)
    qwt(net.text,file="data/network_inference/ComNet2iRef.txt")
    
}

comNet2Co <- function(){
    source("misc.R")
    source("comNet.R")
    
    ## combined the two best performance networks: iRefIndex and co-expression networks
    filename <- "data/network_inference/brainspan_net_top5.txt"; ## brain coexp 5
    netCo <- net_matrix(filename)
    diag(netCo$matrix) <- 0
    netCo$matrix[netCo$matrix > 0] <- 1
    
    filename <- "data/network_inference/glassoBrainNet_ab.txt"; ##
    netLa <- net_matrix(filename)
    diag(netLa$matrix) <- 0
    netLa$matrix <- netLa$matrix + t(netLa$matrix)  
    netLa$matrix[netLa$matrix > 0] <- 1
    
    genes <- union(netCo$node,netLa$node)
    net <- list()
    net$size <- length(genes)
    net$node <- genes
    net$matrix <- matrix(0,net$size,net$size,dimnames=list(genes,genes))
    net$matrix[netCo$node,netCo$node] <- net$matrix[netCo$node,netCo$node] + netCo$matrix[netCo$node,netCo$node]
    net$matrix[netLa$node,netLa$node] <- net$matrix[netLa$node,netLa$node] + netLa$matrix[netLa$node,netLa$node]
    net$matrix[net$matrix > 0] <- 1
    
    ## write the net file
    net$matrix[upper.tri(net$matrix,diag=TRUE)] <- 0
    edges <- which(net$matrix>0,arr.ind=TRUE)
    net.text <- cbind(genes[edges[,1]],genes[edges[,2]],1)
    qwt(net.text,file="data/network_inference/ComNet2Co.txt")
    
}

net_matrix <- function(filename){
    net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
    net.text <- rbind(net.text,net.text[,c(2,1,3)])
    net.text[,3] <- as.numeric(net.text[,3])
    net<- read_net(net.text)
    net
}

net_between <- function(filename,net){
    edgebe <- read.table(filename)
    pedge <- ecdf(edgebe[,4])
    
    edgebe <- edgebe[edgebe[,1] %in% net$node & edgebe[,3] %in% net$node,]
    we <- matrix(0,net$size,net$size, dimnames= list(net$node,net$node))
    subs1 <- match(edgebe[,1],net$node)
    subs2 <- match(edgebe[,3],net$node)
    we[cbind(subs1,subs2)] <- pedge(edgebe[,4])
    we[cbind(subs2,subs1)] <- pedge(edgebe[,4])
    net$we <- we
    
    net
}

PrePPI_net <- function(){
    LRcut <- 600
    idmaps <- read.delim("data/UniprotID_HGNCsy.txt")
    initnet <- read.table("data/PrePPI_HC.txt",header=TRUE,sep="\t")
    
    subs <- initnet[,1] %in% idmaps[,2] & initnet[,2] %in% idmaps[,2]
    initnet <- initnet[subs,]
    initnet[,1] <- idmaps[match(initnet[,1],idmaps[,2]),1]
    initnet[,2] <- idmaps[match(initnet[,2],idmaps[,2]),1]
    initnet[,3] <- initnet[,5]/(initnet[,5] + LRcut)
    initnet <- initnet[,1:3]
    source("~/.Rprofile")
    qwt(initnet,file="data/PrePPI.txt")
}
