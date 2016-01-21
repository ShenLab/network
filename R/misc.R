read_net <- function(net.text)
{
    
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

qwt <- function (x, filer, sep = "\t", flag = 0) 
{
    if (flag == 0) {
        write.table(x, file = filer, quote = FALSE, row.names = FALSE, 
                    col.names = FALSE, sep = sep)
    }
    if (flag == 1) {
        write.table(x, file = filer, quote = FALSE, row.names = TRUE, 
                    col.names = FALSE, sep = sep)
    }
    if (flag == 2) {
        write.table(x, file = filer, quote = FALSE, row.names = FALSE, 
                    col.names = TRUE, sep = sep)
    }
    if (flag == 3) {
        write.table(x, file = filer, quote = FALSE, row.names = TRUE, 
                    col.names = TRUE, sep = sep)
    }
    if (flag == 4) {
        write.table(x, file = filer, quote = TRUE, row.names = TRUE, 
                    col.names = TRUE, sep = sep)
    }
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
