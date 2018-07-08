library(STATegRa)
library(grid)
library(ggplot2)
library(Biobase)
library(mixOmics)
data("STATegRa_S3")
B1 <- createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA,
                               pDataDescr=c("classname"))

# Block2 - miRNA expression data
B2 <- createOmicsExpressionSet(Data=Block2.PCA, pData=ed.PCA, 
                               pDataDescr=c("classname"))
cc <- selectCommonComps(X=Block1.PCA, Y=Block2.PCA, Rmax=3)
#grid.arrange(cc$pssq, cc$pratios, ncol=2)
cc$pratios
PCA.selection(Data=Block1.PCA, fac.sel="single%",
              varthreshold=0.03)$numComps
PCA.selection(Data=Block2.PCA, fac.sel="single%",
              varthreshold=0.03)$numComps
ms <- modelSelection(Input=list(B1,B2), Rmax=4, fac.sel="single%",
                     varthreshold=0.03)
ms
discoRes <- omicsCompAnalysis(Input=list(B1, B2), Names=c("expr", "mirna"),
                              method="DISCOSCA", Rcommon=2, Rspecific=c(2, 2),
                              center=TRUE, scale=TRUE, weight=TRUE)
jiveRes <- omicsCompAnalysis(Input=list(B1, B2), Names=c("expr", "mirna"),
                              method="JIVE", Rcommon=2, Rspecific=c(2, 2),
                              center=TRUE, scale=TRUE, weight=TRUE)
o2plsRes <- omicsCompAnalysis(Input=list(B1, B2),Names=c("expr", "mirna"),
                              method="O2PLS", Rcommon=2, Rspecific=c(2, 2),
                              center=TRUE, scale=TRUE, weight=TRUE)
slotNames(discoRes)
getScores(discoRes, part="common")
getScores(discoRes, part="distinctive", block="1")
getScores(discoRes, part="distinctive", block="2")
getScores(o2plsRes, part="common", block="expr")
getScores(o2plsRes, part="common", block="mirna")
getScores(o2plsRes, part="distinctive", block="1")
getScores(o2plsRes, part="distinctive", block="2")

getVAF(discoRes)
plotVAF(discoRes)
getVAF(jiveRes)
plotVAF(jiveRes)


# Scatterplot of scores variables associated to common components
# DISCO-SCA
plotRes(object=discoRes, comps=c(1, 2), what="scores", type="common",
        combined=FALSE, block="", color="classname", shape=NULL, labels=NULL,
        background=TRUE, palette=NULL, pointSize=2, labelSize=NULL,
        axisSize=NULL, titleSize=NULL) 
# DISCO-SCA scores combined plot for individual components
plotRes(object=discoRes, comps=c(1, 1), what="scores", type="individual",
        combined=TRUE, block="", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=2,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
# JIVE
plotRes(object=jiveRes, comps=c(1, 2), what="scores", type="common",
        combined=FALSE, block="", color="classname", shape=NULL, labels=NULL,
        background=TRUE, palette=NULL, pointSize=2, labelSize=NULL,
        axisSize=NULL, titleSize=NULL) 

# O2PLS 
# Associated to first block
plotRes(object=o2plsRes, comps=c(1, 2), what="scores", type="common",
        combined=FALSE, block="expr", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=2,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Associated to second block
plotRes(object=o2plsRes, comps=c(1, 2), what="scores", type="common",
        combined=FALSE, block="mirna", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=2,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Combined plot of scores variables assocaited to common components
plotRes(object=o2plsRes, comps=c(1, 1), what="scores", type="common",
        combined=TRUE, block="", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=2,
        labelSize=NULL, axisSize=NULL, titleSize=10)


gene<-t(Block1.PCA)
mirna<-t(Block2.PCA)
# Combined plot of scores variables assocaited to common components
plotRes(object=o2plsRes, comps=c(1, 1), what="scores", type="common",
        combined=TRUE, block="", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL) 



# DISCO-SCA scores scatterplot associated to individual components
# Associated to first block
plotRes(object=discoRes, comps=c(1, 2), what="scores", type="individual",
        combined=FALSE, block="expr", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Associated to second block
plotRes(object=discoRes, comps=c(1, 2), what="scores", type="individual",
        combined=FALSE, block="mirna", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)



# DISCO-SCA scores combined plot for individual components
plotRes(object=discoRes, comps=c(1, 1), what="scores", type="individual",
        combined=TRUE, block="", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
# DISCO-SCA combined plot of scores for common and individual components
plotRes(object=discoRes, comps=c(1, 1), what="scores", type="both",
        combined=TRUE, block="expr", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
plotRes(object=discoRes, comps=c(1, 1), what="scores", type="both",
        combined=TRUE, block="mirna", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)

# O2PLS combined plot of scores for common and individual components
plotRes(object=o2plsRes, comps=c(1, 1), what="scores", type="both",
        combined=TRUE, block="expr", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
plotRes(object=o2plsRes, comps=c(1, 1), what="scores", type="both",
        combined=TRUE, block="mirna", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Loadings plot for common components
#Separately for each block
plotRes(object=discoRes, comps=c(1, 2), what="loadings", type="common",
        combined=FALSE, block="expr", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
plotRes(object=discoRes, comps=c(1, 2), what="loadings", type="common",
        combined=FALSE, block="mirna", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)

# Combined plot
plotRes(object=discoRes, comps=c(1, 2), what="loadings", type="common",
        combined=TRUE, block="", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)



# Biplot common part. DISCO-SCA
biplotRes(object=discoRes, type="common", comps=c(1, 2), block="",
          title=NULL, colorCol="classname", sizeValues=c(2, 4),
          shapeValues=c(17, 0), background=TRUE, pointSize=4,
          labelSize=NULL, axisSize=NULL, titleSize=NULL) 
# Biplot common part. O2PLS 
biplotRes(object=o2plsRes, type="common", comps=c(1, 2),
          block="expr", title=NULL, colorCol="classname",
          sizeValues=c(2, 4), shapeValues=c(17, 0),
          background=TRUE, pointSize=4, labelSize=NULL,
          axisSize=NULL, titleSize=NULL)
biplotRes(object=o2plsRes, type="common", comps=c(1, 2),
          block="mirna", title=NULL, colorCol="classname",
          sizeValues=c(2, 4), shapeValues=c(17, 0),
          background=TRUE, pointSize=4, labelSize=NULL,
          axisSize=NULL, titleSize=NULL)
# Biplot distinctive part. O2PLS 
biplotRes(object=discoRes, type="individual", comps=c(1, 2),
          block="expr", title=NULL, colorCol="classname",
          sizeValues=c(2, 4), shapeValues=c(17, 0),
          background=TRUE, pointSize=4, labelSize=NULL,
          axisSize=NULL, titleSize=NULL)
biplotRes(object=discoRes, type="individual", comps=c(1, 2),
          block="mirna", title=NULL, colorCol="classname",
          sizeValues=c(2, 4), shapeValues=c(17, 0),
          background=TRUE, pointSize=4, labelSize=NULL,
          axisSize=NULL, titleSize=NULL)

gene<-t(Block1.PCA)
X<-gene
mirna<-t(Block2.PCA)
Y<-mirna
imgCor(gene, mirna)
class_cca <- rcc(gene, mirna)
cv <- rcc(gene, mirna, ncomp = 3, method = 'shrinkage')
cv$lambda

reg_cca<-rcc(gene,mirna,ncomp = 3,lambda1 = cv$lambda[1],lambda2 = cv$lambda[2])

##circle plot for rcca
plotVar(reg_cca, comp =1:2, cutoff = 0.6, 
        var.names = c(F, F), cex = c(0.5, 0.5),title = 'TCGA, rCCA comp 1-2')

cim(reg_cca, comp = 1:3, xlab = "gene", ylab = "mirna", 
    margins = c(5, 6))

plotIndiv(reg_cca, comp = 1:2, group  = ed.PCA$classname, 
          legend = TRUE, rep.space = "XY-variate",
          title = 'TCGA, rCCA, XY-space')

##3Drcca
plotIndiv(reg_cca, 
          group = ed.PCA$classname, 
          legend =TRUE, 
          ellipse = TRUE, 
          rep.space = 'XY-variate',
          title = '3D-TCGA,rCCA',
          star = TRUE, 
          centroid = TRUE,
          style='3d')


#pls
tumor.pls <- pls(gene, mirna, ncomp = 10, mode = "regression")
tumor.spls <- spls(gene, mirna, ncomp =10, keepX = c(10,10,10), keepY= c(10,10,10), mode = "regression")
plotIndiv(tumor.pls, comp = 1:2, rep.space= 'XY-variate', group = ed.PCA$classname,
          legend = TRUE, title = 'TCGA, PLS comp 1 - 2, XY-space')

plotIndiv(tumor.spls, comp = 1:2, rep.space= 'XY-variate', group = ed.PCA$classname,
          legend = TRUE, title = 'TCGA, sPLS comp 1 - 2, XY-space')
plotVar(tumor.pls, comp =1:2, 
        var.names = list(X.label = TRUE, 
                         Y.label = TRUE), cex = c(4, 5))
plotVar(tumor.spls, comp =1:2, 
        var.names = list(X.label = TRUE, 
                         Y.label = TRUE), cex = c(4, 5),title = 'TCGA, spls comp 1-2')
cim(tumor.pls, comp = 1:3, xlab = "gene", ylab = "mirna", 
    margins = c(5, 6))
cim(tumor.spls, comp = 1:3, xlab = "gene", ylab = "mirna", 
    margins = c(5, 6))

plotVar(tumor.spls, comp =1:2, cutoff = 0.1, 
        var.names = c(F, F), cex = c(1, 1),title = 'TCGA, sPLS comp 1-2')

#3dspls
plotIndiv(tumor.spls, 
          group = ed.PCA$classname, 
          legend =TRUE, 
          ellipse = TRUE, 
          rep.space = 'XY-variate',
          title = '3D-TCGA,sPLS',
          star = TRUE, 
          centroid = TRUE,
          style='3d')

