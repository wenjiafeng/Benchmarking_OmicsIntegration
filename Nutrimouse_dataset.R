library(STATegRa)
library(grid)
library(ggplot2)
library(Biobase)

library(mixOmics)

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
p <- nutrimouse$genotype

nutrimouse.shrink <- rcc(X, Y, ncomp = 3, method = 'shrinkage')
grid1 <- seq(0, 0.2, length = 5) 
grid2 <- seq(0.0001, 0.2, length = 5)

cv <- tune.rcc(X, Y, grid1 = grid1, grid2 = grid2, validation = "loo")
nutrimouse.rcc <- rcc(X,Y, ncomp = 3,  lambda1 = cv$opt.lambda1, 
                      lambda2 = cv$opt.lambda2)
##rcca
plotIndiv(nutrimouse.rcc, comp = 1:2, ind.names = nutrimouse$diet, 
          group = nutrimouse$genotype, rep.space = "XY-variate",
          legend = TRUE, ellipse = TRUE, title = 'Nutrimouse, rCCA XY-space')
##rcca circle
pdf("test.pdf", width = 3, height = 3)
plotVar(nutrimouse.rcc, comp = 1:2, cutoff = 0.5, var.names = c(TRUE, TRUE),
        cex = c(1, 1), title='Nutrimouse, rCCA comp 1 - 2')
dev.off()

##spls
nutrimouse.spls <- spls(X, Y, ncomp =10, keepX = c(10,10,10), keepY= c(10,10,10), mode = "regression")

plotIndiv(nutrimouse.spls, comp = 1:2, rep.space= 'XY-variate', group = nutrimouse$genotype,
          ind.names = nutrimouse$diet,
          legend = TRUE, ellipse = TRUE, title = 'Nutrimouse, sPLS comp 1-2, XY-space')

plotVar(nutrimouse.rcc, comp = 1:2, cutoff = 0.5, var.names = c(TRUE, TRUE),
        cex = c(1, 1), title='Nutrimouse, sPLS comp 1 - 2')

##3Drcca
plotIndiv(nutrimouse.rcc, 
          group = nutrimouse$diet, 
          ind.names = nutrimouse$diet,
          legend =TRUE, 
          ellipse = TRUE, 
          rep.space = 'XY-variate',
          title = '3D-Nutrimouse,rCCA',
          star = TRUE, 
          centroid = TRUE,
          style='3d')
##3Dspls
plotIndiv(nutrimouse.spls, 
          group = nutrimouse$diet, 
          ind.names = nutrimouse$diet,
          legend =TRUE, 
          ellipse = TRUE, 
          rep.space = 'XY-variate',
          title = '3D-Nutrimouse, sPLS',
          star = TRUE, 
          centroid = TRUE,
          style='3d')


##STATegRa

library(STATegRa)
library(grid)
library(ggplot2)
library(Biobase)


#staterga
liqid<-as.matrix(t(X))
gene<-as.matrix(t(Y))
# variable liqid expression data
X1 <- createOmicsExpressionSet(Data=liqid, pData=type,
                               pDataDescr=c("nutrimouse$genotype"))

# variable gene expression data
X2 <- createOmicsExpressionSet(Data=gene, pData=type, 
                               pDataDescr=c("nutrimouse$genotype"))
stat <- selectCommonComps(X=liqid, Y=gene, Rmax=3)
stat$common
stat$pssq
stat$pratios
PCA.selection(Data=liqid, fac.sel="single%",
              varthreshold=0.03)$numComps
#4
PCA.selection(Data=gene, fac.sel="single%",
              varthreshold=0.03)$numComps
#5
ms <- modelSelection(Input=list(X1,X2), Rmax=4, fac.sel="single%",
                     varthreshold=0.03)
ms
#discoRes
discoRes <- omicsCompAnalysis(Input=list(X1, X2), Names=c("liqid", "gene"),
                              method="DISCOSCA", Rcommon=3, Rspecific=c(3, 3),
                              center=TRUE, scale=TRUE, weight=TRUE)
slotNames(discoRes)
getScores(discoRes, part="common")
getScores(discoRes, part="distinctive", block="1")
getScores(discoRes, part="distinctive", block="2")
getVAF(discoRes)
plotVAF(discoRes)
plotRes(object=discoRes, comps=c(1, 2), what="scores", type="common",
        combined=FALSE, block="", shape=NULL, labels=NULL,
        background=TRUE, palette=NULL, pointSize=2, labelSize=NULL,
        axisSize=NULL, titleSize=NULL) 
#DISCO-SCA scores scatterplot associated to individual components
#liqid
plotRes(object=discoRes, comps=c(1, 2), what="scores", type="individual",
        combined=FALSE, block="liqid",  shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=3,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
#gene
plotRes(object=discoRes, comps=c(1, 2), what="scores", type="individual",
        combined=FALSE, block="gene", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
#combined plot for individual components
plotRes(object=discoRes, comps=c(1, 1), what="scores", type="individual",
        combined=TRUE, block="",  shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
#loading
# Combined plot
plotRes(object=discoRes, comps=c(1, 2), what="loadings", type="common",
        combined=TRUE, block="", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
#bipplot
biplotRes(object=discoRes, type="common", comps=c(1, 2), block="",
          title=NULL,sizeValues=c(2, 4),
          shapeValues=c(17, 0), background=TRUE, pointSize=4,
          labelSize=NULL, axisSize=NULL, titleSize=NULL) 



##o2pls
o2plsRes <- omicsCompAnalysis(Input=list(X1, X2),Names=c("liqid", "gene"),
                              method="O2PLS", Rcommon=3, Rspecific=c(3, 3),
                              center=TRUE, scale=TRUE, weight=TRUE)
# Scatterplot of scores variables associated to common components
# Associated to liqid
plotRes(object=o2plsRes, comps=c(1, 2), what="scores", type="common",
        combined=FALSE, block="liqid", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Associated to gene
plotRes(object=o2plsRes, comps=c(1, 2), what="scores", type="common",
        combined=FALSE, block="gene", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Combined plot of scores variables assocaited to common components
plotRes(object=o2plsRes, comps=c(1, 1), what="scores", type="common",
        combined=TRUE, block="", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=3,
        labelSize=NULL, axisSize=NULL, titleSize=NULL) 
#jive scores scatterplot associated to individual components
#liqid
plotRes(object=o2plsRes, comps=c(1, 2), what="scores", type="individual",
        combined=FALSE, block="liqid", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
#gene
plotRes(object=o2plsRes, comps=c(1, 2), what="scores", type="individual",
        combined=FALSE, block="gene",  shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
#combined plot for individual components
plotRes(object=o2plsRes, comps=c(1, 1), what="scores", type="individual",
        combined=TRUE, block="",shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
#bipplot for common part
biplotRes(object=o2plsRes, type="common", comps=c(1, 2),
          block="liqid", title=NULL,
          sizeValues=c(2, 4), shapeValues=c(17, 0),
          background=TRUE, pointSize=4, labelSize=NULL,
          axisSize=NULL, titleSize=NULL)
biplotRes(object=o2plsRes, type="common", comps=c(1, 2),
          block="gene", title=NULL,
          sizeValues=c(2, 4), shapeValues=c(17, 0),
          background=TRUE, pointSize=4, labelSize=NULL,
          axisSize=NULL, titleSize=NULL)


