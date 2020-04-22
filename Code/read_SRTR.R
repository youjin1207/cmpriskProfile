library(sas7bdat)
library(gplots)
library(ColorPalette)
library(RColorBrewer)
BlPalette = colorRampPalette(brewer.pal(9,"PuBu"))(20)
GnPalette = colorRampPalette(brewer.pal(9,"Greens"))(20)
#kidney = read.csv("Data/kidney_20JAN2020.csv", header = TRUE, sep = ",") # read the orignal data 
dat = data.frame(timeto = kidney$X, cause1 = as.integer(kidney$dead), 
                 cause2 = as.integer(kidney$KT), 
                 censor = as.integer(kidney$dead == 0 & kidney$KT == 0),
                 center = as.factor(kidney$CAN_LISTING_CTR_ID), 
                 age = kidney$age50_wl,
                 bmi = kidney$BMI_wl,
                 female = as.integer(kidney$female),
                 race_black = as.integer(kidney$race_Black), 
                 race_hispanic = as.integer(kidney$race_Hispanic), 
                 race_asian = as.integer(kidney$race_Asian),
                 race_other = as.integer(kidney$race_Other),
                 blood_a = as.integer(kidney$blood_a), 
                 blood_ab = as.integer(kidney$blood_ab),
                 blood_b = as.integer(kidney$blood_b),
                 diag_poly = as.integer(kidney$diag_poly),
                 diag_diab = as.integer(kidney$diag_diab),
                 diag_hyper = kidney$diag_hyper,
                 diag_other = kidney$diag_other,
                 copd = kidney$copd,
                 hyperten = kidney$hyperten,
                 malig = kidney$malig,
                 pvd = kidney$pvd,
                 diabetes = kidney$diabetes)

# only consider the centers with >= 25 center size
names(which(table(dat$center) >= 25))
ind = as.character(dat$center) %in% names(which(table(dat$center) >= 25))
partdat = dat[ind,]

# exploratory analysis: correlation between mortality and transplant rates and two covariates (diabetes and age)
prop.death = aggregate(cause1 ~ center, FUN = "mean",data = dat)$cause1
prop.transplant = aggregate(cause2 ~ center, FUN = "mean",data = dat)$cause2
prop.diabetes = aggregate(diabetes ~ center, FUN = "mean",data = dat)$diabetes
prop.female = aggregate(female ~ center, FUN = "mean",data = dat)$female

cor(prop.death, prop.diabetes)
cor(prop.death, prop.female)

cor(prop.transplant, prop.diabetes)
cor(prop.transplant, prop.female)

## calculate population level risk class probability ##
# (1) prognostic score based risk class
dat = partdat
fit1 = coxph(Surv(timeto, cause1) ~ age + bmi + female +     
               race_black + race_hispanic + race_asian +  
               race_other + blood_a + blood_ab +     
               blood_b + diag_poly + diag_diab +     
               diag_hyper + diag_other + copd +        
               hyperten + malig + pvd +       
               diabetes + strata(center), 
             data = dat)
fit2 = coxph(Surv(timeto, cause2) ~ age + bmi + female +     
               race_black + race_hispanic + race_asian +  
               race_other + blood_a + blood_ab +     
               blood_b + diag_poly + diag_diab +     
               diag_hyper + diag_other + copd +        
               hyperten + malig + pvd +       
               diabetes + strata(center), 
             data = dat)

first.scores.whole = predict(fit1, reference = "sample")
second.scores.whole = predict(fit2, reference = "sample")
dat$first.scores.whole = first.scores.whole
dat$second.scores.whole = second.scores.whole

dist.mat = matrix(0, nrow(part1), 2)
score1 = dat$first.scores.whole
score2 = dat$second.scores.whole
dist.mat[,1] = apply(as.matrix(score1), 1, function(x) sum(x > first.scores.raw))
dist.mat[,2] = apply(as.matrix(score2), 1, function(x) sum(x > second.scores.raw))
table.raw = table(dist.mat[,1], dist.mat[,2])
table.raw = as.matrix(table.raw) / nrow(part1)

pdf("Figure/partdat_template_scores.pdf")
par(cex.main=1.8, cex.lab = 4, cex.axis = 4, 
    mar=c(3,2,3,0), tcl = 0.5, oma = c(3, 3, 3, 2), xpd = TRUE)
heatmap.2(as.matrix(table.raw), dendrogram="none",
          Rowv=FALSE, Colv=FALSE,trace= "none",
          ylab = "", cellnote = round(as.matrix(table.raw),2),
          notecol = "black", notecex = 2,
          xlab = "", col = GnPalette,
          #main = "Two-dimensional prognostic scores",
          labCol = c(1, 2, 3, 4, 5),
          labRow = c(1, 2, 3, 4, 5),
          cexRow = 2, cexCol = 2,
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          colsep=1:ncol(table.mat),
          rowsep=1:nrow(table.mat),
          sepcolor = "gray100",
          sepwidth = c(0.005, 0.005),
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.1,4),
          key.title = "", mgp = c(3,1,0),
          srtCol=0,
          key.xtickfun=function(){
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 2, line = -1.1, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.03", 1, cex = 2, line = -1.1, adj = 0.51, outer = TRUE, xpd = NA)
            #mtext(“0.4", 1, cex = 2, line = -1.1, adj = 0.62, outer = TRUE, xpd = NA)
            mtext("0.06", 1, cex = 2, line = -1.1, adj = 0.92, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 0.06, 0.06/20),
          offsetRow = -45, density.info = "none")
mtext("Transplant", 1, cex = 2, outer = TRUE, xpd = NA, line = -5)
mtext("Death", 4, cex = 2, outer = TRUE, xpd = NA, line = 0)
mtext("Two-dimensional prognostic scores", 3, cex = 2, outer = TRUE, xpd = NA, line = -1)
dev.off()

# (2) exact macthing based risk class
matchingdat = partdat
age.cate = ifelse(partdat$age + 50 < 40, 1, 
                  ifelse(partdat$age + 50 < 50, 2, 
                         ifelse(partdat$age + 50 < 60, 3, 4)))
diabetes.cate = ifelse(partdat$diabetes == 0, 1, 2)
bloodtype.cate = ifelse(partdat$blood_a == 1, 1, 
                        ifelse(partdat$blood_b == 1, 2, 
                               ifelse(partdat$blood_ab == 1, 3, 4)))
matchingdat$age.cate = age.cate
matchingdat$diabetes.cate = diabetes.cate
matchingdat$bloodtype.cate = bloodtype.cate

rownames(matchingdat) = as.character(c(1:nrow(matchingdat)))

table.mat = ftable(matchingdat$age.cate, matchingdat$diabetes.cate, matchingdat$bloodtype.cate)
template.mat = as.matrix(table.mat) 
prop.mat = template.mat / nrow(matchingdat)

print(xtable(prop.mat, digits = 3))


pdf("Figure/partdat_template_matching.pdf")
par(cex.main=1.8, cex.lab = 4, cex.axis = 4, 
    mar=c(3,5,3,0), tcl = 0.5, oma = c(3, 3, 3, 2), xpd = FALSE)
heatmap.2(as.matrix(prop.mat), dendrogram="none",
          Rowv=FALSE, Colv=FALSE,trace= "none",
          ylab = "", cellnote = formatC(as.matrix(prop.mat), digits = 3, format = "f"),
          notecol = "black", notecex = 2,
          xlab = "", col = BlPalette, 
          #main = "Exact matching with three key covariates",
          labCol = c("A", "B", "AB", "O"),
          labRow = c("(A1, N)", "(A1, D)", "(A2, N)",
                     "(A2, D)", "(A3, N)", "(A3, D)",
                     "(A4, N)", "(A4, D)"),
          cexRow = 2, cexCol = 2,
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          colsep=1:ncol(prop.mat),
          rowsep=1:nrow(prop.mat),
          sepcolor = "gray100",
          sepwidth = c(0.005, 0.005),
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "", mgp = c(3,1,0),
          srtCol=0,
          key.xtickfun=function(){
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 2, line = -1.1, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.05", 1, cex = 2, line = -1.1, adj = 0.51, outer = TRUE, xpd = NA)
            #mtext(“0.4", 1, cex = 2, line = -1.1, adj = 0.62, outer = TRUE, xpd = NA)
            mtext("0.10", 1, cex = 2, line = -1.1, adj = 0.92, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 0.10, 0.10/20),
          offsetRow = -45, density.info = "none")
mtext("Blood type", 1, cex = 2, outer = TRUE, xpd = NA, line = -5)
mtext("Age X Diabetes", 4, cex = 2, outer = TRUE, xpd = NA, line = 0)
mtext("Exact matching with three key covariates", 3, cex = 2, outer = TRUE, xpd = NA, line = -1)
dev.off()
