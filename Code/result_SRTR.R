source("analysis_SRTR.R")
# read the results from analysis_SRTR.R
# evaluate the center effects at t=1800 days (time.list[18])
noweight.delta1.est = noweight_onetime$delta1[,18]
noweight.delta2.est = noweight_onetime$delta2[,18]
score.delta1.est = score_onetime$delta1[,18]
score.delta2.est = score_onetime$delta2[,18]
matching.delta1.est = matching_onetime$delta1[,18]
matching.delta2.est = matching_onetime$delta2[,18]
# lower/upper bound of the center effects for death
noweight.delta1.low = noweight.delta1.est - 1.96*sqrt((210/212)*noweight_onetime$F1var.if + sum((1/212)^2*noweight_onetime$F1var.if[!is.na(noweight_onetime$F1var.if)]))
noweight.delta1.high = noweight.delta1.est + 1.96*sqrt((210/212)*noweight_onetime$F1var.if + sum((1/212)^2*noweight_onetime$F1var.if[!is.na(noweight_onetime$F1var.if)]))
score.delta1.low = score.delta1.est - 1.96*sqrt((210/212)*score_onetime$F1var.if + sum((1/212)^2*score_onetime$F1var.if[!is.na(noweight_onetime$F1var.if)]))
score.delta1.high = score.delta1.est + 1.96*sqrt((210/212)*score_onetime$F1var.if + sum((1/212)^2*score_onetime$F1var.if[!is.na(noweight_onetime$F1var.if)]))
matching.delta1.low = matching.delta1.est - 1.96*sqrt((210/212)*matching_onetime$F1var.if + sum((1/212)^2*matching_onetime$F1var.if[!is.na(noweight_onetime$F1var.if)]))
matching.delta1.high = matching.delta1.est + 1.96*sqrt((210/212)*matching_onetime$F1var.if + sum((1/212)^2*matching_onetime$F1var.if[!is.na(noweight_onetime$F1var.if)]))
# lower/upper bound of the center effects for transplant
noweight.delta2.low = noweight.delta2.est - 1.96*sqrt((210/212)*noweight_onetime$F2var.if + sum((1/212)^2*noweight_onetime$F2var.if[!is.na(noweight_onetime$F2var.if)]))
noweight.delta2.high = noweight.delta2.est + 1.96*sqrt((210/212)*noweight_onetime$F2var.if + sum((1/212)^2*noweight_onetime$F2var.if[!is.na(noweight_onetime$F2var.if)]))
score.delta2.low = score.delta2.est - 1.96*sqrt((210/212)*score_onetime$F2var.if + sum((1/212)^2*score_onetime$F2var.if[!is.na(noweight_onetime$F2var.if)]))
score.delta2.high = score.delta2.est + 1.96*sqrt((210/212)*score_onetime$F2var.if + sum((1/212)^2*score_onetime$F2var.if[!is.na(noweight_onetime$F2var.if)]))
matching.delta2.low = matching.delta2.est - 1.96*sqrt((210/212)*matching_onetime$F2var.if + sum((1/212)^2*matching_onetime$F2var.if[!is.na(noweight_onetime$F2var.if)]))
matching.delta2.high = matching.delta2.est + 1.96*sqrt((210/212)*matching_onetime$F2var.if + sum((1/212)^2*matching_onetime$F2var.if[!is.na(noweight_onetime$F2var.if)]))

## cross-tabulation with respect to the direction of significance in two center efffect estimates
ind1 = ifelse(score.delta1.low <0 & score.delta1.high <0, 1, 
              ifelse(score.delta1.low <0 & score.delta1.high >0, 2, 
                     ifelse(score.delta1.low >0 & score.delta1.high >0, 3, NA)))
ind2 = ifelse(score.delta2.low <0 & score.delta2.high <0, 1, 
              ifelse(score.delta2.low <0 & score.delta2.high >0, 2, 
                     ifelse(score.delta2.low >0 & score.delta2.high >0, 3, NA)))
table(ind1, ind2)
xtable(table(ind1, ind2))

## estimates from three different methods for the centers with ten largest center effect estimates using prognostic scores
highest.delta1 = matrix(NA, 10, 3)
rownames(highest.delta1) = c(1:10)
colnames(highest.delta1) = c("PS", "Matching", "Unweighted")
for(i in 1:10){
  highest.delta1[i,1] = paste0(formatC(score.delta1.est[score.high.first[i]], digits = 3, format = "f"),
                               " (", formatC(score.delta1.low[score.high.first[i]], digits = 3, format = "f"), ", ", 
                               formatC(score.delta1.high[score.high.first[i]], digits = 3, format = "f"),")")
  highest.delta1[i,2] = paste0(formatC(matching.delta1.est[score.high.first[i]], digits = 3, format = "f"),
                               " (", formatC(matching.delta1.low[score.high.first[i]], digits = 3, format = "f"), ", ", 
                               formatC(matching.delta1.high[score.high.first[i]], digits = 3, format = "f"),")")
  highest.delta1[i,3] = paste0(formatC(noweight.delta1.est[score.high.first[i]], digits = 3, format = "f"),
                               " (", formatC(noweight.delta1.low[score.high.first[i]], digits = 3, format = "f"), ", ", 
                               formatC(noweight.delta1.high[score.high.first[i]], digits = 3, format = "f"),")")
}
print(xtable(highest.delta1))

highest.delta2 = matrix(NA, 10, 3)
rownames(highest.delta2) = c(1:10)
colnames(highest.delta2) = c("PS", "Matching", "Unweighted")
for(i in 1:10){
  highest.delta2[i,1] = paste0(formatC(score.delta2.est[score.high.second[i]], digits = 3, format = "f"),
                               " (", formatC(score.delta2.low[score.high.second[i]], digits = 3, format = "f"), ", ", 
                               formatC(score.delta2.high[score.high.second[i]], digits = 3, format = "f"),")")
  highest.delta2[i,2] = paste0(formatC(matching.delta2.est[score.high.second[i]], digits = 3, format = "f"),
                               " (", formatC(matching.delta2.low[score.high.second[i]], digits = 3, format = "f"), ", ", 
                               formatC(matching.delta2.high[score.high.second[i]], digits = 3, format = "f"),")")
  highest.delta2[i,3] = paste0(formatC(noweight.delta2.est[score.high.second[i]], digits = 3, format = "f"),
                               " (", formatC(noweight.delta2.low[score.high.second[i]], digits = 3, format = "f"), ", ", 
                               formatC(noweight.delta2.high[score.high.second[i]], digits = 3, format = "f"),")")
}
print(xtable(highest.delta2))

## scatter plots to compare three different estimation methods
center.ind = names(table(partdat$center))
size.more200.ind = which(table(partdat$center) > 200) # only consider centers with center size > 200 (189 centers in total)
# prognostic scores vs. unweighted method
pdf("Figure/scatter_ps.pdf",  width = 14, height = 6)
par(mfrow = c(1,2),   mar = c(5,7,5,3),  cex.lab = 1.5, 
    cex.main = 2.0, cex.axis = 1.5, tcl = 0.5, lwd = 2)
plot(score.delta1.est[size.more200.ind], noweight.delta1.est[size.more200.ind],
     ylim = c(-0.4, 0.4), xlim = c(-0.4,0.4), 
     ylab = expression(paste("Unweighted ", hat(delta)[jk])),
     xlab = expression(paste("Prognostic score-weighted ", hat(delta)[jk])),
     main = "Death (k=1)", 
     pch = 19, col = "pink")
lines(x = c(-0.4, 0.4), y = c(-0.4, 0.4))
abline(h = 0, lwd = 0.5, col = "gray")
abline(v = 0, lwd = 0.5, col = "gray")
text(-0.2, 0.3, paste("Cor = ", 
                      formatC(cor(score.delta1.est[size.more200.ind], noweight.delta1.est[size.more200.ind]), 2, format = "f" )),
     cex = 1.5)
plot(score.delta2.est[size.more200.ind], noweight.delta2.est[size.more200.ind],
     ylim = c(-0.4, 0.4), xlim = c(-0.4,0.4), 
     ylab = expression(paste("Unweighted ", hat(delta)[jk])),
     xlab = expression(paste("Prognostic score-weighted ", hat(delta)[jk])),
     main = "Transplant (k=2)", 
     pch = 19, col = "slategray1")
lines(x = c(-0.4, 0.4), y = c(-0.4, 0.4))
abline(h = 0, lwd = 0.5, col = "gray")
abline(v = 0, lwd = 0.5, col = "gray")
text(-0.2, 0.3, paste("Cor = ", 
                      formatC(cor(score.delta2.est[size.more200.ind], noweight.delta2.est[size.more200.ind]), 2, format = "f" )),
     cex = 1.5)
dev.off()

# prognostic scores vs. exact matching method
pdf("Figure/scatter_ps_matching.pdf",  width = 14, height = 6)
par(mfrow = c(1,2),   mar = c(5,7,5,3),  cex.lab = 1.5, 
    cex.main = 2.0, cex.axis = 1.5, tcl = 0.5, lwd = 2)
plot(score.delta1.est[size.more200.ind], matching.delta1.est[size.more200.ind],
     ylim = c(-0.4, 0.4), xlim = c(-0.4,0.4), 
     ylab = expression(paste("Exact matching-weighted ", hat(delta)[jk])),
     xlab = expression(paste("Prognostic score-weighted ", hat(delta)[jk])),
     main = "Death (k=1)", 
     pch = 19, col = "pink")
lines(x = c(-0.4, 0.4), y = c(-0.4, 0.4))
abline(h = 0, lwd = 0.5, col = "gray")
abline(v = 0, lwd = 0.5, col = "gray")
text(-0.2, 0.3, paste("Cor = ", 
                      formatC(cor(score.delta1.est[size.more200.ind], matching.delta1.est[size.more200.ind]), 2, format = "f" )),
     cex = 1.5)
plot(score.delta2.est[size.more200.ind], matching.delta2.est[size.more200.ind],
     ylim = c(-0.4, 0.4), xlim = c(-0.4,0.4), 
     ylab = expression(paste("Exact matching-weighted ", hat(delta)[jk])),
     xlab = expression(paste("Prognostic score-weighted ", hat(delta)[jk])),
     main = "Transplant (k=2)", 
     pch = 19, col = "slategray1")
lines(x = c(-0.4, 0.4), y = c(-0.4, 0.4))
abline(h = 0, lwd = 0.5, col = "gray")
abline(v = 0, lwd = 0.5, col = "gray")
text(-0.2, 0.3, paste("Cor = ", 
                      formatC(cor(score.delta2.est[size.more200.ind], matching.delta2.est[size.more200.ind]), 2, format = "f" )),
     cex = 1.5)
dev.off()

# exact matching vs. unweighted method
pdf("Figure/scatter_matching.pdf",  width = 14, height = 6)
par(mfrow = c(1,2),   mar = c(5,7,5,3),  cex.lab = 1.5, 
    cex.main = 2.0, cex.axis = 1.5, tcl = 0.5, lwd = 2)
plot(matching.delta1.est[size.more200.ind], noweight.delta1.est[size.more200.ind],
     ylim = c(-0.5, 0.5), xlim = c(-0.5,0.5), 
     ylab = expression(paste("Unweighted ", hat(delta)[jk])),
     xlab = expression(paste("Exact matching-weighted ", hat(delta)[jk])),
     main = "Death (k=1)", 
     pch = 19, col = "pink")
points(matching.delta1.est[238], noweight.delta1.est[238], col = "black", pch = 19)
points(matching.delta1.est[142], noweight.delta1.est[142], col = "red", pch = 19)
lines(x = c(-0.5, 0.5), y = c(-0.5, 0.5))
abline(h = 0, lwd = 0.5, col = "gray")
abline(v = 0, lwd = 0.5, col = "gray")
text(-0.3, 0.4, paste("Cor = ", 
                      formatC(cor(matching.delta1.est[size.more200.ind], noweight.delta1.est[size.more200.ind]), 2, format = "f" )),
     cex = 1.5)
plot(matching.delta2.est[size.more200.ind], noweight.delta2.est[size.more200.ind],
     ylim = c(-0.5, 0.5), xlim = c(-0.5,0.5), 
     ylab = expression(paste("Unweighted ", hat(delta)[jk])),
     xlab = expression(paste("Exact matching-weighted ", hat(delta)[jk])),
     main = "Transplant (k=2)", 
     pch = 19, col = "slategray1")
points(matching.delta2.est[238], noweight.delta2.est[238], col = "black", pch = 19)
points(matching.delta2.est[142], noweight.delta2.est[142], col = "red", pch = 19)
lines(x = c(-0.5, 0.5), y = c(-0.4, 0.4))
abline(h = 0, lwd = 0.5, col = "gray")
abline(v = 0, lwd = 0.5, col = "gray")
text(-0.3, 0.4, paste("Cor = ", 
                      formatC(cor(matching.delta2.est[size.more200.ind], noweight.delta2.est[size.more200.ind]), 2, format = "f" )),
     cex = 1.5)
legend("bottomright", c("Center 820", "Center 457"), 
       col = c("black", "red"), pch = 19, cex = 1.3, bty = "n")
dev.off()

# covariates distribution of two centers (center 820 and center 457)
order(matching.delta1.est - noweight.delta1.est)
matching.delta1.est[238]
noweight.delta1.est[238]
matching.delta2.est[238]
noweight.delta2.est[238]

center238 = partdat[partdat$center == centername[238],]
summary(center238$age+50)
summary(partdat$age + 50)
mean(center238$diabetes)
mean(partdat$diabetes)

order(-matching.delta1.est + noweight.delta1.est)
matching.delta1.est[142]
noweight.delta1.est[142]
matching.delta2.est[142]
noweight.delta2.est[142]
center142 = partdat[partdat$center == centername[142],]
summary(center142$age+50)
summary(partdat$age + 50)
mean(center142$diabetes)
mean(partdat$diabetes)

# age distributions
pdf("Figure/cov_age_matching.pdf",  width = 10, height = 6)
par(mfrow = c(1,1),   mar = c(5,7,5,3),  cex.lab = 1.5, 
    cex.main = 2.0, cex.axis = 1.5, tcl = 0.5, lwd = 2)
plot(density(center238$age+50), xlab = "Age", main = "Age distribution",
     xlim = c(10, 80))
polygon(density(center238$age+50), col= rgb(0,0,0.8,alpha=0.4))
lines((density(partdat$age + 50)))
polygon(density(partdat$age+50), col= rgb(0.5,0.5,0.5,alpha=0.3))
lines((density(center142$age + 50)))
polygon(density(center142$age+50), col= rgb(0.8,0,0,alpha=0.4))
legend("topleft", c("Population level", "Center 820", "Center 457"), 
       col = c(rgb(0.5,0.5,0.5,alpha=0.3),rgb(0,0,0.8,alpha=0.4), rgb(0.8,0,0,alpha=0.4)), pch = 19, cex = 1.3, bty = "n")
dev.off()

# diabetic blood pressure distributions
mat = matrix(0, 4, 3)
colnames(mat) = c("Population level", "Center 820", "Center 457")
rownames(mat) = c("Diabetes", "Blood A", "Blood B", "Blood AB")
mat[1,] = c(mean(partdat$diabetes),
            mean(center238$diabetes), # center 820
            mean(center142$diabetes)) # center 475
mat[2,] = c(mean(partdat$blood_a),
            mean(center238$blood_a), # center 820
            mean(center142$blood_a)) # center 475
mat[3,] = c(mean(partdat$blood_b),
            mean(center238$blood_b), # center 820
            mean(center142$blood_b)) # center 475
mat[4,] = c(mean(partdat$blood_ab),
            mean(center238$blood_ab), # center 820
            mean(center142$blood_ab)) # center 475
print(xtable(mat, digits = 3))