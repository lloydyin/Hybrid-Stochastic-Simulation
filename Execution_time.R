windowsFonts(Times=windowsFont("Times New Roman"))
par(family = "Times", mar=c(3,3,2,1), mgp = c(1.6, 0.4, 0),
    cex.main = 1.5,   
    cex.lab = 1.3, 
    cex.axis = 1.1, 
    cex.sub = 1.2) 

lines <- readLines("Execution_time.txt")

DT <- read.table(text = lines[2:8], header = TRUE, sep = "\t")
FRP <- read.table(text = lines[11:17], header = TRUE, sep = "\t")
improvement <- read.table(text = lines[20:24], header = TRUE, sep = "\t")
Ratio <- read.table(text = lines[27:30], header = TRUE, sep = "\t")
Ratio2 <- read.table(text = lines[32:35], header = TRUE, sep = "\t")

par(mfrow=c(1,2))

labels <- paste0("+", improvement$Improvement[1])
p11 = barplot(c(improvement$SSA[1],improvement$New[1]),
        col = adjustcolor(c("black","blue"), alpha.f = 0.5),
        names.arg = c("SSA","New Algorithm"),
        ylab="Execution time (s)",
        ylim = c(0, max(c(improvement$SSA[1],improvement$New[1])) + 5.5))
title(main = "FRP", cex.main = 1.5, font.main = 1)
text(x = p11[2],
     y = improvement$New[1] + 1.5,
     labels = labels,
     cex = 0.9, col = "black", font = 1)

p12 = barplot(c(improvement$SSA[2],improvement$New[2:4]),
        col = adjustcolor(c("black","blue","orange","red"), alpha.f = 0.5),
        names.arg = c("SSA","Random", "Survival", "Uniform"),
        ,ylab="Execution time (s)",
        ylim = c(0, max(c(improvement$SSA[2],improvement$New[2:4])) + 20))
title(main = "DT", cex.main = 1.5, font.main = 1)
labels <- paste0("+", improvement$Improvement[2:4])
text(x = p12[2:4],
     y = improvement$New[2:4] + 10,
     labels = labels,
     cex = 0.9, col = "black", font = 1)


plot(FRP$Scale,FRP$SSA,type="b",col="black",pch=1,ylim=c(-1,max(FRP$SSA)+1),xlab="Number of monomer particles simulated",ylab="Execution time (s)")
title(main = "FRP", cex.main = 1.5, font.main = 1)
lines(FRP$Scale,FRP$New,type="b",col="blue",pch=2)
legend("topleft",c("SSA","New"),col=c("black","blue"),lty=rep(1,2),pch=1:2,cex=0.9)


plot(DT$Scale,DT$SSA,type="b",col="black",pch=16,ylim=c(-10,max(DT$SSA)+10),xlab="Number of monomer particles simulated",ylab="Execution time (s)")
title(main = "DT", cex.main = 1.5, font.main = 1)
lines(DT$Scale,DT$For,type="b",col="blue",pch=17)
lines(DT$Scale,DT$Parallelogram,type="b",col="orange",pch=18)
lines(DT$Scale,DT$Rectangle,type="b",col="red",pch=19)
legend("topleft",c("SSA","Random", "Survival Time", "Uniform"),col=c("black","blue","orange","red"),lty=rep(1,4),pch=16:19,cex=0.9)


Ratio$Ratio <- factor(Ratio$Ratio, levels = c("1000:1", "100:1", "10:1"))
plot(1:3, Ratio$SSA, type = "b", pch = 1, col = "black",
     xaxt = "n", xlab = "Monomer-to-initiator ratio", ylab="Execution time (s)", 
     ylim = range(0,max(Ratio$SSA)))
title(main = "FRP", cex.main = 1.5, font.main = 1)
lines(1:3, Ratio$New, type = "b", pch = 2, col = "blue")
axis(1, at = 1:3, labels = Ratio$Ratio)
legend("bottomright", legend = c("SSA", "New"),
       col = c("black", "blue"), pch = 1:2, lty = 1,cex=0.9)

Ratio2$Ratio <- factor(Ratio2$Ratio, levels = c("1000:1", "100:1", "10:1"))
plot(1:3, Ratio2$SSA, type = "b", pch = 16, col = "black",
     xaxt = "n", xlab = "Monomer-to-initiator ratio", ylab="Execution time (s)", 
     ylim = range(0,max(Ratio2$SSA)))
title(main = "DT", cex.main = 1.5, font.main = 1)
lines(1:3, Ratio2$random, type = "b", pch = 17, col = "blue")
lines(1:3, Ratio2$pxsbx, type = "b", pch = 18, col = "orange")
lines(1:3, Ratio2$jx, type = "b", pch = 19, col = "red")
axis(1, at = 1:3, labels = Ratio2$Ratio)
legend("bottomright", c("SSA","Random", "Survival Time", "Uniform"),
       col = c("black","blue","orange","red"), pch =16:19, lty = 1,cex=0.9)