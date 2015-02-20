#build a model timeline (Fig. 1A)
#4.0in x 2.5in

par(mar=c(0.5,0.5,0.5,0.5))
par(cex=1)

start <- -14
end <- 74

samples <- c(0, 21, 25, 45, 47, 69, 74)

plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='',xlim=c(start-1,end+1),ylim=c(-4,5))

line <- 2


polygon(c(5,5,10,10),c(line+1,line-1,line-1,line+1), col='light blue', border=NA)   # DSS 1 box
polygon(c(26,26,31,31),c(line+1,line-1,line-1,line+1), col='light blue', border=NA) # DSS 2 box
polygon(c(47,47,52,52),c(line+1,line-1,line-1,line+1), col='light blue', border=NA) # DSS 3 box
polygon(c(0.8,0.8,-0.8,-0.8),c(line+1,line-1,line-1,line+1), col='red', border=NA)  # AOM box

segments(start,line,end,line, lwd=4) # timeline

#major ticks
for(i in seq(-10,70,10)){
  segments(i,line+0.8,i,line-0.8, lwd=2)
}

# minor ticks
for(i in seq(-10,70,5)){
  segments(i,line+0.4,i,line-0.5, lwd=1.5)
}

#day label
text(seq(-10,70,10),line-1.5,labels=seq(-10,70,10), cex=0.8, font=2)
text((end+start)/2,line-2.5,'Day', cex=0.8, font=2)

offset <- 2
text(0+offset,line+2.3,'AOM', srt=60, cex=0.6, font=2)
text(7.5+offset,line+2.3,'DSS 1', srt=60, cex=0.6, font=2)
text(28.5+offset,line+2.3,'DSS 2', srt=60, cex=0.6, font=2)
text(49.5+offset,line+2.3,'DSS 3', srt=60, cex=0.6, font=2)


#sampling times
points(x=samples, y=rep(line,length(samples))+0.35, pch=25, bg="black", cex=1)

abx_line <- line-3.5
polygon(c(start, start, end, end),c(abx_line+0.5,abx_line-0.5,abx_line-0.5,abx_line+0.5), col='gray', border=NA)     # AOM box
text((end+start)/2, abx_line, label="Antibiotic treatment", cex=0.7, font=2)

intervention_line <- abx_line-1.2
polygon(c(start, start, 5, 5),c(intervention_line+0.5,intervention_line-0.5,intervention_line-0.5,intervention_line+0.5), col='gray', border=NA)     # AOM box
text((5+start)/2, intervention_line, label="Intervention 1", cex=0.7, font=2)

intervention_line <- intervention_line-1.2
polygon(c(26, 26, end, end),c(intervention_line+0.5,intervention_line-0.5,intervention_line-0.5,intervention_line+0.5), col='gray', border=NA)     # AOM box
text((end+26)/2, intervention_line, label="Intervention 2", cex=0.7, font=2)
