library(qqman)
library(ggplot2)
library(calibrate)

## define function

manhattan1<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
                                                                               "gray60"), chrlabs = NULL,
                      genomewideline = -log10(9.45e-08), highlight1 = NULL, highlight2 = NULL, 
                      highlight3=NULL, highlight4=NULL,
                      highlight5=NULL, highlight6=NULL,
                      highlight7=NULL, highlight8=NULL,
                      highlight9=NULL, highlight10=NULL,
                      highlight11=NULL, highlight12=NULL,
                      logp = TRUE,
                      ...)
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x)))
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x)))
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x)))
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x)))
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]]))
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]]))
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]]))
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]]))
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    options(scipen = 999)
    d$pos = d$BP/1e+06
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index ==
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP +
          lastbase
      }
      ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR ==
                                                           i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
                   las = 1, pch = 21, xlim = c(xmin, xmax), ylim = c(0,
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log10))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos,
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (genomewideline)
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight1)) {
    if (any(!(highlight1 %in% d$SNP)))
      warning("You're trying to highlight1 SNPs that don't exist in your results.")
    d.highlight1 = d[which(d$SNP %in% highlight1), ]
    with(d.highlight1, points(pos, logp, col = "#FA7475", pch = 20,
                              ...))
  }
  if (!is.null(highlight2)) {
    if (any(!(highlight2 %in% d$SNP)))
      warning("You're trying to highlight2 SNPs that don't exist in your results.")
    d.highlight2 = d[which(d$SNP %in% highlight2), ]
    with(d.highlight2, points(pos, logp, col = "#E18A00", pch = 20,
                              ...))
  }
  if (!is.null(highlight3)) {
    if (any(!(highlight3 %in% d$SNP)))
      warning("You're trying to highlight3 SNPs that don't exist in your results.")
    d.highlight3 = d[which(d$SNP %in% highlight3), ]
    with(d.highlight3, points(pos, logp, col = "#BB9D00", pch = 20,
                              ...))
  }
  if (!is.null(highlight4)) {
    if (any(!(highlight4 %in% d$SNP)))
      warning("You're trying to highlight4 SNPs that don't exist in your results.")
    d.highlight4 = d[which(d$SNP %in% highlight4), ]
    with(d.highlight4, points(pos, logp, col = "#81AD00", pch = 20,
                              ...))
  }
  if (!is.null(highlight5)) {
    if (any(!(highlight5 %in% d$SNP)))
      warning("You're trying to highlight5 SNPs that don't exist in your results.")
    d.highlight5 = d[which(d$SNP %in% highlight5), ]
    with(d.highlight5, points(pos, logp, col = "#00B92F", pch = 20,
                              ...))
  }
  if (!is.null(highlight6)) {
    if (any(!(highlight6 %in% d$SNP)))
      warning("You're trying to highlight6 SNPs that don't exist in your results.")
    d.highlight6 = d[which(d$SNP %in% highlight6), ]
    with(d.highlight6, points(pos, logp, col = "#00C087", pch = 20,
                              ...))
  }
  if (!is.null(highlight7)) {
    if (any(!(highlight7 %in% d$SNP)))
      warning("You're trying to highlight7 SNPs that don't exist in your results.")
    d.highlight7 = d[which(d$SNP %in% highlight7), ]
    with(d.highlight7, points(pos, logp, col = "#00BFC2", pch = 20,
                              ...))
  }
  if (!is.null(highlight8)) {
    if (any(!(highlight8 %in% d$SNP)))
      warning("You're trying to highlight8 SNPs that don't exist in your results.")
    d.highlight8 = d[which(d$SNP %in% highlight8), ]
    with(d.highlight8, points(pos, logp, col = "#00B4EE", pch = 20,
                              ...))
  }
  if (!is.null(highlight9)) {
    if (any(!(highlight9 %in% d$SNP)))
      warning("You're trying to highlight9 SNPs that don't exist in your results.")
    d.highlight9 = d[which(d$SNP %in% highlight9), ]
    with(d.highlight9, points(pos, logp, col = "#5B9DFF", pch = 20,
                              ...))
  }
  if (!is.null(highlight10)) {
    if (any(!(highlight10 %in% d$SNP)))
      warning("You're trying to highlight10 SNPs that don't exist in your results.")
    d.highlight10 = d[which(d$SNP %in% highlight10), ]
    with(d.highlight10, points(pos, logp, col = "#C57DFF", pch = 20,
                               ...))
  }
  if (!is.null(highlight11)) {
    if (any(!(highlight11 %in% d$SNP)))
      warning("You're trying to highlight11 SNPs that don't exist in your results.")
    d.highlight11 = d[which(d$SNP %in% highlight11), ]
    with(d.highlight11, points(pos, logp, col = "#F464E4", pch = 20,
                               ...))
  }
  if (!is.null(highlight12)) {
    if (any(!(highlight12 %in% d$SNP)))
      warning("You're trying to highlight11 SNPs that don't exist in your results.")
    d.highlight12 = d[which(d$SNP %in% highlight12), ]
    with(d.highlight12, points(pos, logp, col = "#FF64B1", pch = 20,
                               ...))
  }
}

## a sequence of 8 different colours

gg_color_hue <- function(n) {
  hues = seq(8, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 8
cols = gg_color_hue(n)

## read in results

results <- read.table('manhattan_input.txt',header =T, fill=T)

## Lead SNPs named as gene symbols for annotation

names_to_anno <- c("ADIPOQ","ABO","LEPR","APOB","APOA1","IL6R","FADS1","ADCY3")
new_cols <- c("firebrick1","darkorange","gold","green2","lightskyblue","dodgerblue2",cols[7],cols[8])

## generate manhattan

png("manhattan_mQTL_900dpi.png", width=32, height=12, units='in', res=900)
par(mar=c(14,14,14,14), mgp=c(3,1,0))
manhattan3(results, annotatePval=names_to_anno, highlight1=set1, highlight2=set2,
           highlight3=set3,highlight4=set4,highlight5=set5,
           highlight6=set6,highlight7=set7,highlight8=set8,highlight9=set9,cex.axis=1.5,cex.lab=2,
           chrlabs=as.character(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)))
par(xpd=T)
legend("topright",title=expression(bold("Cardiovascular traits")),c("Adiponectin", "Apolipoprotein A1", "Apolipoprotein B", "Body Mass Index",
               "C-reactive Protein","Cholesterol","Interleukin-6","Low Density Lipoprotein"), pch= c(20,20,20,20,
                                                                                                     20,20,20,20), col=new_cols, cex=2)
text(27000000,14,labels = c('CELSR2'), cex=2) ## adjust placement of CELSR2 slightly
mtext('22', side=1,cex=1.5,line=1,at=2862081700)
#text(0.5,0.5,labels = c('22'), cex=1.5)
dev.off()
