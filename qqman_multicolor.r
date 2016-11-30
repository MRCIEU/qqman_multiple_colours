manhattan1<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
"gray60"), chrlabs = NULL,
genomewideline = -log10(1e-07), highlight1 = NULL, highlight2 = NULL, 
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
las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
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
with(d.highlight1, points(pos, logp, col = "#F8766D", pch = 20,
...))
}
if (!is.null(highlight2)) {
if (any(!(highlight2 %in% d$SNP)))
warning("You're trying to highlight2 SNPs that don't exist in your results.")
d.highlight2 = d[which(d$SNP %in% highlight2), ]
with(d.highlight2, points(pos, logp, col = "#DE8C00", pch = 20,
...))
}
if (!is.null(highlight3)) {
if (any(!(highlight3 %in% d$SNP)))
warning("You're trying to highlight3 SNPs that don't exist in your results.")
d.highlight3 = d[which(d$SNP %in% highlight3), ]
with(d.highlight3, points(pos, logp, col = "#B79F00", pch = 20,
...))
}
if (!is.null(highlight4)) {
  if (any(!(highlight4 %in% d$SNP)))
    warning("You're trying to highlight4 SNPs that don't exist in your results.")
  d.highlight4 = d[which(d$SNP %in% highlight4), ]
  with(d.highlight4, points(pos, logp, col = "#7CAE00", pch = 20,
                            ...))
}
if (!is.null(highlight5)) {
  if (any(!(highlight5 %in% d$SNP)))
    warning("You're trying to highlight5 SNPs that don't exist in your results.")
  d.highlight5 = d[which(d$SNP %in% highlight5), ]
  with(d.highlight5, points(pos, logp, col = "#00BA38", pch = 20,
                            ...))
}
if (!is.null(highlight6)) {
  if (any(!(highlight6 %in% d$SNP)))
    warning("You're trying to highlight6 SNPs that don't exist in your results.")
  d.highlight6 = d[which(d$SNP %in% highlight6), ]
  with(d.highlight6, points(pos, logp, col = "#00C08B", pch = 20,
                            ...))
}
if (!is.null(highlight7)) {
  if (any(!(highlight7 %in% d$SNP)))
    warning("You're trying to highlight7 SNPs that don't exist in your results.")
  d.highlight7 = d[which(d$SNP %in% highlight7), ]
  with(d.highlight7, points(pos, logp, col = "#00BFC4", pch = 20,
                            ...))
}
if (!is.null(highlight8)) {
  if (any(!(highlight8 %in% d$SNP)))
    warning("You're trying to highlight8 SNPs that don't exist in your results.")
  d.highlight8 = d[which(d$SNP %in% highlight8), ]
  with(d.highlight8, points(pos, logp, col = "#00B4F0", pch = 20,
                            ...))
}
if (!is.null(highlight9)) {
  if (any(!(highlight9 %in% d$SNP)))
    warning("You're trying to highlight9 SNPs that don't exist in your results.")
  d.highlight9 = d[which(d$SNP %in% highlight9), ]
  with(d.highlight9, points(pos, logp, col = "#619CFF", pch = 20,
                            ...))
}
if (!is.null(highlight10)) {
  if (any(!(highlight10 %in% d$SNP)))
    warning("You're trying to highlight10 SNPs that don't exist in your results.")
  d.highlight10 = d[which(d$SNP %in% highlight10), ]
  with(d.highlight10, points(pos, logp, col = "#C77CFF", pch = 20,
                            ...))
}
if (!is.null(highlight11)) {
  if (any(!(highlight11 %in% d$SNP)))
    warning("You're trying to highlight11 SNPs that don't exist in your results.")
  d.highlight11 = d[which(d$SNP %in% highlight11), ]
  with(d.highlight11, points(pos, logp, col = "#F564E3", pch = 20,
                            ...))
}
if (!is.null(highlight12)) {
if (any(!(highlight12 %in% d$SNP)))
warning("You're trying to highlight11 SNPs that don't exist in your results.")
d.highlight12 = d[which(d$SNP %in% highlight12), ]
with(d.highlight12, points(pos, logp, col = "#FF64B0", pch = 20,
...))
}
}