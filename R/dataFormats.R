
# example sample data
Ind <- c("HG00096", "HG00097", "HG00099", "NA20827")
SEX <- c("male", "female", "female", "male")
POP <- c("GBR", "GBR", "GBR", "TSI")
sample_example <- data.frame(Ind, SEX, POP)

# example allele counts data
chr <- c(22, 22, 22)
start <- c(135031, 135031, 135031)
end <- c(135032, 135032, 135032)
ref <- c("G", "G", "G")
alt <- c("A", "A", "A")
Ind <- c("HG00276", "HG00282", "NA11831")
refCount <- c(19, 12, 10)
altCount <- c(0, 0, 0)
ASE_example <- data.frame(chr, start, end, ref, alt, Ind, refCount, altCount)

# example legend data
ID <- c("10:60515:C:T", "rs148087467:60523:T:G", "rs147855157:61372:A:C")
pos <- c(60515, 60523, 61372)
allele0 <- c("C", "T", "A")
allele1 <- c("T", "G", "C")
LEG_example <- data.frame(ID, pos, allele0, allele1)
