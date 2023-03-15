#!/usr/bin/Rscript

## read in

data = read.csv("data_binaries_2.0.0.csv", header = F)
energy = data$V2
rpA = data$V9
rsA = data$V8
r_sigma = data$V55
abs_rpA_minus_rpB = data$V23
ipA_minus_eaB_div_rpA = data$V47
abs_rsA_minus_rsB = data$V22
eaB_minus_ipB_div_rpA_sq = data$V52
ipA_minus_eaA_div_rpA = data$V44
ipA_minus_eaB_div_rsA = data$V43
hA_minus_lA = data$V34

# Simliarity (plots)
# 01,02
# 03,04,06
# 05,08,09
# 07
# 10

## It is expected that the energy is the cause of all of them
pair01 = data.frame(energy, rpA)
pair02 = data.frame(energy, rsA)
pair03 = data.frame(energy, r_sigma)
pair04 = data.frame(energy, abs_rpA_minus_rpB)
pair05 = data.frame(energy, ipA_minus_eaB_div_rpA)
pair06 = data.frame(energy, abs_rsA_minus_rsB)
pair07 = data.frame(energy, eaB_minus_ipB_div_rpA_sq)
pair08 = data.frame(energy, ipA_minus_eaA_div_rpA)
pair09 = data.frame(energy, ipA_minus_eaB_div_rsA)
pair10 = data.frame(energy, hA_minus_lA)

pair11 = data.frame(energy, rpA, rsA, r_sigma, abs_rpA_minus_rpB, ipA_minus_eaB_div_rpA, abs_rsA_minus_rsB, eaB_minus_ipB_div_rpA_sq, ipA_minus_eaA_div_rpA, ipA_minus_eaB_div_rsA, hA_minus_lA)
pair12 = data.frame(energy, rpA, r_sigma, ipA_minus_eaB_div_rpA, eaB_minus_ipB_div_rpA_sq, hA_minus_lA)

## write files
write.table(pair01, file="pair01.tab", quote=F, row.names=F, col.names=F, sep="\t")
write.table(pair02, file="pair02.tab", quote=F, row.names=F, col.names=F, sep="\t")
write.table(pair03, file="pair03.tab", quote=F, row.names=F, col.names=F, sep="\t")
write.table(pair04, file="pair04.tab", quote=F, row.names=F, col.names=F, sep="\t")
write.table(pair05, file="pair05.tab", quote=F, row.names=F, col.names=F, sep="\t")
write.table(pair06, file="pair06.tab", quote=F, row.names=F, col.names=F, sep="\t")
write.table(pair07, file="pair07.tab", quote=F, row.names=F, col.names=F, sep="\t")
write.table(pair08, file="pair08.tab", quote=F, row.names=F, col.names=F, sep="\t")
write.table(pair09, file="pair09.tab", quote=F, row.names=F, col.names=F, sep="\t")
write.table(pair10, file="pair10.tab", quote=F, row.names=F, col.names=F, sep="\t")
write.table(pair11, file="pair10.tab", quote=F, row.names=F, col.names=F, sep="\t")
write.table(pair12, file="pair10.tab", quote=F, row.names=F, col.names=F, sep="\t")
