PredictDiag <- function(Hbvarinats) {
  WT <- HbAB$`sp|P68871|HBB_HUMAN`
  re <- list()
  for (i in 1:length(names(Hbvariants))){
    MT <- matrix(Hbvariants[[i]], byrow = TRUE)
    z <- MT==WT
    re[[i]]<- data.frame(variant=names(Hbvariants)[i], mut.site_N = which(z==FALSE)-1,mut.site_C=148-which(z==FALSE),
                         WT.AA= WT[which(z==FALSE)],MT.AA= MT[which(z==FALSE)])
  }
  re1 <- do.call(rbind,re)
  re1$Mutation <- paste(re1$WT.AA, re1$mut.site_N,re1$MT.AA)
  IDs <- re1$mut.site_N
  mutMass <- function(residue)
  { if (residue == "A")
    mass = 71.03711
  if (residue == "R")
    mass = 156.10111
  if (residue == "N")
    mass = 114.04293
  if (residue == "D")
    mass = 115.02694
  if (residue == "E")
    mass = 129.04259
  if (residue == "Q")
    mass = 128.05858
  if (residue == "G")
    mass = 57.02146
  if (residue == "H")
    mass = 137.05891
  if (residue == "I")
    mass = 113.08406
  if (residue == "L")
    mass = 113.08406
  if (residue == "K")
    mass = 128.09496
  if (residue == "M")
    mass = 131.04049
  if (residue == "F")
    mass = 147.06841
  if (residue == "P")
    mass = 97.05276
  if (residue == "S")
    mass = 87.03203
  if (residue == "T")
    mass = 101.04768
  if (residue == "W")
    mass = 186.07931
  if (residue == "Y")
    mass = 163.06333
  if (residue == "V")
    mass = 99.06841
  if (residue == "C")
    mass = 103.00919
  return(mass)
  }
  re1$mutMass <- round(sapply(re1$MT.AA,mutMass)-sapply(re1$WT.AA,mutMass), 5)

  df <- dplyr::left_join(re1, diag, by = "mut.site_N")
  df2 <- df[, c(1,2,3,6,7,9,10,11,12)]

  #2
  #bc ions
  HbA.bc1 <- subset(HbA.B, Ion.type=="b"|Ion.type=="c")
  HbA.bc <-as.data.table(HbA.bc1)
  #yz ions
  HbA.yz1 <- subset(HbA.B, Ion.type=="y"|Ion.type=="z")
  HbA.yz <-as.data.table(HbA.yz1)
  #bc ions
  m <- length(df2$diag.N.term.start)
  re <- list()
  for (i in 1:m){
    s <- as.numeric(df2$diag.N.term.start[i])
    e <- as.numeric(df2$diag.N.term.end[i])
    re[[i]] <- cbind(variant=df2$variant[i],HbA.bc[Ion.num %between% c(s, e)], Mutation = df2$Mutation[i], MutMass.Da= df2$mutMass[i])
  }
  re2 <- do.call(rbind,re)
  #yz ions
  n <- length(df2$diag.C.term.start)
  re3 <- list()
  for (j in 1:n){
    s <- as.numeric(df2$diag.C.term.start[j])
    e <- as.numeric(df2$diag.C.term.end[j])
    re3[[j]] <- cbind(variant=df2$variant[j],HbA.yz[Ion.num %between% c(s, e)],Mutation = df2$Mutation[j], MutMass.Da= df2$mutMass[j])
  }
  re4 <- do.call(rbind,re3)

  re5 <- rbind(re2, re4) %>% filter(!is.na(variant))
  re6 <- re5[with(re5, order(variant))]

  ID2 <- names(Hbvariants)
  out2 <- list()
  for (j in 1:length(ID2)){
    Hbseq1 <- getSequence(Hbvariants, as.string = TRUE)
    Hbseq2 <- sub("M","",Hbseq1[[j]])
    input <- monomz(Hbseq2, fragments = "by")

    input2 <- monomz(Hbseq2, fragments = "cz")

    out1 <- rbind(input,input2)
    out1$varints <- ID2[j]
    out2[[j]] <- out1
  }
  all.frags <- do.call(rbind,out2)

  ID3 <- unique(re6$variant)
  out <- list()
  for (i in 1:length(ID3)){
    sub <- subset(re6, variant==ID3[i])
    sub %>% mutate_if(is.numeric, as.character) -> sub
    sub2 <- subset(all.frags, varints==ID3[i])
    sub2 %>% mutate_if(is.numeric, as.character) -> sub2
    out[[i]] <-inner_join(sub2, sub, by = c("Ion_type"="Ion.type", "Ion_num"="Ion.num"))
  }
  relist <- do.call(rbind,out)
  relist2 <- relist %>% separate(variant,c("x","Variant"), sep = 3)
  names(relist2)[5] <- "Ref_Mass"
  relist2$Name <- paste(relist2$Ion,relist2$Variant, sep = "_" )
  relist3 <- relist2[, c(12,2,3,4,5,8,10,11)]
return(relist3)
}
