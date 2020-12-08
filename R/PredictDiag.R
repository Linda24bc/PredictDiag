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
  find.diag <- list()
  for (j in 1:length(IDs)){
    find <- diag[diag$mut.site == IDs[j], ]
    find.diag[[j]] <- cbind(variant= re1$variant[j],mut.site_C=re1$mut.site_C[j], find, Mutation= re1$Mutation[j] )
  }
  find.diag1 <- do.call(rbind,find.diag) %>% cbind(MT.AA=re1$MT.AA)
  find.diag2 <- find.diag1[, c(1,2,3,4,10,5,6,7,8,9)]
  #2
  #bc ions
  HbA.bc1 <- subset(HbA.B, Ion.type=="b"|Ion.type=="c")
  HbA.bc <-as.data.table(HbA.bc1)
  #yz ions
  HbA.yz1 <- subset(HbA.B, Ion.type=="y"|Ion.type=="z")
  HbA.yz <-as.data.table(HbA.yz1)
  #bc ions
  m <- length(find.diag2$diag.N.term.start)
  re <- list()
  for (i in 1:m){
    s <- as.numeric(find.diag2$diag.N.term.start[i])
    e <- as.numeric(find.diag2$diag.N.term.end[i])
    re[[i]] <- cbind(variant=find.diag2$variant[i],HbA.bc[Ion.num %between% c(s, e)], Mutation = find.diag2$Mutation[i])
  }
  re2 <- do.call(rbind,re)
  #yz ions
  n <- length(find.diag2$diag.C.term.start)
  re3 <- list()
  for (j in 1:n){
    s <- as.numeric(find.diag2$diag.C.term.start[j])
    e <- as.numeric(find.diag2$diag.C.term.end[j])
    re3[[j]] <- cbind(variant=find.diag2$variant[j],HbA.yz[Ion.num %between% c(s, e)],Mutation = find.diag2$Mutation[j])
  }
  re4 <- do.call(rbind,re3)

  re5 <- rbind(re2, re4) %>% filter(!is.na(variant))
  re6 <- re5[with(re5, order(variant))]

  ID2 <- names(Hbvariants)
  out2 <- list()
  for (j in 1:length(ID2)){
    Hbseq1 <- getSequence(Hbvariants, as.string = TRUE)
    Hbseq2 <- sub("M","",Hbseq1[[j]])
    input <- FragmentPeptide(Hbseq2, fragments = "by", IAA = FALSE,N15 = FALSE)
    sub.input <- subset(input, select=c(ms2type, ms2mz))
    sub.input2 <- sub.input %>% separate(ms2type,c("ion", "CS"), "]") %>%
      separate(ion, c("type", "ion_num"), sep = 2) %>%
      separate(type, c("X", "ion_type"), sep = 1)
    sub.input3 <- subset(sub.input2, select=c("ion_type","ion_num","CS","ms2mz" ), CS=="1+")
    input2 <- FragmentPeptide(Hbseq2, fragments = "cz", IAA = FALSE,N15 = FALSE)
    sub.input4 <- subset(input2, select=c(ms2type, ms2mz))
    sub.input5 <- sub.input4 %>% separate(ms2type,c("ion", "CS"), "]") %>%
      separate(ion, c("type", "ion_num"), sep = 2) %>%
      separate(type, c("X", "ion_type"), sep = 1)
    sub.input6 <- subset(sub.input5, select=c("ion_type","ion_num","CS","ms2mz" ), CS=="1+")
    out1 <- rbind(sub.input3,sub.input6)
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
    out[[i]] <-inner_join(sub2, sub, by = c("ion_type"="Ion.type", "ion_num"="Ion.num"))
  }
  relist <- do.call(rbind,out)
  relist2 <- relist %>% separate(variant,c("x","Variant"), sep = 3)
  names(relist2)[4] <- "Ref_Mass"
  relist2$Name <- paste(relist2$ion_type,relist2$Variant, sep = "_" )
  relist3 <- relist2[, c(10,8,1,2,4,7,9)]
}
