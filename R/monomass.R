monomz <- function (sequence, fragments = "by")
{
  results_list <- vector("list")
  for (sequence_number in 1:length(sequence)) {
    peptide_vector <- strsplit(sequence[sequence_number], split = "")[[1]]
    peptide_length <- length(peptide_vector)
    if (peptide_length < 2)
      stop("sequence must contain two or more residues")
    C <- 12
    H <- 1.007825035
    O <- 15.99491463
    S <- 31.9720707
    P <- 30.973762
    N <- 14.0030740
    proton <- 1.0072764668
    residueMass <- function(residue) {
      if (residue == "A")
        mass = C * 3 + H * 5 + N + O
      if (residue == "R")
        mass = C * 6 + H * 12 + N * 4 + O
      if (residue == "N")
        mass = C * 4 + H * 6 + N * 2 + O * 2
      if (residue == "D")
        mass = C * 4 + H * 5 + N + O * 3
      if (residue == "E")
        mass = C * 5 + H * 7 + N + O * 3
      if (residue == "Q")
        mass = C * 5 + H * 8 + N * 2 + O * 2
      if (residue == "G")
        mass = C * 2 + H * 3 + N + O
      if (residue == "H")
        mass = C * 6 + H * 7 + N * 3 + O
      if (residue == "I")
        mass = C * 6 + H * 11 + N + O
      if (residue == "L")
        mass = C * 6 + H * 11 + N + O
      if (residue == "K")
        mass = C * 6 + H * 12 + N * 2 + O
      if (residue == "M")
        mass = C * 5 + H * 9 + N + O + S
      if (residue == "F")
        mass = C * 9 + H * 9 + N + O
      if (residue == "P")
        mass = C * 5 + H * 7 + N + O
      if (residue == "S")
        mass = C * 3 + H * 5 + N + O * 2
      if (residue == "T")
        mass = C * 4 + H * 7 + N + O * 2
      if (residue == "W")
        mass = C * 11 + H * 10 + N * 2 + O
      if (residue == "Y")
        mass = C * 9 + H * 9 + N + O * 2
      if (residue == "V")
        mass = C * 5 + H * 9 + N + O
      if (residue == "C")
        mass = C * 3 + H * 5 + N + O + S

      return(mass)
    }
    masses <- sapply(peptide_vector, residueMass)
    pm <- sum(masses)
    p1 <- round(pm + H * 2 + O + proton, digits = 5)
    if (fragments == "by") {
      b1 <- vector(mode = "numeric", length = 0)
      bi <- vector(mode = "integer", length = 0)
      y1 <- vector(mode = "numeric", length = 0)
      yi <- vector(mode = "integer", length = 0)
      for (i in 1:(peptide_length - 1)) {
        mass <- sum(masses[1:i])
        b1[i] <- round(mass + proton, digits = 5)
        bi[i] <- i
      }
      for (j in 2:peptide_length) {
        mass <- sum(masses[j:peptide_length])
        y1[j - 1] <- round(mass + H * 2 + O + proton, digits = 5)
        yi[j - 1] <- peptide_length - j + 1
      }
      ms1z1 <- rep(p1, times = (length(bi) + length(yi)))
      b1.type <- paste("b", bi, sep = "")
      y1.type <- paste("y", yi, sep = "")
      ms2type <- c(b1.type, y1.type)
      MH <- c(b1,y1)
    }
    if (fragments == "cz") {
      c1 <- vector(mode = "numeric", length = 0)
      ci <- vector(mode = "integer", length = 0)
      z1 <- vector(mode = "numeric", length = 0)
      zi <- vector(mode = "integer", length = 0)
      for (i in 1:(peptide_length - 1)) {
        mass <- sum(masses[1:i])
        c1[i] <- round(mass + 3 * H + N + proton, digits = 5)
        ci[i] <- i
      }
      for (j in 2:peptide_length) {
        mass <- sum(masses[j:peptide_length])
        z1[j - 1] <- round(mass + O - N + proton, digits = 5)
        zi[j - 1] <- peptide_length - j + 1
      }
      ms1z1 <- rep(p1, times = (length(ci) +length(zi)))
      c1.type <- paste("c", ci, sep = "")
      z1.type <- paste("z", zi, sep = "")
      ms2type <- c(c1.type, z1.type)
      MH <- c(c1, z1)
    }
    results_list[[sequence_number]] <- data.frame(ms1z1,Ion = ms2type,ms2type, MH)%>% separate(ms2type, c("Ion_type", "Ion_num"), sep = 1)
  }
  return(as.data.frame(do.call("rbind", results_list)))
}
