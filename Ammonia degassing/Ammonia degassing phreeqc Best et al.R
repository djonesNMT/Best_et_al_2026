#setwd()

library(phreeqc)

phrLoadDatabaseString(phreeqc.dat)
phrSetOutputStringsOn(TRUE)

phrLoadDatabaseString(minteq.v4.dat)

################
#ok, there's a better way but just creating individual files for now...

# PC21-3a
PC21_3a <- readLines("PC21_3a_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(PC21_3a)
PC1 <- phrGetOutputStrings()
PC1

PC21_3a_v2 <- readLines("PC21_3a_NH3_option2sulfide.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(PC21_3a_v2)
PC1v2 <- phrGetOutputStrings()
PC1v2

# PC21-3b
PC21_3b <- readLines("PC21_3b_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(PC21_3b)
PC2 <- phrGetOutputStrings()
PC2

# PC22-3a
PC22_3a <- readLines("PC22_3a_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(PC22_3a)
PC3 <- phrGetOutputStrings()
PC3

# PC22-3b
PC22_3b <- readLines("PC22_3b_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(PC22_3b)
PC4 <- phrGetOutputStrings()
PC4

# PC22-6
PC22_6 <- readLines("PC22_6_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(PC22_6)
PC5 <- phrGetOutputStrings()
PC5

# PC22-10
PC22_10 <- readLines("PC22_10_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(PC22_10)
PC6 <- phrGetOutputStrings()
PC6

# PC21-10a
PC21_10a <- readLines("PC21_10a_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(PC21_10a)
PC7 <- phrGetOutputStrings()
PC7

# PC21-10b
PC21_10b <- readLines("PC21_10b_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(PC21_10b)
PC8 <- phrGetOutputStrings()
PC8

# GB21-3
GB21_3 <- readLines("GB21_3_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(GB21_3)
GB1 <- phrGetOutputStrings()
GB1

# GB22-3
GB22_3 <- readLines("GB22_3_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(GB22_3)
GB2 <- phrGetOutputStrings()
GB2

# RS22-2
RS22_2 <- readLines("RS22_2_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(RS22_2)
RS1 <- phrGetOutputStrings()
RS1

# RS21-3A
RS21_3A <- readLines("RS21_3a_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(RS21_3A)
RS2 <- phrGetOutputStrings()
RS2

# RS21-3B
RS21_3B <- readLines("RS21_3b_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(RS21_3B)
RS3 <- phrGetOutputStrings()
RS3

# RS21-3C
RS21_3C <- readLines("RS21_3c_NH3.txt")
phrLoadDatabaseString(phreeqc.dat)
phrRunString(RS21_3C)
RS4 <- phrGetOutputStrings()
RS4


