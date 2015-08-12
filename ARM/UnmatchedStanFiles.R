# Find .stan files without R sample file.
# This uses a very simple search, and should be used a rough guide only
library(stringr)
stanFiles = data.frame(
  fullName = dir(pattern = "*.stan", recursive = TRUE, full.names = TRUE),
  stringsAsFactors = FALSE )
stanFiles$shortName = str_sub(basename(stanFiles$fullName),1,-6)
stanFiles = stanFiles[order(stanFiles$shortName),]

rFiles = dir(pattern = "*.R", recursive = TRUE, full.names = TRUE)

stanInRFiles = sapply(rFiles, function(file){
  f = readLines(file)
  mat = na.omit(str_match(f, 'stan\\(file *= *\\\'(.*)\\.stan')[,2])
  mat[!is.na(mat)]
})
stanInRFiles = sort(unique(unlist(stanInRFiles)))

# Stan files without sample
notThere = stanFiles[stanFiles$shortName %in% setdiff(stanFiles$shortName, stanInRFiles),]

cat(sort(notThere$fullName), sep="\n")
