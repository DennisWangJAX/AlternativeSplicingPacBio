populateMatrix = function(query, subject){
  info = findOverlaps(query, subject)
  ind = 1:length(subject)
  return((ind %in% subjectHits(info)))
}
