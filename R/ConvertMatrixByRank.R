ConvertMatrixByRank <- function(mtrx, rank, ind1, ind2) {
  rank[is.na(rank)] <- 0;

  if (length(unique(sign(rank)))==1)
    rank <- sign(rowMeans(mtrx[, ind2, drop=FALSE])-rowMeans(mtrx[, ind1, drop=FALSE]))*rank;

  qn <- abs(qnorm(((1:length(rank))-0.5)/(2*length(rank))));
  rk <- sign(rank)*qn[(length(rank)+1)-abs(rank(abs(rank)))];
  names(rk) <- names(rank);

  mx <- mtrx;
  mx[, ind1] <- apply(mtrx[, ind1, drop=FALSE], 2, function(x) x-rowMeans(mtrx[, ind1, drop=FALSE]));
  mx[, ind2] <- apply(mtrx[, ind2, drop=FALSE], 2, function(x) x-rowMeans(mtrx[, ind2, drop=FALSE]));

  mx;
}
