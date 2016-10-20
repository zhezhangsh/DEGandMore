DeGPseq <- function(mtrx, grps, paired=FALSE) {
  require(DEGandMore);
  require(GPseq); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by GPSeq; performing unpaired test instead.\n");
  
  u1 <- summary(rowMeans(mtrx[rowSums(mtrx)>0, grps[[1]]]))[[5]];
  u2 <- summary(rowMeans(mtrx[rowSums(mtrx)>0, grps[[2]]]))[[5]];
  w  <- as.vector(u2)/as.vector(u1); 
  
  chi <- apply(mtrx, 1, function(d) {
    dx <- d[grps[[1]]];
    dy <- d[grps[[2]]];
    ox <- generalized.poisson.likelihood(dx);
    oy <- generalized.poisson.likelihood(dy);
    likelihood_ratio_tissue_generalized_poisson(dx, ox$lambda, ox$theta, dy, oy$lambda, oy$theta, w=w)[[2]]; 
  }); 
  
  pv <- pchisq(chi, 1, lower.tail = FALSE, log.p = TRUE); 
  pv <- exp(pv); 
  qv <- p.adjust(pv, method='BH');

  m1 <- rowMeans(mtrx[, grps[[1]]]); 
  m2 <- rowMeans(mtrx[, grps[[2]]])/w;
  l2 <- log2(pmax(m2, min(m2[m2>0])/2)) - log2(pmax(m1, min(m1[m1>0])/2)); 
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps);
}


generalized.poisson.likelihood <- function (y) {
  n = length(y)
  y_bar = mean(y)
  var_y = var(y)
  ly = max(y)
  if (var_y == 0 || ly < 2) {
    return(list(mark = 0, theta = -1, lambda = -1, y_bar = y_bar, 
                length = n))
  } else {
    l = 1 - sqrt(y_bar/var_y)
    obs = rep(0, ly + 1)
    for (i in 1:n) {
      obs[y[i] + 1] = obs[y[i] + 1] + 1
    }
  }
  mark = 0
  flag = 1
  for (its in 1:10000) {
    hl_i = rep(0, ly + 1)
    il_i = rep(0, ly + 1)
    tmp = rep(0, ly + 1)
    tmp = y_bar + ((1:ly) - y_bar) * l
    if (sum(tmp[2:ly] == 0) == 0) {
      hl_i = (1:(ly - 1)) / tmp[2:ly] * (2:ly) * obs[3:(ly + 1)]; 
      il_i = (1:(ly - 1)) / (tmp[2:ly]^2) * (2:ly) * obs[3:(ly + 1)] * ((2:ly) - y_bar); 
    } else {
      flag = 0
      break
    }
    hl = sum(hl_i)
    il = sum(il_i)
    hl = hl - n * y_bar
    if (flag == 1 && il != 0) {
      delta_l = hl/il
      z1 = max(abs(l), 0.01)
      z2 = abs(delta_l)/z1
      l = l + delta_l
      if (z2 <= 1e-06 && l < 1) {
        mark = 1
        break
      }
    } else {
      break
    }
  }
  theta = y_bar * (1 - l)
  if (l < 0) {
    tm = floor(-1 * theta/l)
    if (tm < 4) {
      mark = 0
    }
    z1 = -1 * theta/tm
    if (z1 < -1 && l < -1) {
      mark = 0
    }
  }
  output = c(mark, theta, l, y_bar)
  return(list(mark = mark, theta = theta, lambda = l, y_bar = y_bar, length = n))
}