# Optimize clustering by finding the best k
OptimizeClustering <- function(mtrx, k0=3, k1=round(sqrt(ncol(mtrx))), method=c('hierarchical', 'kmeans', 'som'), 
                               evaluation=c('silhouette', 'dunn', 'davies')) {
  
  method <- tolower(method[1]);
  evaluation <- tolower(evaluation[1]);
  
  if (method == 'hierarchical') {
    hc <- hclust(as.dist(1-cor(mtrx))); 
    dt <- dist(t(mtrx)); 
    cl <- lapply(k0:k1, function(k) cutree(hc, k=k));
  } else if (method == 'kmeans') {
    cl <- lapply(k0:k1, function(k) kmeans(t(mtrx), centers = k)$cluster); 
  } else if (method == 'som') {
    require(som);
    cl <- lapply(k0:k1, function(k) { 
      k <- as.numeric(k); 
      sm <- som(t(mtrx), 1, k); 
      cl <- paste(sm$visual$x, sm$visual$y, sep='_'); 
      cl <- as.integer(factor(cl)); 
      names(cl) <- colnames(mtrx);
      cl;
    });
  } else stop('Error: unknown clustering method', method);
  names(cl) <- k0:k1;
  
  if (evaluation == 'silhouette') {
    require(cluster);
    dt <- dist(t(mtrx)); 
    sc <- lapply(cl, function(c) silhouette(c, dt));
    sl <- sapply(sc, function(s) mean(s[, 3]));
    ii <- which(sl==max(sl))[1]; 
    
    sl <- sc[[ii]]; 
    sl <- cbind(sl[, 1], sl[, 2], sl[, 3]);
    colnames(sl) <- c('cluster', 'neighbor', 'width'); 
  } else if (evaluation == 'dunn') {
    require(clValid);
    dt <- dist(t(mtrx));
    dn <- sapply(cl, function(cl) dunn(dt, clusters=cl));
    ii <- which(dn==max(dn)); 
    
    sl <- silhouette(cl[[ii]], dt); 
    sl <- cbind(sl[, 1], sl[, 2], sl[, 3]);
    colnames(sl) <- c('cluster', 'neighbor', 'width'); 
  } else if (evaluation == 'davies') {
    require(clusterSim);
    db <- sapply(cl, function(cl) index.DB(t(mtrx), cl)$DB);
    ii <- which(db==min(db)); 
    
    sl <- silhouette(cl[[ii]], dist(t(mtrx))); 
    sl <- cbind(sl[, 1], sl[, 2], sl[, 3]);
    colnames(sl) <- c('cluster', 'neighbor', 'width'); 
  };
  
  n <- sapply(split(sl[, 1], sl[, 1]), length);
  b <- sapply(split(sl[, 2], sl[, 1]), mean);
  w <- sapply(split(sl[, 3], sl[, 1]), mean);
  s <- cbind(N=n, Mean_Neighbor=b, Mean_Width=w); 
  
  list(parameter=c(min=k0, max=k1, optimal=(k0:k1)[ii], method=method, evaluation=evaluation), cluster=cl[[ii]],
       summary=s, silhouette=sl, all=cl); 
}