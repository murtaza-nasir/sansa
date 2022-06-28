#' Title
#'
#' @param x Input predictor as a dataframe
#' @param y Target variable as factor
#' @param lambda Lambda parameter to select distribution of synthetic variables
#' @param ksel K parameter to choose how many neighbors are used in calculations
#'
#' @return A list with two elements: x contains predictors with synthetic data, y contains target with synthetic data.
#' @export
#'
#' @examples
#'
#' library(sansa)
#' library(ggplot2)
#' minority = data.frame(x1 = rnorm(10, 10, 3),
#'                       x2 = rnorm(10, 25, 10),
#'                       target = "true")
#' majority = data.frame(x1 = rnorm(100, 4, 2),
#'                       x2 = rnorm(100, 30, 10),
#'                       target = "false")
#'
#' dataset = rbind(minority, majority)
#'
#' ggplot(dataset) + geom_point(aes(x1, x2, color = target))
#' sansaobject = sansa(x = dataset[,1:2], y = dataset$target, lambda = 1, ksel = 3)
#'
#' balanced <- sansaobject$x
#' balanced$target = sansaobject$y
#'
#' ggplot(balanced) + geom_point(aes(x1, x2, color = target))
#'
#' @importFrom data.table .SD
#' @import ggplot2
sansa <-function (x, y, lambda = 0, ksel = 3) {

  set.seed(1001)
  dat <- if (is.data.frame(x)) x else as.data.frame(x)
  duplicates <- duplicated(dat)
  dat <- dat[!duplicates,]
  y <- y[!duplicates]
  maj <- max(table(y))
  min <- min(table(y))
  dif <- maj - min

  if (.75 * min < ksel) ksel = floor(.75 * min)

  #nclust <- round(sqrt(min))

  mf1 = lambda #lambda

  constant = 30

  dat$y <- y
  minclass <- names(table(y)[which.min(table(y))])
  omins <- dat[dat$y==minclass,]
  omaxs <- dat[dat$y!=minclass,]
  ominorig <- omins
  cols <- ncol(x)

  ominorig$dists <- rowSums(FNN::knn.dist(ominorig[,-ncol(ominorig)], k = ksel)) #find euclidean distance of 10 nearest neighbors and sum it. add as column
  ominorig <- ominorig[order(ominorig$dists, decreasing = T),] #sort descending by distance

  dists <- FNN::knn.index(ominorig[,1:cols], k = ksel) #get indexes of 10 nearest neighbors for each row
  dimdist <- ominorig[0,1:cols] #create empty table

  for(i in 1:nrow(ominorig)) #loop over all rows for minority class
  {
    dimdist[i,1:cols] <- ominorig[i,1:cols] #copy row from original
    for(j in 1:ksel)
    {
      dimdist[i,1:cols] <- dimdist[i,1:cols] + ominorig[i,1:cols]- ominorig[dists[i,j],1:cols] #sum the difference between dimensions for each row and its 10 closest neighbors
    }
    dimdist[i,1:cols] <- dimdist[i,1:cols] / constant #average the above. Basically find the average difference in each dimension for 10 nearest neighbors
  }

  replist <- (ominorig$dists^mf1) #distort distances based on parameter to create replication index
  replist <- round(replist / sum(replist) * dif) #round replication index

  minexp <- ominorig[rep(seq_len(nrow(ominorig)), c(replist)),1:ncol(x)] #replicate original rows based on replication index
  distexp <- dimdist[rep(seq_len(nrow(dimdist)), c(replist)),1:ncol(x)] #replicate average distance rows based on replication index

  comvar <- as.matrix(data.table::setDT(dat)[, lapply(.SD, function(x) max(x) - min(x)), by= y][,-"y"])
  comvar <- data.table::as.data.table(sweep(comvar, 2, comvar[1,],"/"))

  distexp <- distexp*matrix(stats::rnorm(nrow(distexp)*ncol(x), mean = 0, sd = 1), ncol = ncol(x))

  minexp <- minexp + distexp #mean(sapply(x, sd))*.05

  minexp$y <- minclass

  dat <- rbind(omins, minexp, omaxs[,1:(ncol(x)+1)])
  dat <- dat[sample(nrow(dat),nrow(dat)),]

  colu <- -ncol(dat)

  list(x = dat[, colu],
       y = dat$y)
}


# Use like this:
# recursk(x, y, lambda = lambda, ksel = k) <- x is a dataframe with predictors, y is a vector of the target variable
