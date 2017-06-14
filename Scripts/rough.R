## Packages needed
library(sdsBAT)
library(signal)
library(tidyverse)

## Data
raw <- sdsBAT::calls[ , -c(1, 2)]
load("~/R/Data/clean_aligned_details.Rdata")
## Raw Spectrograms
x <- raw$chirp1[[1]]
Fs <- 500000L # recording sampling rate
fftn <- 2**8
step <- fftn/8

spg <- specgram(x, n = fftn, Fs = Fs, overlap = fftn - step)
mag <- abs(spg$S)
image(t(20*log10(mag)), axes = FALSE, col = gray(0:255 / 255))

get_spectrogram <- function(x, n = min(256, length(x)), Fs = 2, window = hanning(n),
  overlap = ceiling(length(window)/2)){
  if(!is.na(x)){spg <- specgram(x[[1]], n, Fs, window, overlap)}else(spg <- NA)
  return(spg)
}

## All spectrograms
spg <- raw %>%
  apply(c(1,2), get_spectrogram, n = fftn, Fs = Fs, overlap = fftn - step) %>%
  as.data.frame()

## Isolate the spectrogram matrix
isolate_var <- function(x, var = 1){
  if(!is.na(x)){S <- x[[1]][[var]]}else(S <- NA)
  return(S)
}

S <- spg %>%
  apply( c(1,2), isolate_var) %>%
  apply(c(1,2), function(x) abs(x[[1]])) %>%
  apply(c(1,2), function(x) x[[1]] / max(x[[1]])) %>%
  apply(c(1,2), function(x) log10(x[[1]]) ) %>%
  as.data.frame()

image(t(S$chirp1[[1]]), axes = F)
spg[1,1]

test <- smooth.spline(spg$chirp1[[1]])

image(t(A[[1]]), axes = F)

f <- seq(0, 250, length.out = nrow(A[[1]]))
f_keep <- f[f>9 & f<212]

spectrograms <- array(dim = c(dim(A[[1]]), length(A)))
for(i in seq.int(length(A))){
  spectrograms[,,i] <- A[[i]]
}
dim(spectrograms)

spectrograms <- spectrograms[f_keep,,]

image(t(spectrograms[,,4]))

vec_spec <- apply(spectrograms, 3, function(x) c(x))

PCA <- prcomp(t(vec_spec))
var <- PCA$sdev**2 / sum(PCA$sdev**2)

plot(cumsum(var[1:10]))

components <- array(dim = c(dim(spectrograms[,,1]), 10))
for(i in seq.int(10)){
  components[,,i] <- PCA$rotation[,i]
}
image(components[,,9])

sp <- apply(vec_spec, 1, function(x) tapply(x, D$V2, mean))

mod <- lm(sp[1,] ~ PCA$rotation[,1:5])
length(sp[1,])
dim(PCA$rotation[,1:5])

components <- PCA$rotation[,1:5] 
matplot(components, type = 'l')

plot(PCA$sdev[1]*components[,1], type = 'l')
components <- t(apply(components, 1, '*', PCA$sdev[1:5]))
matplot(components, type = 'l')

get_weights <- function(response, predictors){
  mod <- lm(response~predictors)
  return(mod$coefficients)
}

weights <- apply(sp, 1, get_weights, predictors = components)

hyp <- apply(weights, 1, sdsBAT:::pou_type2mle, phylogeny, logl_function = pou_logl_fast, 
  optim_function = 'uobyqa', upper_initialisation = c(1,1,1), n_restarts = 5)
hyp

internal_weights <- predictive_distribution(weights, hyp)

components <- cbind(rep(1, nrow(components)), components)

map_est <- internal_weights$predictive_mean %*% t(components) 

map_est <- apply(map_est, 1, '+', PCA$center)

map_spec <- array(dim = c(208, 50, 21))
for(i in 1:21){
  map_spec[,,i] <- map_est[,i]
}

image(t(map_spec[,,1]), axes =  F)

anc <- map_spec[,,1]

sum(f <= 9)
sum(f >= 212)
anc <- rbind(matrix(0, ncol = 50, nrow = sum(f <= 9)),
  anc,
  matrix(0, ncol = 50, nrow = sum(f >= 212)))
write.csv(anc, 'ancestral_spectrogram.csv')
