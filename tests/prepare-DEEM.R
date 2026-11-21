library(MR.DEEM)

summary_data <- data.frame(
  chr = 1L, pos = c(110L, 210L, 310L, 410L),
  gamma_s = c(.20, .12, .01, .21), sigma_s = rep(.05, 4),
  gamma_h = c(.18, .11, .02, .19), sigma_x = c(.04, .05, .06, .055),
  Gamma_h = c(.09, .07, .015, .095), sigma_y = c(.05, .06, .07, .065)
)
blocks <- data.frame(chr = rep("chr1", 2), start = c(100, 400), stop = c(400, 500))
reference <- matrix(c(1, .20, .10, .20, 1, .25, .10, .25, 1), 3L, 3L)

prepare <- function(dat, blocks, reference) {
  prepare.DEEM(dat, paste(dat$chr, dat$pos, sep = ":"),
               list(reference), blocks, ns = 1000, ne = 1000, no = 1000)
}

many <- prepare(summary_data[1:3, ], blocks[1, ], list(reference))
stopifnot(many[[1L]] == 1, identical(lengths(many[[4L]]), 2L),
          identical(dim(many[[5L]][[1L]]), c(2L, 2L)))

mixed <- prepare(summary_data, blocks, list(reference, matrix(1, 1L, 1L)))
stopifnot(mixed[[1L]] == 2, identical(lengths(mixed[[4L]]), c(2L, 1L)),
          identical(dim(mixed[[5L]][[2L]]), c(1L, 1L)),
          isTRUE(all.equal(as.numeric(mixed[[5L]][[2L]]), .065^2 + .095^2 / 998)))

filtered <- summary_data[1:3, ]
filtered$gamma_s[2L] <- .01
reduced <- prepare(filtered, blocks[1, ], list(reference))
stopifnot(reduced[[1L]] == 1, identical(lengths(reduced[[4L]]), 1L),
          identical(dim(reduced[[5L]][[1L]]), c(1L, 1L)))
