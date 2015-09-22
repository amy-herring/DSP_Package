
context("'Mcmc-UpdateNaSex.R' fcns")


# Test 'tfMat' -----------------------------------------------------------------

tfMat <- list( nIsOne = matrix(c(FALSE, 
                                 TRUE), nrow=2),
               nIsTwo = matrix(c(FALSE, FALSE,
                                 TRUE,  FALSE,
                                 FALSE, TRUE,
                                 TRUE,  TRUE), byrow=TRUE, nrow=4) )

test_that("'getTfPerms' tests", {
  expect_identical(getTfPerms(1), tfMat$nIsOne)
  expect_identical(getTfPerms(2), tfMat$nIsTwo)
})




# Test 'getPowSet' -------------------------------------------------------------

out0 <- list(numeric(0))
out1 <- list(numeric(0), 1)
out2 <- list(numeric(0), 1, 2, c(1, 2))
out3 <- list(numeric(0), 1, 2, c(1, 2), 3, c(1, 3), c(2, 3), c(1, 2, 3))

test_that("'getPowSet' tests", {
  expect_identical(getPowSet(numeric(0)), out0)
  expect_identical(getPowSet(1), out1)
  expect_identical(getPowSet(c(1,2)), out2)
  expect_identical(getPowSet(c(1,2,3)), out3)
})




# Test 'getCycSexPartit' -------------------------------------------------------

cycIdx <- list(c(1, 2, 3), 4, c(5, 6), c(7, 8, 9))
# missing: 1, 2, 4, 9
naSexBool <- c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE)

cycSexPartit <- getCycSexPartit(cycIdx, naSexBool)
out <- list( known = list(3, numeric(0), c(5, 6), c(7, 8)),
             miss  = list(c(1, 2), 4, numeric(0), 9),
             perms = list( list(numeric(0), 1, 2, c(1, 2)),
                           list(numeric(0), 4),
                           list(numeric(0)),
                           list(numeric(0), 9) ) )

test_that("'getCycSexPartit' tests", {
  expect_identical(cycSexPartit$known, out$known)
  expect_identical(cycSexPartit$miss, out$miss)
  expect_identical(cycSexPartit$perms, out$perms)
})




# Test 'getCycPermsIdx' ----------------------------------------------------------

cycPermsIdx <- getCycPermsIdx(cycIdx, naSexBool)
out <- list( list( 3, c(1, 3), c(2, 3), c(1, 2, 3) ),
             list( numeric(0), 4 ),
             list( c(5, 6) ),
             list( c(7, 8), c(7, 8, 9) ) )

test_that("'getCycPermsIdx' tests", {
  expect_identical(cycPermsIdx, out)
})




# Test 'getCycSexPriors' -------------------------------------------------------

cycSexPriors <- getCycSexPriors(cycIdx, which(naSexBool), c(0.2, 0.3, 0.6, 0.25))
out <- list(c(0.2, 0.3), 0.6, numeric(0), 0.25)

test_that("'getCycSexPriors' tests", {
  expect_identical(cycSexPriors, out)
})





# Test 'getSexPermsLik' --------------------------------------------------------


sexPriorLik <- getSexPriorLik(cycSexPriors)
out <- list( c((1 - 0.2) * (1 - 0.3), 0.2 * (1 - 0.3), (1 - 0.2) * 0.3, 0.2 * 0.3),
             c(1 - 0.6, 0.6),
             c(1),
             c(1 - 0.25, 0.25) )

test_that("'getSexPriorLik' tests", {
  expect_identical(sexPriorLik, out)
})




# Test 'getPostSexProbs' -------------------------------------------------------

pregCycBool <- c(TRUE, FALSE, TRUE, FALSE)
uBeta <- c(-1, 0, 1, -1, 0, 1, -1, 0, 1)
xiDay <- rep(1, length(uBeta))

postSexProbs <- getPostSexProbs(pregCycBool, uBeta, xiDay, cycPermsIdx, sexPriorLik)
out <- list( (1 - exp( -c(exp(1), 
                          exp(-1) + exp(1), 
                          exp(0) + exp(1), 
                          exp(-1) + exp(0) + exp(1)) )) * sexPriorLik[[1]],
             c(1, exp(-exp(-1))) * sexPriorLik[[2]],
             1,
             exp( -c(exp(-1) + exp(0),
                     exp(-1) + exp(0) + exp(1)) ) * sexPriorLik[[4]] )

# If a pregnancy occurs then posterior probability of no sex is 0
noSexForPreg <- getPostSexProbs(TRUE, uBeta, xiDay, 
                                list(list(numeric(0), 1)), list(c(0.75, 0.25)))

test_that("'getPostSexProbs' tests", {
  expect_identical(postSexProbs, out)
  expect_identical(noSexForPreg[1], 0)
})




# Test 'sampNaSex' -------------------------------------------------------

updateCycIdx <- sampCycIdx(pregCycBool, uBeta, xiDay, cycPermsIdx, sexPriorLik)
checkInList <- function(x, theList) (TRUE %in% sapply(theList, function(y) identical(x, y)))

test_that("'sampNaSex' tests", {
  expect_identical(checkInList(updateCycIdx[[1]], cycPermsIdx[[1]]), TRUE)
  expect_identical(checkInList(updateCycIdx[[2]], cycPermsIdx[[2]]), TRUE)
  expect_identical(checkInList(updateCycIdx[[3]], cycPermsIdx[[3]]), TRUE)
  expect_identical(checkInList(updateCycIdx[[4]], cycPermsIdx[[4]]), TRUE)
})




# Test 'getIdx' ----------------------------------------------------------------

# x <- list(1:3, 4, 5:6)
# sexBool <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
# out <- list(c(1L, 3L), 5L)
# 
# test_that("'getIdx' tests", {
#   expect_identical(getIdx(x, sexBool), out)
# })

                  



