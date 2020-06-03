context("Testing Omixer Functions")

testingData <- sampleList <- 
    tibble(sampleId=str_pad(1:48, 4, pad="0"),
    sex=as_factor(sample(c("m", "f"), 48, replace=TRUE)), 
    age=round(rnorm(48, mean=30, sd=8), 0), 
    smoke=as_factor(sample(c("yes", "ex", "never"), 48, replace=TRUE)),
    date=sample(seq(as.Date('2008/01/01'), as.Date('2016/01/01'), 
        by="day"), 48))

randVars <- c("sex", "age", "smoke", "date")

compareDataRand <- omixerRand(testingData, iterNum=10, wells=48, div="row", 
    randVars=randVars)

compareDataSpec <- omixerSpecific(testingData, seed=compareDataRand$seed[1], wells=48, div="row", 
    randVars=randVars)

test_that("All observations are kept after randomization.", {
  expect_true(all(table(compareDataRand$smoke) == table(testingData$smoke)))
})

test_that("omixerSpecific and omixerRand match with matching seeds.", {
  expect_true(all(compareDataSpec == compareDataRand))
})



