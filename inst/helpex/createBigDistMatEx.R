library(SSN)
# NOT RUN
# mf04 <- importSSN(system.file("lsndata/MiddleFork04.ssn",
#	        package = "SSN"), o.write = TRUE)
# use SpatialStreamNetwork object mf04 that was already created
data(mf04)
# for examples, copy MiddleFork04.ssn directory to R's temporary directory
copyLSN2temp()
#make sure mf04p has the correct path, will vary for each users installation
mf04 <- updatePath(mf04, paste0(tempdir(),'/MiddleFork04.ssn'))

# create distance matrix among observed data points
createBigDistMat(mf04, o.write = TRUE, no.cores = 1)

# create distance matrix among observed data points
#     and between observed and prediction points
data(mf04p)
mf04p <- updatePath(mf04p, paste0(tempdir(),'/MiddleFork04.ssn'))
createBigDistMat(mf04p, predpts = "pred1km", o.write = TRUE, no.cores = 1)

# NOT RUN include prediction to prediction site distances
      # createBigDistMat(mf04p, predpts = "pred1km", o.write = TRUE,
      #     amongpreds = TRUE, no.cores = 4)
