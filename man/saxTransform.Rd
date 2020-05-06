\name{saxTransform}
\alias{saxTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{saxTransform
}
\description{
This function is added to the package to enable users converting numeric data to discrete data. This is due to the fact that Ghost designed for discrete data and this function discretize numeric data and prepare them for the ghost algorithm.
}
\usage{
saxTransform(data_frame, buckets,skipColumnVec,constraint_row)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_frame}{
   A data frame with numeric values.
}
  \item{buckets}{
  The Input data range is divided to this number.
}
 \item{skipColumnVec}{
  Column number that is not used in the algorithm.
}
 \item{constraint_row}{
  Column number that is considered for constant column.
}
}

\references{
1- Rawassizadeh, Reza, Hamidreza Keshavarz, and Michael Pazzani. "Ghost Imputation: Accurately Reconstructing Missing Data of the Off Period." IEEE Transactions on Knowledge and Data Engineering (to appear).

2- Lin, J., Keogh, E., Lonardi, S., & Chiu, B. (2003). A symbolic representation of time series, with implications for streaming algorithms. In Proceedings of the 8th ACM SIGMOD workshop on Research issues in data mining and knowledge discovery (pp. 2-11). ACM.
}
\author{
Siyavash Shabani, s.shabani.aut@gmail.com, Reza Rawassizadeh, rrawassizadeh@acm.org
}


\examples{

data(sax_test)

#### Input dataframe-----------------------
#   S0  S1  S2 S3
#1   1   2  54 65
#2   1  NA  21 54
#3   2  34  32 87
#4   1  23  58 52
#5   1  43  75 56
#6   2  12  20 95
#7   1  54  14 87
#8   3  -6  NA 30
#9   2   5 -60 32
#10  1 -85  58 25
#11  2  78  95 45
#12  3  52  52 62
#13  2  20  NA 58
#14  3  NA -62 78
#15  1  20 -10 96
#16  1  30  -6 NA
#17  1  12 -85 45
#18  1  NA  78 20
#19  1  23  95 NA

saxTransform(sax_test,buckets =10,skipColumnVec=1, constraint_row=1)

### Output data----------------------------------------------
#     S0  S1  S2   S3
# [1,] "1" "2" "54" "65"
# [2,] "1" ""  "f"  "h"
# [3,] "2" "g" "g"  "j"
# [4,] "1" "g" "h"  "h"
# [5,] "1" "h" "i"  "h"
# [6,] "2" "f" "f"  "k"
# [7,] "1" "h" "f"  "j"
# [8,] "3" "e" ""   "g"
# [9,] "2" "f" "b"  "g"
#[10,] "1" "a" "h"  "g"
#[11,] "2" "j" "k"  "h"
#[12,] "3" "h" "h"  "i"
#[13,] "2" "f" ""   "h"
#[14,] "3" ""  "b"  "j"
#[15,] "1" "f" "e"  "k"
#[16,] "1" "g" "e"  ""
#[17,] "1" "f" "a"  "h"
#[18,] "1" ""  "j"  "f"
#[19,] "1" "g" "k"  ""
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
