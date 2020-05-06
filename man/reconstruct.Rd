\name{reconstruct}
\alias{reconstruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{reconstruct: Missing Data Segments Imputation in Multivariate Streams

}
\description{
Ghost is an accurate imputation algorithm for reconstructing the missing segment in multi-variate data streams. Inspired by single-shot learning, it reconstructs the missing segment by identifying the first similar segment in the stream. Nevertheless, there should be one column of data available, i.e. a constraint column. The values of columns can be characters ( A, B, C,etc.). The result of the imputed dataset will be returned a .csv file.
}
\usage{
reconstruct(data_frame, constraintCol, wSize, direction_save,epsilon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_frame}{
     A data frame with missing values.
}
\item{constraintCol}{
    The column number that all of its fields have data (without missing values). This column is considered as a constraint and it can always produce data, even if the system is shut down.
}
\item{wSize}{
   Length of a window that is used for data reconstruction, before and after the missing row(s).
}
\item{direction_save}{
   A direction for saving the output .CSV file(see details).
}
\item{epsilon}{
   A similarity coefficient that is used for searching in the algorithm (see details).
}
}

\details{
More information about operation of algorithm is prepared in algorithm's article:
https://www.researchgate.net/publication/332779980_Ghost_Imputation_Accurately_Reconstructing_Missing_Data_of_the_Off_Period

Epsilon:
The algorithm searches the data for the closest similar segment. As the first step, the algorithm determines prior and posterior segments missing part (the size of the segment will be given by wSize). As the second step, the algorithm starts to find the similar segment that passes the segment size and constraint similiarty condition. Sometimes, finding windows with exact similarity is impossible in a dataset. To mitigate this issue, and finding windows with approximate similarities the user can define the minimum percentage of similarity for searching the dataset with Epsilon coefficient.

direction_save:
If the user inserts the Direction_save, the output file will be saved in the specified folder. Contrarily, if the user does not insert the Direction_save, the output file will be saved in the Environment R.
}

\references{
Rawassizadeh, Reza, Hamidreza Keshavarz, and Michael Pazzani. "Ghost Imputation: Accurately Reconstructing Missing Data of the Off Period." IEEE Transactions on Knowledge and Data Engineering (to appear).
}
\author{
 Siyavash Shabani,s.shabani.aut@gmail.com
 Reza Rawassizadeh,rrawassizadeh@acm.org

}

\examples{
#An example of the operation of the Algorithm.

data(test_ghost_csv)

## sample dataset----------------------------------
#   S0 S1 S2 S3
#1   5  F  G  H
#2   5  B  N  T
#3   4     P  O
#4   1  X  C  B
#5   1  N     X
#6   1  R  R  R
#7   1  W     W
#8   1  W  W  W
#9   2
#10  2
#11  1  O  K  O
#12  1  B     O
#13  2     S  D
#14  1  W  W
#15  1  W  S  W
#16  2  P  I  M
#17  2  R  U
#18  1  O  K  O
#19  1  B     O
#20  1  R  R  R
#21  5  F  G  H
#22  5  B  N  T
#23  4
#24  1  X  C  B
#25  1  N     X


reconstruct(test_ghost_csv,1,2,epsilon=0.4)

### output---------------------------------------------------
#   S0 S1 S2 S3
#1   5  F  G  H
#2   5  B  N  T
#3   4     P  O
#4   1  X  C  B
#5   1  N     X
#6   1  R  R  R
#7   1  W     W
#8   1  W  W  W
#9   2  P  I  M
#10  2  R  U
#11  1  O  K  O
#12  1  B     O
#13  2     S  D
#14  1  W  W
#15  1  W  S  W
#16  2  P  I  M
#17  2  R  U
#18  1  O  K  O
#19  1  B     O
#20  1  R  R  R
#21  5  F  G  H
#22  5  B  N  T
#23  4     P  O
#24  1  X  C  B
#25  1  N     X

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
