suppressPackageStartupMessages(library(SingleR))
args <- commandArgs(trailingOnly=TRUE)
capture.output({
  print(formals(trainSingleR)); print(formals(classifySingleR));
  print(formals(pruneScores)); print(formals(aggregateReference));
  cat('\nPRUNING DOCUMENTATION\n')
  cat(tools::Rd2txt(utils:::.getHelpFile(help('pruneScores', package='SingleR'))))
}, file=args[[1]])
