# clean .len file
library(stringr)
# read .len file
p_len <- "data/intermediate/popabc/18_loci.len"
len <- readLines(p_len)
# get loci names
loci <-
  grep(x = len, pattern = ".nex$", value = TRUE) |>
  stringr::str_extract("([A-Za-z0-9_]+)(\\.nex$)", group = 1)
# remove commented lines
len_no_comments <- len[!grepl(len, pattern = "^#")]
ss <- grep(x = len_no_comments, pattern = "^s$")
ssx <- sort(c(ss, (ss - 1), (ss - 2)))
len_no_comments2 <- len_no_comments[-ssx]
# create header
header_len <-
  c(paste("#",  date()),
  "#Population 1 - 'kinabalu'",
  "#Population 2 - 'tambuyukon'",
  "#Population 3 - 'trusmadi'",
  paste("#Locus", loci),
  "",
  "3",
  length(loci),
  paste(rep("s", length(loci)), collapse = " "),
  ""
  )

# concat new header
len_head <- c(header_len, len_no_comments2)

# remove 2 consecutive line breaks
lb <- which(len_head == "")
lb_simple_lb <- len_head[-lb[(lb[-1] - lb) == 1]]

# write clean len
p_out <- "data/intermediate/popabc/18_loci-clean.len"
if(file.exists(p_out))
  file.remove(p_out)
writeLines(lb_simple_lb,
           con = p_out)
