# make draft to post & build

library(fs)
library(blogdown)

fs::file_move("content/post/",
              "content/english/post/")
blogdown::build_site()