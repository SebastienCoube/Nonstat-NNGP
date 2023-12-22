complete = seq(3840) [sapply(sapply(seq(3840), function(i)paste("res", i, "complete.RDS", sep = "")), function(x)x%in%list.files())]
saveRDS(list(started = complete, complete = complete), "started_complete.RDS")
