##
## Repair path on other platforms
##
attr(xs, "filepaths") <- file.path(.find.package("mzmatch.Data"),
                                 gsub("/", .Platform$file.sep,
                                      attr(xs, "filepaths")))
