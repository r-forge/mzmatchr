"unsafeFindInterval3bis" <- function (x, vec, rightmost.closed = FALSE, all.inside = FALSE, left.open = FALSE) {
            .Internal(findInterval(vec=vec, x=x,
                                   rightmost.closed=rightmost.closed,
                                   all.inside=all.inside, left.open = left.open))

    }