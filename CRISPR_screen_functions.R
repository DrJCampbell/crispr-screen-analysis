#
#
#


remove_controls <- function(x, pos="PLK1", neg="Olfr", id_colname="shrna.id"){
	ids <- x[,id_colname]
	pos_rows <-grep(pos,ids)
	neg_rows <-grep(neg,ids)
	pos_neg_rows <- c(pos_rows, neg_rows)
	return(x[-pos_neg_rows,])
}


get_pptm <- function(x){
	x_total_count_sum <- sum(
		x$total.hits
		)
	pptm <- x$total.hits / (x_total_count_sum / 10^7)
	pptm_psuedo <- make_pseudo_counts(pptm)
	x_with_pptm <- cbind(
		x,
		pptm_psuedo
		)
	return(x_with_pptm)
}

make_pseudo_counts <- function(x){
	pseudo_counts <- x + 0.5
	return(pseudo_counts)
}


# log2
# find ratio
# center and scale

znorm <- function(x0, x1, response_colname="pptm_psuedo", id_colname="shrna.id", min_counts=50, original_counts_colname="total.hits"){

	# find rows where the original counts for the
	# t0 group are below a threshold. exclude these
	# from both the t0 and t1 columns	
	rows_to_drop <- which(
		x0[,original_counts_colname] < min_counts
		)
	if(length(rows_to_drop) == 0){
		x0_cleaned <- x0
		x1_cleaned <- x1
	}else{
		x0_cleaned <- x0[-rows_to_drop,]
		x1_cleaned <- x1[-rows_to_drop,]
	}
	# log2 x0 and x1
	x0_log <- log2(x0_cleaned[,response_colname])
	x1_log <- log2(x1_cleaned[,response_colname])
	
	# find ratio
	x1_x0 <- x1_log / x0_log
	
	# center and scale
	x1_x0_med <- median(x1_x0)
	x1_x0_mad <- mad(x1_x0)
	x1_x0_zscore <- ( x1_x0 - x1_x0_med ) / x1_x0_mad
	names(x1_x0_zscore) <- x0_cleaned[,id_colname]
	tables_with_zscores <- cbind(
		x0_cleaned,
		x1_cleaned,
		x1_x0_zscore
		)
	return(tables_with_zscores)
}

make_scatter_plot <- function(
	x0,
	x1,
	gene="PARP1",
	filename="plot.pdf",
	x0name="library",
	x1name="sample",
	response_colname="pptm_psuedo",
	main=""
	){
	
	pdf(filename, width=5, height=5)

	plot(
		x0[,response_colname],
		x1[,response_colname],
		log="xy",
		xlab=x0name,
		ylab=x1name,
		main=main
		)

	rows_to_mark <- grep(gene, x0$shrna.id)

	for(row in rows_to_mark){
		points(
			x0[row,response_colname],
			x1[row,response_colname],
			pch=19,
			col="red"
			)
	}

	dev.off()
}



