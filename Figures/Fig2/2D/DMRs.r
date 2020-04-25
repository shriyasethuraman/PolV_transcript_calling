#Find clusters of DNA methylation
#takes: matrix, methylation cutoff, maximum distance between methyl marks, minimum number of methyl marks
#returns a matrix with positions and average intensities of methylation clusters
DMRs <- function(matrix1, cutoff, distance, number)
{
	output <- matrix(0, nrow=nrow(matrix1), ncol=ncol(matrix1))	#define output array
	DMR_list <- c(NA,NA,NA,NA)
	DMR_counter <- 0;
	for (count_row in 1:nrow(matrix1))								#go over every row of the matrix
	{
		#cat("Row ", count_row, "\n")
		count_col <- 1												#reset column counter
		DMR <- 0													#declare flag if within a possible DMR
		current_distance <- 0										#declare variable to measure current distance from previous methyl mark

		repeat																										#
		{																											#
			#cat("Looking for first meC ", count_col, " ", matrix1[count_row, count_col], "\n")
			if(matrix1[count_row, count_col] > cutoff && is.finite(matrix1[count_row, count_col])) break()			#go to first methylated position 
			if(count_col >= ncol(matrix1)) break()																	#if no methylation, break at the end
			count_col <- count_col + 1																				#
		}																											#

		while(count_col < ncol(matrix1))
		{
			if(matrix1[count_row, count_col] > cutoff && is.finite(matrix1[count_row, count_col]))					#if current position is methylated
			{
				if(DMR == 0)											#if not within a DMR
				{
					
					start <- count_col									#define start of a new DMR
					end <- count_col									#define end of a new DMR
					meC_sum <- matrix1[count_row, count_col]			#start adding methylation intensities
					DMR <- 1											#flag that within a DMR
					meC_counter <- 1									#start counting methylated Cs within DMR
					current_distance <- 0								#start measuring distance from DMR
					#cat("1\t", count_row, "\t", count_col, "\t", matrix1[count_row, count_col], "\t", DMR, "\t", start, "\t", end, "\t", meC_counter, "\t", current_distance, "\n")
				}else													#if within a DMR
				{
					current_distance <- 0								#start measuring distance from DMR
					end <- count_col									#define end of the DMR
					meC_sum <- meC_sum + matrix1[count_row, count_col]	#add methylation intensity
					meC_counter <- meC_counter + 1						#count methylated Cs within DMR
					#cat("2\t", count_row, "\t", count_col, "\t", matrix1[count_row, count_col], "\t", DMR, "\t", start, "\t", end, "\t", meC_counter, "\t", current_distance, "\n")
				}
			}else														#if current position not methylated
			{
				current_distance <- current_distance + 1				#increment distance
				if(current_distance > distance || count_col + current_distance >= ncol(matrix1))							#if maximum distance is exceeded
				{
					if(meC_counter < number && DMR == 1)							#if number of meC lower than the minimum - reset DMR
					{
						
						count_col <- start +1									#go back to the start of the unsuccesful DMR
						DMR <- 0											#set DMR flag at 0
						meC_counter <- 0									#reset counter
						start <- 0											#reset start
						end <- 0											#reset end
						#cat("3\t", count_row, "\t", count_col, "\t", matrix1[count_row, count_col], "\t", DMR, "\t", start, "\t", end, "\t", meC_counter, "\t", current_distance, "\n")
					}
					if(meC_counter >= number && DMR == 1)							#if number of meC higher or equal to the minimum - save DMR
					{
						for(write_position in start:end)
						{
							output[count_row, write_position] <- meC_sum / meC_counter
						}
						DMR_counter <- DMR_counter + 1
						DMR_list <- cbind(DMR_list, c(count_row, start, end, meC_sum / meC_counter))
						#cat("CALLED DMR\t", count_row, "\t", count_col, "\t", matrix1[count_row, count_col], "\t", DMR, "\t", start, "\t", end, "\t", meC_counter, "\t", current_distance, "\t", meC_sum / meC_counter, "\n")
						DMR <- 0											#set DMR flag at 0
						meC_counter <- 0									#reset counter
						count_col <- end
						start <- 0											#reset start
						end <- 0											#reset end
					}
				}
			
			}
			count_col <- count_col + 1								#increment column counter		
		}
	}
	cat("Produced ", DMR_counter, " DMRs\n")
	return(output)

}