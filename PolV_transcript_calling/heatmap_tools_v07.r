#tools for generating heatmaps of Pol V transcripts
#load to R using source('heatmap_tools_v07.r')
#generate input files using perl convert_2d.pl input column mode > output
#input - file from IntersectBed
#column - number of column with data to be extracted staring at 0
#mode 0 - extract data, 1 - extract lengths
#output - CSV file to use with calculate_values() in R

#version 0.7 - essential parts with randomization-based statistics for generating figures
#version 0.61- essential parts only - for Gudrun
#version 0.6 - organize into essential and additional parts
#version 0.5 - added ability to compare whole transcripts
#version 0.4 - added ability to work with DNA methylation data
#version 0.3 - added ability to draw vertical lines which delimit ends of transcripts
#version 0.2 - added ability to draw unscaled 3' ends

#install required libraries
if (!require("gplots")) {
 install.packages("gplots", dependencies = TRUE)
 library(gplots)
 }
if (!require("RColorBrewer")) {
 install.packages("RColorBrewer", dependencies = TRUE)
 library(RColorBrewer)
 }

#calculate_values calculates values from two tables
#values will be used for generating a heatmap
#arguments: 
#file1 - value file 1 (experimental), has to be a CSV file
#file2 - value file 2 (control), has to be a CSV file
#RPM1 - value for RPM normalization of file 1 (1M/read count)
#RPM2 - value for RPM normalization of file 2 (1M/read count)
#mode 1 - perform subtraction (use RPM 1 to keep raw read counts), 2 - perform division
#returns matrix with calculated data
calculate_values <- function(file1, file2, RPM1, RPM2, mode)		#
{
	table1 <- read.csv(file1,header = F)									#load CSV files
	table2 <- read.csv(file2,header = F)
	matrix1 <- data.matrix(table1[,2:ncol(table1)])				#convert files to matrixes without rowname column
	matrix2 <- data.matrix(table2[,2:ncol(table2)])
	rownames1 <- data.matrix(table1[,1])						#load row names
	rownames2 <- data.matrix(table2[,1])
	if(length(setdiff(rownames1, rownames2))!=0)				#check if rownames are identical between input files
	{
		cat("WARNING: ROW NAMES OF INPUT DATA FILES DO NOT MATCH\n")
	}
	
	if(mode == "1"){											#mode 1 - perform subtraction of RPM values
		result <- (matrix1 / RPM1) - (matrix2 / RPM2)
	}
	if(mode == "2"){											#mode 2 - perform divisions of RPM values
																#add equivalent of 0.5 reads in File 1 RPM to both values
		result <- (((matrix1 / RPM1) + (0.5 / RPM1)) / ((matrix2 / RPM2) + (0.5 / RPM1)))
	}
	rownames(result) <- rownames1								#give proper row names to the output
	return(result)
}

#load_sizes - loads file with sequence lengths
#checks if it matches the heatmap matrix
#arguments: 
#sizes - file name of the sizes file and 
#data_matrix - data to be used for heatmap generation
#returns a list of sequence lengths
load_sizes <- function(sizes, data_matrix)
{
	table1 <- read.csv(sizes)									#load CSV file
	matrix1 <- data.matrix(table1[,2])							#convert files to matrixes without rowname column
	rownames1 <- data.matrix(table1[,1])						#load row names
	rownames(matrix1) <- rownames1								#give proper row names
	if(length(setdiff(rownames(matrix1), rownames(data_matrix)))!=0)	#check if rownames are identical between input files
	{
		cat("WARNING: ROW NAMES OF SIZES DO NOT MATCH\n")
	}
	return(matrix1)
}

#percent_trim - trims certain percent of sequences from both ends of the spectrum
#arguments: 
#matrix - input matrix; 
#percent - what percentage is trimmed from top and bottom (eg. 10 = 10 percent)
#returns trimmed data matrix or list
percent_trim <- function(matrix, percent)
{
	start <- round(nrow(matrix) * percent / 100, digits = 0)					#calculate first position to keep
	end <- nrow(matrix) - round(nrow(matrix) * percent / 100, digits = 0)		#calculate last position to keep
	output_matrix <- matrix[start:end,]											#produce trimmed matrix
	return(output_matrix)
}

#scale_heatmap - generates heatmap with transcribed regions scaled to the same lengths
#arguments: 
#input_matrix - a matrix with 300 bp upstream of TSS and at least 300 bp downstream of TTS of the longest transcript
#sizes - list of sizes of all transcripts
#returns a matrix with scaled heatmap data
scaler <- function(input_matrix, sizes)
{
	scaled_length <- 500															#length of scaled region
	output <- matrix(NA, nrow=nrow(input_matrix), ncol=scaled_length + 600)			#create output scaled array
	for (count_row in 1:nrow(output)){												#go over every row of the array
		#cat("row ")
		#cat(count_row)																#print progress
		#cat("\n")
		for(count_col in 1:ncol(output)){											#go over every colunm in the row
			if(count_col <= 300){													#for first 300nt and after end of transcript no scaling
				output[count_row, count_col] <- input_matrix[count_row, count_col]	#copy with no scaling
			}
			if(count_col > scaled_length + 300){									#for first 300nt and after end of transcript no scaling
				output[count_row, count_col] <- input_matrix[count_row, (count_col-scaled_length+sizes[count_row])]	#copy with no scaling
			}		
			if(count_col >= 300 && count_col <= scaled_length + 300)	{			#within transcript perform scaling
			
				n_begin <- 300 + floor(sizes[count_row] * (count_col - 1 - 300) / scaled_length)	#calculate starting point in input array
				n_end <- 300 + ceiling(sizes[count_row] * (count_col - 300) / scaled_length)		#calculate end point in input array
			
				#cat("row ",count_row,", length", sizes[count_row],", column ", count_col, ", n_begin ", n_begin, ", n_end", n_end, "\n")
			
				average <- c()														#declare empty vector
				for(counter in n_begin:n_end){										#go over cells scalled
					average <- c(average, input_matrix[count_row, counter])			#add value to vector
				}
				output[count_row, count_col] <- mean(average)						#record average to scaled matrix
			}
		}	
	}
	return(output)
}

#profiler - generates heatmap without scaling around 5' end
#arguments:
#input_matrix - a matrix with 300 bp upstream of TSS and at least 300 bp downstream of TTS of the longest transcript
#sizes - list of sizes of all transcripts
#length - length of transcriped sequence that will be included
#end - 5 - 5' end, 3 - 3' end
#returns a matrix with regions outside of transcripts replaced with NA
profiler <- function(input_matrix, sizes, length, end)
{
	#length <- 300																		#lenght included sequences downstream of end
	if(end == "5")																		#profile around 5' end
	{
		output <- input_matrix[,1:(length + 300)]										#create output array
		for (count_row in 1:nrow(output)){												#go over every row of the input matrix
			#cat("row", count_row, "\n")
			if(sizes[count_row] <= length)
			{
				output[count_row, (300 + sizes[count_row]):ncol(output)] <- NA			#make every cell beyond 3'end NA
			
			}
		}
	}
	if(end == "3")																		#profile around 3' end
	{
		output <- matrix(NA, nrow=nrow(input_matrix), ncol=length + 300)	
		for (count_row in 1:nrow(output)){												#go over every row of the input matrix
			#cat("row", count_row, "\n")
			if(sizes[count_row] < length)
			{
				#cat("Smaller, row ", count_row, "Length ", sizes[count_row], "\n")
				output[count_row, (length - sizes[count_row]):(length + 300)] <- input_matrix[count_row, 300:(300+sizes[count_row]+300)] 
			}
			if(sizes[count_row] >= length)			
			{
				#cat("Larger, row ", count_row, "Length ", sizes[count_row], "length+300", (length + 300), "300+sizes[count_row]-length", (300+sizes[count_row]-length), "300+sizes[count_row]+300", (300+sizes[count_row]+300), "\n")
				output[count_row, 1:(length + 300)] <- input_matrix[count_row, (300+sizes[count_row]-length+1):(300+sizes[count_row]+300)] 
			}
		}
	}	
	return(output)
}

#performs heatmap analysis and saves results to files
#arguments:
#input1 - CSV file with 1st set of data eg. "data.csv"
#input2 - CSV file with 2nd set of data
#sizes - CSV file with lengths
#RPM1 and 2 - RPM value for sequence 1 and 2
#USAGE: heatmap_pipeline("file1.csv", "file2.csv", "sizes.cvs", 30, 20)  
heatmap_pipeline <- function(input1, input2, sizes, RPM1, RPM2)
{
	values <- calculate_values(input1, input2, RPM1, RPM2, 2)	#1-subtraction, 2-division
	cat("Values loaded\n")
	sizes <- load_sizes(sizes, values)
	cat("Sizes loaded and row names checked\n")
	values_trim <- percent_trim(values, 5)
	sizes_trim <- percent_trim(sizes, 5)
	cat("Data trimmed\n")
	profile_unscaled5 <- profiler(values_trim, sizes_trim, 500, 5)
	profile_unscaled3 <- profiler(values_trim, sizes_trim, 500, 3)
	cat("Unscaled profiles calculated\n")
	profile_scaled <- scaler(values_trim, sizes_trim)
	cat("Scaled profile calculated\n")
	png(paste(input1,input2,"unscaled5.png"), width=8.5, height=11, units="in", res=300)
	figure(profile_unscaled5, 300, 300)
	dev.off()
	return(profile_unscaled5)
	cat("Unscaled profile 5 saved\n")
	png(paste(input1,input2,"unscaled3.png"), width=8.5, height=11, units="in", res=300)
	figure(profile_unscaled3, 500, 500)
	dev.off()
	cat("Unscaled profile 3 saved\n")	
	png(paste(input1,input2,"scaled.png"), width=8.5, height=11, units="in", res=300)
	figure(profile_scaled, 300, 800)
	dev.off()
	cat("Scaled profile saved\n")	
	sig5 <- significance(profile_unscaled5, 100, 350, 100, 1000)
	cat("significance on 5 end: ")
	sig5	
}

#calculates data for a scaled profile
calculate_scaled_data <- function(input1, input2, sizes, RPM1, RPM2){
	values <- calculate_values(input1, input2, RPM1, RPM2, 2)	#1-subtraction, 2-division
	cat("Values loaded\n")
	sizes <- load_sizes(sizes, values)
	cat("Sizes loaded and row names checked\n")
	values_trim <- percent_trim(values, 5)
	sizes_trim <- percent_trim(sizes, 5)
	cat("Data trimmed\n")
	profile_scaled <- scaler(values_trim, sizes_trim)
	cat("Scaled profile calculated\n")
	return(profile_scaled)
}

#calculates data for a scaled profile
calculate_unscaled_data <- function(input1, input2, sizes, RPM1, RPM2){
	values <- calculate_values(input1, input2, RPM1, RPM2, 2)	#1-subtraction, 2-division
	cat("Values loaded\n")
	sizes <- load_sizes(sizes, values)
	cat("Sizes loaded and row names checked\n")
	profile_unscaled <- profiler(values, sizes, 500, 5)  			#use last option to change ends: 5 for 5', 3 for 3'
	cat("Scaled profile calculated\n")
	return(profile_unscaled)
}


#figure - generates a figure with a plot and a heatmap
#arguments 
#input_matrix - a matrix with heatmap data
#break1 and 2 - positions of vertical lines
#plots a plot using linear values and a heatmap using log values
figure <- function(input_matrix, break1, break2)
{
	breaks <- matrix(NA, nrow=nrow(input_matrix), ncol=2)
	breaks[,1] <- break1																	#create matrix with column breaks
	breaks[,2] <- break2
	#reduced <- log2(input_matrix[seq(1, nrow(input_matrix), 5),])							#keeps every 5th sequence, log2 transformation of values
	reduced <- input_matrix[seq(1, nrow(input_matrix), 1),]								#keeps every 10th sequence
	#my_palette <- colorRampPalette(c("green", "yellow", "red"))(n=299)						#define color palette - green - red
	my_palette <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)				#define color palette - shades of red
	col_breaks = c(seq(quantile(reduced, 0.10, na.rm=TRUE),quantile(reduced, 0.90, na.rm=TRUE),length=300))	#define color breaks using percentiles
	#col_breaks = c(seq(0,5,length=300))																	#define color breaks using specific vlaues
	heatmap.2(reduced, col=my_palette, trace="none", breaks=col_breaks,dendrogram="none", Colv="NA", 		#generate heatmap
	Rowv="NA", key=FALSE, labRow=FALSE, labCol=FALSE, lmat = rbind(c(4,3,0),c(2,1,3)), lwid=c(1,5,0.3), lhei=c(4.5,7), colsep=breaks)
	
  	par(fig=c(0.035,0.955,0.58,1), new=TRUE)																#define figure arrangement
	plot(colMeans(input_matrix, na.rm=TRUE), type="l", ylab="", xlab="", cex.axis=1)				#generate profile
	abline(v=break1, col="grey")
	abline(v=break2, col="grey")
}

#calculate pvalue of enrichment using a randomization test
significance <- function(input_matrix, reg1, reg2, region_length, repeats)
{
  me <- function(x) mean(x, na.rm=TRUE)                   #define me mean function able to deal with NAs
  sub_matrix1 <- input_matrix[,reg1:(reg1+region_length)] #extract sub matrix1
  list1 <- as.numeric(apply(sub_matrix1, 1, me))          #calculate row means of sub matrix1
  sub_matrix2 <- input_matrix[,reg1:(reg2+region_length)] #the same for sub matrix2
  list2 <- as.numeric(apply(sub_matrix2, 1, me))
  averages <- c(list1, list2)                             #combine lists of averages
  categories <- c(rep(0, times=length(list1)), rep(1, times=length(list2))) #populate the categories vector (list1-0, list2-1)
                                                          #calculate permutations
  dist <- replicate(repeats, diff( by(averages, sample(categories, length(categories), FALSE), mean) ) )
  #hist(dist, col = "black", breaks = 100)                 #draw histogram
  #abline(v = diff(by(averages, categories, mean)), col = "blue", lwd = 2)
  pval <- sum(dist > diff(by(averages, categories, mean)))/repeats #calculate p value
  if(pval == 0) {pval = 1/repeats}
  means <- by(averages, categories, mean)
  result <- c()
  result$pval <- pval
  result$ratio <- means[2] / means[1]
  return(result)
}





