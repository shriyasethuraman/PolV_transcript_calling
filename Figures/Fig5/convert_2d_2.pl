#script used to convert single bp data from an intersectBed output
#generates CSV table for use in R
#usage:
#column numbers start with 0 
#tasks: 0-generate table, 1-generate lengths
#does not properly handle non-existing columns!
#perl convert_2d.pl filename column task > output.csv

use autodie;

$input_file = $ARGV[0];			#First argument - name of the input file
$column = $ARGV[1];				#Second argument - column in the input file which contains output data
$task = $ARGV[2];				#Third argument - task; 0-generate table, 1-generate lengths
$sequence_length = 9295;		#Length of the sequences in the input file; has to be equal for all sequences! (normally 3180)
$column_position = 7;			#Column which contains position in the sequence
$columnID = 9;					#Column which contains sequence ID
$column_length = 6;				#Column with sequence length

open INPUT, $input_file;
if($task == 0){
	while(<INPUT>)					#go through every line of the file
	{
		chomp($_);
		@line = split(/\t/, $_);			#split line of the file into array
		if($line[$column_position] == 1)	#first line of a new sequence
		{
			print "$line[$columnID]_$line[$column_length]";		#print sequence name ID_length
			print ",";
			print $line[$column];			#print requested value
			print ",";
		}
		elsif($line[$column_position] == $sequence_length)	#last line of a sequence
		{		
			print $line[$column];			#print requested value
			print "\n";						#end of line	
		}
		elsif($line[$column_position] > $sequence_length)
		{
			die "Error: sequence longer than defined in the script";
		}
		else
		{
			print $line[$column];			#print requested value
			print ",";
		}
	}
}
if($task == 1){
	$previousID = -1;
	while(<INPUT>)					#go through every line of the file
	{
		chomp($_);
		@line = split(/\t/, $_);								#split line of the file into array
		if($line[$columnID] != $previousID){					#checks if a new sequence is found
			print "$line[$columnID]_$line[$column_length]";		#print sequence name ID_length
			print ",";
			print $line[$column_length];						#print sequence length
			print "\n";
			$previousID = $line[$columnID];						#saves column ID
		}
	}
}
