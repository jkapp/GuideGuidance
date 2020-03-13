# findCommonMers
The first module of a greater program that will automate large scale guide oligo creation.

As guides are meant to target unique genomic regions it's beneficial to know how many times a specific stretch of bases is represented in a genome. If a stretch of bases is represented more than once then unintended cutting may occur. findCommonMers parses the forward and reverse direction of a reference genome and counts how many times all present 17mers are represented. 

findCommonMers takes in a reference genome in fasta format (-r) and a user defined count threshold (-t) to generate a set of common kmers. However, on it's own findCommonMers prints out all 17mers that are >= the user defined threshold. 

# Limitations
This shit here will DEFINITELY crash your mac air if you run it on a reletively large genome but if you have a server then it's fine and fast. There are over 17 billion possible 17mers which are populated in a dictionary as strings. Might want to store things as 1's and 0's, write temp files to the disk, and learn C. 
