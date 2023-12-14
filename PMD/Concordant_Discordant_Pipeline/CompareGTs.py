import sys

#Input files (file1 = Bcftools vcf file, file2 = Atlas vcf file)
file1_lines = open(sys.argv[1], 'r').readlines()
file2_lines = open(sys.argv[2], 'r').readlines()

#Output files
conc_list = open('conc_list.txt', 'w')
bcfdisc_list = open('bcfdisc_list.txt', 'w')
atlasdisc_list=open('atlasdisc_list.txt', 'w')

#Creates a dictionary for the Bcftools file where the keys are tuples (col1=Chromosome, col2=SNP Position) and values are col3 (GT)
file1_dict = {}
for line in file1_lines:
    col1, col2, col3 = line.strip().split('\t')
    file1_dict[(col1, col2)] = col3

#Creates a dictionary for the Atlas file the same as Bcftools
file2_dict = {}
for line in file2_lines:
    col1, col2, col3 = line.strip().split('\t')
    file2_dict[(col1, col2)] = col3

#Compares entries between the two dictionaries
for key in file1_dict:
    #Checks if a key (SNP) exists in both dictionaries and its values are not "./." or "."
    if key in file2_dict and file1_dict[key] not in ["./.", "."] and file2_dict[key] not in ["./.", "."]:
        #If values (GTs) for a key (SNP) in both dictionaries are the same, write to the conc_list file
        if file1_dict[key] == file2_dict[key]:
            conc_list.write(key[0] + "\t" + key[1] + "\n")
        #If values (GTs) differ, write to the bcfdisc_list and atlasdisc_list files
        else:
            bcfdisc_list.write(key[0] + "\t" + key[1] + "\n")
            atlasdisc_list.write(key[0] + "\t" + key[1] + "\n")
    #If the key (SNP) is only present in the first dictionary and is not "./." or ".", write to the bcfdisc_list file
    elif file1_dict[key] not in ["./.", "."]:
        bcfdisc_list.write(key[0] + "\t" + key[1] + "\n")

#Checks for keys (SNPs) that exist in the second dictionary but not in the first
for key in file2_dict:
    #If a key (SNP) only exists in the second dictionary and its value is not "./." or ".", write to the atlasdisc_list file
    if key not in file1_dict:
        if file2_dict[key] not in ["./.", "."]:
            atlasdisc_list.write(key[0] + "\t" + key[1] + "\n")
    #If the key (SNP) exists in both, but only the second dictionary has a meaningful value, write to the atlasdisc_list file
    elif (file1_dict[key] in ["./.", "."]) and (file2_dict[key] not in ["./.", "."]):
        atlasdisc_list.write(key[0] + "\t" + key[1] + "\n")

#Closes the output files
conc_list.close()
bcfdisc_list.close()
atlasdisc_list.close()

