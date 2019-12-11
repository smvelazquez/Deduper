import argparse
import re

forward_reads = {}
reverse_reads =  {}
real_umis = []


def get_args():
     parser = argparse.ArgumentParser()
     parser.add_argument("-sam", "--FILEPATH_SAM", help="the absolute FILEPATH for your SAM file", required = True)
     parser.add_argument("-umis", "--FILEPATH_UMIS", help="the absolute FILEPATH with your known UMIs", required = True)
     return parser.parse_args()
args = get_args()

file = open(args.FILEPATH_SAM, "r")
umis = open(args.FILEPATH_UMIS, "r")
mapped_reads = open("matched_deduped.sam", "wt")
duplicate_reads = open("duplicates.sam", "wt")
#read in your UMIs
for line in umis:
    um = line.strip()
    real_umis.append(um)
print(real_umis)

def find_UMI(string):
    '''Extracts the UMI from the header string that is passed into it. NOTE: The UMIs must be 8 bp long'''
    UMI = re.findall(r'[ATGC]{8}', string)
    return(UMI)

def sum_of_cigar(cigar_string):
    '''Adds all of the numbers in the cigar string. NOTE: Should only be used in adjust_rev_position function.'''
    integers = re.findall(r'(\d+)', cigar_string)
    integers = [int(x) for x in integers]
    return(sum(integers))

def adjust_rev_position(cigar_string, position):
    '''This function will correct the position for the reverse reads. It will take the length of
    cigar string and add it to the reported left-most position from the SAM file to produce the
    right-most position for the reverse read. From here, it will adjust for soft clipping at the end
    of the cigar string and add that number to the position. If it encounters a D or an N it will also
    add these numbers to the position.'''
    sum_D = 0
    sum_N = 0
    sum_S = 0
    total = sum_of_cigar(cigar_string)
    right_pos = int(position) + total
    if 'S' in cigar_string:
        matchS = re.findall(r'(\d+)S', cigar_string)
        intlistS = [int(x) for x in matchS]
        sum_S = intlistS[-1]
    elif 'N' in cigar_string:
        matchN = re.findall(r'(\d+)N', cigar_string)
        intlistN = [int(x) for x in matchN]
        sum_N = [sum(x) for x in intlistN]
    elif 'D' in cigar_string:
        matchD = re.findall(r'(\d+)D', cigar_string)
        intlistD = [int(x) for x in matchD]
        sum_D = [sum(x) for x in intlistD]
    true_pos = right_pos - sum_S + sum_N + sum_D
    return(true_pos)

def adjust_forward_position(cigar_string, position):
    '''This function will correct the position of the forward reads if there is an S, otherwise the position will
    not be touched.'''
    sum_S = 0
    true_pos = 0
    if 'S' in cigar_string:
        matchS = re.findall(r'(\d+)S', cigar_string)
        intlistS = [int(x) for x in matchS]
        sum_S =[sum(x) for x in intlistS]
        true_pos = int(position) + sum_S
    else:
        true_pos = int(position)
    return(true_pos)

def is_this_a_duplicate(curr_UMI_string, curr_new_position, curr_chro, previous_reads):
    '''
    curr_read is UMI_string, new_position, chro
    previous_reads is this:
        keys: a two-tuple, first member is chromomsome/scaffold, second member is adjusted position
        vals: a list of UMIs (e.g. 'AACGCCAT')
    '''
    for location, umi_list in previous_reads.items():
        #print("is this a duplicate")
        chrom, adjusted_pos  = location
        if curr_UMI_string in umi_list and curr_chro == chrom and curr_new_position == adjusted_pos:
            #print("duplicate true")
            return True
        else:
            #print("duplicate false")
            return False

def deduper(file):
    '''This will check the orientation of your read and place into the appropriate dictionary (either forward or reverse).
    The key will be the chromosome and adjusted position, respectively; the value will be the UMI. From here, all other
    functions will be integrated so duplicates will be ignored and unique reads will be written out'''
    forward_count = 0
    reverse_count = 0
    #read in your file line by line
    for line in file:
        read = line # store your read
        line.strip("\n")

        #print("reading through lines")
        #break the sam file up into its important pieces
        if line[0] != '@':
            parts = line.split("\t")
            #return(parts)
            flags = parts[1]
            real_flags = int(flags)
            header = parts[0]
            chro = parts[2]
            position = parts[3]
            cigar_string = parts[5]

            #check the flags
            #if reverse do the following
            if((real_flags & 16)==16):
                print("made it to reverse")
                # reverse strand
                UMI_string = find_UMI(header) #pull out the UMI
                if UMI_string not in real_umis:
                    duplicate_reads.write(read)
                new_position = adjust_rev_position(cigar_string, position) #adjust the given position
                reverse_reads.setdefault((chro, new_position), UMI_string) #set the default of your empty dictionary
                #to be your first read; dictionary will be in the following format: {(chro, adj_pos): "UMI"}

                #use the is_is_this_a_duplicate to determine if the next read is a duplicate or not, this section will write out all unique reads and ignore duplicates
                if is_this_a_duplicate(UMI_string, forward_position, chro, forward_reads) == False:
                    #print("done")
                    reverse_reads[chro, new_position]= UMI_string
                    #print(read)
                    #print("no duplicates found")
                    mapped_reads.write(read)
                    reverse_count +=1
                print(reverse_reads)
            #now check for the forward reads, repeat the above for forward
            elif ((real_flags & 16) != 16):
                #print("made it to forward")
                UMI_string = find_UMI(header) #pull out the UMI
                print(UMI_string)
                if UMI_string not in real_umis:
                    duplicate_reads.write(read)
                forward_position = adjust_forward_position(cigar_string, position)
                forward_reads.setdefault((chro, forward_position), UMI_string)
                #print(forward_reads)
                #now check to see if the next reads are duplicates of the first, it'll keep any reads that aren't and ignore those who are
                if is_this_a_duplicate(UMI_string, forward_position, chro, forward_reads) == False:
                    forward_reads[chro, forward_position]= UMI_string #if read not found, add it to the dictionary and write it out
                    #print(read)
                    #print("no duplicates found")
                    mapped_reads.write(read)
                    forward_count +=1
    #print(forward_reads)
    mapped_reads.close()
    return("plus strand count:", forward_count, "minus strand count:", reverse_count)

deduper(file)
n):
    '''This function will correct the position for the reverse reads. It will take the length of
    cigar string and add it to the reported left-most position from the SAM file to produce the
    right-most position for the reverse read. From here, it will adjust for soft clipping at the end
    of the cigar string and add that number to the position. If it encounters a D or an N it will also
    add these numbers to the position.'''
    sum_D = 0
    sum_N = 0
    sum_S = 0
    total = sum_of_cigar(cigar_string)
    right_pos = int(position) + total
    if 'S' in cigar_string:
        matchS = re.findall(r'(\d+)S', cigar_string)
        intlistS = [int(x) for x in matchS]
        sum_S = intlistS[-1]
    elif 'N' in cigar_string:
        matchN = re.findall(r'(\d+)N', cigar_string)
        intlistN = [int(x) for x in matchN]
        sum_N = [sum(x) for x in intlistN]
    elif 'D' in cigar_string:
        matchD = re.findall(r'(\d+)D', cigar_string)
        intlistD = [int(x) for x in matchD]
        sum_D = [sum(x) for x in intlistD]
    true_pos = right_pos - sum_S + sum_N + sum_D
    return(true_pos)

def adjust_forward_position(cigar_string, position):
    '''This function will correct the position of the forward reads if there is an S, otherwise the position will
    not be touched.'''
    sum_S = 0
    true_pos = 0
    if 'S' in cigar_string:
        matchS = re.findall(r'(\d+)S', cigar_string)
        intlistS = [int(x) for x in matchS]
        sum_S =[sum(x) for x in intlistS]
        true_pos = int(position) + sum_S
    else:
        true_pos = int(position)
    return(true_pos)

def is_this_a_duplicate(curr_UMI_string, curr_new_position, curr_chro, previous_reads):
    '''
    curr_read is UMI_string, new_position, chro
    previous_reads is this:
        keys: a two-tuple, first member is chromomsome/scaffold, second member is adjusted position
        vals: a list of UMIs (e.g. 'AACGCCAT')
    '''
    for location, umi_list in previous_reads.items():
        #print("is this a duplicate")
        chrom, adjusted_pos  = location
        if curr_UMI_string in umi_list and curr_chro == chrom and curr_new_position == adjusted_pos:
            #print("duplicate true")
            return True
        else:
            #print("duplicate false")
            return False

def deduper(file):
    '''This will check the orientation of your read and place into the appropriate dictionary (either forward or reverse).
    The key will be the chromosome and adjusted position, respectively; the value will be the UMI. From here, all other
    functions will be integrated so duplicates will be ignored and unique reads will be written out'''

    #read in your file line by line
    for line in file:
        read = line # store your read
        line.strip("\n")

        #break the sam file up into its important pieces
        if line[0] != '@':
            parts = line.split("\t")
            #return(parts)
            flags = parts[1]
            real_flags = int(flags)
            header = parts[0]
            chro = parts[2]
            position = parts[3]
            cigar_string = parts[5]

            #check the flags
            if((real_flags & 256) == 256):
                continue

            #if reverse do the following
            if((real_flags & 16)==16):
                # reverse strand
                UMI_string = find_UMI(header) #pull out the UMI
                if UMI_string not in real_umis:
                    continue
                new_position = adjust_rev_position(cigar_string, position) #adjust the given position
                reverse_reads.setdefault((chro, new_position), UMI_string) #set the default of your empty dictionary
                #to be your first read; dictionary will be in the following format: {(chro, adj_pos): "UMI"}

                #use the is_is_this_a_duplicate to determine if the next read is a duplicate or not, this section will write out all unique reads and ignore duplicates
                if is_this_a_duplicate(UMI_string, forward_position, chro, forward_reads) == False:
                    #print("done")
                    reverse_reads[chro, new_position]= UMI_string
                    #print(read)
                    #print("no duplicates found")
                    mapped_reads.write(read + "\n")
                    reverse_count +=1

            #now check for the forward reads, repeat the above for forward
            elif ((real_flags & 16) != 16):
                UMI_string = find_UMI(header) #pull out the UMI
                if UMI_string not in real_umis:
                    continue
                forward_position = adjust_forward_position(cigar_string, position)
                forward_reads.setdefault((chro, forward_position), UMI_string)

                #now check to see if the next reads are duplicates of the first, it'll keep any reads that aren't and ignore those who are
                if is_this_a_duplicate(UMI_string, forward_position, chro, forward_reads) == False:
                    forward_reads[chro, forward_position]= UMI_string #if read not found, add it to the dictionary and write it out
                    #print(read)
                    #print("no duplicates found")
                    mapped_reads.write(read + "\n")
                    forward_count +=1
    mapped_reads.close()
    return("plus strand count:", forward_count, "minus strand count:", reverse_reads)
    
    deduper(file)
