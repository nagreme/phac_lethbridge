#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# February 24, 2017
# Purpose: Parse .clstr files write gene names clusters to a file.
#          Contains objects to facilitate interacting with .clstr files
#   sys.argv[1] = .clstr file
#----------------------------------------------------------------------

import sys
import re

#---------------------------------Main---------------------------------

def main():

    #open input .clstr file
    inFile = open(sys.argv[1], 'rU')
    #open output file
    outFile = open("{}.csv".format(sys.argv[1][0:sys.argv[1].find(".")]), 'w')

    curr_clust = None

    #parse through
    for line in inFile:
        #if line starts with a >,
        if (">Cluster" in line):
            #create a cluster object, if curr cluster wasn't null, print all of its members to the output file
            if (curr_clust):
                outFile.write(str(curr_clust))

            #and then update the curr cluster pointer
            curr_clust = Cluster(line)

        #if line doesn't start with >
        else:
            #add it to the current cluster as a cluster member
            curr_clust.add_member(line)

    #print the last cluster
    outFile.write(str(curr_clust))

    #close the files
    inFile.close()
    outFile.close()


class ClustMember:
    """A member of a cluster
    Will parse a line from a clstr file to create the object"""
    def __init__(self, line):
        match = re.search("[0-9]+\t([0-9]+)nt, >(.*?) (.*)\.\.\. (\*)?(at [+-]/)?([0-9]+\.[0-9][0-9])?%?",line)
        self.length = int(match.group(1)) #Sequence length
        self.seq_id = match.group(2) #Sequence ID
        self.gene_name = match.group(3) #Sequence gene name
        self.representative = True if match.group(4) else False #Is this sequence the cluster representative?
        self.identity = None if match.group(4) else float(match.group(6)) #If not, what is its percent identity to the representative?

    def __str__(self):
        return self.gene_name

    def __repr__(self):
        return self.gene_name


class Cluster:
    """Contains a list of cluster members"""
    def __init__(self,line):
        match = re.search(">.* ([0-9]*)",line)
        self.cluster_id = match.group(1) #Cluster ID
        self.members = [] #list of ClustMember objects

    def __str__(self):
        #build a comma-separated list of the members
        mems = ""
        if (len(self.members) > 0):
            mems = str(self.members[0])
            for i in range(1,len(self.members)):
                mems += "\t{}".format(self.members[i])

        return "Cluster {}\t{}\n".format(self.cluster_id, mems)

    #modification March 1st to ensure that the representative is the first member in the list
    def add_member(self,member):
        if (type(member) == str):
            member = ClustMember(member)
        elif (type(member) != ClustMember):
            member = None
            print("That wasn't a valid member of any sort\n")

        if (member):
            if (member.representative):
                self.members.insert(0,member)
            else:
                self.members.append(member)

    def print_cluster(self):
        print("Cluster_",end="")
        print(self.cluster_id,end=",")
        self.print_members()

    def print_members(self):
        if (len(self.members) > 0):
            print(self.members[0],end="")

            for i in range(1,len(self.members)):
                print(",",end="")
                print(self.members[i],end="")

            print()




if __name__ == '__main__':
    main()
