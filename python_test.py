#!/usr/bin/env python

from __future__ import division
import sys
import re

# def empty(seq):
#     try:
#         return all(empty(x) for x in seq)
#     except TypeError:
#         return False



def main():

    # curr_clust = None
    #
    # for line in [">Cluster 13\n","0	3939nt, >ERR348862_00288 group_1460... *\n",">Cluster 14\n","0	3933nt, >0_meta_02819 predicted membrane-associated_2... *\n","1	3933nt, >ERR343065_01918 predicted membrane-associated... at +/99.34%\n"]:
    #
    #     # print(line)
    #
    #     #if line starts with a >,
    #     if (">Cluster" in line):
    #         #create a cluster object, if curr cluster wasn't null, print all of its members to the output file
    #         if (curr_clust):
    #             print(curr_clust)
    #
    #         #and then update the curr cluster pointer
    #         curr_clust = Cluster(line)
    #
    #     #if line doesn't start with >
    #     else:
    #         #add it to the current cluster as a cluster member
    #         curr_clust.add_member(line)
    #
    # #print the last cluster
    # print(curr_clust)



    # seq_id_dict = {}
    # seq_count = 1
    # line = "##sequence-region 874761_01096 1 387"
    # match = re.search("(##sequence-region )([0-9]*_[0-9]*)( .*)",line)
    # seq_id_dict[match.group(2)] = "0_metagenome_{}".format(seq_count)
    # seq_count += 1
    # print("{first}{id}{last_part}".format(first=match.group(1), id=seq_id_dict[match.group(2)], last_part=match.group(3)))

    # print 'Hello world!', sys.argv[1]
    # print 3/4
    # print 3//4
    #
    #
    # str = "a,bdas lshjgdva sldahvls"
    #
    # if "das" in str:
    #     print "foo"
    #
    # match = re.findall(r'a..', str)
    #
    # for m in match:
    #     print m
    #
    #
    # f = open('test.txt', 'w')
    # f.write("test output file")
    # f.close()
    #
    # start = (re.search('(^\w+\.[0-9])(\[)',"AL111168.1[1001573:1001837](+)").group(1))
    #
    # print(start)
    # print("\n")
    #
    # var = float("52.78")
    #
    # var += 1
    #
    # print(var)

    # arr = []
    # arr.append("sadf")
    # arr.append("asdfas")
    # arr.append("rew")
    # arr.append("ssdfkjfsd")
    # arr.append("s")
    # # arr = ["sadf","asdfas","rew","ssdfkjfsd","s"]
    #
    # print(arr)
    #
    # print(sorted(arr,key = lambda x : len(x)))

    #
    # temp = re.search("(5 +5)", "4565     567")
    # print(temp.group(1))

    # print (len(sys.argv))
    # if (sys.argv[2]):
    #     print("foo")

    # print(sys.argv)
    # print(len(sys.argv))
    # for index in range(2,len(sys.argv)):
    #     print(index)



    # arr = [[],[[1]],[]]
    # arre = [[],[[]],[]]
    #
    # print(empty(arr))
    # print(empty(arre))

#     test() #prints 3
#     test(x=2) #prints 2


#
#
#
# def test(x=3):
#         print(x)

def test_funct():
    print("foo")

# class ClustMember:
#     """A member of a cluster
#     Will parse a line from a clstr file to create the object"""
#     def __init__(self, line):
#         match = re.search("[0-9]+\t([0-9]+)nt, >(.*?) (.*)\.\.\. (\*)?(at [+-]/)?([0-9]+\.[0-9][0-9])?%?",line)
#         self.length = int(match.group(1)) #Sequence length
#         self.seq_id = match.group(2) #Sequence ID
#         self.gene_name = match.group(3) #Sequence gene name
#         self.representative = True if match.group(4) else False #Is this sequence the cluster representative?
#         self.identity = None if match.group(4) else float(match.group(6)) #If not, what is its percent identity to the representative?
#
#     def __str__(self):
#         return self.seq_id
#
#     def __repr__(self):
#         return self.seq_id
#
#
# class Cluster:
#     """Contains a list of cluster members"""
#     def __init__(self,line):
#         match = re.search(">.* ([0-9]*)",line)
#         self.cluster_id = match.group(1) #Cluster ID
#         self.members = [] #list of ClustMember objects
#
#     def __str__(self):
#         mems = ""
#         if (len(self.members) > 0):
#             mems = str(self.members[0])
#             for i in range(1,len(self.members)):
#                 mems += ",{}".format(self.members[i])
#
#         return "Cluster {},{}".format(self.cluster_id, mems)
#
#     def add_member(self,member):
#         if (type(member) == ClustMember):
#             self.members.append(member)
#         elif (type(member) == str):
#             self.members.append(ClustMember(member))
#         else:
#             print("That wasn't a valid member of any sort\n")
#
#     def print_cluster(self):
#         print("Cluster_",end="")
#         print(self.cluster_id,end=",")
#         self.print_members()
#
#     def print_members(self):
#         if (len(self.members) > 0):
#             print(self.members[0],end="")
#
#             for i in range(1,len(self.members)):
#                 print(",",end="")
#                 print(self.members[i],end="")
#
#             print()






if __name__ == '__main__':
    main()
