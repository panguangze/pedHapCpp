import os
import sys

lst_file = open(sys.argv[1])
lst_out = open(sys.argv[2], "w")
pos = int(sys.argv[3])
# poses = sys.argv[3].split(",")
# poses_list = range(int(poses[0]), int(poses[1]) + 1)
keep = False
for line in lst_file.readlines():
    line_split = line.split(" ")
    lst_count = int(line_split[0])
    # print(list(range(2, 2+2*lst_count, 2)))
    for i in range(2, 2+2*lst_count, 2):
        if pos in range(int(line_split[i]), int(line_split[i]) + len(line_split[i+1])):
            lst_out.write(line_split[i]+"\t"+line_split[i+1]+"\n")
            # keep = True
            break
    # if keep:
    #     lst_out.write(line)
    #     keep = False
