import sys
import os

in_file = open(sys.argv[1])
out_file = open(sys.argv[2],"w")
min_support_snp = 3

for line in in_file.readlines():
    line = line.strip()
    l_array = line.split("\t")
    current_count = 0
    start = ""
    end = ""
    is_nahr = False
    is_bnd = False
    if l_array[0] == "chrom":
        out_file.write(line+"\n")
        continue
    if len(l_array) == 1:
        recom_event = l_array[0].split("_")
        if "NAHR" in l_array[0]:
            out_file.write(recom_event[0]+"\t"+recom_event[0]+"\t"+str(len(l_array))+"\tNAHR\n")
        if "BND" in l_array[0]:
            out_file.write(recom_event[0]+"\t"+recom_event[0]+"\t"+str(len(l_array))+"\tBND\n")
    else:
        if len(l_array) <= min_support_snp and "BND" not in line:
            continue
        # print(l_array)
        start = l_array[0].split("_")[0]
        end = l_array[-1].split("_")[0]
        out_file.write(start+"\t"+end+"\t"+str(len(l_array))+"\tHR\n")
    # for item in l_array:
        #     item_array = item.split("_")
        #     out_file
    # if l_array[0] == "#":
    #     if current_count >=2 or is_nahr:
    #         out_file.write(start)
    #         current_count = 0
    #         current_content = ""
    # elif l_array[4] == "NAHR":
    #     is_nahr = True
    #     current_content = l_array[0] + "\t" + "NAHR"
    # elif "BND" in l_array[4]:
    #     is_bnd = True
    #     current_content = l_array[0] + "\t" + "BND"
    # else:
