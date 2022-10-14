import sys
import os

in_file = open(sys.argv[1])
out_file = open(sys.argv[2],"w")
bed_out = open(sys.argv[3],"w")
min_support_snp = 10
current_chr = ""
for line in in_file.readlines():
    line = line.strip()
    l_array = line.split("\t")
    current_count = 0
    start = ""
    end = ""
    is_nahr = False
    is_bnd = False
    if l_array[0] == "chrom":
        current_chr = l_array[1]
        out_file.write(line+"\n")
        continue
    if len(l_array) == 1:
        recom_event = l_array[0].split("_")
        if "DENOVO" in l_array[0]:
            out_file.write(recom_event[0]+"\t"+recom_event[0]+"\t"+str(len(l_array))+"\tDENOVO\n")
        elif "SIMPLER" in l_array[0]:
            out_file.write(recom_event[0]+"\t"+recom_event[0]+"\t"+str(len(l_array))+"\t_SIMPLERUNNAHR\n")
            bed_out.write(current_chr+"\t"+str(recom_event[5])+"\t"+str(recom_event[6])+"\n")
        elif "_NAHR" in l_array[0]:
            out_file.write(recom_event[0]+"\t"+recom_event[0]+"\t"+str(len(l_array))+"\t_NAHR\n")
        # if "BND" in l_array[0]:
        #     out_file.write(recom_event[0]+"\t"+recom_event[0]+"\t"+str(len(l_array))+"\tBND\n")
    else:
        # print(len(l_array),min_support_snp)
        # if len(l_array) <= min_support_snp:
        #     continue
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
