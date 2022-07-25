import os
import sys

# def flip(s):
#     r = ""
#     for i in s:
#         if i == '1':
#             r = r+ '0'
#         elif i == '0':
#             r = r + '1'
#         else:
#             r = i
#     return r
# def flip_s(i):
#     r = ""
#     if i == '1':
#         r = r+ '0'
#     elif i == '0':
#         r = r + '1'
#     else:
#         r = i
#     return r
#
# def compare(hap1,hap2):
#     score1 = 0
#     score2 = 0
#     for i in range(0,len(hap1)):
#         if hap1[i] == hap2[i]:
#             score1 = score1 + 1
#         if hap1[i] == flip_s(hap2[i]):
#             score2 = score2 + 1
#     return
#     # s = min(s1,s2)
#     # e = min(e1,e2)
#     # for i in range(0,)
#     # kmax = max(k1[0], k2[0])
#     # kmin = min(k1[0], k2[0])
#     # c = kmax - kmin
#     # for i in range(c-1)
# lst_file = open(sys.argv[1])
# lst_out = open(sys.argv[2], "w")
# # pos = int(sys.argv[3])
# # spos =
# poses = sys.argv[3].split(",")
# sub_hap = sys.argv[4]
# # poses_list = range(int(poses[0]), int(poses[1]) + 1)
# keep = False
# for line in lst_file.readlines():
#     line_split = line.split(" ")
#     lst_count = int(line_split[0])
#     # print(list(range(2, 2+2*lst_count, 2)))
#     for i in range(2, 2+2*lst_count, 2):
#         for p in poses:
#             if p in range(int(line_split[i]), int(line_split[i]) + len(line_split[i+1])):
#                 lst_out.write(line_split[i]+"\t"+line_split[i+1]+"\t"+str(int(line_split[i]) + len(line_split[i+1]))+"\n")
#                 # keep = True
#                 break
    # if keep:
    #     lst_out.write(line)
    #     keep = False
# def parse_vcf(vcf_file, lst1, lst2):
#     r1 = [[]]
#     r2 = [[]]
#     vcf_stream = open(vcf_file)
#     idx1 = 0
#     idx2 = 0
#     i = 468011
#     current_1 = ""
#     current_2 = ""
#     for line in vcf_stream.readlines():
#         if "#" in line:
#             continue
#         else:
#             if i == 469457:
#                 kk = 33
#                 pass
#             line_split = line.split("\t")
#             # if idx1 < len(lst1):
#             #     if i in lst1[idx1]:
#             #         if current_1 == "":
#             #             current_1 = line_split[9][0]
#             #         else:
#             #             current_1 = current_1 + line_split[9][0]
#             #         if len(current_1)!=1 and lst1[idx1][-1] == i:
#             #             larray = lst1[idx1]
#             #             rs = ""
#             #             for idx, item in enumerate(larray):
#             #                 c_i = larray.index(item)
#             #                 rs = rs+current_1[c_i]
#             #             r1[-1].append(rs)
#             #     else:
#             #         # if len(current_1) != 1:
#             #         #     r1[-1].append(current_1)
#             #         # else:
#             #         #     if len(current_1) == 4:
#             #         #         r1.append([current_1])
#             #         current_1 = ""
#             #         if lst1[idx1][-1] < i:
#             #             idx1 = idx1 + 1
#
#             if idx2 < len(lst2):
#                 if i in lst2[idx2]:
#                     if current_2 == "":
#                         current_2 = line_split[9][2]
#                     else:
#                         current_2 = current_2 + line_split[9][2]
#                     if len(current_2)!=1 and lst2[idx2][-1] == i:
#                         larray = lst2[idx2]
#                         rs = ""
#                         prev_item = 0
#                         c_i = 0
#                         for idx, item in enumerate(larray):
#                             if item <= prev_item:
#                                 c_i = larray.index(item)
#                             else:
#                                 prev_item = item
#                         for idx in range(0,4):
#                             rs = rs + current_2[idx]
#                         for idx in range(c_i,len(current_2)):
#                             rs = rs + current_2[idx]
#                         # for idx,item in enumerate(current_2):
#                         #     rs = rs + current_2[idx]
#                         #     if c_i > 3:
#                         #         rs = rs+current_2[c_i]
#                         r1[-1].append(rs)
#                         # r2[-1].append(current_2)
#                 else:
#                     current_2 = ""
#                     # else:
#                     #     if len(current_2) == 4:
#                     #         r2.append([current_2])
#                     if lst2[idx2][-1] < i:
#                         idx2 = idx2 + 1
#         i = i+1
#     print(r2)
#     return r1,r2

txt_f = sys.argv[1]
vcf_f = sys.argv[2]
lst_f = sys.argv[3]

# lst_out = open(sys.argv[3], "w")
# pos = int(sys.argv[3])

def parse_txt(txt_file):
    result = {}
    txt_stream = open(txt_file)
    current_sample = ""
    for line in txt_stream.readlines():
        line_split = line.strip().split("\t")
        if "chrom" in line:
            result[line_split[2]] = []
            current_sample = line_split[2]
        else:
            if len(line_split) < 2:
                continue
            # r = [i for i range(int(line_split[0].split("_")[0])]
            s1 = int(line_split[0].split("_")[0])
            s2 = int(line_split[-2].split("_")[0])
            result[current_sample].append([s1-2,s1-1,s1,s1+1,s2-1,s2,s2+1,s2+2])
    return result

def parse_vcf(vcf_file, lst1, lst2):
    r1 = []
    r2 = []
    vcf_r = []
    vcf_stream = open(vcf_file)
    for line in vcf_stream.readlines():
        if "#" in line:
            continue
        else:
            line_split = line.split("\t")
            vcf_r.append(line_split[9][0:3])
    for item in lst1:
        tmp_s = str(item[0])+"_"+str(item[-1]) + "_"

        for i in item:
            if vcf_r[i - 468012][1] == "/":
                tmp_s = tmp_s + "x"
            else:
                tmp_s = tmp_s + vcf_r[i - 468012][0]
        r1.append(tmp_s)
    for item in lst2:
        tmp_s = str(item[0])+"_"+str(item[-1]) + "_"
        for i in item:
            if vcf_r[i - 468012][1] == "/":
                tmp_s = tmp_s + "x"
            else:
                tmp_s = tmp_s + vcf_r[i - 468012][0]
        r2.append(tmp_s)
    # print(r1)
    return r1,r2


supports = {}

def get_range(pos,end,lst_str,lst):
    for item in lst:
        item_split = item.split("_")
        if int(item_split[1]) < pos or int(item_split[0]) > end:
            continue
        else:
            if item in supports.keys():
                supports[item].append(str(pos)+"_"+lst_str)
            else:
                supports[item] = [str(pos)+"_"+lst_str]

def parser_lst(lst_f, lst1):
    keep = False
    for line in open(lst_f).readlines():
        line_split = line.split(" ")
        lst_count = int(line_split[0])
        if int(line_split[2]) < 468726:
            continue
        if int(line_split[2]) > 504012:
            break
        # print(list(range(2, 2+2*lst_count, 2)))
        for i in range(2, 2+2*lst_count, 2):
            get_range(int(line_split[i]),int(line_split[i]) + len(line_split[i+1]) - 1, line_split[i+1], lst1)
            # for p in range():
            #     if p in range(int(line_split[i]), int(line_split[i]) + len(line_split[i+1])):


def reverse(it_str):
    it_split = it_str.split("_")
    r = it_split[0]+"_"
    for i in it_split[1]:
        if i == '0':
            r = r+'1'
        elif i == '1':
            r = r+'0'
        else:
            r = r+i
    return r

def parse_support(supports):
    result={}
    for key,value in supports.items():
        key_split = key.split("_")
        k1 = str(int(key_split[0])+1)+"_"+key_split[2][1:3]
        k2 = str(int(key_split[1])+1)+"_"+key_split[2][5:7]
        for item in value:
            # item_split = item.split("_")
            r_item = reverse(item)
            if k1 in item or k1 in r_item or k2 in item or k2 in r_item:
                if key in result:
                    result[key] = result[key] + 1
                else:
                    result[key] = 1
    return result


result = parse_txt(txt_f)
print(result)
r1,r2= parse_vcf(vcf_f, result["HG004"], result["HG003"])
parser_lst(lst_f, r1)
rrr = parse_support(supports)
print(rrr)
