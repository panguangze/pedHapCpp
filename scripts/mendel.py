import vcf
import sys
vcf_in = sys.argv[1]
vcf_out = open(sys.argv[2],"w")

def mendel_conflic(v1,v2,v3):
    g1 = v1.split(":")[0]
    g2 = v2.split(":")[0]
    g3 = v3.split(":")[0]

    if ((g1[0] == g2[0] or g1[0] == g2[2]) and (g1[2] == g3[0] or g1[0] == g3[2])) or ((g1[2] == g2[0] or g1[2] == g2[2]) and (g1[0] == g3[0] or g1[0] == g3[2])):
        return False
    return True

def check(v):
    s_v = v.split(":")
    result = []
    if s_v[0] == "./.":
        if len(s_v[5].split(",")) != 3:
            return v
        if int(s_v[2]) <5:
            return v
        pl =  [int(x) for x in s_v[5].split(",")]
        # if len(pl) != 3:
        l = [i for i,val in enumerate(pl) if val==0]
        # maxv = max(pl)
        # maxv_idx = pl.index(maxv)
        # zero_idx = pl.index(0)
        real_gt=""
        for i in range(0,len(pl)):
            if i in l:
            # if i == maxv_idx or i == zero_idx:
            #     continue
                if i == 0:
                    real_gt = "0/0"
                if i == 1:
                    real_gt = "0/1"
                if i == 2:
                    real_gt = "1/1"
                s_v[0] = real_gt
                result.append(":".join(s_v))
    else:
        result.append(":".join(s_v))
    return result

for line in open(vcf_in):
    if "#" in line:
        vcf_out.write(line)
        continue
    splited_v = line.split("\t")

    c = check(splited_v[9])
    f = check(splited_v[10])
    m = check(splited_v[11])

    if len([i for i,val in enumerate([len(c), len(f), len(m)]) if val==2]) >=2:
        vcf_out.write(line)

    if len(c) == 2:
        if (!mendel_conflic(c[0], f[0], m[0]) and !mendel_conflic(c[1], f[0], m[0])) or (mendel_conflic(c[0], f[0], m[0]) and mendel_conflic(c[1], f[0], m[0])):
            vcf_out.write(line)
            # pass
        elif mendel_conflic(c[0], f[0], m[0]):
            splited_v[9] = c[0]
            splited_v[10] = f[0]
            splited_v[11] = m[0]
            line = "\t".join(splited_v)
            vcf_out.write(line)
        elif mendel_conflic(c[1], f[0], m[0]):
            splited_v[9] = c[1]
            splited_v[10] = f[0]
            splited_v[11] = m[0]
            line = "\t".join(splited_v)
            vcf_out.write(line)
        else:
            vcf_out.write(line)

    if len(f) == 2:
        if (!mendel_conflic(c[0], f[0], m[0]) and !mendel_conflic(c[0], f[1], m[0])) or (mendel_conflic(c[0], f[0], m[0]) and mendel_conflic(c[0], f[1], m[0])):
            vcf_out.write(line)
            # pass
        elif mendel_conflic(c[0], f[0], m[0]):
            splited_v[9] = c[0]
            splited_v[10] = f[0]
            splited_v[11] = m[0]
            line = "\t".join(splited_v)
            vcf_out.write(line)
        elif mendel_conflic(c[0], f[1], m[0]):
            splited_v[9] = c[0]
            splited_v[10] = f[1]
            splited_v[11] = m[0]
            line = "\t".join(splited_v)
            vcf_out.write(line)
        else:
            vcf_out.write(line)

    if len(m) == 2:
        if (!mendel_conflic(c[0], f[0], m[0]) and !mendel_conflic(c[0], f[0], m[1])) or (mendel_conflic(c[0], f[0], m[0]) and mendel_conflic(c[0], f[0], m[1])):
            vcf_out.write(line)
            # pass
        elif mendel_conflic(c[0], f[0], m[0]):
            splited_v[9] = c[0]
            splited_v[10] = f[0]
            splited_v[11] = m[0]
            line = "\t".join(splited_v)
            vcf_out.write(line)
        elif mendel_conflic(c[0], f[0], m[1]):
            splited_v[9] = c[0]
            splited_v[10] = f[0]
            splited_v[11] = m[1]
            line = "\t".join(splited_v)
            vcf_out.write(line)
        else:
            vcf_out.write(line)

    # splited_v[9] = check(splited_v[9])
    # splited_v[10] = check(splited_v[10])
    # splited_v[11] = check(splited_v[11])

    # if mendel_conflic(splited_v[9], splited_v[10], splited_v[11]):
    #     line = "\t".join(splited_v)
    #     vcf_out.write(line)
    # else:
    #     # line = "\t".join(splited_v)
    #     vcf_out.write(line)