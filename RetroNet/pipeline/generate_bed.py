###In retrov3 version
import numpy as np
import os
import sys

def read_call(file_name):
    dic={}
    with open(file_name,"r") as file:
        for line in file:
            columns = line.strip().split()
            column_name=columns[0].split("_")[0]
            if column_name in dic:
                dic[column_name].append(columns[1:3]+columns[4:8]) 
            else:
                dic[column_name] = [columns[1:3]+columns[4:8]]
    return dic

## read the parameter
outpath=sys.argv[1]
sub=sys.argv[2]
TEclass=sys.argv[3]
ver=sys.argv[4]
masterpath=sys.argv[5]
cutoff=sys.argv[6]
hg=sys.argv[7]
##path of the calls 

di0=read_call(f"{outpath}/{sub}/retro_v{ver}_0/{TEclass}/{sub}.{TEclass}.SR.PE.calls")
di1=read_call(f"{outpath}/{sub}/retro_v{ver}_1/{TEclass}/{sub}.{TEclass}.SR.PE.calls")


##gain from bed

uniq_dict={}
with open ("{}/{}_Inspected_{}_cut{}.txt".format(outpath + "/" + sub + "/RetroNet",TEclass,sub,cutoff),"r") as file:
    for line in file:
        columns = line.strip().split()
        if float(columns[1]) > float(cutoff) and columns[2] == '0' and columns[3] == 'PASS' and "green" not in columns[-1]:
            ##read the position
            read_pos_all=columns[0].split("_")
            read_pos=tuple([read_pos_all[-3],int(read_pos_all[-2]),int(read_pos_all[-2])+1,read_pos_all[-5]])
            if tuple(read_pos) in uniq_dict:
                uniq_dict[read_pos].append(float(columns[1])) 
            else:
                uniq_dict[read_pos] = [float(columns[1])]
#print(uniq_dict)


if TEclass == "LINE":
    TEfamily = "L1HS"
elif TEclass == "ALU":
    TEfamily = "AluYa5"
elif TEclass == "SVA":
    TEfamily = "SVA_E"


##mkdir folder 
os.system(f"mkdir {outpath}/{sub}/retro_v{ver}")
with open(f'{outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed.temp',"w") as file2:
    for key,value in uniq_dict.items():
        if key[3] == "strand0":
            for i in di0[key[0]]:
                if int(i[0]) == int(key[1]):
                    if len(i[-1].split(","))>3:
                        dire=i[-1].split(",")[3]
                    else:
                        dire=" "
                    #if int(i[2]) < 16:
                    os.system(f"{masterpath}/visual/RetroVis.sh -i {sub} -t {TEclass} -f {TEfamily} -c {key[0]} -d {i[0]} -e {i[1]} -r {ver} -s 0 -p {outpath} -m {masterpath} -g {hg}")
                    file2.write(key[0]+ "\t"+ '\t'.join(i[:-1])+ "\t" + key[3]+ "\t"+ dire+ "\t" + f'{max(value)}/{np.median(value)}' +"\n")
        if key[3] == "strand1":
            for i in di1[key[0]]:
                if int(i[0]) == int(key[1]):
                    if len(i[-1].split(","))>3:
                        dire=i[-1].split(",")[3]
                    else:
                        dire=" "
                    #if int(i[2]) < 16:
                    os.system(f"{masterpath}/visual/RetroVis.sh -i {sub} -t {TEclass} -f {TEfamily} -c {key[0]} -d {i[0]} -e {i[1]} -r {ver} -s 1 -p {outpath} -m {masterpath} -g {hg}")
                    file2.write(key[0]+ "\t"+ '\t'.join(i[:-1])+ "\t" + key[3]+ "\t"+ dire+ "\t" + f'{max(value)}/{np.median(value)}' +"\n")
        
with open (f'{outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.head',"w") as file3:
    file3.write("## In column insert_strand, strand 1 means that MEI occurs in plus strand; strand 0 means that MEI occurs in minus strand;"+'\n')
    file3.write("## In column support_read strand, + represent that supporting read comes from upstream in the human plus strand, - represent that supporting read comes from downstream in the human plus strand;"+'\n')
    file3.write("## In column probability, the first number is the highest probability and the second number is the median probability." + "\n")
    file3.write("\t".join(["chrom","start_pos","stop_pos","num_support_read","num_split_read","num_pair_read","insert_strand","support_read_direction","Probability"])+"\n")

os.system(f'sort -k1,1V -k2,2n -o {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed.temp {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed.temp')
os.system(f'cat {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.head {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed.temp > {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed')
os.system(f'rm {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed.temp')
os.system(f'rm {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.head')
