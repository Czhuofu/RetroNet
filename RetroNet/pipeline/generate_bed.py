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

def is_position_in_bed(chrom, start, end, bed_file):
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            bed_chrom = fields[0]
            bed_start = int(fields[1])
            bed_end = int(fields[2])

            # check the position with the bed file
            if chrom == bed_chrom:
                if set(range(int(start),int(end))) & set(range(bed_start-200,bed_end+200)):
                    print(f'{chrom} {start} {end} is in {bed_file}')
                    return True
    return False

def read_position_in_other_strand(key,di):
    for i in di[key[0]]:
        if set(range(int(key[1]),int(key[2]))) & set(range(int(i[0])-50,int(i[1])+50)):
            if len(i[-1].split(","))>3:
                ###calculate the direction and the number of reads from different direction
                dire=i[-1].split(",")[3]
                plus=0 if '+' not in i[-1] else i[-1].count('+')-1
                minus=0 if '-' not in i[-1] else i[-1].count('-')-1
                return [i[2],i[3],i[4],dire,plus,minus]
            else:
                dire=" "
                return [i[2],i[3],i[4]]
    return []

def detect_sv(total,dire,plus,minus):
    if int(total)<11:
        return " "
    elif not ('+' in dire and '-' in dire):
        return "May be a SV"
    elif plus<minus*0.1 or minus<plus*0.1:
        return "May be a SV"
    else:
        return " "

def detect_unbalance(this_strand,other_strand):
    #print(this_strand,other_strand)
    if len(other_strand)<4:
        return [this_strand[0],this_strand[1],this_strand[2],this_strand[3],this_strand[4],this_strand[5],this_strand[-1]]
    elif int(other_strand[0])<2:
        return [this_strand[0],this_strand[1],this_strand[2],this_strand[3],this_strand[4],this_strand[5],this_strand[-1]]
    else:
        dire_both = this_strand[5]+other_strand[3]
        if '+' in dire_both and '-' in dire_both:
            dire_both = '+/-'
        else:
            dire_both = '+' if '+' in dire_both else '-'
        all_pair_read=int(this_strand[4])+ int(other_strand[2])
        return [this_strand[0],this_strand[1],this_strand[2],this_strand[3],this_strand[4],dire_both,detect_sv(all_pair_read, dire_both,int(other_strand[-2])+int(this_strand[-3]),int(other_strand[-1])+int(this_strand[-2]))]
        ##need to finish the output
    
def check_l1_activity(p1,p2):
    if p1[0] and p2[1] in ["A","N"]:
        if p1[1] in ["C","N"] and p2[0] in ["T","N"]:
            if p1[2] and p2[2] in ["A","G","N"]:
                return True
    return False
    
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
            if TEclass == "LINE":
                if not check_l1_activity(columns[4],columns[5]):
                    continue
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
                    if not (is_position_in_bed(key[0],i[0],i[1],f'{masterpath}/1kgenome/{hg}_1kgenome/{TEclass}.1kgenome.{hg}.bed') or any(un in key[0] for un in ['Un','random','alt','GL'])):
                        print(f"{masterpath}/visual/RetroVis.sh -i {sub} -t {TEclass} -f {TEfamily} -c {key[0]} -d {i[0]} -e {i[1]} -r {ver} -s 0 -p {outpath} -m {masterpath} -g {hg}")
                        os.system(f"{masterpath}/visual/RetroVis.sh -i {sub} -t {TEclass} -f {TEfamily} -c {key[0]} -d {i[0]} -e {i[1]} -r {ver} -s 0 -p {outpath} -m {masterpath} -g {hg}")
                        plus_this_strand=0 if '+' not in i[-1] else i[-1].count('+')-1
                        minus_this_strand=0 if '-' not in i[-1] else i[-1].count('-')-1
                        this_strand=[i[0],i[1],i[2],i[3],i[4],dire,plus_this_strand,minus_this_strand,detect_sv(i[4],dire,plus_this_strand,minus_this_strand)]
                        output=detect_unbalance(this_strand,read_position_in_other_strand(key,di1))
                        file2.write(key[0]+ "\t"+ '\t'.join(output[:5])+ "\t" + key[3]+ "\t"+ output[5]+ "\t" + f'{max(value)}/{np.median(value)}' + "\t" + output[6] + "\n")
        if key[3] == "strand1":
            for i in di1[key[0]]:
                if int(i[0]) == int(key[1]):
                    if len(i[-1].split(","))>3:
                        dire=i[-1].split(",")[3]
                    else:
                        dire=" "
                    #if int(i[2]) < 16:
                    if not (is_position_in_bed(key[0],i[0],i[1],f'{masterpath}/1kgenome/{hg}_1kgenome/{TEclass}.1kgenome.{hg}.bed') or any(un in key[0] for un in ['Un','random','alt','GL'])):
                        os.system(f"{masterpath}/visual/RetroVis.sh -i {sub} -t {TEclass} -f {TEfamily} -c {key[0]} -d {i[0]} -e {i[1]} -r {ver} -s 1 -p {outpath} -m {masterpath} -g {hg}")
                        print(f"{masterpath}/visual/RetroVis.sh -i {sub} -t {TEclass} -f {TEfamily} -c {key[0]} -d {i[0]} -e {i[1]} -r {ver} -s 1 -p {outpath} -m {masterpath} -g {hg}")
                        plus_this_strand=0 if '+' not in i[-1] else i[-1].count('+')-1
                        minus_this_strand=0 if '-' not in i[-1] else i[-1].count('-')-1
                        this_strand=[i[0],i[1],i[2],i[3],i[4],dire,plus_this_strand,minus_this_strand,detect_sv(i[4],dire,plus_this_strand,minus_this_strand)]
                        output=detect_unbalance(this_strand,read_position_in_other_strand(key,di0))
                        file2.write(key[0]+ "\t"+ '\t'.join(output[:5])+ "\t" + key[3]+ "\t"+ output[5]+ "\t" + f'{max(value)}/{np.median(value)}' + "\t" + output[6] + "\n")
                        
with open (f'{outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.head',"w") as file3:
    file3.write("## In column insert_strand, strand 1 means that MEI occurs in plus strand; strand 0 means that MEI occurs in minus strand;"+'\n')
    file3.write("## In column support_read strand, + represent that supporting read comes from upstream in the human plus strand, - represent that supporting read comes from downstream in the human plus strand;"+'\n')
    file3.write("## In column probability, the first number is the highest probability and the second number is the median probability." + "\n")
    file3.write("\t".join(["chrom","start_pos","stop_pos","num_support_read","num_split_read","num_pair_read","insert_strand","support_read_direction","Probability","annotation"])+"\n")

os.system(f'sort -k1,1V -k2,2n -o {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed.temp {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed.temp')
os.system(f'bedtools merge -i {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed.temp -c 4,5,6,7,8,9,10 -o sum,sum,sum,collapse,collapse,collapse,collapse -delim "|" > {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed.temp1')
os.system(f'cat {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.head {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed.temp1 > {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed')
os.system(f'rm {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.bed.temp*')
os.system(f'rm {outpath}/{sub}/retro_v{ver}/{sub}.{TEclass}.head')
