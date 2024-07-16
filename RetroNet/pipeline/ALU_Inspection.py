import os
import sys
import torch
import torchvision
from PIL import Image

Totensor=torchvision.transforms.Compose([torchvision.transforms.ToTensor()])

# Initializing #
DNA_base = "ACTG"
cutoff = 0.5
outpath = "."
sub = "test"
ver = 4
TEclass = "Alu"

# Create upstream_arrow #
upstream_arrow = torch.zeros(3,5,3)
upstream_arrow[2,:,:] = 1
upstream_arrow[0:2,0,1:3] = 1
upstream_arrow[0:2,1,2] = 1
upstream_arrow[0:2,4,1:3] = 1
upstream_arrow[0:2,3,2] = 1
# Create downstream_arrow
downstream_arrow = torch.zeros(3,5,3)
downstream_arrow[2,:,:] = 1
downstream_arrow[0:2,0,0:2] = 1
downstream_arrow[0:2,1,0] = 1
downstream_arrow[0:2,4,0:2] = 1
downstream_arrow[0:2,3,0] = 1

img = sys.argv[1]
ver = sys.argv[2]

# Read Basic Infomation #
temp_tensor=Totensor(Image.open(img).convert('RGB'))
Detect_blue = torch.sum(temp_tensor[1:3,:,:],dim=0)
Detect_red = torch.sum(temp_tensor[0:2,:,:],dim=0)
Detect_black = torch.sum(temp_tensor,dim=0)

## Detect Arrow ##
reads_direction = []
arrow_row = []
arrow_min = []
arrow_max = []
clip_map = []
clip_min = []
clip_max = []
for j in [10,38,66,94]:
    arrow_detect = Detect_blue[j:j+5,:].eq(1).nonzero(as_tuple=True)[1]
    if arrow_detect.numel() != 0:
        arrow_min.append(min(arrow_detect))
        arrow_max.append(max(arrow_detect))
        arrow_row.append(j)
        for k in range(min(arrow_detect),max(arrow_detect)+1):
            if temp_tensor[:,j:j+5,k:k+3].equal(downstream_arrow):
                reads_direction.append("downstream")
                break
            if temp_tensor[:,j:j+5,k:k+3].equal(upstream_arrow):
                reads_direction.append("upstream")
                break
        clip_detect = Detect_black[j:j+5,:].eq(0).nonzero(as_tuple=True)[1]
        if clip_detect.numel() != 0:
            clip_min.append(min(clip_detect))
            clip_max.append(max(clip_detect))
            clip_map_detect = Detect_red[j:j+5,:].eq(1).nonzero(as_tuple=True)[1]
            if clip_map_detect.numel() != 0:
                clip_map.append("mapped")
            else:
                clip_map.append("N")
        else:
            clip_map.append("N")
            clip_min.append("N")
            clip_max.append("N")

if reads_direction.count("N") == 1:
    N_index = reads_direction.index("N")
    arrow_min[N_index] = 9999 - 150
    arrow_max[N_index] = 9999

align_min = []
align_max = []
align_color = []
align_mapability = []
transduction3_color = []
transduction3_min = []
transduction3_max = []
transduction5_color = []
transduction5_min = []
transduction5_max = []
#
for k in [17,45,73,101]:
    align_min_tmp = []
    align_max_tmp = []
    align_color_tmp = []
    align_mapability_tmp = []
    transduction3_color_tmp = []
    transduction3_min_tmp = []
    transduction3_max_tmp = []
    transduction5_color_tmp = []
    transduction5_min_tmp = []
    transduction5_max_tmp = []
    for m in [0,1,2,3]:
        red_detect = Detect_black[k+m*5:k+m*5+4,:].eq(1).nonzero(as_tuple=True)[1]
        magenta_detect = Detect_black[k+m*5:k+m*5+4,:].eq(2).nonzero(as_tuple=True)[1]
        if red_detect.numel() != 0:
            black_count = Detect_black[k+m*5:k+m*5+4,min(red_detect):max(red_detect)].eq(0).nonzero(as_tuple=True)[1].numel()
            align_mapability_tmp.append((1 - black_count/red_detect.numel()))
            align_color_tmp.append("red")
            align_min_tmp.append(min(red_detect))
            align_max_tmp.append(max(red_detect))
            green_detect = Detect_black[k+m*5-1:k+m*5,max(red_detect)+1:].eq(1).nonzero(as_tuple=True)[1]
            black_detect = Detect_black[k+m*5-1:k+m*5,max(red_detect)+1:].eq(0).nonzero(as_tuple=True)[1]
            if green_detect.numel() != 0:
                transduction3_color_tmp.append("green")
                transduction3_min_tmp.append(min(green_detect) + max(red_detect))
                transduction3_max_tmp.append(max(green_detect) + max(red_detect))
            elif black_detect.numel() != 0:
                transduction3_color_tmp.append("black")
                transduction3_min_tmp.append(min(black_detect) + max(red_detect))
                transduction3_max_tmp.append(max(black_detect) + max(red_detect))
            else:
                transduction3_color_tmp.append("N")
                transduction3_min_tmp.append("N")
                transduction3_max_tmp.append("N")
            green_detect = Detect_black[k+m*5-1:k+m*5,:min(red_detect)+1].eq(1).nonzero(as_tuple=True)[1]
            black_detect = Detect_black[k+m*5-1:k+m*5,:min(red_detect)+1].eq(0).nonzero(as_tuple=True)[1]
            if green_detect.numel() != 0:
                transduction5_color_tmp.append("green")
                transduction5_min_tmp.append(min(green_detect))
                transduction5_max_tmp.append(max(green_detect))
            elif black_detect.numel() != 0:
                transduction5_color_tmp.append("black")
                transduction5_min_tmp.append(min(black_detect))
                transduction5_max_tmp.append(max(black_detect))
            else:
                transduction5_color_tmp.append("N")
                transduction5_min_tmp.append("N")
                transduction5_max_tmp.append("N")
        elif magenta_detect.numel() != 0:
            black_count = Detect_black[k+m*5:k+m*5+4,min(magenta_detect):max(magenta_detect)].eq(0).nonzero(as_tuple=True)[1].numel()
            align_mapability_tmp.append((1 - black_count/magenta_detect.numel()))
            align_color_tmp.append("magenta")
            align_min_tmp.append(min(magenta_detect))
            align_max_tmp.append(max(magenta_detect))
            green_detect = Detect_black[k+m*5-1:k+m*5,max(magenta_detect)+1:].eq(1).nonzero(as_tuple=True)[1]
            black_detect = Detect_black[k+m*5-1:k+m*5,max(magenta_detect)+1:].eq(0).nonzero(as_tuple=True)[1]
            if green_detect.numel() != 0:
                transduction3_color_tmp.append("green")
                transduction3_min_tmp.append(min(green_detect) + max(magenta_detect))
                transduction3_max_tmp.append(max(green_detect) + max(magenta_detect))
            elif black_detect.numel() != 0:
                transduction3_color_tmp.append("black")
                transduction3_min_tmp.append(min(black_detect) + max(magenta_detect))
                transduction3_max_tmp.append(max(black_detect) + max(magenta_detect))
            else:
                transduction3_color_tmp.append("N")
                transduction3_min_tmp.append("N")
                transduction3_max_tmp.append("N")
            green_detect = Detect_black[k+m*5-1:k+m*5,:min(magenta_detect)+1].eq(1).nonzero(as_tuple=True)[1]
            black_detect = Detect_black[k+m*5-1:k+m*5,:min(magenta_detect)+1].eq(0).nonzero(as_tuple=True)[1]
            if green_detect.numel() != 0:
                transduction5_color_tmp.append("green")
                transduction5_min_tmp.append(min(green_detect))
                transduction5_max_tmp.append(max(green_detect))
            elif black_detect.numel() != 0:
                transduction5_color_tmp.append("black")
                transduction5_min_tmp.append(min(black_detect))
                transduction5_max_tmp.append(max(black_detect))
            else:
                transduction5_color_tmp.append("N")
                transduction5_min_tmp.append("N")
                transduction5_max_tmp.append("N")
        else:
            align_mapability_tmp.append("N")
            align_color_tmp.append("N")
            align_min_tmp.append("N")
            align_max_tmp.append("N")
            transduction3_color_tmp.append("N")
            transduction3_min_tmp.append("N")
            transduction3_max_tmp.append("N")
            transduction5_color_tmp.append("N")
            transduction5_min_tmp.append("N")
            transduction5_max_tmp.append("N")
    align_mapability.append(align_mapability_tmp)
    align_color.append(align_color_tmp)
    align_min.append(align_min_tmp)
    align_max.append(align_max_tmp)
    transduction3_color.append(transduction3_color_tmp)
    transduction3_min.append(transduction3_min_tmp)
    transduction3_max.append(transduction3_max_tmp)
    transduction5_color.append(transduction5_color_tmp)
    transduction5_min.append(transduction5_min_tmp)
    transduction5_max.append(transduction5_max_tmp)

Inspection = ""
## Manual Inspection ##
# 1 # Upstream + Downstream
if len(reads_direction) != 2:
    Inspection = Inspection + "miss arrow; "
else:
    if reads_direction[0] != reads_direction[1]:
        # 1.1 ##  ->     ##
        #     ##      <- ##
        if reads_direction[0] == "upstream":
            if int(arrow_min[1]) - int(arrow_max[0]) > 600:
                Inspection = Inspection + "Gap between two arrow is larger than 600bp; "
            if clip_map[0] == "mapped":
                # 1.1.1 ## ->-     ##
                #       ##     -<- ##
                if clip_map[1] == "mapped":
                    # (1) Arrow Conflict #
                    if int(clip_max[0]) < int(arrow_max[0]) or int(clip_min[1]) > int(arrow_min[1]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    if (int(arrow_min[0]) < int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_min[1])) or (int(arrow_min[0]) > int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_max[1])):
                        if (int(arrow_max[0]) - int(arrow_min[1])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[0]) > int(arrow_max[1]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict #
                    if align_color[1].count("N") != 4:
                        # Two Clip-pair #
                        if align_color[3].count("N") != 4:
                            if (align_max[3][0] == "N" or int(align_max[3][0]) != 441) and (align_max[3][1] == "N" or int(align_max[3][1]) != 448) and (align_max[3][2] == "N" or int(align_max[3][2]) != 441) and (align_max[3][3] == "N" or int(align_max[3][3]) != 441):
                                Inspection = Inspection + "No PolyA; "
                            junc_5 = 0
                            junc_3 = 0
                            for e in [0,1,2,3]:
                                if align_min[0][e] != "N" and align_min[3][e] != "N" and int(align_min[0][e]) - int(align_min[3][e]) > junc_5:
                                    junc_5 = int(align_min[0][e]) - int(align_min[3][e])
                                if align_max[0][e] != "N" and align_max[3][e] != "N" and int(align_max[0][e]) - int(align_max[3][e]) > junc_3:
                                    junc_3 = int(align_max[0][e]) - int(align_max[3][e])
                            if junc_5 > 5 or junc_3 > 5:
                                Inspection = Inspection + "Junction site difference; "
                            if transduction5_color[2].count("N") != 4:
                                strange_align_count = 0
                                for e in [0,1,2,3]:
                                    if transduction5_min[2][e] != "N" and align_min[1][e] != "N" and abs(int(transduction5_min[2][e]) - int(transduction5_max[2][e])) > 10 and (int(transduction5_max[2][e]) - 10) > int(align_min[1][e]):
                                        strange_align_count = strange_align_count + 1
                                if strange_align_count >= 2:
                                    Inspection = Inspection + "Unknown clip within Alu Insertion; "
                            if transduction3_color[1].count("N") != 4:
                                strange_align_count = 0
                                for e in [0,1,2,3]:
                                    if transduction3_min[1][e] != "N" and align_max[2][e] != "N" and abs(int(transduction3_min[1][e]) - int(transduction3_max[1][e])) > 10 and (int(transduction3_min[1][e]) + 10) < int(align_max[2][e]):
                                        strange_align_count = strange_align_count + 1
                                if strange_align_count >= 2:
                                    Inspection = Inspection + "Unknown clip within Alu Insertion; "
                        # Clip-pair + Split #
                        else:
                            if (align_max[2][0] == "N" or int(align_max[2][0]) != 441) and (align_max[2][1] == "N" or int(align_max[2][1]) != 448) and (align_max[2][2] == "N" or int(align_max[2][2]) != 441) and (align_max[2][3] == "N" or int(align_max[2][3]) != 441):
                                Inspection = Inspection + "No PolyA; "
                            junc_5 = 0
                            junc_3 = 0
                            for e in [0,1,2,3]:
                                if align_min[0][e] != "N" and align_min[2][e] != "N" and int(align_min[0][e]) - int(align_min[2][e]) > junc_5:
                                    junc_5 = int(align_min[0][e]) - int(align_min[2][e])
                                if align_max[0][e] != "N" and align_max[2][e] != "N" and int(align_max[0][e]) - int(align_max[2][e]) > junc_3:
                                    junc_3 = int(align_max[0][e]) - int(align_max[2][e])
                            if junc_5 > 5 or junc_3 > 5:
                                Inspection = Inspection + "Junction site difference; "
                    else:
                        # Split + Clip-pair #
                        if align_color[3].count("N") != 4:
                            if (align_max[3][0] == "N" or int(align_max[3][0]) != 441) and (align_max[3][1] == "N" or int(align_max[3][1]) != 448) and (align_max[3][2] == "N" or int(align_max[3][2]) != 441) and (align_max[3][3] == "N" or int(align_max[3][3]) != 441):
                                Inspection = Inspection + "No PolyA; "
                            junc_5 = 0
                            junc_3 = 0
                            for e in [0,1,2,3]:
                                if align_min[0][e] != "N" and align_min[3][e] != "N" and int(align_min[0][e]) - int(align_min[3][e]) > junc_5:
                                    junc_5 = int(align_min[0][e]) - int(align_min[3][e])
                                if align_max[0][e] != "N" and align_max[3][e] != "N" and int(align_max[0][e]) - int(align_max[3][e]) > junc_3:
                                    junc_3 = int(align_max[0][e]) - int(align_max[3][e])
                            if junc_5 > 5 or junc_3 > 5:
                                Inspection = Inspection + "Junction site difference; "
                        # Split + Split #
                        else:
                            if (align_max[2][0] == "N" or int(align_max[2][0]) != 441) and (align_max[2][1] == "N" or int(align_max[2][1]) != 448) and (align_max[2][2] == "N" or int(align_max[2][2]) != 441) and (align_max[2][3] == "N" or int(align_max[2][3]) != 441):
                                Inspection = Inspection + "No PolyA; "
                            junc_5 = 0
                            junc_3 = 0
                            for e in [0,1,2,3]:
                                if align_min[0][e] != "N" and align_min[2][e] != "N" and int(align_min[0][e]) - int(align_min[2][e]) > junc_5:
                                    junc_5 = int(align_min[0][e]) - int(align_min[2][e])
                                if align_max[0][e] != "N" and align_max[2][e] != "N" and int(align_max[0][e]) - int(align_max[2][e]) > junc_3:
                                    junc_3 = int(align_max[0][e]) - int(align_max[2][e])
                            if junc_5 > 5 or junc_3 > 5:
                                Inspection = Inspection + "Junction site difference; "
                # 1.1.2 ## ->-     ##
                #       ##     <-- ##
                else:
                    # (1) Arrow Conflict #
                    if int(clip_max[0]) < int(arrow_max[0]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    if (int(arrow_min[0]) < int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_min[1])) or (int(arrow_min[0]) > int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_max[1])):
                        if (int(arrow_max[0]) - int(arrow_min[1])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[0]) > int(arrow_max[1]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict #
                    junc_5 = 0
                    for e in [0,1,2,3]:
                        if align_min[0][e] != "N" and align_min[2][e] != "N" and int(align_min[0][e]) - int(align_min[2][e]) > junc_5:
                            junc_5 = int(align_min[0][e]) - int(align_min[2][e])
                    if junc_5 > 5:
                        Inspection = Inspection + "Junction site difference; "
                    # Clip-Pair + Pair
                    if align_color[1].count("N") != 4:
                        if transduction5_color[2].count("N") != 4:
                            strange_align_count = 0
                            for e in [0,1,2,3]:
                                if transduction5_min[2][e] != "N" and align_min[1][e] != "N" and abs(int(transduction5_min[2][e]) - int(transduction5_max[2][e])) > 10 and (int(transduction5_max[2][e]) - 10) > int(align_min[1][e]):
                                    strange_align_count = strange_align_count + 1
                            if strange_align_count >= 2:
                                Inspection = Inspection + "Unknown clip within Alu Insertion; "
                        if transduction3_color[1].count("N") != 4:
                            strange_align_count = 0
                            for e in [0,1,2,3]:
                                if transduction3_min[1][e] != "N" and align_max[2][e] != "N" and abs(int(transduction3_min[1][e]) - int(transduction3_max[1][e])) > 10 and (int(transduction3_min[1][e]) + 10) < int(align_max[2][e]):
                                    strange_align_count = strange_align_count + 1
                            if strange_align_count >= 2:
                                Inspection = Inspection + "Unknown clip within Alu Insertion; "
            else:
                # 1.1.3 ## -->     ##
                #       ##     -<- ##
                if clip_map[1] == "mapped":
                    # (1) Arrow Conflict #
                    if int(clip_min[1]) > int(arrow_min[1]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    if (int(arrow_min[0]) < int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_min[1])) or (int(arrow_min[0]) > int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_max[1])):
                        if (int(arrow_max[0]) - int(arrow_min[1])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[0]) > int(arrow_max[1]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict #
                    if align_color[3].count("N") != 4:
                        # Pair + Clip-pair #
                        if (align_max[3][0] == "N" or int(align_max[3][0]) != 441) and (align_max[3][1] == "N" or int(align_max[3][1]) != 448) and (align_max[3][2] == "N" or int(align_max[3][2]) != 441) and (align_max[3][3] == "N" or int(align_max[3][3]) != 441):
                            Inspection = Inspection + "No PolyA; "
                        junc_3 = 0
                        for e in [0,1,2,3]:
                            if align_max[1][e] != "N" and align_max[3][e] != "N" and int(align_max[1][e]) - int(align_max[3][e]) > junc_3:
                                junc_3 = int(align_max[1][e]) - int(align_max[3][e])
                        if junc_3 > 5:
                            Inspection = Inspection + "Junction site difference; "
                        if transduction5_color[2].count("N") != 4:
                            strange_align_count = 0
                            for e in [0,1,2,3]:
                                if transduction5_min[2][e] != "N" and align_min[1][e] != "N" and abs(int(transduction5_min[2][e]) - int(transduction5_max[2][e])) > 10 and (int(transduction5_max[2][e]) - 10) > int(align_min[1][e]):
                                    strange_align_count = strange_align_count + 1
                            if strange_align_count >= 2:
                                Inspection = Inspection + "Unknown clip within Alu Insertion; "
                        if transduction3_color[1].count("N") != 4:
                            strange_align_count = 0
                            for e in [0,1,2,3]:
                                if transduction3_min[1][e] != "N" and align_max[2][e] != "N" and abs(int(transduction3_min[1][e]) - int(transduction3_max[1][e])) > 10 and (int(transduction3_min[1][e]) + 10) < int(align_max[2][e]):
                                    strange_align_count = strange_align_count + 1
                            if strange_align_count >= 2:
                                Inspection = Inspection + "Unknown clip within Alu Insertion; "
                    else:
                        # Pair + Split #
                        if (align_max[2][0] == "N" or int(align_max[2][0]) != 441) and (align_max[2][1] == "N" or int(align_max[2][1]) != 448) and (align_max[2][2] == "N" or int(align_max[2][2]) != 441) and (align_max[2][3] == "N" or int(align_max[2][3]) != 441):
                            Inspection = Inspection + "No PolyA; "
                        junc_3 = 0
                        for e in [0,1,2,3]:
                            if align_max[1][e] != "N" and align_max[2][e] != "N" and int(align_max[1][e]) - int(align_max[2][e]) > junc_3:
                                junc_3 = int(align_max[1][e]) - int(align_max[2][e])
                        if junc_3 > 5:
                            Inspection = Inspection + "Junction site difference; "
                # 1.1.4 ## -->     ##
                #       ##     <-- ##
                        ## Two Pair #
                else:
                    # (1) Arrow Conflict #
                    if (int(arrow_min[0]) < int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_min[1])) or (int(arrow_min[0]) > int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_max[1])):
                        if (int(arrow_max[0]) - int(arrow_min[1])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[0]) > int(arrow_max[1]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict
                    if transduction5_color[2].count("N") != 4:
                        strange_align_count = 0
                        for e in [0,1,2,3]:
                            if transduction5_min[2][e] != "N" and align_min[1][e] != "N" and abs(int(transduction5_min[2][e]) - int(transduction5_max[2][e])) > 10 and (int(transduction5_max[2][e]) - 10) > int(align_min[1][e]):
                                strange_align_count = strange_align_count + 1
                        if strange_align_count >= 2:
                            Inspection = Inspection + "Unknown clip within Alu Insertion; "
                    if transduction3_color[1].count("N") != 4:
                        strange_align_count = 0
                        for e in [0,1,2,3]:
                            if transduction3_min[1][e] != "N" and align_max[2][e] != "N" and abs(int(transduction3_min[1][e]) - int(transduction3_max[1][e])) > 10 and (int(transduction3_min[1][e]) + 10) < int(align_max[2][e]):
                                strange_align_count = strange_align_count + 1
                        if strange_align_count >= 2:
                            Inspection = Inspection + "Unknown clip within Alu Insertion; "
        # 1.2 ##      <- ##
        #     ##  ->     ##
        else:
            if int(arrow_min[1]) - int(arrow_max[0]) > 600:
                Inspection = Inspection + "Gap between two arrow is larger than 600bp; "
            if clip_map[0] == "mapped":
                # 1.2.1 ##     -<- ##
                #       ## ->-     ##
                if clip_map[1] == "mapped":
                    if int(clip_max[1]) < int(arrow_max[1]) or int(clip_min[0]) > int(arrow_min[0]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    if (int(arrow_min[1]) < int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_min[0])) or (int(arrow_min[1]) > int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_max[0])):
                        if (int(arrow_max[1]) - int(arrow_min[0])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[1]) > int(arrow_max[0]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict #
                    if align_color[3].count("N") != 4:
                        # Two Clip-pair #
                        if align_color[1].count("N") != 4:
                            if (align_max[1][0] == "N" or int(align_max[1][0]) != 441) and (align_max[1][1] == "N" or int(align_max[1][1]) != 448) and (align_max[1][2] == "N" or int(align_max[1][2]) != 441) and (align_max[1][3] == "N" or int(align_max[1][3]) != 441):
                                Inspection = Inspection + "No PolyA; "
                            junc_5 = 0
                            junc_3 = 0
                            for e in [0,1,2,3]:
                                if align_min[2][e] != "N" and align_min[1][e] != "N" and int(align_min[2][e]) - int(align_min[1][e]) > junc_5:
                                    junc_5 = int(align_min[2][e]) - int(align_min[1][e])
                                if align_max[2][e] != "N" and align_max[1][e] != "N" and int(align_max[2][e]) - int(align_max[1][e]) > junc_3:
                                    junc_3 = int(align_max[2][e]) - int(align_max[1][e])
                            if junc_5 > 5 or junc_3 > 5:
                                Inspection = Inspection + "Junction site difference; "
                            if transduction5_color[0].count("N") != 4:
                                strange_align_count = 0
                                for e in [0,1,2,3]:
                                    if transduction5_min[0][e] != "N" and align_min[3][e] != "N" and abs(int(transduction5_min[0][e]) - int(transduction5_max[0][e])) > 10 and (int(transduction5_max[0][e]) - 10) > int(align_min[3][e]):
                                        strange_align_count = strange_align_count + 1
                                if strange_align_count >= 2:
                                    Inspection = Inspection + "Unknown clip within Alu Insertion; "
                            if transduction3_color[3].count("N") != 4:
                                strange_align_count = 0
                                for e in [0,1,2,3]:
                                    if transduction3_min[3][e] != "N" and align_max[0][e] != "N" and abs(int(transduction3_min[3][e]) - int(transduction3_max[3][e])) > 10 and (int(transduction3_min[3][e]) + 10) < int(align_max[0][e]):
                                        strange_align_count = strange_align_count + 1
                                if strange_align_count >= 2:
                                    Inspection = Inspection + "Unknown clip within Alu Insertion; "
                        # Clip-pair + Split #
                        else:
                            if (align_max[0][0] == "N" or int(align_max[0][0]) != 441) and (align_max[0][1] == "N" or int(align_max[0][1]) != 448) and (align_max[0][2] == "N" or int(align_max[0][2]) != 441) and (align_max[0][3] == "N" or int(align_max[0][3]) != 441):
                                Inspection = Inspection + "No PolyA; "
                            junc_5 = 0
                            junc_3 = 0
                            for e in [0,1,2,3]:
                                if align_min[2][e] != "N" and align_min[0][e] != "N" and int(align_min[2][e]) - int(align_min[0][e]) > junc_5:
                                    junc_5 = int(align_min[2][e]) - int(align_min[0][e])
                                if align_max[2][e] != "N" and align_max[0][e] != "N" and int(align_max[2][e]) - int(align_max[0][e]) > junc_3:
                                    junc_3 = int(align_max[2][e]) - int(align_max[0][e])
                            if junc_5 > 5 or junc_3 > 5:
                                Inspection = Inspection + "Junction site difference; "
                    else:
                        # Split + Clip-pair #
                        if align_color[1].count("N") != 4:
                            if (align_max[1][0] == "N" or int(align_max[1][0]) != 441) and (align_max[1][1] == "N" or int(align_max[1][1]) != 448) and (align_max[1][2] == "N" or int(align_max[1][2]) != 441) and (align_max[1][3] == "N" or int(align_max[1][3]) != 441):
                                Inspection = Inspection + "No PolyA; "
                            junc_5 = 0
                            junc_3 = 0
                            for e in [0,1,2,3]:
                                if align_min[2][e] != "N" and align_min[1][e] != "N" and int(align_min[2][e]) - int(align_min[1][e]) > junc_5:
                                    junc_5 = int(align_min[2][e]) - int(align_min[1][e])
                                if align_max[2][e] != "N" and align_max[1][e] != "N" and int(align_max[2][e]) - int(align_max[1][e]) > junc_3:
                                    junc_3 = int(align_max[2][e]) - int(align_max[1][e])
                            if junc_5 > 5 or junc_3 > 5:
                                Inspection = Inspection + "Junction site difference; "
                        # Split + Split #
                        else:
                            if (align_max[0][0] == "N" or int(align_max[0][0]) != 441) and (align_max[0][1] == "N" or int(align_max[0][1]) != 448) and (align_max[0][2] == "N" or int(align_max[0][2]) != 441) and (align_max[0][3] == "N" or int(align_max[2][3]) != 441):
                                Inspection = Inspection + "No PolyA; "
                            junc_5 = 0
                            junc_3 = 0
                            for e in [0,1,2,3]:
                                if align_min[2][e] != "N" and align_min[0][e] != "N" and int(align_min[2][e]) - int(align_min[0][e]) > junc_5:
                                    junc_5 = int(align_min[2][e]) - int(align_min[0][e])
                                if align_max[2][e] != "N" and align_max[0][e] != "N" and int(align_max[2][e]) - int(align_max[0][e]) > junc_3:
                                    junc_3 = int(align_max[2][e]) - int(align_max[0][e])
                            if junc_5 > 5 or junc_3 > 5:
                                Inspection = Inspection + "Junction site difference; "
                # 1.2.2 ##     -<- ##
                #       ## -->     ##
                else:
                    if int(clip_min[0]) > int(arrow_min[0]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    if (int(arrow_min[1]) < int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_min[0])) or (int(arrow_min[1]) > int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_max[0])):
                        if (int(arrow_max[1]) - int(arrow_min[0])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[1]) > int(arrow_max[0]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict #
                    if align_color[1].count("N") != 4:
                        # Pair + Clip-pair #
                        if (align_max[1][0] == "N" or int(align_max[1][0]) != 441) and (align_max[1][1] == "N" or int(align_max[1][1]) != 448) and (align_max[1][2] == "N" or int(align_max[1][2]) != 441) and (align_max[1][3] == "N" or int(align_max[1][3]) != 441):
                            Inspection = Inspection + "No PolyA; "
                        junc_5 = 0
                        for e in [0,1,2,3]:
                            if align_min[3][e] != "N" and align_min[1][e] != "N" and int(align_min[3][e]) - int(align_min[1][e]) > junc_5:
                                junc_5 = int(align_min[3][e]) - int(align_min[1][e])
                        if junc_5 > 5:
                            Inspection = Inspection + "Junction site difference; "
                        if transduction5_color[0].count("N") != 4:
                            strange_align_count = 0
                            for e in [0,1,2,3]:
                                if transduction5_min[0][e] != "N" and align_min[3][e] != "N" and abs(int(transduction5_min[0][e]) - int(transduction5_max[0][e])) > 10 and (int(transduction5_max[0][e]) - 10) > int(align_min[3][e]):
                                    strange_align_count = strange_align_count + 1
                            if strange_align_count >= 2:
                                Inspection = Inspection + "Unknown clip within Alu Insertion; "
                        if transduction3_color[3].count("N") != 4:
                            strange_align_count = 0
                            for e in [0,1,2,3]:
                                if transduction3_min[3][e] != "N" and align_min[0][e] != "N" and abs(int(transduction3_min[3][e]) - int(transduction3_max[3][e])) > 10 and (int(transduction3_min[3][e]) + 10) < int(align_max[0][e]):
                                    strange_align_count = strange_align_count + 1
                            if strange_align_count >= 2:
                                Inspection = Inspection + "Unknown clip within Alu Insertion; "
                    else:
                        # Pair + Split #
                        if (align_max[0][0] == "N" or int(align_max[0][0]) != 441) and (align_max[0][1] == "N" or int(align_max[0][1]) != 448) and (align_max[0][2] == "N" or int(align_max[0][2]) != 441) and (align_max[0][3] == "N" or int(align_max[0][3]) != 441):
                            Inspection = Inspection + "No PolyA; "
                        junc_5 = 0
                        for e in [0,1,2,3]:
                            if align_min[3][e] != "N" and align_min[0][e] != "N" and int(align_min[3][e]) - int(align_min[0][e]) > junc_5:
                                junc_5 = int(align_min[3][e]) - int(align_min[0][e])
                        if junc_5 > 5:
                            Inspection = Inspection + "Junction site difference; "
            else:
                # 1.2.3 ##     <-- ##
                #       ## ->-     ##
                if clip_map[1] == "mapped":
                    # (1) Arrow Conflict #
                    if int(clip_max[1]) < int(arrow_max[1]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    if (int(arrow_min[1]) < int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_min[0])) or (int(arrow_min[1]) > int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_max[0])):
                        if (int(arrow_max[1]) - int(arrow_min[0])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[1]) > int(arrow_max[0]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict
                    junc_3 = 0
                    for e in [0,1,2,3]:
                        if align_max[2][e] != "N" and align_max[0][e] != "N" and int(align_max[2][e]) - int(align_max[0][e]) > junc_3:
                            junc_3 = int(align_max[2][e]) - int(align_max[0][e])
                    if junc_3 > 5:
                        Inspection = Inspection + "Junction site difference; "
                    if align_color[3].count("N") != 4:
                        if transduction5_color[0].count("N") != 4:
                            strange_align_count = 0
                            for e in [0,1,2,3]:
                                if transduction5_min[0][e] != "N" and align_min[3][e] != "N" and abs(int(transduction5_min[0][e]) - int(transduction5_max[0][e])) > 10 and (int(transduction5_max[0][e]) - 10) > int(align_min[3][e]):
                                    strange_align_count = strange_align_count + 1
                            if strange_align_count >= 2:
                                Inspection = Inspection + "Unknown clip within Alu Insertion; "
                        if transduction3_color[3].count("N") != 4:
                            strange_align_count = 0
                            for e in [0,1,2,3]:
                                if transduction3_min[3][e] != "N" and align_min[0][e] != "N" and abs(int(transduction3_min[3][e]) - int(transduction3_max[3][e])) > 10 and (int(transduction3_min[3][e]) + 10) < int(align_max[0][e]):
                                    strange_align_count = strange_align_count + 1
                            if strange_align_count >= 2:
                                Inspection = Inspection + "Unknown clip within Alu Insertion; "
                # 1.2.4 ##     <-- ##
                #       ## -->     ##
                        ## Two Pair #
                else:
                    # (1) Arrow Conflict #
                    if (int(arrow_min[1]) < int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_min[0])) or (int(arrow_min[1]) > int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_max[0])):
                        if (int(arrow_max[1]) - int(arrow_min[0])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[1]) > int(arrow_max[0]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict
                    if transduction5_color[0].count("N") != 4:
                        strange_align_count = 0
                        for e in [0,1,2,3]:
                            if transduction5_min[0][e] != "N" and align_min[3][e] != "N" and abs(int(transduction5_min[0][e]) - int(transduction5_max[0][e])) > 10 and (int(transduction5_max[0][e]) - 10) > int(align_min[3][e]):
                                strange_align_count = strange_align_count + 1
                        if strange_align_count >= 2:
                            Inspection = Inspection + "Unknown clip within Alu Insertion; "
                    if transduction3_color[3].count("N") != 4:
                        strange_align_count = 0
                        for e in [0,1,2,3]:
                            if transduction3_min[3][e] != "N" and align_min[0][e] != "N" and abs(int(transduction3_min[3][e]) - int(transduction3_max[3][e])) > 10 and (int(transduction3_min[3][e]) + 10) < int(align_max[0][e]):
                                strange_align_count = strange_align_count + 1
                        if strange_align_count >= 2:
                            Inspection = Inspection + "Unknown clip within Alu Insertion; "
    # 2 # Only Upstream or Downstream
    else:
        # 2.1 ## -->     ##
        #     ##   -->   ##
        if reads_direction[0] == "upstream":
            if abs(int(arrow_max[0]) - int(arrow_max[1])) > 500:
                Inspection = Inspection + "Upstream arrows gap lager than 500bp"
            if clip_map[0] == "mapped":
                # 2.1.1 ## ->-     ##
                #       ## ->-     ##
                if clip_map[1] == "mapped":
                    # (1) Arrow Conflict #
                    if clip_min[0] != "N" and clip_min[1] != "N":
                        if abs(int(clip_min[0]) - int(clip_min[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    if int(clip_max[1]) < int(arrow_max[1]) or int(clip_max[0]) < int(arrow_max[0]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    # (2) Alignment Conflict #
                    abs_junction = 300
                    for e in [0,1,2,3]:
                        if align_min[0][e] != "N" and align_min[2][e] !="N" and abs(int(align_min[0][e]) - int(align_min[2][e])) < abs_junction:
                            abs_junction = abs(int(align_min[0][e]) - int(align_min[2][e]))
                    if abs_junction > 5:
                        Inspection = Inspection + "Junction site difference; "
                # 2.1.2 ##   ->-   ##
                #       ## -->     ##
                else:
                    # (1) Arrow Conflict #
                    if clip_min[0] != "N" and clip_min[1] != "N":
                        if abs(int(clip_min[0]) - int(clip_min[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    if int(clip_max[0]) < int(arrow_max[0]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    # (2) Alignment Conflict #
                    junc_5 = 0
                    for e in [0,1,2,3]:
                        if align_min[0][e] != "N" and align_min[3][e] != "N" and int(align_min[0][e]) - int(align_min[3][e]) > junc_5:
                            junc_5 = int(align_min[0][e]) - int(align_min[3][e])
                    if junc_5 > 5:
                        Inspection = Inspection + "Junction site difference; "
            else:
                # 2.1.3 ## -->     ##
                #       ##   ->-   ##
                if clip_map[1] == "mapped":
                    # (1) Arrow Conflict #
                    if clip_min[0] != "N" and clip_min[1] != "N":
                        if abs(int(clip_min[0]) - int(clip_min[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    if int(clip_max[1]) < int(arrow_max[1]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    junc_5 = 0
                    for e in [0,1,2,3]:
                        if align_min[2][e] != "N" and align_min[1][e] != "N" and int(align_min[2][e]) - int(align_min[1][e]) > junc_5:
                            junc_5 = int(align_min[2][e]) - int(align_min[1][e])
                    if junc_5 > 5:
                        Inspection = Inspection + "Junction site difference; "
                # 2.1.4 ## -->     ##
                #       ##   -->   ##
                else:
                    # (1) Arrow Conflict #
                    if clip_min[0] != "N" and clip_min[1] != "N":
                        if abs(int(clip_min[0]) - int(clip_min[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
        # 2.2 ##     <-- ##
        #     ##   <--   ##
        else:
            if abs(int(arrow_min[0]) - int(arrow_min[1])) > 500:
                Inspection = Inspection + "Downstream arrows gap lager than 500bp"
            if clip_map[0] == "mapped":
                # 2.2.1 ##     -<- ##
                #       ##     -<- ##
                if clip_map[1] == "mapped":
                    # (1) Arrow Conflict #
                    if clip_max[0] != "N" and clip_max[1] != "N":
                        if abs(int(clip_max[0]) - int(clip_max[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    if int(clip_min[0]) > int(arrow_min[0]) or int(clip_min[1]) > int(arrow_min[1]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    # (2) Alignment Conflict #
                    if align_color[1].count("N") != 4:
                        if align_color[3].count("N") != 4:
                            ## Two Clip-Pair ##
                            if ((align_max[1][0] == "N" or int(align_max[1][0]) != 441) and (align_max[1][1] == "N" or int(align_max[1][1]) != 448) and (align_max[1][2] == "N" or int(align_max[1][2]) != 441) and (align_max[1][3] == "N" or int(align_max[1][3]) != 441)) or ((align_max[3][0] == "N" or int(align_max[3][0]) != 441) and (align_max[3][1] == "N" or int(align_max[3][1]) != 448) and (align_max[3][2] == "N" or int(align_max[3][2]) != 441) and (align_max[3][3] == "N" or int(align_max[3][3]) != 441)):
                                Inspection = Inspection + "No PolyA; "
                            abs_junction = 300
                            for e in [0,1,2,3]:
                                if align_max[1][e] != "N" and align_max[3][e] !="N" and abs(int(align_max[1][e]) - int(align_max[3][e])) < abs_junction:
                                    abs_junction = abs(int(align_max[1][e]) - int(align_max[3][e]))
                            if abs_junction > 5:
                                Inspection = Inspection + "Junction site difference; "
                        else:
                            ## Clip-Pair + Split ##
                            if ((align_max[1][0] == "N" or int(align_max[1][0]) != 441) and (align_max[1][1] == "N" or int(align_max[1][1]) != 448) and (align_max[1][2] == "N" or int(align_max[1][2]) != 441) and (align_max[1][3] == "N" or int(align_max[1][3]) != 441)) or ((align_max[2][0] == "N" or int(align_max[2][0]) != 441) and (align_max[2][1] == "N" or int(align_max[2][1]) != 448) and (align_max[2][2] == "N" or int(align_max[2][2]) != 441) and (align_max[2][3] == "N" or int(align_max[2][3]) != 441)):
                                Inspection = Inspection + "No PolyA; "
                            abs_junction = 300
                            for e in [0,1,2,3]:
                                if align_max[1][e] != "N" and align_max[2][e] !="N" and abs(int(align_max[1][e]) - int(align_max[2][e])) < abs_junction:
                                    abs_junction = abs(int(align_max[1][e]) - int(align_max[2][e]))
                            if abs_junction > 5:
                                Inspection = Inspection + "Junction site difference; "
                    else:
                        ## Split + Clip-Pair ##
                        if align_color[3].count("N") != 4:
                            if ((align_max[0][0] == "N" or int(align_max[0][0]) != 441) and (align_max[0][1] == "N" or int(align_max[0][1]) != 448) and (align_max[0][2] == "N" or int(align_max[0][2]) != 441) and (align_max[0][3] == "N" or int(align_max[0][3]) != 441)) or ((align_max[3][0] == "N" or int(align_max[3][0]) != 441) and (align_max[3][1] == "N" or int(align_max[3][1]) != 448) and (align_max[3][2] == "N" or int(align_max[3][2]) != 441) and (align_max[3][3] == "N" or int(align_max[3][3]) != 441)):
                                Inspection = Inspection + "No PolyA; "
                            abs_junction = 300
                            for e in [0,1,2,3]:
                                if align_max[0][e] != "N" and align_max[3][e] !="N" and abs(int(align_max[0][e]) - int(align_max[3][e])) < abs_junction:
                                    abs_junction = abs(int(align_max[0][e]) - int(align_max[3][e]))
                            if abs_junction > 5:
                                Inspection = Inspection + "Junction site difference; "
                        ## Two Split ##
                        else:
                            if ((align_max[0][0] == "N" or int(align_max[0][0]) != 441) and (align_max[0][1] == "N" or int(align_max[0][1]) != 448) and (align_max[0][2] == "N" or int(align_max[0][2]) != 441) and (align_max[0][3] == "N" or int(align_max[0][3]) != 441)) or ((align_max[2][0] == "N" or int(align_max[2][0]) != 441) and (align_max[2][1] == "N" or int(align_max[2][1]) != 448) and (align_max[2][2] == "N" or int(align_max[2][2]) != 441) and (align_max[2][3] == "N" or int(align_max[2][3]) != 441)):
                                Inspection = Inspection + "No PolyA; "
                            abs_junction = 300
                            for e in [0,1,2,3]:
                                if align_max[0][e] != "N" and align_max[2][e] !="N" and abs(int(align_max[0][e]) - int(align_max[2][e])) < abs_junction:
                                    abs_junction = abs(int(align_max[0][e]) - int(align_max[2][e]))
                            if abs_junction > 5:
                                Inspection = Inspection + "Junction site difference; "
                # 2.2.2 ##   -<-   ##
                #       ##     <-- ##
                else:
                    # (1) Arrow Conflict #
                    if clip_max[0] != "N" and clip_max[1] != "N":
                        if abs(int(clip_max[0]) - int(clip_max[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    if int(clip_min[0]) > int(arrow_min[0]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    # (2) Alignment Conflict #
                    if align_color[1].count("N") != 4:
                        ## Clip-Pair + Pair
                        if (align_max[1][0] == "N" or int(align_max[1][0]) != 441) and (align_max[1][1] == "N" or int(align_max[1][1]) != 448) and (align_max[1][2] == "N" or int(align_max[1][2]) != 441) and (align_max[1][3] == "N" or int(align_max[1][3]) != 441):
                            Inspection = Inspection + "No PolyA; "
                        junc_3 = 0
                        for e in [0,1,2,3]:
                            if align_max[2][e] != "N" and align_max[1][e] != "N" and int(align_max[2][e]) - int(align_max[1][e]) > junc_3:
                                junc_3 = int(align_max[2][e]) - int(align_max[1][e])
                        if junc_3 > 5:
                            Inspection = Inspection + "Junction site difference; "
                    else:
                        ## Split + Pair ##
                        if (align_max[0][0] == "N" or int(align_max[0][0]) != 441) and (align_max[0][1] == "N" or int(align_max[0][1]) != 448) and (align_max[0][2] == "N" or int(align_max[0][2]) != 441) and (align_max[0][3] == "N" or int(align_max[0][3]) != 441):
                            Inspection = Inspection + "No PolyA; "
                        junc_3 = 0
                        for e in [0,1,2,3]:
                            if align_max[2][e] != "N" and align_max[0][e] != "N" and int(align_max[2][e]) - int(align_max[0][e]) > junc_3:
                                junc_3 = int(align_max[2][e]) - int(align_max[0][e])
                        if junc_3 > 5:
                            Inspection = Inspection + "Junction site difference; "
            else:
                # 2.1.3 ##     <-- ##
                #       ##   -<-   ##
                if clip_map[1] == "mapped":
                    # (1) Arrow Conflict #
                    if clip_max[0] != "N" and clip_max[1] != "N":
                        if abs(int(clip_max[0]) - int(clip_max[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    if int(clip_min[1]) > int(arrow_min[1]):
                        Inspection = Inspection + "Arrow with Alu tail; "
                    # (2) Alignment Conflict #
                    if align_color[3] != "N":
                        ##  Pair + Clip-Pair ##
                        if (align_max[3][0] == "N" or int(align_max[3][0]) != 441) and (align_max[3][1] == "N" or int(align_max[3][1]) != 448) and (align_max[3][2] == "N" or int(align_max[3][2]) != 441) and (align_max[3][3] == "N" or int(align_max[3][3]) != 441):
                            Inspection = Inspection + "No PolyA; "
                        junc_3 = 0
                        for e in [0,1,2,3]:
                            if align_max[0][e] != "N" and align_max[3][e] != "N" and int(align_max[0][e]) - int(align_max[3][e]) > junc_3:
                                junc_3 = int(align_max[0][e]) - int(align_max[3][e])
                        if junc_3 > 5:
                            Inspection = Inspection + "Junction site difference; "
                    else:
                        ## Pair + Split ##
                        if (align_max[2][0] == "N" or int(align_max[2][0]) != 441) and (align_max[2][1] == "N" or int(align_max[2][1]) != 448) and (align_max[2][2] == "N" or int(align_max[2][2]) != 441) and (align_max[2][3] == "N" or int(align_max[2][3]) != 441):
                            Inspection = Inspection + "No PolyA; "
                        junc_3 = 0
                        for e in [0,1,2,3]:
                            if align_max[0][e] != "N" and align_max[2][e] != "N" and int(align_max[0][e]) - int(align_max[2][e]) > junc_3:
                                junc_3 = int(align_max[0][e]) - int(align_max[2][e])
                        if junc_3 > 5:
                            Inspection = Inspection + "Junction site difference; "
                # 2.1.4 ##     <-- ##
                #       ##    <--  ##
                else:
                    ## Two Pair-End ##
                    # (1) Arrow Conflict #
                    if clip_max[0] != "N" and clip_max[1] != "N":
                        if abs(int(clip_max[0]) - int(clip_max[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
    ## PCR duplication #
    PCR_count = -1
    if reads_direction[0] == reads_direction[1]:
        if align_color[0] == align_color[2] and align_color[1] == align_color[3]:
            if abs(int(arrow_min[0]) - int(arrow_min[1])) < 5:
                PCR_count = PCR_count + 1
            if abs(int(arrow_max[0]) - int(arrow_max[1])) < 5:
                PCR_count = PCR_count + 1
            pcr_align = -1
            for e in [0,1,2,3]:
                pcr_tmp_count = 0
                if align_color[0][e] != "N":
                    if abs(int(align_min[0][e]) - int(align_min[2][e])) < 5:
                        pcr_tmp_count = pcr_tmp_count + 1
                    if abs(int(align_max[0][e]) - int(align_max[2][e])) < 5:
                        pcr_tmp_count = pcr_tmp_count + 1
                if align_color[1][e] != "N":
                    if abs(int(align_min[1][e]) - int(align_min[3][e])) < 5:
                        pcr_tmp_count = pcr_tmp_count + 1
                    if abs(int(align_max[1][e]) - int(align_max[3][e])) < 5:
                        pcr_tmp_count = pcr_tmp_count + 1
                if pcr_align ==  -1 or pcr_align > pcr_tmp_count:
                    pcr_align = pcr_tmp_count
            if pcr_align != -1:
                PCR_count = PCR_count + pcr_align
    if PCR_count >= 3:
        Inspection = Inspection + "PCR duplicate? "

    ## Overlapping alignment ##
    if align_color[0].count("N") != 4:
        if align_color[2].count("N") !=4:
            overlap_mappability = -1
            for e in [0,1,2,3]:
                # Overlapping #
                if align_color[0][e] != "N" and align_color[2][e] != "N":
                    if int(align_max[0][e]) > int(align_min[2][e]) and int(align_min[0][e]) < int(align_max[2][e]) and abs(max(int(align_min[0][e]),int(align_min[2][e])) - min(int(align_max[0][e]),int(align_max[2][e]))) > 20:
                        overlap_notEqual_count = 0
                        for k in range(max(int(align_min[0][e]),int(align_min[2][e])),min(int(align_max[0][e]),int(align_max[2][e]))+1):
                            if temp_tensor[:,17+e*5:21+e*5,k].equal(temp_tensor[:,73+e*5:77+e*5,k]):
                                overlap_notEqual_count = overlap_notEqual_count + 1
                        map_tmp = overlap_notEqual_count/(min(int(align_max[0][e]),int(align_max[2][e])) - max(int(align_min[0][e]),int(align_min[2][e])) + 1)
                        if map_tmp > overlap_mappability:
                            overlap_mappability = map_tmp
            if overlap_mappability != -1 and overlap_mappability < 0.9:
                Inspection = Inspection + "Overlap alignment is too low; "
        if align_color[3].count("N") !=4:
            overlap_mappability = -1
            for e in [0,1,2,3]:
                # Overlapping #
                if align_color[0][e] != "N" and align_color[3][e] != "N":
                    if int(align_max[0][e]) > int(align_min[3][e]) and int(align_min[0][e]) < int(align_max[3][e]) and abs(max(int(align_min[0][e]),int(align_min[3][e])) - min(int(align_max[0][e]),int(align_max[3][e]))) > 20:
                        overlap_notEqual_count = 0
                        for k in range(max(int(align_min[0][e]),int(align_min[3][e])),min(int(align_max[0][e]),int(align_max[3][e]))+1):
                            if temp_tensor[:,17+e*5:21+e*5,k].equal(temp_tensor[:,101+e*5:105+e*5,k]):
                                overlap_notEqual_count = overlap_notEqual_count + 1
                        map_tmp = overlap_notEqual_count/(min(int(align_max[0][e]),int(align_max[3][e])) - max(int(align_min[0][e]),int(align_min[3][e])) + 1)
                        if map_tmp > overlap_mappability:
                            overlap_mappability = map_tmp
            if overlap_mappability != -1 and overlap_mappability < 0.9:
                Inspection = Inspection + "Overlap alignment is too low; "
    if align_color[1].count("N") != 4:
        if align_color[2].count("N") !=4:
            overlap_mappability = -1
            for e in [0,1,2,3]:
                # Overlapping #
                if align_color[1][e] != "N" and align_color[2][e] != "N":
                    if int(align_max[1][e]) > int(align_min[2][e]) and int(align_min[1][e]) < int(align_max[2][e]) and abs(max(int(align_min[1][e]),int(align_min[2][e])) - min(int(align_max[1][e]),int(align_max[2][e]))) > 20:
                        overlap_notEqual_count = 0
                        for k in range(max(int(align_min[1][e]),int(align_min[2][e])),min(int(align_max[1][e]),int(align_max[2][e]))+1):
                            if temp_tensor[:,45+e*5:49+e*5,k].equal(temp_tensor[:,73+e*5:77+e*5,k]):
                                overlap_notEqual_count = overlap_notEqual_count + 1
                        map_tmp = overlap_notEqual_count/(min(int(align_max[1][e]),int(align_max[2][e])) - max(int(align_min[1][e]),int(align_min[2][e])) + 1)
                        if map_tmp > overlap_mappability:
                            overlap_mappability = map_tmp
            if overlap_mappability != -1 and overlap_mappability < 0.9:
                Inspection = Inspection + "Overlap alignment is too low; "
        if align_color[3].count("N") !=4:
            overlap_mappability = -1
            for e in [0,1,2,3]:
                # Overlapping #
                if align_color[1][e] != "N" and align_color[3][e] != "N":
                    if int(align_max[1][e]) > int(align_min[3][e]) and int(align_min[1][e]) < int(align_max[3][e]) and abs(max(int(align_min[1][e]),int(align_min[3][e])) - min(int(align_max[1][e]),int(align_max[3][e]))) > 20:
                        overlap_notEqual_count = 0
                        for k in range(max(int(align_min[1][e]),int(align_min[3][e])),min(int(align_max[1][e]),int(align_max[3][e]))+1):
                            if temp_tensor[:,45+e*5:49+e*5,k].equal(temp_tensor[:,101+e*5:105+e*5,k]):
                                overlap_notEqual_count = overlap_notEqual_count + 1
                        map_tmp = overlap_notEqual_count/(min(int(align_max[1][e]),int(align_max[3][e])) - max(int(align_min[1][e]),int(align_min[3][e])) + 1)
                        if map_tmp > overlap_mappability:
                            overlap_mappability = map_tmp
            if overlap_mappability != -1 and overlap_mappability < 0.9:
                Inspection = Inspection + "Overlap alignment is too low; "

    ## read1 and read2 mapability ##
    read_align = -1
    if align_color[0].count("N") != 4:
        for e in [0,1,2,3]:
            if align_color[0][e] != "N" and abs(int(align_max[0][e]) - int(align_min[0][e])) > 20:
                if float(align_mapability[0][e]) > read_align:
                    read_align = float(align_mapability[0][e])
        if read_align != -1 and read_align < 0.9:
            Inspection = Inspection + "read1 mappability is lower than 90%; "
    read_align = -1
    if align_color[1].count("N") != 4:
        for e in [0,1,2,3]:
            if align_color[1][e] != "N" and abs(int(align_max[1][e]) - int(align_min[1][e])) > 20:
                if float(align_mapability[1][e]) > read_align:
                    read_align = float(align_mapability[1][e])
        if read_align != -1 and read_align < 0.9:
            Inspection = Inspection + "read1 mappability is lower than 90%; "
    read_align = -1
    if align_color[2].count("N") != 4:
        for e in [0,1,2,3]:
            if align_color[2][e] != "N" and abs(int(align_max[2][e]) - int(align_min[2][e])) > 20:
                if float(align_mapability[2][e]) > read_align:
                    read_align = float(align_mapability[2][e])
        if read_align != -1 and read_align < 0.9:
            Inspection = Inspection + "read2 mappability is lower than 90%; "
    read_align = -1
    if align_color[3].count("N") != 4:
        for e in [0,1,2,3]:
            if align_color[3][e] != "N" and abs(int(align_max[3][e]) - int(align_min[3][e])) > 20:
                if float(align_mapability[3][e]) > read_align:
                    read_align = float(align_mapability[3][e])
        if read_align != -1 and read_align < 0.9:
            Inspection = Inspection + "read2 mappability is lower than 90%; "

if Inspection == "":
    Inspection = "PASS"

## Mark transductions ##
transduction3_mark = ""
if transduction3_color[0].count("N") != 4:
    if transduction3_color[0].count("green") != 0:
        transduction_len = 0
        for e in [0,1,2,3]:
            if transduction3_min[0][e] != "N" and transduction3_max[0][e] != "N" and abs(int(transduction3_min[0][e]) - int(transduction3_max[0][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[0][e]) - int(transduction3_max[0][e]))
        if transduction_len >= 10:
            transduction3_mark = transduction3_mark + "read1 3'green " + str(transduction_len) + "bp; "
if transduction3_color[1].count("N") != 4:
    if transduction3_color[1].count("green") != 0:
        transduction_len = 0
        for e in [0,1,2,3]:
            if transduction3_min[1][e] != "N" and transduction3_max[1][e] != "N" and abs(int(transduction3_min[1][e]) - int(transduction3_max[1][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[1][e]) - int(transduction3_max[1][e]))
        if transduction_len >= 10:
            transduction3_mark = transduction3_mark + "read1 3'green " + str(transduction_len) + "bp; "
if transduction3_color[2].count("N") != 4:
    if transduction3_color[2].count("green") != 0:
        transduction_len = 0
        for e in [0,1,2,3]:
            if transduction3_min[2][e] != "N" and transduction3_max[2][e] != "N" and abs(int(transduction3_min[2][e]) - int(transduction3_max[2][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[2][e]) - int(transduction3_max[2][e]))
        if transduction_len >= 10:
            transduction3_mark = transduction3_mark + "read2 3'green " + str(transduction_len) + "bp; "
if transduction3_color[3].count("N") != 4:
    if transduction3_color[3].count("green") != 0:
        transduction_len = 0
        for e in [0,1,2,3]:
            if transduction3_min[3][e] != "N" and transduction3_max[3][e] != "N" and abs(int(transduction3_min[3][e]) - int(transduction3_max[3][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[3][e]) - int(transduction3_max[3][e]))
        if transduction_len >= 10:
            transduction3_mark = transduction3_mark + "read2 3'green " + str(transduction_len) + "bp; "

transduction5_mark = ""
if transduction5_color[0].count("N") != 4:
    if transduction5_color[0].count("green") != 0:
        transduction_len = 0
        for e in [0,1,2,3]:
            if transduction5_min[0][e] != "N" and transduction5_max[0][e] != "N" and abs(int(transduction5_min[0][e]) - int(transduction5_max[0][e])) > transduction_len:
                transduction_len = abs(int(transduction5_min[0][e]) - int(transduction5_max[0][e]))
        if transduction_len >= 10:
            transduction5_mark = transduction5_mark + "read1 5'green " + str(transduction_len) + "bp; "
if transduction5_color[1].count("N") != 4:
    if transduction5_color[1].count("green") != 0:
        transduction_len = 0
        for e in [0,1,2,3]:
            if transduction5_min[1][e] != "N" and transduction5_max[1][e] != "N" and abs(int(transduction5_min[1][e]) - int(transduction5_max[1][e])) > transduction_len:
                transduction_len = abs(int(transduction5_min[1][e]) - int(transduction5_max[1][e]))
        if transduction_len >= 10:
            transduction5_mark = transduction5_mark + "read1 5'green " + str(transduction_len) + "bp; "
if transduction5_color[2].count("N") != 4:
    if transduction5_color[2].count("green") != 0:
        transduction_len = 0
        for e in [0,1,2,3]:
            if transduction5_min[2][e] != "N" and transduction5_max[2][e] != "N" and abs(int(transduction5_min[2][e]) - int(transduction5_max[2][e])) > transduction_len:
                transduction_len = abs(int(transduction5_min[2][e]) - int(transduction5_max[2][e]))
        if transduction_len >= 10:
            transduction5_mark = transduction5_mark + "read2 5'green " + str(transduction_len) + "bp; "
if transduction5_color[3].count("N") != 4:
    if transduction5_color[3].count("green") != 0:
        transduction_len = 0
        for e in [0,1,2,3]:
            if transduction5_min[3][e] != "N" and transduction5_max[3][e] != "N" and abs(int(transduction5_min[3][e]) - int(transduction5_max[3][e])) > transduction_len:
                transduction_len = abs(int(transduction5_min[3][e]) - int(transduction5_max[3][e]))
        if transduction_len >= 10:
            transduction5_mark = transduction5_mark + "read2 5'green " + str(transduction_len) + "bp; "

f = open("Inspection_tmp_ALU_{}.txt".format(ver),"w")
f.write("{}\t{}\t{}\t".format(Inspection,transduction3_mark,transduction5_mark))
f.close()
