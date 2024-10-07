import os
import sys
import numpy
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
TEclass = "SVA"

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
for j in [10,28,46,64]:
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

for k in [17,35,53,71]:
    align_min_tmp = []
    align_max_tmp = []
    align_color_tmp = []
    align_mapability_tmp = []
    transduction3_color_tmp = []
    transduction3_min_tmp = []
    transduction3_max_tmp = []
    for m in [0,1]:
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
                transduction3_min_tmp.append(min(green_detect))
                transduction3_max_tmp.append(max(green_detect))
            elif black_detect.numel() != 0:
                transduction3_color_tmp.append("black")
                transduction3_min_tmp.append(min(black_detect))
                transduction3_max_tmp.append(max(black_detect))
            else:
                transduction3_color_tmp.append("N")
                transduction3_min_tmp.append("N")
                transduction3_max_tmp.append("N")
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
                transduction3_min_tmp.append(min(green_detect))
                transduction3_max_tmp.append(max(green_detect))
            elif black_detect.numel() != 0:
                transduction3_color_tmp.append("black")
                transduction3_min_tmp.append(min(black_detect))
                transduction3_max_tmp.append(max(black_detect))
            else:
                transduction3_color_tmp.append("N")
                transduction3_min_tmp.append("N")
                transduction3_max_tmp.append("N")
        else:
            align_mapability_tmp.append("N")
            align_color_tmp.append("N")
            align_min_tmp.append("N")
            align_max_tmp.append("N")
            transduction3_color_tmp.append("N")
            transduction3_min_tmp.append("N")
            transduction3_max_tmp.append("N")
    align_mapability.append(align_mapability_tmp)
    align_color.append(align_color_tmp)
    align_min.append(align_min_tmp)
    align_max.append(align_max_tmp)
    transduction3_color.append(transduction3_color_tmp)
    transduction3_min.append(transduction3_min_tmp)
    transduction3_max.append(transduction3_max_tmp)

for e in [0,1,2,3]:
    if align_color[e].count("N") == 1:
        if align_color[e].index("N") == 0:
            align_min[e][0] = align_min[e][1] - 5
            align_max[e][0] = align_max[e][1] - 5
        else:
            align_min[e][1] = align_min[e][0] + 5
            align_max[e][1] = align_max[e][0] + 5

align_count_read1 = 0
align_count_read2 = 0
for e in [0,1]:
    if align_color[e].count("N") != 2:
        align_count_read1 = align_count_read1 + 1
    if align_color[e+2].count("N") !=2:
        align_count_read2 = align_count_read2 + 1
    

Inspection = ""
PolyA_len = ""
## Manual Inspection ##
# 1 # Upstream + Downstream
if align_count_read1 == 0 or align_count_read2 == 0:
    Inspection = "miss alignment! "
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
                        Inspection = Inspection + "Arrow with SVA tail; "
                    if (int(arrow_min[0]) < int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_min[1])) or (int(arrow_min[0]) > int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_max[1])):
                        if (int(arrow_max[0]) - int(arrow_min[1])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[0]) > int(arrow_max[1]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict #
                    if align_color[1].count("N") != 2:
                        # Two Clip-pair #
                        if align_color[3].count("N") != 2:
                            cross_over_len1 = int(align_min[1][0]) - int(align_max[2][0])
                            cross_over_len2 = int(align_min[1][1]) - int(align_max[2][1])
                            if cross_over_len1 > 0 and cross_over_len2 > 0 and min(cross_over_len1,cross_over_len2) > 600:
                                Inspection = Inspection + "Alignments cross over; "
                            if int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159):
                                Inspection = Inspection + "all mapping in polyA; "
                            if int(align_max[3][0]) < (1360 + 159) and int(align_max[3][1]) < (1355 + 159):
                                Inspection = Inspection + "No polyA; "
                            else:
                                if int(align_max[3][0]) < (1370 + 159) and int(align_max[3][1]) < (1365 + 159):
                                    Inspection = Inspection + "PolyA too short; "
                                tmp_ployA_len = max((int(align_max[3][0]) - 1359 - 159),(int(align_max[3][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                            if (int(align_min[0][0]) - int(align_min[3][0]) > 5 or int(align_max[0][0]) - int(align_max[3][0]) > 5) and (int(align_min[0][1]) - int(align_min[3][1]) > 5 or int(align_max[0][1]) - int(align_max[3][1]) > 5):
                                Inspection = Inspection + "Junction site difference; "
                        # Clip-pair + Split #
                        else:
                            if int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159):
                                Inspection = Inspection + "all mapping in polyA; "
                            if int(align_max[2][0]) < (1360 + 159) and int(align_max[2][1]) < (1355 + 159):
                                Inspection = Inspection + "No polyA; "
                            else:
                                if int(align_max[2][0]) < (1370 + 159) and int(align_max[2][1]) < (1365 + 159):
                                    Inspection = Inspection + "PolyA too short; "
                                tmp_ployA_len = max((int(align_max[2][0]) - 1359 - 159),(int(align_max[2][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                            if (int(align_min[0][0]) - int(align_min[2][0]) > 5 or int(align_max[0][0]) - int(align_max[2][0]) > 5) and (int(align_min[0][1]) - int(align_min[2][1]) > 5 or int(align_max[0][1]) - int(align_max[2][1]) > 5):
                                Inspection = Inspection + "Junction site difference; "
                    else:
                        # Split + Clip-pair #
                        if int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159):
                                Inspection = Inspection + "all mapping in polyA; "
                        if align_color[3].count("N") != 2:
                            if int(align_max[3][0]) < (1360 + 159) and int(align_max[3][1]) < (1355 + 159):
                                Inspection = Inspection + "No polyA; "
                            else:
                                if int(align_max[3][0]) < (1370 + 159) and int(align_max[3][1]) < (1365 + 159):
                                    Inspection = Inspection + "PolyA too short; "
                                tmp_ployA_len = max((int(align_max[3][0]) - 1359 - 159),(int(align_max[3][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                            if (int(align_min[0][0]) - int(align_min[3][0]) > 5 or int(align_max[0][0]) - int(align_max[3][0]) > 5) and (int(align_min[0][1]) - int(align_min[3][1]) > 5 or int(align_max[0][1]) - int(align_max[3][1]) > 5):
                                Inspection = Inspection + "Junction site difference; "
                        # Split + Split #
                        else:
                            if int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159):
                                Inspection = Inspection + "all mapping in polyA; "
                            if int(align_max[2][0]) < (1360 + 159) and int(align_max[2][1]) < (1355 + 159):
                                Inspection = Inspection + "No polyA; "
                            else:
                                if int(align_max[2][0]) < (1370 + 159) and int(align_max[2][1]) < (1365 + 159):
                                    Inspection = Inspection + "PolyA too short; "
                                tmp_ployA_len = max((int(align_max[2][0]) - 1359 - 159),(int(align_max[2][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                            if (int(align_min[0][0]) - int(align_min[2][0]) > 5 or int(align_max[0][0]) - int(align_max[2][0]) > 5) and (int(align_min[0][1]) - int(align_min[2][1]) > 5 or int(align_max[0][1]) - int(align_max[2][1]) > 5):
                                Inspection = Inspection + "Junction site difference; "
                # 1.1.2 ## ->-     ##
                #       ##     <-- ##
                else:
                    # (1) Arrow Conflict #
                    if int(clip_max[0]) < int(arrow_max[0]):
                        Inspection = Inspection + "Arrow with SVA tail; "
                    if (int(arrow_min[0]) < int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_min[1])) or (int(arrow_min[0]) > int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_max[1])):
                        if (int(arrow_max[0]) - int(arrow_min[1])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[0]) > int(arrow_max[1]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict #
                    if align_color[1].count("N") != 2:
                        # Clip-pair + Pair #
                        cross_over_len1 = int(align_min[1][0]) - int(align_max[2][0])
                        cross_over_len2 = int(align_min[1][1]) - int(align_max[2][1])
                        if cross_over_len1 > 0 and cross_over_len2 > 0 and min(cross_over_len1,cross_over_len2) > 600:
                            Inspection = Inspection + "Alignments cross over; "
                        if int(align_max[2][0]) < (1381 + 159 - 400) and int(align_max[2][1]) < (1374 + 159 - 400):
                            Inspection = Inspection + "Large gap to PolyA; "
                        if int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if (int(align_min[0][0]) - int(align_min[2][0]) > 5) and (int(align_min[0][1]) - int(align_min[2][1]) > 5):
                            Inspection = Inspection + "Junction site difference; "
                    else:
                        # Split + Pair #
                        if int(align_max[2][0]) < (1381 + 159 - 400) and int(align_max[2][1]) < (1374 + 159 - 400):
                            Inspection = Inspection + "Large gap to PolyA; "
                        if int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159):
                            Inspection = Inspection + "all mapping in polyA"
                        if (int(align_min[0][0]) - int(align_min[2][0]) > 5) and (int(align_min[0][1]) - int(align_min[2][1]) > 5):
                            Inspection = Inspection + "Junction site difference; "
            else:
                # 1.1.3 ## -->     ##
                #       ##     -<- ##
                if clip_map[1] == "mapped":
                    # (1) Arrow Conflict #
                    if int(clip_min[1]) > int(arrow_min[1]):
                        Inspection = Inspection + "Arrow with SVA tail; "
                    if (int(arrow_min[0]) < int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_min[1])) or (int(arrow_min[0]) > int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_max[1])):
                        if (int(arrow_max[0]) - int(arrow_min[1])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[0]) > int(arrow_max[1]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict #
                    if align_color[3].count("N") != 2:
                        # Pair + Clip-pair #
                        cross_over_len1 = int(align_min[1][0]) - int(align_max[2][0])
                        cross_over_len2 = int(align_min[1][1]) - int(align_max[2][1])
                        if cross_over_len1 > 0 and cross_over_len2 > 0 and min(cross_over_len1,cross_over_len2) > 600:
                            Inspection = Inspection + "Alignments cross over; "
                        if int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[3][0]) < (1360 + 159) and int(align_max[3][1]) < (1355 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[3][0]) < (1370 + 159) and int(align_max[3][1]) < (1365 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            tmp_ployA_len = max((int(align_max[3][0]) - 1359 - 159),(int(align_max[3][1]) - 1355 - 159))
                            PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                        if (int(align_min[1][0]) - int(align_min[3][0]) > 5) and (int(align_min[1][1]) - int(align_min[3][1]) > 5):
                            Inspection = Inspection + "Junction site difference; "             
                    else:
                        # Pair + Split #
                        if int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[2][0]) < (1360 + 159) and int(align_max[2][1]) < (1355 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[2][0]) < (1370 + 159) and int(align_max[2][1]) < (1365 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            tmp_ployA_len = max((int(align_max[2][0]) - 1359 - 159),(int(align_max[2][1]) - 1355 - 159))
                            PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                        if (int(align_min[1][0]) - int(align_min[2][0]) > 5) and (int(align_min[1][1]) - int(align_min[2][1]) > 5):
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
                    # (2) Alignment Conflict #
                    cross_over_len1 = int(align_min[1][0]) - int(align_max[2][0])
                    cross_over_len2 = int(align_min[1][1]) - int(align_max[2][1])
                    if cross_over_len1 > 0 and cross_over_len2 > 0 and min(cross_over_len1,cross_over_len2) > 600:
                        Inspection = Inspection + "Alignments cross over; "
                    if int(align_max[2][0]) < (1381 + 159 - 400) and int(align_max[2][1]) < (1374 + 159 - 400):
                        Inspection = Inspection + "Large gap to PolyA; "
                    if int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159):
                        Inspection = Inspection + "all mapping in polyA; "
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
                        Inspection = Inspection + "Arrow with SVA tail; "
                    if (int(arrow_min[1]) < int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_min[0])) or (int(arrow_min[1]) > int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_max[0])):
                        if (int(arrow_max[1]) - int(arrow_min[0])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[1]) > int(arrow_max[0]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict #
                    if align_color[3].count("N") != 2:
                        # Two Clip-pair #
                        if align_color[1].count("N") != 2:
                            cross_over_len1 = int(align_min[3][0]) - int(align_max[0][0])
                            cross_over_len2 = int(align_min[3][1]) - int(align_max[0][1])
                            if cross_over_len1 > 0 and cross_over_len2 > 0 and min(cross_over_len1,cross_over_len2) > 600:
                                Inspection = Inspection + "Alignments cross over; "
                            if int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159):
                                Inspection = Inspection + "all mapping in polyA; "
                            if int(align_max[1][0]) < (1360 + 159) and int(align_max[1][1]) < (1355 + 159):
                                Inspection = Inspection + "No polyA; "
                            else:
                                if int(align_max[1][0]) < (1370 + 159) and int(align_max[1][1]) < (1365 + 159):
                                    Inspection = Inspection + "PolyA too short; "
                                tmp_ployA_len = max((int(align_max[1][0]) - 1359 - 159),(int(align_max[1][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                            if (int(align_min[2][0]) - int(align_min[1][0]) > 5 or int(align_max[2][0]) - int(align_max[1][0]) > 5) and (int(align_min[2][1]) - int(align_min[1][1]) > 5 or int(align_max[2][1]) - int(align_max[1][1]) > 5):
                                Inspection = Inspection + "Junction site difference; "
                        # Clip-pair + Split #
                        else:
                            if int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159):
                                Inspection = Inspection + "all mapping in polyA; "
                            if int(align_max[0][0]) < (1360 + 159) and int(align_max[0][1]) < (1355 + 159):
                                Inspection = Inspection + "No polyA; "
                            else:
                                if int(align_max[0][0]) < (1370 + 159) and int(align_max[0][1]) < (1365 + 159):
                                    Inspection = Inspection + "PolyA too short; "
                                tmp_ployA_len = max((int(align_max[0][0]) - 1359 - 159),(int(align_max[0][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                            if (int(align_min[2][0]) - int(align_min[0][0]) > 5 or int(align_max[2][0]) - int(align_max[0][0]) > 5) and (int(align_min[2][1]) - int(align_min[0][1]) > 5 or int(align_max[2][1]) - int(align_max[0][1]) > 5):
                                Inspection = Inspection + "Junction site difference; "
                    else:
                        # Split + Clip-pair #
                        if align_color[1].count("N") != 2:
                            if int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159):
                                Inspection = Inspection + "all mapping in polyA; "
                            if int(align_max[1][0]) < (1360 + 159) and int(align_max[1][1]) < (1355 + 159):
                                Inspection = Inspection + "No polyA; "
                            else:
                                if int(align_max[1][0]) < (1370 + 159) and int(align_max[1][1]) < (1365 + 159):
                                    Inspection = Inspection + "PolyA too short; "
                                tmp_ployA_len = max((int(align_max[1][0]) - 1359 - 159),(int(align_max[1][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                            if (int(align_min[2][0]) - int(align_min[1][0]) > 5 or int(align_max[2][0]) - int(align_max[1][0]) > 5) and (int(align_min[2][1]) - int(align_min[1][1]) > 5 or int(align_max[2][1]) - int(align_max[1][1]) > 5):
                                Inspection = Inspection + "Junction site difference; "
                        # Split + Split #
                        else:
                            if int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159):
                                Inspection = Inspection + "all mapping in polyA; "
                            if int(align_max[0][0]) < (1360 + 159) and int(align_max[0][1]) < (1355 + 159):
                                Inspection = Inspection + "No polyA; "
                            else:
                                if int(align_max[0][0]) < (1370 + 159) and int(align_max[0][1]) < (1365 + 159):
                                    Inspection = Inspection + "PolyA too short; "
                                tmp_ployA_len = max((int(align_max[0][0]) - 1359 - 159),(int(align_max[0][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                            if (int(align_min[2][0]) - int(align_min[0][0]) > 5 or int(align_max[2][0]) - int(align_max[0][0]) > 5) and (int(align_min[2][1]) - int(align_min[0][1]) > 5 or int(align_max[2][1]) - int(align_max[0][1]) > 5):
                                Inspection = Inspection + "Junction site difference; "
                # 1.2.2 ##     -<- ##
                #       ## -->     ##
                else:
                    if int(clip_min[0]) > int(arrow_min[0]):
                        Inspection = Inspection + "Arrow with SVA tail; "
                    if (int(arrow_min[1]) < int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_min[0])) or (int(arrow_min[1]) > int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_max[0])):
                        if (int(arrow_max[1]) - int(arrow_min[0])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[1]) > int(arrow_max[0]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict #
                    if align_color[1].count("N") != 2:
                        # Pair + Clip-pair #
                        cross_over_len1 = int(align_min[3][0]) - int(align_max[0][0])
                        cross_over_len2 = int(align_min[3][1]) - int(align_max[0][1])
                        if cross_over_len1 > 0 and cross_over_len2 > 0 and min(cross_over_len1,cross_over_len2) > 600:
                            Inspection = Inspection + "Alignments cross over; "
                        if int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[1][0]) < (1360 + 159) and int(align_max[1][1]) < (1355 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[1][0]) < (1370 + 159) and int(align_max[1][1]) < (1365 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            tmp_ployA_len = max((int(align_max[1][0]) - 1359 - 159),(int(align_max[1][1]) - 1355 - 159))
                            PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                        if (int(align_max[3][0]) - int(align_max[1][0]) > 5) and (int(align_max[3][1]) - int(align_max[1][1]) > 5):
                            Inspection = Inspection + "Junction site difference; "              
                    else:
                        # Pair + Split #
                        if int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[0][0]) < (1360 + 159) and int(align_max[0][1]) < (1355 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[0][0]) < (1370 + 159) and int(align_max[0][1]) < (1365 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            tmp_ployA_len = max((int(align_max[0][0]) - 1359 - 159),(int(align_max[0][1]) - 1355 - 159))
                            PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                        if (int(align_max[3][0]) - int(align_max[0][0]) > 5) and (int(align_max[3][1]) - int(align_max[0][1]) > 5):
                            Inspection = Inspection + "Junction site difference; " 
            else:
                # 1.2.3 ##     <-- ##
                #       ## ->-     ##
                if clip_map[1] == "mapped":
                    # (1) Arrow Conflict #
                    if int(clip_max[1]) < int(arrow_max[1]):
                        Inspection = Inspection + "Arrow with SVA tail; "
                    if (int(arrow_min[1]) < int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_min[0])) or (int(arrow_min[1]) > int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_max[0])):
                        if (int(arrow_max[1]) - int(arrow_min[0])) > 25:
                            Inspection = Inspection + "TSD is longer than 25bp; "
                    if int(arrow_min[1]) > int(arrow_max[0]):
                        Inspection = Inspection + "Arrow Conflict; "
                    # (2) Alignment Conflict #
                    if align_color[3].count("N") != 2:
                        # Clip-pair + Pair #
                        cross_over_len1 = int(align_min[3][0]) - int(align_max[0][0])
                        cross_over_len2 = int(align_min[3][1]) - int(align_max[0][1])
                        if cross_over_len1 > 0 and cross_over_len2 > 0 and min(cross_over_len1,cross_over_len2) > 600:
                            Inspection = Inspection + "Alignments cross over; "
                        if int(align_max[0][0]) < (1381 + 159 - 400) and int(align_max[0][1]) < (1374 + 159 - 400):
                            Inspection = Inspection + "Large gap to PolyA; "
                        if int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if (int(align_min[2][0]) - int(align_min[0][0]) > 5) and (int(align_min[2][1]) - int(align_min[0][1]) > 5):
                            Inspection = Inspection + "Junction site difference; "
                    else:
                        # Split + Pair #
                        if int(align_max[0][0]) < (1381 + 159 - 400) and int(align_max[0][1]) < (1374 + 159 - 400):
                            Inspection = Inspection + "Large gap to PolyA; "
                        if int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159):
                            Inspection = Inspection + "all mapping in polyA"
                        if (int(align_min[2][0]) - int(align_min[0][0]) > 5) and (int(align_min[2][1]) - int(align_min[0][1]) > 5):
                            Inspection = Inspection + "Junction site difference; "
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
                    # (2) Alignment Conflict #
                    cross_over_len1 = int(align_min[3][0]) - int(align_max[0][0])
                    cross_over_len2 = int(align_min[3][1]) - int(align_max[0][1])
                    if cross_over_len1 > 0 and cross_over_len2 > 0 and min(cross_over_len1,cross_over_len2) > 600:
                        Inspection = Inspection + "Alignments cross over; "
                    if int(align_max[0][0]) < (1381 + 159 - 400) and int(align_max[0][1]) < (1374 + 159 - 400):
                        Inspection = Inspection + "Large gap to PolyA; "
                    if int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159):
                        Inspection = Inspection + "all mapping in polyA; "
    # 2 # Only Upstream or Downstream
    else:
        # 2.1 ## -->     ##
        #     ##   -->   ##
        if reads_direction[0] == "upstream":
            if clip_map[0] == "mapped":
                # 2.1.1 ## ->-     ##
                #       ## ->-     ##
                if clip_map[1] == "mapped":
                    # (1) Arrow Conflict #
                    if clip_min[0] != "N" and clip_min[1] != "N":
                        if abs(int(clip_min[0]) - int(clip_min[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    if int(clip_max[1]) < int(arrow_max[1]) or int(clip_max[0]) < int(arrow_max[0]):
                        Inspection = Inspection + "Arrow with SVA tail; "
                    # (2) Alignment Conflict #
                    try:
                        if abs(int(align_min[0][0]) - int(align_min[2][0])) > 5 and abs(int(align_min[0][1]) - int(align_min[2][1])) > 5:
                            Inspection = Inspection + "Junction site difference; "
                        if (align_color[1].count("N") != 2 and int(align_min[1][0]) > (1350 + 159) and int(align_min[1][1]) > (1345 +159)) or (align_color[3].count("N") != 2 and int(align_min[3][0]) > (1350 + 159) and int(align_min[3][1]) > (1345 +159)):
                            Inspection = Inspection + "all mapping in polyA; "
                    except:
                        print(img)
                # 2.1.2 ##   ->-   ##
                #       ## -->     ##
                else:
                    # (1) Arrow Conflict #
                    if clip_min[0] != "N" and clip_min[1] != "N":
                        if abs(int(clip_min[0]) - int(clip_min[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    if int(clip_max[0]) < int(arrow_max[0]):
                        Inspection = Inspection + "Arrow with SVA tail; "
                    # (2) Alignment Conflict #
                    if align_color[1].count("N") != 2 and abs(int(align_min[1][0]) - int(align_min[3][0])) > 500 and abs(int(align_min[1][1]) - int(align_min[3][1])) > 500:
                        Inspection = Inspection + "Distance between alignments too long; "
                    if (int(align_min[0][0]) - int(align_min[3][0]) > 5) and (int(align_min[0][1]) - int(align_min[3][1]) > 5):
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
                        Inspection = Inspection + "Arrow with SVA tail; "
                    # (2) Alignment Conflict #
                    if align_color[3].count("N") != 2 and abs(int(align_min[1][0]) - int(align_min[3][0])) > 500 and abs(int(align_min[1][1]) - int(align_min[3][1])) > 500:
                        Inspection = Inspection + "Distance between alignments too long; "
                    if (int(align_min[1][0]) > (1350 + 159) and int(align_min[1][1]) > (1345 +159)) or (align_color[3].count("N") != 2 and int(align_min[3][0]) > (1350 + 159) and int(align_min[3][1]) > (1345 +159)):
                        Inspection = Inspection + "all mapping in polyA; "
                    if (int(align_min[2][0]) - int(align_min[1][0]) > 5) and (int(align_min[2][1]) - int(align_min[1][1]) > 5):
                        Inspection = Inspection + "Junction site difference; "
                # 2.1.4 ## -->     ##
                #       ##   -->   ##
                else:
                    # (1) Arrow Conflict #
                    if clip_min[0] != "N" and clip_min[1] != "N":
                        if abs(int(clip_min[0]) - int(clip_min[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    # (2) Alignment Conflict #
                    if abs(int(align_min[1][0]) - int(align_min[3][0])) > 500 and abs(int(align_min[1][1]) - int(align_min[3][1])) > 500:
                        Inspection = Inspection + "Distance between alignments too long; "
                    if (int(align_min[1][0]) > (1350 + 159) and int(align_min[1][1]) > (1345 +159)) or (int(align_min[3][0]) > (1350 + 159) and int(align_min[3][1]) > (1345 +159)):
                        Inspection = Inspection + "all mapping in polyA; "
        # 2.2 ##     <-- ##
        #     ##   <--   ##
        else:
            if clip_map[0] == "mapped":
                # 2.2.1 ##     -<- ##
                #       ##     -<- ##
                if clip_map[1] == "mapped":
                    # (1) Arrow Conflict #
                    if clip_max[0] != "N" and clip_max[1] != "N":
                        if abs(int(clip_max[0]) - int(clip_max[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    if int(clip_min[0]) > int(arrow_min[0]) or int(clip_min[1]) > int(arrow_min[1]):
                        Inspection = Inspection + "Arrow with SVA tail; "
                    # (2) Alignment Conflict #
                    if align_color[1].count("N") != 2:
                        if align_color[3].count("N") != 2:
                            ## Two Clip-Pair ##
                            if (int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159)) or (int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159)):
                                Inspection = Inspection + "all mapping in polyA; "
                            if (int(align_max[1][0]) < (1360 + 159) and int(align_max[1][1]) < (1355 + 159)) or (int(align_max[3][0]) < (1360 + 159) and int(align_max[3][1]) < (1355 + 159)):
                                Inspection = Inspection + "No polyA; "
                            else:
                                if (int(align_max[1][0]) < (1370 + 159) and int(align_max[1][1]) < (1365 + 159)) or (int(align_max[3][0]) < (1370 + 159) and int(align_max[3][1]) < (1365 + 159)):
                                    Inspection = Inspection + "PolyA too short; "
                                tmp_ployA_len = max((int(align_max[1][0]) - 1359 - 159),(int(align_max[1][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                                tmp_ployA_len = max((int(align_max[3][0]) - 1359 - 159),(int(align_max[3][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                            if abs(int(align_max[1][0]) - int(align_max[3][0])) > 5 and abs(int(align_max[1][1]) - int(align_max[3][1])) > 5:
                                Inspection = Inspection + "Junction site difference; "
                        else:
                            ## Clip-Pair + Split ##
                            if (int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159)) or (int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159)):
                                Inspection = Inspection + "all mapping in polyA; "
                            if (int(align_max[1][0]) < (1360 + 159) and int(align_max[1][1]) < (1355 + 159)) or (int(align_max[2][0]) < (1360 + 159) and int(align_max[2][1]) < (1355 + 159)):
                                Inspection = Inspection + "No polyA; "
                            else:
                                if (int(align_max[1][0]) < (1370 + 159) and int(align_max[1][1]) < (1365 + 159)) or (int(align_max[2][0]) < (1370 + 159) and int(align_max[2][1]) < (1365 + 159)):
                                    Inspection = Inspection + "PolyA too short; "
                                tmp_ployA_len = max((int(align_max[1][0]) - 1359 - 159),(int(align_max[1][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                                tmp_ployA_len = max((int(align_max[2][0]) - 1359 - 159),(int(align_max[2][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                            if abs(int(align_max[1][0]) - int(align_max[2][0])) > 5 and abs(int(align_max[1][1]) - int(align_max[2][1])) > 5:
                                Inspection = Inspection + "Junction site difference; "
                    else:
                        ## Split + Clip-Pair ##
                        if align_color[3].count("N") != 2:
                            if (int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159)) or (int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159)):
                                Inspection = Inspection + "all mapping in polyA; "
                            if (int(align_max[0][0]) < (1360 + 159) and int(align_max[0][1]) < (1355 + 159)) or (int(align_max[3][0]) < (1360 + 159) and int(align_max[3][1]) < (1355 + 159)):
                                Inspection = Inspection + "No polyA; "
                            else:
                                if (int(align_max[0][0]) < (1370 + 159) and int(align_max[0][1]) < (1365 + 159)) or (int(align_max[3][0]) < (1370 + 159) and int(align_max[3][1]) < (1365 + 159)):
                                    Inspection = Inspection + "PolyA too short; "
                                tmp_ployA_len = max((int(align_max[0][0]) - 1359 - 159),(int(align_max[0][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                                tmp_ployA_len = max((int(align_max[3][0]) - 1359 - 159),(int(align_max[3][1]) - 1355 - 159))
                                PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                            if abs(int(align_max[0][0]) - int(align_max[3][0])) > 5 and abs(int(align_max[0][1]) - int(align_max[3][1])) > 5:
                                Inspection = Inspection + "Junction site difference; "
                        ## Two Split ##
                        else:
                            try:
                                if (int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159)) or (int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159)):
                                    Inspection = Inspection + "all mapping in polyA; "
                                if (int(align_max[0][0]) < (1360 + 159) and int(align_max[0][1]) < (1355 + 159)) or (int(align_max[2][0]) < (1360 + 159) and int(align_max[2][1]) < (1355 + 159)):
                                    Inspection = Inspection + "No polyA; "
                                else:
                                    if (int(align_max[0][0]) < (1370 + 159) and int(align_max[0][1]) < (1365 + 159)) or (int(align_max[2][0]) < (1370 + 159) and int(align_max[2][1]) < (1365 + 159)):
                                        Inspection = Inspection + "PolyA too short; "
                                    tmp_ployA_len = max((int(align_max[0][0]) - 1359 - 159),(int(align_max[0][1]) - 1355 - 159))
                                    PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                                    tmp_ployA_len = max((int(align_max[2][0]) - 1359 - 159),(int(align_max[2][1]) - 1355 - 159))
                                    PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                                if abs(int(align_max[0][0]) - int(align_max[2][0])) > 5 and abs(int(align_max[0][1]) - int(align_max[2][1])) > 5:
                                    Inspection = Inspection + "Junction site difference; "
                            except:
                                print(img)
                # 2.2.2 ##   -<-   ##
                #       ##     <-- ##
                else:
                    # (1) Arrow Conflict #
                    if clip_max[0] != "N" and clip_max[1] != "N":
                        if abs(int(clip_max[0]) - int(clip_max[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    if int(clip_min[0]) > int(arrow_min[0]):
                        Inspection = Inspection + "Arrow with SVA tail; "
                    # (2) Alignment Conflict #
                    if align_color[1].count("N") != 2:
                        ## Clip-Pair + Pair
                        if abs(int(align_max[0][0]) - int(align_max[2][0])) > 500 and abs(int(align_max[0][1]) - int(align_max[2][1])) > 500:
                            Inspection = Inspection + "Distance between alignments too long; "
                        if (int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159)) or (int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159)):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[1][0]) < (1360 + 159) and int(align_max[1][1]) < (1355 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[1][0]) < (1370 + 159) and int(align_max[1][1]) < (1365 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            tmp_ployA_len = max((int(align_max[1][0]) - 1359 - 159),(int(align_max[1][1]) - 1355 - 159))
                            PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                        if (int(align_max[2][0]) - int(align_max[1][0]) > 5) and (int(align_max[2][1]) - int(align_max[1][1]) > 5):
                            Inspection = Inspection + "Junction site difference; "
                    else:
                        ## Split + Pair ##
                        if (int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159)) or (int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159)):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[0][0]) < (1360 + 159) and int(align_max[0][1]) < (1355 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[0][0]) < (1370 + 159) and int(align_max[0][1]) < (1365 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            tmp_ployA_len = max((int(align_max[0][0]) - 1359 - 159),(int(align_max[0][1]) - 1355 - 159))
                            PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                        if int(align_max[2][0]) < (1381 + 159 - 400) and int(align_max[2][1]) < (1374 + 159 - 400):
                            Inspection = Inspection + "Large gap to PolyA; "
                        if (int(align_max[2][0]) - int(align_max[0][0]) > 5) and (int(align_max[2][1]) - int(align_max[0][1]) > 5):
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
                        Inspection = Inspection + "Arrow with SVA tail; "
                    # (2) Alignment Conflict #
                    if align_color[3].count("N") != 2:
                        ##  Pair + Clip-Pair ##
                        if abs(int(align_max[0][0]) - int(align_max[2][0])) > 500 and abs(int(align_max[0][1]) - int(align_max[2][1])) > 500:
                            Inspection = Inspection + "Distance between alignments too long; "
                        if (int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159)) or (int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159)):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[3][0]) < (1360 + 159) and int(align_max[3][1]) < (1355 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[3][0]) < (1370 + 159) and int(align_max[3][1]) < (1365 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            tmp_ployA_len = max((int(align_max[3][0]) - 1359 - 159),(int(align_max[3][1]) - 1355 - 159))
                            PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                        if (int(align_max[0][0]) - int(align_max[3][0]) > 5) and (int(align_max[0][1]) - int(align_max[3][1]) > 5):
                            Inspection = Inspection + "Junction site difference; "
                    else:
                        ## Pair + Split ##
                        if (int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159)) or (int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159)):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[2][0]) < (1360 + 159) and int(align_max[2][1]) < (1355 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[2][0]) < (1370 + 159) and int(align_max[2][1]) < (1365 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            tmp_ployA_len = max((int(align_max[2][0]) - 1359 - 159),(int(align_max[2][1]) - 1355 - 159))
                            PolyA_len = PolyA_len + "PolyA " + str(tmp_ployA_len) + "bp; "
                        if int(align_max[0][0]) < (1381 + 159 - 400) and int(align_max[0][1]) < (1374 + 159 - 400):
                            Inspection = Inspection + "Large gap to PolyA; "
                        if (int(align_max[0][0]) - int(align_max[2][0]) > 5) and (int(align_max[0][1]) - int(align_max[2][1]) > 5):
                            Inspection = Inspection + "Junction site difference; "
                # 2.1.4 ##     <-- ##
                #       ##    <--  ##
                else:
                    ## Two Pair-End ##
                    # (1) Arrow Conflict #
                    if clip_max[0] != "N" and clip_max[1] != "N":
                        if abs(int(clip_max[0]) - int(clip_max[1])) > 5:
                            Inspection = Inspection + "Clip site difference"
                    # (2) Alignment Conflict #
                    if (int(align_min[0][0]) > (1350 + 159) and int(align_min[0][1]) > (1345 +159)) or (int(align_min[2][0]) > (1350 + 159) and int(align_min[2][1]) > (1345 +159)):
                        Inspection = Inspection + "all mapping in polyA; "
                    if (int(align_max[0][0]) < (1381 + 159 - 400) and int(align_max[0][1]) < (1374 + 159 - 400)) or (int(align_max[2][0]) < (1381 + 159 - 400) and int(align_max[2][1]) < (1374 + 159 - 400)):
                        Inspection = Inspection + "Large gap to PolyA; "
                    if abs(int(align_max[0][0]) - int(align_max[2][0])) > 500 and abs(int(align_max[0][1]) - int(align_max[2][1])) > 500:
                        Inspection = Inspection + "Distance between alignments too long; "

    PCR_count = 0
    if reads_direction[0] == reads_direction[1]:
        if align_color[0] == align_color[2] and align_color[1] == align_color[3]:
            if abs(int(arrow_min[0]) - int(arrow_min[1])) < 5:
                PCR_count = PCR_count + 1
            if abs(int(arrow_max[0]) - int(arrow_max[1])) < 5:
                PCR_count = PCR_count + 1
            if align_color[0].count("N") != 2:
                if (abs(int(align_min[0][0]) - int(align_min[2][0])) < 5) and (abs(int(align_min[0][1]) - int(align_min[2][1])) < 5):
                    PCR_count = PCR_count + 1
                if (abs(int(align_max[0][0]) - int(align_max[2][0])) < 5) and (abs(int(align_max[0][1]) - int(align_max[2][1])) < 5):
                    PCR_count = PCR_count + 1
            if align_color[1].count("N") != 2:
                if (abs(int(align_min[1][0]) - int(align_min[3][0])) < 5) and (abs(int(align_min[1][1]) - int(align_min[3][1])) < 5):
                    PCR_count = PCR_count + 1
                if (abs(int(align_max[1][0]) - int(align_max[3][0])) < 5) and (abs(int(align_max[1][1]) - int(align_max[3][1])) < 5):
                    PCR_count = PCR_count + 1
    if PCR_count >= 3:
        Inspection = Inspection + "PCR duplicate? "

## Overlapping alignment ##
#[17,35,53,71]
#[17,45,73,101]
if align_color[0].count("N") != 2:
    if align_color[2].count("N") != 2:
        overlap_mappability = -1
        for e in [0,1]:
            # Overlapping #
            if align_color[0][e] != "N" and align_color[2][e] != "N":
                if int(align_max[0][e]) > int(align_min[2][e]) and int(align_min[0][e]) < int(align_max[2][e]) and abs(max(int(align_min[0][e]),int(align_min[2][e])) - min(int(align_max[0][e]),int(align_max[2][e]))) > 20 and max(int(align_min[0][e]),int(align_min[2][e])) > 210:
                    overlap_notEqual_count = 0
                    for k in range(max(int(align_min[0][e]),int(align_min[2][e])),min(int(align_max[0][e]),int(align_max[2][e]))+1):
                        if temp_tensor[:,17+e*5:21+e*5,k].equal(temp_tensor[:,53+e*5:57+e*5,k]):
                            overlap_notEqual_count = overlap_notEqual_count + 1
                    map_tmp = overlap_notEqual_count/(min(int(align_max[0][e]),int(align_max[2][e])) - max(int(align_min[0][e]),int(align_min[2][e])) + 1)
                    if map_tmp > overlap_mappability:
                        overlap_mappability = map_tmp
        if overlap_mappability != -1 and overlap_mappability < 0.9:
            Inspection = Inspection + "Overlap alignment is too low; "
    if align_color[3].count("N") != 2:
        overlap_mappability = -1
        for e in [0,1]:
            # Overlapping #
            if align_color[0][e] != "N" and align_color[3][e] != "N":
                if int(align_max[0][e]) > int(align_min[3][e]) and int(align_min[0][e]) < int(align_max[3][e]) and abs(max(int(align_min[0][e]),int(align_min[3][e])) - min(int(align_max[0][e]),int(align_max[3][e]))) > 20 and max(int(align_min[0][e]),int(align_min[3][e])) > 210:
                    overlap_notEqual_count = 0
                    for k in range(max(int(align_min[0][e]),int(align_min[3][e])),min(int(align_max[0][e]),int(align_max[3][e]))+1):
                        if temp_tensor[:,17+e*5:21+e*5,k].equal(temp_tensor[:,71+e*5:75+e*5,k]):
                            overlap_notEqual_count = overlap_notEqual_count + 1
                    map_tmp = overlap_notEqual_count/(min(int(align_max[0][e]),int(align_max[3][e])) - max(int(align_min[0][e]),int(align_min[3][e])) + 1)
                    if map_tmp > overlap_mappability:
                        overlap_mappability = map_tmp
        if overlap_mappability != -1 and overlap_mappability < 0.9:
            Inspection = Inspection + "Overlap alignment is too low; "
if align_color[1].count("N") != 2:
    if align_color[2].count("N") != 2:
        overlap_mappability = -1
        for e in [0,1]:
            # Overlapping #
            if align_color[1][e] != "N" and align_color[2][e] != "N":
                if int(align_max[1][e]) > int(align_min[2][e]) and int(align_min[1][e]) < int(align_max[2][e]) and abs(max(int(align_min[1][e]),int(align_min[2][e])) - min(int(align_max[1][e]),int(align_max[2][e]))) > 20 and max(int(align_min[1][e]),int(align_min[2][e])) > 210:
                    overlap_notEqual_count = 0
                    for k in range(max(int(align_min[1][e]),int(align_min[2][e])),min(int(align_max[1][e]),int(align_max[2][e]))+1):
                        if temp_tensor[:,35+e*5:39+e*5,k].equal(temp_tensor[:,53+e*5:57+e*5,k]):
                            overlap_notEqual_count = overlap_notEqual_count + 1
                    map_tmp = overlap_notEqual_count/(min(int(align_max[1][e]),int(align_max[2][e])) - max(int(align_min[1][e]),int(align_min[2][e])) + 1)
                    if map_tmp > overlap_mappability:
                        overlap_mappability = map_tmp
        if overlap_mappability != -1 and overlap_mappability < 0.9:
            Inspection = Inspection + "Overlap alignment is too low; "
    if align_color[3].count("N") != 2:
        overlap_mappability = -1
        for e in [0,1]:
            # Overlapping #
            if align_color[1][e] != "N" and align_color[3][e] != "N":
                if int(align_max[1][e]) > int(align_min[3][e]) and int(align_min[1][e]) < int(align_max[3][e]) and abs(max(int(align_min[1][e]),int(align_min[3][e])) - min(int(align_max[1][e]),int(align_max[3][e]))) > 20 and max(int(align_min[1][e]),int(align_min[3][e])) > 210:
                    overlap_notEqual_count = 0
                    for k in range(max(int(align_min[1][e]),int(align_min[3][e])),min(int(align_max[1][e]),int(align_max[3][e]))+1):
                        if temp_tensor[:,35+e*5:39+e*5,k].equal(temp_tensor[:,71+e*5:75+e*5,k]):
                            overlap_notEqual_count = overlap_notEqual_count + 1
                    map_tmp = overlap_notEqual_count/(min(int(align_max[1][e]),int(align_max[3][e])) - max(int(align_min[1][e]),int(align_min[3][e])) + 1)
                    if map_tmp > overlap_mappability:
                        overlap_mappability = map_tmp
        if overlap_mappability != -1 and overlap_mappability < 0.9:
            Inspection = Inspection + "Overlap alignment is too low; "

## read1 and read2 mapability ##
read_align = -1
if align_color[0].count("N") != 2:
    for e in [0,1]:
        if align_color[0][e] != "N" and abs(int(align_max[0][e]) - int(align_min[0][e])) > 20 and align_min[0][e] > 210:
            if float(align_mapability[0][e]) > read_align:
                read_align = float(align_mapability[0][e])
    if read_align != -1 and read_align < 0.9:
        Inspection = Inspection + "read1 mappability is lower than 90%; "
read_align = -1
if align_color[1].count("N") != 2:
    for e in [0,1]:
        if align_color[1][e] != "N" and abs(int(align_max[1][e]) - int(align_min[1][e])) > 20 and align_min[1][e] > 210:
            if float(align_mapability[1][e]) > read_align:
                read_align = float(align_mapability[1][e])
    if read_align != -1 and read_align < 0.9:
        Inspection = Inspection + "read1 mappability is lower than 90%; "
read_align = -1
if align_color[2].count("N") != 2:
    for e in [0,1]:
        if align_color[2][e] != "N" and abs(int(align_max[2][e]) - int(align_min[2][e])) > 20 and align_min[2][e] > 210:
            if float(align_mapability[2][e]) > read_align:
                read_align = float(align_mapability[2][e])
    if read_align != -1 and read_align < 0.9:
        Inspection = Inspection + "read2 mappability is lower than 90%; "
read_align = -1
if align_color[3].count("N") != 2:
    for e in [0,1]:
        if align_color[3][e] != "N" and abs(int(align_max[3][e]) - int(align_min[3][e])) > 20 and align_min[3][e] > 210:
            if float(align_mapability[3][e]) > read_align:
                read_align = float(align_mapability[3][e])
    if read_align != -1 and read_align < 0.9:
        Inspection = Inspection + "read2 mappability is lower than 90%; "

if Inspection == "":
    Inspection = "PASS"

## Mark transductions ##
transduction3_mark = ""
if transduction3_color[0].count("N") != 2:
    clip_from = max(align_max[0])
    if transduction3_color[0].count("green") != 0:
        transduction_len = 0
        for e in [0,1]:
            if transduction3_min[0][e] != "N" and transduction3_max[0][e] != "N" and abs(int(transduction3_min[0][e]) - int(transduction3_max[0][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[0][e]) - int(transduction3_max[0][e]))
        transduction3_mark = transduction3_mark + "from:" + str(int(clip_from - 159)) + " " + reads_direction[0] + " " + " read1 3'green " + str(transduction_len) + "bp; "
    elif transduction3_color[0].count("black") != 0:
        transduction_len = 0
        for e in [0,1]:
            if transduction3_min[0][e] != "N" and transduction3_max[0][e] != "N" and abs(int(transduction3_min[0][e]) - int(transduction3_max[0][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[0][e]) - int(transduction3_max[0][e]))
        transduction3_mark = transduction3_mark + "from:" + str(int(clip_from - 159)) + " " + reads_direction[0] + " " + "read1 3'black clip " + str(transduction_len) + "bp; "
if transduction3_color[1].count("N") != 2:
    clip_from = max(align_max[1])
    if transduction3_color[1].count("green") != 0:
        transduction_len = 0
        for e in [0,1]:
            if transduction3_min[1][e] != "N" and transduction3_max[1][e] != "N" and abs(int(transduction3_min[1][e]) - int(transduction3_max[1][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[1][e]) - int(transduction3_max[1][e]))
        transduction3_mark = transduction3_mark + "from:" + str(int(clip_from - 159)) + " " + reads_direction[0] + " " + "read1 3'green " + str(transduction_len) + "bp; "
    elif transduction3_color[1].count("black") != 0:
        transduction_len = 0
        for e in [0,1]:
            if transduction3_min[1][e] != "N" and transduction3_max[1][e] != "N" and abs(int(transduction3_min[1][e]) - int(transduction3_max[1][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[1][e]) - int(transduction3_max[1][e]))
        transduction3_mark = transduction3_mark + "from:" + str(int(clip_from - 159)) + " " + reads_direction[0] + " " + "read1 3'black clip " + str(transduction_len) + "bp; "
if transduction3_color[2].count("N") != 2:
    clip_from = max(align_max[2])
    if transduction3_color[2].count("green") != 0:
        transduction_len = 0
        for e in [0,1]:
            if transduction3_min[2][e] != "N" and transduction3_max[2][e] != "N" and abs(int(transduction3_min[2][e]) - int(transduction3_max[2][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[2][e]) - int(transduction3_max[2][e]))
        transduction3_mark = transduction3_mark + "from:" + str(int(clip_from - 159)) + " " + reads_direction[1] + " " + "read2 3'green " + str(transduction_len) + "bp; "
    elif transduction3_color[2].count("black") != 0:
        transduction_len = 0
        for e in [0,1]:
            if transduction3_min[2][e] != "N" and transduction3_max[2][e] != "N" and abs(int(transduction3_min[2][e]) - int(transduction3_max[2][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[2][e]) - int(transduction3_max[2][e]))
        transduction3_mark = transduction3_mark + "from:" + str(int(clip_from - 159)) + " " + reads_direction[1] + " " + "read2 3'black clip " + str(transduction_len) + "bp; "
if transduction3_color[3].count("N") != 2:
    clip_from = max(align_max[3])
    if transduction3_color[3].count("green") != 0:
        transduction_len = 0
        for e in [0,1]:
            if transduction3_min[3][e] != "N" and transduction3_max[3][e] != "N" and abs(int(transduction3_min[3][e]) - int(transduction3_max[3][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[3][e]) - int(transduction3_max[3][e]))
        transduction3_mark = transduction3_mark + "from:" + str(int(clip_from - 159)) + " " + reads_direction[1] + " " + "read2 3'green " + str(transduction_len) + "bp; "
    elif transduction3_color[2].count("black") != 0:
        transduction_len = 0
        for e in [0,1]:
            if transduction3_min[3][e] != "N" and transduction3_max[3][e] != "N" and abs(int(transduction3_min[3][e]) - int(transduction3_max[3][e])) > transduction_len:
                transduction_len = abs(int(transduction3_min[3][e]) - int(transduction3_max[3][e]))
        transduction3_mark = transduction3_mark + "from:" + str(int(clip_from - 159)) + " " + reads_direction[1] + " " + "read2 3'black clip " + str(transduction_len) + "bp; "

f = open("Inspection_tmp_SVA_{}.txt".format(ver),"w")
f.write("{}\t{}\t{}\t".format(Inspection,PolyA_len,transduction3_mark))
f.close()
