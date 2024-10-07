import os
import sys
import numpy
import torch
import torchvision
from PIL import Image

Totensor=torchvision.transforms.Compose([torchvision.transforms.ToTensor()])

# Initializing #
LINE_ref = "GGGAGGAGGAGCCAAGATGGCCGAATAGGAACAGCTCCGGTCTACAGCTCCCAGCGTGAGCGACGCAGAAGACGGGTGATTTCTGCATTTCCATCTGAGGTACCGGGTTCATCTCACTAGGGAGTGCCAGACAGTGGGCGCAGGCCAGTGTGTGTGCGCACCGTGCGCGAGCCGAAGCAGGGCGAGGCATTGCCTCACCTGGGAAGCGCAAGGGGTCAGGGAGTTCCCTTTCCGAGTCAAAGAAAGGGGTGACGGACGCACCTGGAAAATCGGGTCACTCCCACCCGAATATTGCGCTTTTCAGACCGGCTTAAGAAACGGCGCACCACGAGACTATATCCCACACCTGGCTCAGAGGGTCCTACGCCCACGGAATCTCGCTGATTGCTAGCACAGCAGTCTGAGATCAAACTGCAAGGCGGCAACGAGGCTGGGGGAGGGGCGCCCGCCATTGCCCAGGCTTGCTTAGGTAAACAAAGCAGCCGGGAAGCTCGAACTGGGTGGAGCCCACCACAGCTCAAGGAGGCCTGCCTGCCTCTGTAGGCTCCACCTCTGGGGGCAGGGCACAGACAAACAAAAAGACAGCAGTAACCTCTGCAGACTTAAGTGTCCCTGTCTGACAGCTTTGAAGAGAGCAGTGGTTCTCCCAGCACGCAGCTGGAGATCTGAGAACGGGCAGACTGCCTCCTCAAGTGGGTCCCTGACCCCTGACCCCCGAGCAGCCTAACTGGGAGGCACCCCCCAGCAGGGGCACACTGACACCTCACACGGCAGGGTATTCCAACAGACCTGCAGCTGAGGGTCCTGTCTGTTAGAAGGAAAACTAACAACCAGAAAGGACATCTACACCGAAAACCCATCTGTACATCACCATCATCAAAGACCAAAAGTAGATAAAACCACAAAGATGGGGAAAAAACAGAACAGAAAAACTGGAAACTCTAAAACGCAGAGCGCCTCTCCTCCTCCAAAGGAACGCAGTTCCTCACCAGCAACAGAACAAAGCTGGATGGAGAATGATTTTGACGAGCTGAGAGAAGAAGGCTTCAGACGATCAAATTACTCTGAGCTACGGGAGGACATTCAAACCAAAGGCAAAGAAGTTGAAAACTTTGAAAAAAATTTAGAAGAATGTATAACTAGAATAACCAATACAGAGAAGTGCTTAAAGGAGCTGATGGAGCTGAAAACCAAGGCTCGAGAACTACGTGAAGAATGCAGAAGCCTCAGGAGCCGATGCGATCAACTGGAAGAAAGGGTATCAGCAATGGAAGATGAAATGAATGAAATGAAGCGAGAAGGGAAGTTTAGAGAAAAAAGAATAAAAAGAAATGAGCAAAGCCTCCAAGAAATATGGGACTATGTGAAAAGACCAAATCTACGTCTGATTGGTGTACCTGAAAGTGATGTGGAGAATGGAACCAAGTTGGAAAACACTCTGCAGGATATTATCCAGGAGAACTTCCCCAATCTAGCAAGGCAGGCCAACGTTCAGATTCAGGAAATACAGAGAACGCCACAAAGATACTCCTCGAGAAGAGCAACTCCAAGACACATAATTGTCAGATTCACCAAAGTTGAAATGAAGGAAAAAATGTTAAGGGCAGCCAGAGAGAAAGGTCGGGTTACCCTCAAAGGAAAGCCCATCAGACTAACAGCGGATCTCTCGGCAGAAACCCTACAAGCCAGAAGAGAGTGGGGGCCAATATTCAACATTCTTAAAGAAAAGAATTTTCAACCCAGAATTTCATATCCAGCCAAACTAAGCTTCATAAGTGAAGGAGAAATAAAATACTTTATAGACAAGCAAATGCTGAGAGATTTTGTCACCACCAGGCCTGCCCTAAAAGAGCTCCTGAAGGAAGCGCTAAACATGGAAAGGAACAACCGGTACCAGCCGCTGCAAAATCATGCCAAAATGTAAAGACCATCGAGACTAGGAAGAAACTGCATCAACTAATGAGCAAAATCACCAGCTAACATCATAATGACAGGATCAAATTCACACATAACAATATTAACTTTAAATATAAATGGACTAAATTCTGCAATTAAAAGACACAGACTGGCAAGTTGGATAAAGAGTCAAGACCCATCAGTGTGCTGTATTCAGGAAACCCATCTCACGTGCAGAGACACACATAGGCTCAAAATAAAAGGATGGAGGAAGATCTACCAAGCCAATGGAAAACAAAAAAAGGCAGGGGTTGCAATCCTAGTCTCTGATAAAACAGACTTTAAACCAACAAAGATCAAAAGAGACAAAGAAGGCCATTACATAATGGTAAAGGGATCAATTCAACAAGAGGAGCTAACTATCCTAAATATTTATGCACCCAATACAGGAGCACCCAGATTCATAAAGCAAGTCCTCAGTGACCTACAAAGAGACTTAGACTCCCACACATTAATAATGGGAGACTTTAACACCCCACTGTCAACATTAGACAGATCAACGAGACAGAAAGTCAACAAGGATACCCAGGAATTGAACTCAGCTCTGCACCAAGCAGACCTAATAGACATCTACAGAACTCTCCACCCCAAATCAACAGAATATACATTTTTTTCAGCACCACACCACACCTATTCCAAAATTGACCACATAGTTGGAAGTAAAGCTCTCCTCAGCAAATGTAAAAGAACAGAAATTATAACAAACTATCTCTCAGACCACAGTGCAATCAAACTAGAACTCAGGATTAAGAATCTCACTCAAAGCCGCTCAACTACATGGAAACTGAACAACCTGCTCCTGAATGACTACTGGGTACATAACGAAATGAAGGCAGAAATAAAGATGTTCTTTGAAACCAACGAGAACAAAGACACCACATACCAGAATCTCTGGGACGCATTCAAAGCAGTGTGTAGAGGGAAATTTATAGCACTAAATGCCTACAAGAGAAAGCAGGAAAGATCCAAAATTGACACCCTAACATCACAATTAAAAGAACTAGAAAAGCAAGAGCAAACACATTCAAAAGCTAGCAGAAGGCAAGAAATAACTAAAATCAGAGCAGAACTGAAGGAAATAGAGACACAAAAAACCCTTCAAAAAATCAATGAATCCAGGAGCTGGTTTTTTGAAAGGATCAACAAAATTGATAGACCGCTAGCAAGACTAATAAAGAAAAAAAGAGAGAAGAATCAAATAGACACAATAAAAAATGATAAAGGGGATATCACCACCGATCCCACAGAAATACAAACTACCATCAGAGAATACTACAAACACCTCTACGCAAATAAACTAGAAAATCTAGAAGAAATGGATACATTCCTCGACACATACACTCTCCCAAGACTAAACCAGGAAGAAGTTGAATCTCTGAATAGACCAATAACAGGCTCTGAAATTGTGGCAATAATCAATAGTTTACCAACCAAAAAGAGTCCAGGACCAGATGGATTCACAGCCGAATTCTACCAGAGGTACAAGGAGGAACTGGTACCATTCCTTCTGAAACTATTCCAATCAATAGAAAAAGAGGGAATCCTCCCTAACTCATTTTATGAGGCCAGCATCATTCTGATACCAAAGCCGGGCAGAGACACAACCAAAAAAGAGAATTTTAGACCAATATCCTTGATGAACATTGATGCAAAAATCCTCAATAAAATACTGGCAAACCGAATCCAGCAGCACATCAAAAAGCTTATCCACCATGATCAAGTGGGCTTCATCCCTGGGATGCAAGGCTGGTTCAATATACGCAAATCAATAAATGTAATCCAGCATATAAACAGAGCCAAAGACAAAAACCACATGATTATCTCAATAGATGCAGAAAAAGCCTTTGACAAAATTCAACAACCCTTCATGCTAAAAACTCTCAATAAATTAGGTATTGATGGGACGTATTTCAAAATAATAAGAGCTATCTATGACAAACCCACAGCCAATATCATACTGAATGGGCAAAAACTGGAAGCATTCCCTTTGAAAACTGGCACAAGACAGGGATGCCCTCTCTCACCGCTCCTATTCAACATAGTGTTGGAAGTTCTGGCCAGGGCAATCAGGCAGGAGAAGGAAATAAAGGGTATTCAATTAGGAAAAGAGGAAGTCAAATTGTCCCTGTTTGCAGACGACATGATTGTTTATCTAGAAAACCCCATCGTCTCAGCCCAAAATCTCCTTAAGCTGATAAGCAACTTCAGCAAAGTCTCAGGATACAAAATCAATGTACAAAAATCACAAGCATTCTTATACACCAACAACAGACAAACAGAGAGCCAAATCATGGGTGAACTCCCATTCACAATTGCTTCAAAGAGAATAAAATACCTAGGAATCCAACTTACAAGGGATGTGAAGGACCTCTTCAAGGAGAACTACAAACCACTGCTCAAGGAAATAAAAGAGGACACAAACAAATGGAAGAACATTCCATGCTCATGGGTAGGAAGAATCAATATCGTGAAAATGGCCATACTGCCCAAGGTAATTTACAGATTCAATGCCATCCCCATCAAGCTACCAATGACTTTCTTCACAGAATTGGAAAAAACTACTTTAAAGTTCATATGGAACCAAAAAAGAGCCCGCATCGCCAAGTCAATCCTAAGCCAAAAGAACAAAGCTGGAGGCATCACACTACCTGACTTCAAACTATACTACAAGGCTACAGTAACCAAAACAGCATGGTACTGGTACCAAAACAGAGATATAGATCAATGGAACAGAACAGAGCCCTCAGAAATAATGCCGCATATCTACAACTATCTGATCTTTGACAAACCTGAGAAAAACAAGCAATGGGGAAAGGATTCCCTATTTAATAAATGGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTATACAAAAATCAATTCAAGATGGATTAAAGATTTAAACGTTAGACCTAAAACCATAAAAACCCTAGAAGAAAACCTAGGCATTACCATTCAGGACATAGGCGTGGGCAAGGACTTCATGTCCAAAACACCAAAAGCAATGGCAACAAAAGCCAAAATTGACAAATGGGATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTACCATCAGAGTGAACAGGCAACCTACAACATGGGAGAAAATTTTCGCAACCTACTCATCTGACAAAGGGCTAATATCCAGAATCTACAATGAACTCAAACAAATTTACAAGAAAAAAACAAACAACCCCATCAAAAAGTGGGCGAAGGACATGAACAGACACTTCTCAAAAGAAGACATTTATGCAGCCAAAAAACACATGAAGAAATGCTCATCATCACTGGCCATCAGAGAAATGCAAATCAAAACCACTATGAGATATCATCTCACACCAGTTAGAATGGCAATCATTAAAAAGTCAGGAAACAACAGGTGCTGGAGAGGATGTGGAGAAATAGGAACACTTTTACACTGTTGGTGGGACTGTAAACTAGTTCAACCATTGTGGAAGTCAGTGTGGCGATTCCTCAGGGATCTAGAACTAGAAATACCATTTGACCCAGCCATCCCATTACTGGGTATATACCCAAAGGACTATAAATCATGCTGCTATAAAGACACATGCACACGTATGTTTATTGCGGCACTATTCACAATAGCAAAGACTTGGAACCAACCCAAATGTCCAACAATGATAGACTGGATTAAGAAAATGTGGCACATATACACCATGGAATACTATGCAGCCATAAAAAATGATGAGTTCATATCCTTTGTAGGGACATGGATGAAATTGGAAACCATCATTCTCAGTAAACTATCGCAAGAACAAAAAACCAAACACCGCATATTCTCACTCATAGGTGGGAATTGAACAATGAGATCACATGGACACAGGAAGGGGAATATCACACTCTGGGGACTGTGGTGGGGTCGGGGGAGGGGGGAGGGATAGCATTGGGAGATATACCTAATGCTAGATGACACGTTAGTGGGTGCAGCGCACCAGCATGGCACATGTATACATATGTAACTAACCTGCACAATGTGCACATGTACCCTAAAACTTAGAGTATAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAATAAAAA"
DNA_base = "ACTG"
cutoff = 0.5
outpath = "."
sub = "test"
ver = 4
TEclass = "LINE"

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
for j in [10,23,36,49]:
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
## Detect Alignment ##
align_min = []
align_max = []
align_color = []
align_mapability = []
#
for k in [17,30,43,56]:
    red_detect = Detect_black[k:k+4,:].eq(1).nonzero(as_tuple=True)[1]
    magenta_detect = Detect_black[k:k+4,:].eq(2).nonzero(as_tuple=True)[1]
    if red_detect.numel() != 0:
        black_count = Detect_black[k:k+4,min(red_detect):max(red_detect)].eq(0).nonzero(as_tuple=True)[1].numel()
        align_mapability.append((1 - black_count/red_detect.numel()))
        align_color.append("red")
        align_min.append(min(red_detect))
        align_max.append(max(red_detect))
    elif magenta_detect.numel() != 0:
        black_count = Detect_black[k:k+4,min(magenta_detect):max(magenta_detect)].eq(0).nonzero(as_tuple=True)[1].numel()
        align_mapability.append((1 - black_count/magenta_detect.numel()))
        align_color.append("magenta")
        align_min.append(min(magenta_detect))
        align_max.append(max(magenta_detect))
    else:
        align_mapability.append("N")
        align_color.append("N")
        align_min.append("N")
        align_max.append("N")

Inspection = ""
PolyA_len = ""
## Manual Inspection ##
# 1 # Upstream + Downstream
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
                    Inspection = Inspection + "Arrow with L1 tail; "
                if (int(arrow_min[0]) < int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_min[1])) or (int(arrow_min[0]) > int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_max[1])):
                    if (int(arrow_max[0]) - int(arrow_min[1])) > 25:
                        Inspection = Inspection + "TSD is longer than 25bp; "
                if int(arrow_min[0]) > int(arrow_max[1]):
                    Inspection = Inspection + "Arrow Conflict; "
                # (2) Alignment Conflict #
                if align_color[1] != "N":
                    # Two Clip-pair #
                    if align_color[3] != "N":
                        if int(align_min[1]) - int(align_max[2]) > 600:
                            Inspection = Inspection + "Alignments cross over; "
                        if int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[3]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[3]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[3]) - 6024 -159) + "bp; "
                        if int(align_min[0]) - int(align_min[3]) > 5 or int(align_max[0]) - int(align_max[3]) > 5:
                            Inspection = Inspection + "Junction site difference; "
                    # Clip-pair + Split #
                    else:
                        if int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[2]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[2]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[2]) - 6024 -159) + "bp; "
                        if int(align_min[0]) - int(align_min[2]) > 5 or int(align_max[0]) - int(align_max[2]) > 5:
                            Inspection = Inspection + "Junction site difference; "
                else:
                    # Split + Clip-pair #
                    if int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                    if align_color[3] != "N":
                        if int(align_max[3]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[3]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[3]) - 6024 -159) + "bp; "
                        if int(align_min[0]) - int(align_min[3]) > 5 or int(align_max[0]) - int(align_max[3]) > 5:
                            Inspection = Inspection + "Junction site difference; "
                    # Split + Split #
                    else:
                        if int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA"
                        if int(align_max[2]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[2]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[2]) - 6024 -159) + "bp; "
                        if int(align_min[0]) - int(align_min[2]) > 5 or int(align_max[0]) - int(align_max[2]) > 5:
                            Inspection = Inspection + "Junction site difference; "
            # 1.1.2 ## ->-     ##
            #       ##     <-- ##
            else:
                # (1) Arrow Conflict #
                if int(clip_max[0]) < int(arrow_max[0]):
                    Inspection = Inspection + "Arrow with L1 tail; "
                if (int(arrow_min[0]) < int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_min[1])) or (int(arrow_min[0]) > int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_max[1])):
                    if (int(arrow_max[0]) - int(arrow_min[1])) > 25:
                        Inspection = Inspection + "TSD is longer than 25bp; "
                if int(arrow_min[0]) > int(arrow_max[1]):
                    Inspection = Inspection + "Arrow Conflict; "
                # (2) Alignment Conflict #
                if align_color[1] != "N":
                    # Clip-pair + Pair #
                    if int(align_min[1]) - int(align_max[2]) > 600:
                        Inspection = Inspection + "Alignments cross over; "
                    if int(align_max[2]) < (6064 + 159 - 400):
                        Inspection = Inspection + "Large gap to PolyA; "
                    if int(align_min[2]) > (6012 + 159):
                        Inspection = Inspection + "all mapping in polyA; "
                    if int(align_min[0]) - int(align_min[2]) > 5:
                        Inspection = Inspection + "Junction site difference; "
                else:
                    # Split + Pair #
                    if int(align_max[2]) < (6064 + 159 - 400):
                        Inspection = Inspection + "Large gap to PolyA; "
                    if int(align_min[2]) > (6012 + 159):
                        Inspection = Inspection + "all mapping in polyA; "
                    if int(align_min[0]) - int(align_min[2]) > 5:
                        Inspection = Inspection + "Junction site difference; "
        else:
            # 1.1.3 ## -->     ##
            #       ##     -<- ##
            if clip_map[1] == "mapped":
                # (1) Arrow Conflict #
                if int(clip_min[1]) > int(arrow_min[1]):
                    Inspection = Inspection + "Arrow with L1 tail; "
                if (int(arrow_min[0]) < int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_min[1])) or (int(arrow_min[0]) > int(arrow_min[1]) and int(arrow_max[0]) > int(arrow_max[1])):
                    if (int(arrow_max[0]) - int(arrow_min[1])) > 25:
                        Inspection = Inspection + "TSD is longer than 25bp; "
                if int(arrow_min[0]) > int(arrow_max[1]):
                    Inspection = Inspection + "Arrow Conflict; "
                # (2) Alignment Conflict #
                if align_color[3] != "N":
                    # Pair + Clip-pair #
                    if int(align_min[1]) - int(align_max[2]) > 600:
                        Inspection = Inspection + "Alignments cross over; "
                    if int(align_min[2]) > (6012 + 159):
                        Inspection = Inspection + "all mapping in polyA; "
                    if int(align_max[3]) < (6024 + 159):
                        Inspection = Inspection + "No polyA; "
                    else:
                        if int(align_max[3]) < (6040 + 159):
                            Inspection = Inspection + "PolyA too short; "
                        PolyA_len = PolyA_len + "PolyA " + str(int(align_max[3]) - 6024 -159) + "bp; "
                    if int(align_max[1]) - int(align_max[3]) > 5:
                        Inspection = Inspection + "Junction site difference; "              
                else:
                    # Pair + Split #
                    if int(align_min[2]) > (6012 + 159):
                        Inspection = Inspection + "all mapping in polyA; "
                    if int(align_max[2]) < (6024 + 159):
                        Inspection = Inspection + "No polyA; "
                    else:
                        if int(align_max[2]) < (6040 + 159):
                            Inspection = Inspection + "PolyA too short; "
                        PolyA_len = PolyA_len + "PolyA " + str(int(align_max[2]) - 6024 -159) + "bp; "
                    if int(align_max[1]) - int(align_max[2]) > 5:
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
                if int(align_min[1]) - int(align_max[2]) > 600:
                    Inspection = Inspection + "Alignments cross over; "
                if int(align_max[2]) < (6064 + 159 - 400):
                    Inspection = Inspection + "Large gap to PolyA; "
                if int(align_min[2]) > (6012 + 159):
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
                    Inspection = Inspection + "Arrow with L1 tail; "
                if (int(arrow_min[1]) < int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_min[0])) or (int(arrow_min[1]) > int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_max[0])):
                    if (int(arrow_max[1]) - int(arrow_min[0])) > 25:
                        Inspection = Inspection + "TSD is longer than 25bp; "
                if int(arrow_min[1]) > int(arrow_max[0]):
                    Inspection = Inspection + "Arrow Conflict; "
                # (2) Alignment Conflict #
                if align_color[3] != "N":
                    # Two Clip-pair #
                    if align_color[1] != "N":
                        if int(align_min[3]) - int(align_max[0]) > 500:
                            Inspection = Inspection + "Alignments cross over; "
                        if int(align_min[0]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[1]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[1]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[1]) - 6024 -159) + "bp; "
                        if int(align_min[2]) - int(align_min[1]) > 5 or int(align_max[2]) - int(align_max[1]) > 5:
                            Inspection = Inspection + "Junction site difference; "
                    # Clip-pair + Split #
                    else:
                        if int(align_min[0]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[0]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[0]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[0]) - 6024 -159) + "bp; "
                        if int(align_min[2]) - int(align_min[0]) > 5 or int(align_max[2]) - int(align_max[0]) > 5:
                            Inspection = Inspection + "Junction site difference; "
                else:
                    # Split + Clip-pair #
                    if int(align_min[0]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA"
                    if align_color[1] != "N":
                        if int(align_max[1]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[1]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[1]) - 6024 -159) + "bp; "
                        if int(align_min[2]) - int(align_min[1]) > 5 or int(align_max[2]) - int(align_max[1]) > 5:
                            Inspection = Inspection + "Junction site difference; "
                    # Split + Split #
                    else:
                        if int(align_min[0]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[0]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[0]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[0]) - 6024 -159) + "bp; "
                        if int(align_min[2]) - int(align_min[0]) > 5 or int(align_max[2]) - int(align_max[0]) > 5:
                            Inspection = Inspection + "Junction site difference; "
            # 1.2.2 ##     -<- ##
            #       ## -->     ##
            else:
                if int(clip_min[0]) > int(arrow_min[0]):
                    Inspection = Inspection + "Arrow with L1 tail; "
                if (int(arrow_min[1]) < int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_min[0])) or (int(arrow_min[1]) > int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_max[0])):
                    if (int(arrow_max[1]) - int(arrow_min[0])) > 25:
                        Inspection = Inspection + "TSD is longer than 25bp; "
                if int(arrow_min[1]) > int(arrow_max[0]):
                    Inspection = Inspection + "Arrow Conflict; "
                # (2) Alignment Conflict #
                if align_color[1] != "N":
                    # Pair + Clip-pair #
                    if int(align_min[3]) - int(align_max[0]) > 600:
                        Inspection = Inspection + "Alignments cross over; "
                    if int(align_min[0]) > (6012 + 159):
                        Inspection = Inspection + "all mapping in polyA; "
                    if int(align_max[1]) < (6024 + 159):
                        Inspection = Inspection + "No polyA; "
                    else:
                        if int(align_max[1]) < (6040 + 159):
                            Inspection = Inspection + "PolyA too short; "
                        PolyA_len = PolyA_len + "PolyA " + str(int(align_max[1]) - 6024 -159) + "bp; "
                    if int(align_max[3]) - int(align_max[1]) > 5:
                        Inspection = Inspection + "Junction site difference; "                  
                else:
                    # Pair + Split #
                    if int(align_min[0]) > (6012 + 159):
                        Inspection = Inspection + "all mapping in polyA; "
                    if int(align_max[0]) < (6024 + 159):
                        Inspection = Inspection + "No polyA; "
                    else:
                        if int(align_max[0]) < (6040 + 159):
                            Inspection = Inspection + "PolyA too short; "
                        PolyA_len = PolyA_len + "PolyA " + str(int(align_max[0]) - 6024 -159) + "bp; "
                    if int(align_max[3]) - int(align_max[0]) > 5:
                        Inspection = Inspection + "Junction site difference; " 
        else:
            # 1.2.3 ##     <-- ##
            #       ## ->-     ##
            if clip_map[1] == "mapped":
                # (1) Arrow Conflict #
                if int(clip_max[1]) < int(arrow_max[1]):
                    Inspection = Inspection + "Arrow with L1 tail; "
                if (int(arrow_min[1]) < int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_min[0])) or (int(arrow_min[1]) > int(arrow_min[0]) and int(arrow_max[1]) > int(arrow_max[0])):
                    if (int(arrow_max[1]) - int(arrow_min[0])) > 25:
                        Inspection = Inspection + "TSD is longer than 25bp; "
                if int(arrow_min[1]) > int(arrow_max[0]):
                    Inspection = Inspection + "Arrow Conflict; "
                # (2) Alignment Conflict #
                if align_color[3] != "N":
                    # Clip-pair + Pair #
                    if int(align_min[3]) - int(align_max[0]) > 600:
                        Inspection = Inspection + "Alignments cross over; "
                    if int(align_max[0]) < (6064 + 159 - 400):
                        Inspection = Inspection + "Large gap to PolyA; "
                    if int(align_min[0]) > (6012 + 159):
                        Inspection = Inspection + "all mapping in polyA; "
                    if int(align_min[2]) - int(align_min[0]) > 5:
                        Inspection = Inspection + "Junction site difference; "
                else:
                    # Split + Pair #
                    if int(align_max[0]) < (6064 + 159 - 400):
                        Inspection = Inspection + "Large gap to PolyA; "
                    if int(align_min[0]) > (6012 + 159):
                        Inspection = Inspection + "all mapping in polyA; "
                    if int(align_min[2]) - int(align_min[0]) > 5:
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
                if int(align_min[3]) - int(align_max[0]) > 600:
                    Inspection = Inspection + "Alignments cross over; "
                if int(align_max[0]) < (6064 + 159 - 400):
                    Inspection = Inspection + "Large gap to PolyA; "
                if int(align_min[0]) > (6012 + 159):
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
                    Inspection = Inspection + "Arrow with L1 tail; "
                # (2) Alignment Conflict #
                if abs(int(align_min[0]) - int(align_min[2])) > 5:
                    Inspection = Inspection + "Junction site difference; "
                if (align_color[1] !="N" and int(align_min[1]) > (6012 + 159)) or (align_color[3] !="N" and int(align_min[3]) > (6012 + 159)):
                    Inspection = Inspection + "all mapping in polyA"
            # 2.1.2 ##   ->-   ##
            #       ## -->     ##
            else:
                # (1) Arrow Conflict #
                if clip_min[0] != "N" and clip_min[1] != "N":
                    if abs(int(clip_min[0]) - int(clip_min[1])) > 5:
                        Inspection = Inspection + "Clip site difference"
                if int(clip_max[0]) < int(arrow_max[0]):
                    Inspection = Inspection + "Arrow with L1 tail; "
                # (2) Alignment Conflict #
                if align_color[1] != "N" and abs(int(align_min[1]) - int(align_min[3])) > 500:
                    Inspection = Inspection + "Distance between alignments too long; "
                if int(align_min[0]) - int(align_min[3]) > 5:
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
                    Inspection = Inspection + "Arrow with L1 tail; "
                # (2) Alignment Conflict #
                if align_color[3] != "N" and abs(int(align_min[1]) - int(align_min[3])) > 500:
                    Inspection = Inspection + "Distance between alignments too long; "
                if int(align_min[1]) > (6012 + 159) or (align_color[3] !="N" and int(align_min[3]) > (6012 + 159)):
                    Inspection = Inspection + "all mapping in polyA; "
                if int(align_min[2]) - int(align_min[1]) > 5:
                    Inspection = Inspection + "Junction site difference; "
            # 2.1.4 ## -->     ##
            #       ##   -->   ##
            else:
                # (1) Arrow Conflict #
                if clip_min[0] != "N" and clip_min[1] != "N":
                    if abs(int(clip_min[0]) - int(clip_min[1])) > 5:
                        Inspection = Inspection + "Clip site difference"
                # (2) Alignment Conflict #
                if abs(int(align_min[1]) - int(align_min[3])) > 500:
                    Inspection = Inspection + "Distance between alignments too long; "
                if int(align_min[1]) > (6012 + 159) or int(align_min[3]) > (6012 + 159):
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
                    Inspection = Inspection + "Arrow with L1 tail; "
                # (2) Alignment Conflict #
                if align_color[1] != "N":
                    if align_color[3] != "N":
                        ## Two Clip-Pair ##
                        if int(align_min[0]) > (6012 + 159) or int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[1]) < (6024 + 159) or int(align_max[3]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[1]) < (6040 + 159) or int(align_max[3]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[1]) - 6024 -159) + "bp; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[3]) - 6024 -159) + "bp; "
                        if abs(int(align_max[1]) - int(align_max[3])) > 5:
                            Inspection = Inspection + "Junction site difference; "
                    else:
                        ## Clip-Pair + Split ##
                        if int(align_min[0]) > (6012 + 159) or int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[1]) < (6024 + 159) or int(align_max[2]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[1]) < (6040 + 159) or int(align_max[2]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[1]) - 6024 -159) + "bp; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[2]) - 6024 -159) + "bp; "
                        if abs(int(align_max[1]) - int(align_max[2])) > 5:
                            Inspection = Inspection + "Junction site difference; "
                else:
                    ## Split + Clip-Pair ##
                    if align_color[3] != "N":
                        if int(align_min[0]) > (6012 + 159) or int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[0]) < (6024 + 159) or int(align_max[3]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[0]) < (6040 + 159) or int(align_max[3]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[0]) - 6024 -159) + "bp; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[3]) - 6024 -159) + "bp; "
                        if abs(int(align_max[0]) - int(align_max[3])) > 5:
                            Inspection = Inspection + "Junction site difference; "
                    ## Two Split ##
                    else:
                        if int(align_min[0]) > (6012 + 159) or int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                        if int(align_max[0]) < (6024 + 159) or int(align_max[2]) < (6024 + 159):
                            Inspection = Inspection + "No polyA; "
                        else:
                            if int(align_max[0]) < (6040 + 159) or int(align_max[2]) < (6040 + 159):
                                Inspection = Inspection + "PolyA too short; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[0]) - 6024 -159) + "bp; "
                            PolyA_len = PolyA_len + "PolyA " + str(int(align_max[2]) - 6024 -159) + "bp; "
                        if abs(int(align_max[0]) - int(align_max[2])) > 5:
                            Inspection = Inspection + "Junction site difference; "
            # 2.2.2 ##   -<-   ##
            #       ##     <-- ##
            else:
                # (1) Arrow Conflict #
                if clip_max[0] != "N" and clip_max[1] != "N":
                    if abs(int(clip_max[0]) - int(clip_max[1])) > 5:
                        Inspection = Inspection + "Clip site difference"
                if int(clip_min[0]) > int(arrow_min[0]):
                    Inspection = Inspection + "Arrow with L1 tail; "
                # (2) Alignment Conflict #
                if align_color[1] != "N":
                    ## Clip-Pair + Pair
                    if abs(int(align_max[0]) - int(align_max[2])) > 500:
                        Inspection = Inspection + "Distance between alignments too long; "
                    if int(align_min[0]) > (6012 + 159) or int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                    if int(align_max[1]) < (6024 + 159):
                        Inspection = Inspection + "No polyA; "
                    else:
                        if int(align_max[1]) < (6040 + 159):
                            Inspection = Inspection + "PolyA too short; "
                        PolyA_len = PolyA_len + "PolyA " + str(int(align_max[1]) - 6024 -159) + "bp; "
                    if int(align_max[2]) - int(align_max[1]) > 5:
                        Inspection = Inspection + "Junction site difference; "
                else:
                    ## Split + Pair ##
                    if int(align_min[0]) > (6012 + 159) or int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                    if int(align_max[0]) < (6024 + 159):
                        Inspection = Inspection + "No polyA; "
                    else:
                        if int(align_max[0]) < (6040 + 159):
                            Inspection = Inspection + "PolyA too short; "
                        PolyA_len = PolyA_len + "PolyA " + str(int(align_max[0]) - 6024 -159) + "bp; "
                    if int(align_max[2]) < (6064 + 159 - 400):
                        Inspection = Inspection + "Large gap to PolyA; "
                    if int(align_max[2]) - int(align_max[0]) > 5:
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
                    Inspection = Inspection + "Arrow with L1 tail; "
                # (2) Alignment Conflict #
                if align_color[3] != "N":
                    ##  Pair + Clip-Pair ##
                    if abs(int(align_max[0]) - int(align_max[2])) > 500:
                        Inspection = Inspection + "Distance between alignments too long; "
                    if int(align_min[0]) > (6012 + 159) or int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                    if int(align_max[3]) < (6024 + 159):
                        Inspection = Inspection + "No polyA; "
                    else:
                        if int(align_max[3]) < (6040 + 159):
                            Inspection = Inspection + "PolyA too short; "
                        PolyA_len = PolyA_len + "PolyA " + str(int(align_max[3]) - 6024 -159) + "bp; "
                    if int(align_max[0]) - int(align_max[3]) > 5:
                        Inspection = Inspection + "Junction site difference; "
                else:
                    ## Pair + Split ##
                    if int(align_min[0]) > (6012 + 159) or int(align_min[2]) > (6012 + 159):
                            Inspection = Inspection + "all mapping in polyA; "
                    if int(align_max[2]) < (6024 + 159):
                        Inspection = Inspection + "No polyA; "
                    else:
                        if int(align_max[2]) < (6040 + 159):
                            Inspection = Inspection + "PolyA too short; "
                        PolyA_len = PolyA_len + "PolyA " + str(int(align_max[2]) - 6024 -159) + "bp; "
                    if int(align_max[0]) < (6064 + 159 - 400):
                        Inspection = Inspection + "Large gap to PolyA; "
                    if int(align_max[0]) - int(align_max[2]) > 5:
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
                if int(align_max[0]) < (6064 + 159 - 400) or int(align_max[2]) < (6064 + 159 - 400):
                    Inspection = Inspection + "Large gap to PolyA; "
                if abs(int(align_min[0]) - int(align_min[2])) > 500:
                    Inspection = Inspection + "Distance between alignments too long; "
                if int(align_min[0]) > (6012 + 159) or int(align_min[2]) > (6012 + 159):
                    Inspection = Inspection + "all mapping in polyA; "

## PCR duplication #
PCR_count = 0
if reads_direction[0] == reads_direction[1]:
    if align_color[0] == align_color[2] and align_color[1] == align_color[3]:
        if abs(int(arrow_min[0]) - int(arrow_min[1])) < 5:
            PCR_count = PCR_count + 1
        if abs(int(arrow_max[0]) - int(arrow_max[1])) < 5:
            PCR_count = PCR_count + 1
        if align_color[0] != "N":
            if abs(int(align_min[0]) - int(align_min[2])) < 5:
                PCR_count = PCR_count + 1
            if abs(int(align_max[0]) - int(align_max[2])) < 5:
                PCR_count = PCR_count + 1
        if align_color[1] != "N":
            if abs(int(align_min[1]) - int(align_min[3])) < 5:
                PCR_count = PCR_count + 1
            if abs(int(align_max[1]) - int(align_max[3])) < 5:
                PCR_count = PCR_count + 1
if PCR_count >= 3:
    Inspection = Inspection + "PCR duplicate? "

## Overlapping alignment ##
if align_color[0] != "N" and int(align_max[0]) < (6020 + 159):
    if align_color[2] != "N" and int(align_max[2]) < (6020 + 159):
        # Overlapping #
        if int(align_max[0]) > int(align_min[2]) and int(align_min[0]) < int(align_max[2]) and abs(max(int(align_min[0]),int(align_min[2])) - min(int(align_max[0]),int(align_max[2]))) > 20:
            overlap_notEqual_count = 0
            for k in range(max(int(align_min[0]),int(align_min[2])),min(int(align_max[0]),int(align_max[2]))+1):
                if temp_tensor[:,17:21,k].equal(temp_tensor[:,43:47,k]):
                    overlap_notEqual_count = overlap_notEqual_count + 1
            if (1 - overlap_notEqual_count/(min(int(align_max[0]),int(align_max[2])) - max(int(align_min[0]),int(align_min[2])) + 1)) > 0.10:
                Inspection = Inspection + "Overlap alignment is too low; "
    if align_color[3] != "N" and int(align_max[3]) < (6020 + 159):
        if int(align_max[0]) > int(align_min[3]) and int(align_min[0]) < int(align_max[3]) and abs(max(int(align_min[0]),int(align_min[3])) - min(int(align_max[0]),int(align_max[3]))) > 20:
            overlap_notEqual_count = 0
            for k in range(max(int(align_min[0]),int(align_min[3])),min(int(align_max[0]),int(align_max[3]))+1):
                if temp_tensor[:,17:21,k].equal(temp_tensor[:,56:60,k]):
                    overlap_notEqual_count = overlap_notEqual_count + 1
            if (1 - overlap_notEqual_count/(min(int(align_max[0]),int(align_max[3])) - max(int(align_min[0]),int(align_min[3])) + 1)) > 0.10:
                Inspection = Inspection + "Overlap alignment is too low; "
if align_color[1] != "N" and int(align_max[1]) < (6020 + 159):
    if align_color[2] != "N" and int(align_max[2]) < (6020 + 159):
        # Overlapping #
        if int(align_max[1]) > int(align_min[2]) and int(align_min[1]) < int(align_max[2]) and abs(max(int(align_min[1]),int(align_min[2])) - min(int(align_max[1]),int(align_max[2]))) > 20:
            overlap_notEqual_count = 0
            for k in range(max(int(align_min[1]),int(align_min[2])),min(int(align_max[1]),int(align_max[2]))+1):
                if temp_tensor[:,30:34,k].equal(temp_tensor[:,43:47,k]):
                    overlap_notEqual_count = overlap_notEqual_count + 1
            if (1 - overlap_notEqual_count/(min(int(align_max[1]),int(align_max[2])) - max(int(align_min[1]),int(align_min[2])) + 1)) > 0.10:
                Inspection = Inspection + "Overlap alignment is too low; "
    if align_color[3] != "N" and int(align_max[3]) < (6020 + 159):
        if int(align_max[1]) > int(align_min[3]) and int(align_min[1]) < int(align_max[3]) and abs(max(int(align_min[1]),int(align_min[3])) - min(int(align_max[1]),int(align_max[3]))) > 20:
            overlap_notEqual_count = 0
            for k in range(max(int(align_min[1]),int(align_min[3])),min(int(align_max[1]),int(align_max[3]))+1):
                if temp_tensor[:,30:34,k].equal(temp_tensor[:,56:60,k]):
                    overlap_notEqual_count = overlap_notEqual_count + 1
            if (1 - overlap_notEqual_count/(min(int(align_max[1]),int(align_max[3])) - max(int(align_min[1]),int(align_min[3])) + 1)) > 0.10:
                Inspection = Inspection + "Overlap alignment is too low; "

## read1 and read2 mapability ##
if align_mapability[0] != "N" and abs(int(align_max[0]) - int(align_min[0])) > 20 and int (align_max[0]) < (6020 + 159):
    if float(align_mapability[0]) < 0.90:
        Inspection = Inspection + "read1 mappability is lower than 90%; "
if align_mapability[1] != "N" and abs(int(align_max[1]) - int(align_min[1])) > 20 and int (align_max[1]) < (6020 + 159):
    if float(align_mapability[1]) < 0.90:
        Inspection = Inspection + "read1 mappability is lower than 90%; "
if align_mapability[2] != "N" and abs(int(align_max[2]) - int(align_min[2])) > 20 and int (align_max[2]) < (6020 + 159):
    if float(align_mapability[2]) < 0.90:
        Inspection = Inspection + "read2 mappability is lower than 90%; "
if align_mapability[3] != "N" and abs(int(align_max[3]) - int(align_min[3])) > 20 and int (align_max[3]) < (6020 + 159):
    if float(align_mapability[3]) < 0.90:
        Inspection = Inspection + "read2 mappability is lower than 90%; "


if Inspection == "":
    Inspection = "PASS"

## ACA & TAG ##
align1 = "5927-5929: NNN 6010-6012: NNN"
if Detect_black[17:21,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1].numel() != 0:
    for m in range(Detect_black[17:21,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1].numel()):
        tmp1 = Detect_black[17:21,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1][m].item()
        tmp2 = Detect_black[17:21,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[0][m].item()
        align1 = align1[0:tmp1+11] + DNA_base[tmp2] + align1[tmp1+12:]
if Detect_black[17:21,(6010+159):(6012+159)].eq(1).nonzero(as_tuple=True)[1].numel() != 0:
    for m in range(Detect_black[17:21,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[1].numel()):
        tmp1 = Detect_black[17:21,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[1][m].item()
        tmp2 = Detect_black[17:21,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[0][m].item()
        align1 = align1[0:tmp1+26] + DNA_base[tmp2] + align1[tmp1+27:]
# align2 #
align2 = "5927-5929: NNN 6010-6012: NNN"
if Detect_black[30:34,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1].numel() != 0:
    for m in range(Detect_black[30:34,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1].numel()):
        tmp1 = Detect_black[30:34,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1][m].item()
        tmp2 = Detect_black[30:34,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[0][m].item()
        align2 = align2[0:tmp1+11] + DNA_base[tmp2] + align2[tmp1+12:]
if Detect_black[30:34,(6010+159):(6012+159)].eq(1).nonzero(as_tuple=True)[1].numel() != 0:
    for m in range(Detect_black[30:34,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[1].numel()):
        tmp1 = Detect_black[30:34,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[1][m].item()
        tmp2 = Detect_black[30:34,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[0][m].item()
        align2 = align2[0:tmp1+26] + DNA_base[tmp2] + align2[tmp1+27:]
# align3 #
align3 = "5927-5929: NNN 6010-6012: NNN"
if Detect_black[43:47,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1].numel() != 0:
    for m in range(Detect_black[43:47,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1].numel()):
        tmp1 = Detect_black[43:47,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1][m].item()
        tmp2 = Detect_black[43:47,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[0][m].item()
        align3 = align3[0:tmp1+11] + DNA_base[tmp2] + align3[tmp1+12:]
if Detect_black[43:47,(6010+159):(6012+159)].eq(1).nonzero(as_tuple=True)[1].numel() != 0:
    for m in range(Detect_black[43:47,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[1].numel()):
        tmp1 = Detect_black[43:47,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[1][m].item()
        tmp2 = Detect_black[43:47,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[0][m].item()
        align3 = align3[0:tmp1+26] + DNA_base[tmp2] + align3[tmp1+27:]
# align4 #
align4 = "5927-5929: NNN 6010-6012: NNN"
if Detect_black[56:60,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1].numel() != 0:
    for m in range(Detect_black[56:60,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1].numel()):
        tmp1 = Detect_black[56:60,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[1][m].item()
        tmp2 = Detect_black[56:60,(5927+159):(5930+159)].eq(1).nonzero(as_tuple=True)[0][m].item()
        align4 = align4[0:tmp1+11] + DNA_base[tmp2] + align4[tmp1+12:]
if Detect_black[56:60,(6010+159):(6012+159)].eq(1).nonzero(as_tuple=True)[1].numel() != 0:
    for m in range(Detect_black[56:60,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[1].numel()):
        tmp1 = Detect_black[56:60,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[1][m].item()
        tmp2 = Detect_black[56:60,(6010+159):(6013+159)].eq(1).nonzero(as_tuple=True)[0][m].item()
        align4 = align4[0:tmp1+26] + DNA_base[tmp2] + align4[tmp1+27:]

ACA1 = align1.split()[1]
ACA2 = align2.split()[1]
ACA3 = align3.split()[1]
ACA4 = align4.split()[1]
TAG1 = align1.split()[3]
TAG2 = align2.split()[3]
TAG3 = align3.split()[3]
TAG4 = align4.split()[3]
ACA = ""
TAG = ""
ACAset = list(set([ACA1,ACA2,ACA3,ACA4]))
TAGset = list(set([TAG1,TAG2,TAG3,TAG4]))
for tmp_ACA in ACAset:
    if tmp_ACA == "ACA" or tmp_ACA == "ACG":
        ACA = tmp_ACA
        break
    elif tmp_ACA != "NNN":
            if ACA == "":
                ACA = tmp_ACA
            else:
                ACA = ACA + "/" + tmp_ACA
if ACA == "":
    ACA = "NNN"
for tmp_TAG in TAGset:
    if tmp_TAG == "TAG":
        TAG = tmp_TAG
        break
    elif tmp_TAG != "NNN":
            if TAG == "":
                TAG = tmp_TAG
            else:
                TAG = TAG + "/" + tmp_TAG
if TAG == "":
    TAG = "NNN"

f = open("Inspection_tmp_LINE_{}.txt".format(ver),"w")
f.write("{}\t{}\t{}\t{}\t".format(Inspection,ACA,TAG,PolyA_len))

## Microhomology ##
if reads_direction[0] == "upstream":
    if clip_map[0] == "mapped":
        f.write("microhomology: " + LINE_ref[(int(align_min[0])-170):(int(align_min[0])-150)] + "; ")
    else:
        if reads_direction[1] == "upstream":
            if clip_map[1] == "mapped":
                f.write("microhomology: " + LINE_ref[(int(align_min[2])-170):(int(align_min[2])-150)] + " ")
else:
    if reads_direction[1] == "upstream":
            if clip_map[1] == "mapped":
                f.write("microhomology: " + LINE_ref[(int(align_min[2])-170):(int(align_min[2])-150)] + " ")
f.close()

